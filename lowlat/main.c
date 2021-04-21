#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <fftw3.h>
#include "pipe.h"
#include <portaudio.h>
#include <time.h>

clock_t start, diff;


/*
 * we are working on multiple blocks all with size N
 *
 */

pthread_t ptAccumulator;
pthread_t *ptConvolution;

pipe_t *pAccQueue;

int iBlockSize=65536; // size in frames of the output accumulator buffer blocks 
int ABlockSize;	// raw byte size of a block.
int iNBlocks=16; // number of accumulator buffer blocks

struct CFormatStruct {
	double A;
	double B;
//	double C;
//	double D;
};

typedef struct CFormatStruct CFormat;

struct qOutElementStruct {
	int iBlockTarget; // where will this element be stored in the output buffer?
	int nBlocks; // how many blocks does the element span? (min=2) even a 1 block input will cover 2 blocks due OLAP-ADD
	CFormat* buffer; // this buffer must be allocated prior to referencing, and will be free'ed as soon as the block has been transferred.
};

typedef struct qOutElementStruct qOutElement;
CFormat* cfAccumulatorBuffer;
CFormat* cfInputBuffer;
CFormat* cfFilterBuffer;


void 		TransferToAccumulator(int startingBlock, int nBlocks, CFormat* buffer);
CFormat* 	getDelayTap(int iDelayTap) ;
void 		addInputBlock() ;

//
// this thread function pushes data into the accumulator buffer from the queue.
// other threads push blocks into the queue from FFT and Direct Convolution
// routines.
//
void*
fAccumulator(void* pPipeConsumer) {
	int rc;
	qOutElement qOE;

	while (1) {
		rc=pipe_pop((pipe_consumer_t*)pPipeConsumer,&qOE,1);
		TransferToAccumulator(qOE.iBlockTarget,qOE.nBlocks,qOE.buffer);
		//printf("Pulled %u block from the queue for segment %u\n",rc,qOE.iBlockTarget);
		//free(qOE.buffer);
	}
	pthread_exit(NULL);
	return 0;
}

struct sConvoBuilder {
	int nBlocks; // how many blocks will this convolution processor handle?
	int iBlockTarget; // where will this convolution processor put its data?
	pipe_producer_t* ppHandle; // handle to pipe Producer structure
};

typedef struct sConvoBuilder sConvoBuilder;

void*
fConvolution(void* sCB) {
	sConvoBuilder* sConvoBuilderPtr=(sConvoBuilder*)sCB;
	qOutElement qOE;
	qOE.buffer=getDelayTap(sConvoBuilderPtr->iBlockTarget);
	memset(qOE.buffer,sConvoBuilderPtr->iBlockTarget,ABlockSize*sConvoBuilderPtr->nBlocks);

	qOE.iBlockTarget=sConvoBuilderPtr->iBlockTarget;
	qOE.nBlocks=sConvoBuilderPtr->nBlocks;

	pipe_push(sConvoBuilderPtr->ppHandle,&qOE,1);	

	//printf("Convolution pushed something into Queue for segment %u\n",sConvoBuilderPtr->iBlockTarget);
	usleep(10000);
	pthread_exit(NULL);
	return 0;
}

void
TransferToAccumulator(int startingBlock, int nBlocks, CFormat* buffer)
{
	int offset1=0;	
	int offset2=0;	
	int BlockOffset1=0;
	int BlockOffset2=0;

	int startingBlockOffset=startingBlock*iBlockSize;

	for (int i=0;i<nBlocks;i++) {
		BlockOffset2=i*iBlockSize;
		BlockOffset1=startingBlockOffset+BlockOffset2;
		for (int j=0;j<iBlockSize;j++) {
			offset1=BlockOffset1+j;
			offset2=BlockOffset2+j;
			(cfAccumulatorBuffer+offset1)->A+=(buffer+offset2)->A;
			(cfAccumulatorBuffer+offset1)->B+=(buffer+offset2)->B;
			(cfAccumulatorBuffer+offset1)->C+=(buffer+offset2)->C;
			(cfAccumulatorBuffer+offset1)->D+=(buffer+offset2)->D;
			
		}
	}
}

/*
 * final consumer of the blocks.
 * pulls a block and plays it.
 * 
 * Once the block is consumed all the memory is freed.
 *
 */
void*
fGetBlock(void* pPipe) {
	int rc=0;
	
	while (1) {
		pipe_consumer_t *cons=pipe_consumer_new((pipe_t*)pPipe);
		//pipe_pop(cons)
		//printf("Pulled block %u from the accumulator and playing it.\n",iBlk);
		pipe_consumer_free(cons);
	}
}


/*
 * This gets called from the convolution threads
 * 
 * gotta put a mutex or something around iBlockCounter.
 *
 */
int iBlockCounter=0;
CFormat** inputBlocks;
int* iInputBlockFlags;
pthread_mutex_t	pmtBlockCounter;

/*
 * getDelayTap returns a block of newly allocated memory
 * containing input audio. iDelayTap=0..n refers to the 
 * requested block starting from the most recent.
 *
 * When a block is 'consumed' then a flag is set allowing it
 * to be replaced in sequence.
 */
CFormat* getDelayTap(int iDelayTap) {
	//printf("inside getDelayTap with iDelayTap=%u\n",iDelayTap);

	pthread_mutex_lock(&pmtBlockCounter);
	int iBlockID=(iDelayTap + iBlockCounter) % iNBlocks;
	pthread_mutex_unlock(&pmtBlockCounter);
	//printf("inside getDelayTap with iBlockID=%u\n",iBlockID);
	CFormat* transferBuffer=(CFormat*)malloc(ABlockSize);
	memcpy(transferBuffer,inputBlocks[iBlockID],ABlockSize);
	iInputBlockFlags[iBlockID]=0;
	int j=0;
	for (int i=0;i<iNBlocks;i++) {
		j+=iInputBlockFlags[i];
	}
	if (j>0) addInputBlock();
	return transferBuffer;
	
}

CFormat* getBlock() {
	return (CFormat*)malloc(ABlockSize);
}

/*
 * Copy a new block of audio into the input buffer
 * and advance the current block pointer.
 */
 
void addInputBlock() {
	pthread_mutex_lock(&pmtBlockCounter);
	
	//printf("--------------------------------------------------------adding a new block at %u\n",iBlockCounter);
	
	// copy the block into the input buffer
	// advance the block pointer around the ring
	CFormat *block=getBlock(); // pull a block off the filesystem.
	memcpy(inputBlocks[iBlockCounter],block,ABlockSize);
	free(block);

	// wrap the block counter around the ring if necessary
	if (++iBlockCounter>=iNBlocks)
		iBlockCounter=0;

	
	// reset all the block flags
	for (int i=0;i<iNBlocks;i++) {
		iInputBlockFlags[i]=1;
	}
	pthread_mutex_unlock(&pmtBlockCounter);
}

/*
 * 
 */
int 
main(int argc,char **argv) {

	pthread_mutex_init(&pmtBlockCounter,NULL);

	// This is the array of pointers that holds the incoming 
	// data blocks.
	//
	// we clear the datablocks to Zero so that we only return
	// useful data.
	inputBlocks=malloc(sizeof(CFormat*)*iNBlocks);
	memset((void*)inputBlocks,0,sizeof(CFormat*)*iNBlocks);

	// allocate thread handle space for convolution threads.
	ptConvolution=malloc(sizeof(pthread_t)*iNBlocks);

	ABlockSize=sizeof(CFormat)*iBlockSize;

	cfInputBuffer=(CFormat*)malloc(ABlockSize*iNBlocks);
	cfFilterBuffer=(CFormat*)malloc(ABlockSize*iNBlocks);

	cfAccumulatorBuffer=(CFormat*)malloc(ABlockSize*iNBlocks);

	/*
	 * iInputBlockFlags handles the input block ring buffer and
	 * advances the block pointer as blocks are consumed.
	 */
	iInputBlockFlags=(int*)malloc(sizeof(int)*iNBlocks);

	for (int i=0;i<iNBlocks;i++) {
		inputBlocks[i]=cfInputBuffer+(i*iBlockSize);
		iInputBlockFlags[i]=1;
	}
	
	// seed the random number thingy
	srand(time(NULL));

	printf("low latency convolution processor\n");

	pAccQueue = pipe_new(sizeof(qOutElement),0);

	int rc=0;
	
	pipe_consumer_t* pPipeConsumer=pipe_consumer_new(pAccQueue);
	rc=pthread_create (&ptAccumulator, NULL, fAccumulator,(void*)pPipeConsumer);

	if (rc) {
		printf("ERROR creating accumulator %d\n",rc);
		exit(1);
	}

	/*
	 * instantiate a ConvoBuilder structure
	 * and build a bunch of  convolution threads.
	 *
	 *

	 * TODO: put a loop around these commands and then use a pthread_join to wait for each one to finish.
	 * this may help with structure and thread management.
	 */
	sConvoBuilder sCB[iNBlocks];

	for (int i=0;i<iNBlocks;i++) {
		sCB[i].nBlocks=1; // 1 block
		sCB[i].iBlockTarget=i; // where will it put the data?
		sCB[i].ppHandle=pipe_producer_new(pAccQueue);
	}

	/*
	 * clear the pipe structure as we already have
	 * the producer and consumer structures.
	 */
	pipe_free (pAccQueue);

	/*
	 * Just keep doing this until....
	 */
	while (1) {

		/*
		 * Create a block of convolution threads and then wait for them to complete
		 */
		for (int i=0;i<iNBlocks;i++) {
			rc=pthread_create (ptConvolution+i, NULL, fConvolution,(void*)&sCB[i]);
		}

		/*
		 * wait for all the convolution threads to finish before starting another block of them.
		 */
		for (int i=0;i<iNBlocks;i++) {
			pthread_join(*(ptConvolution+i),NULL);
		}
		//printf("done convolution loop****************************\n");
		diff = clock() - start;
		start=clock();

		int msec = diff * 1000 / CLOCKS_PER_SEC;
		float rt=(float)iBlockSize/48000.0;
		float ratio=(rt*1000.0)/(float)msec;
		printf("Time taken %d seconds %d milliseconds to process %f (s) audio = %f times real time\n", msec/1000, msec%1000,rt,ratio);
		
	}
	

	if (rc) {
		printf("ERROR creating accumulator %d\n",rc);
		exit(1);
	}



	/*
	 * Wait on completed thread.
	 */
	printf("threads created and running..... \n");
	pthread_join(ptAccumulator,NULL);

	/*
	 * Cleanup and free memory
	 */

	pthread_mutex_destroy(&pmtBlockCounter);

	free (cfAccumulatorBuffer);
	free (cfInputBuffer);
	free (cfFilterBuffer);

	exit(0);
}
