#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "pipe.h"
#include <pthread.h>

/*
 * we are working on multiple blocks all with size N
 *
 *
 *
 */

pthread_t ptAccumulator;
pthread_t ptProducer;

pipe_t *pAccQueue;

int iBlockSize; // size in frames of the output accumulator buffer blocks 
int ABlockSize;// raw byte size of a block.
int iNBlocks; // number of accumulator buffer blocks

struct CFormatStruct {
	double A;
	double B;
	double C;
	double D;
};

typedef struct CFormatStruct CFormat;

struct qOutElementStruct {
	int iBlockTarget; // where will this element be stored in the output buffer?
	int iBlockNumbers; // how many blocks does the element span? (min=2) even a 1 block input will cover 2 blocks due OLAP-ADD
	CFormat* buffer; // this buffer must be allocated prior to referencing, and will be free'ed as soon as the block has been transferred.
};

typedef struct qOutElementStruct qOutElement;
CFormat* cfAccumulatorBuffer;

//
// this thread function pushes data into the accumulator buffer from the queue.
// other threads push blocks into the queue from FFT and Direct Convolution
// routines.
//
void*
fAccumulator(void* pPipe) {
	int buf;
	int rc;
	while (1) {
		pipe_consumer_t* cons = pipe_consumer_new((pipe_t*)pPipe);
		rc=pipe_pop(cons,&buf,1);
		printf("got %u from the queue %u\n",rc,buf);
		pipe_consumer_free(cons);
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
	int rc=0;
	while (1) {
		int r=rand();
		pipe_push(sConvoBuilderPtr->ppHandle,&r,1);	
		printf("pushed %u\n",r);
		usleep(100000);
	}
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
	free(buffer);
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

int main(int argc,char **argv) {

	iBlockSize=2048;
	iNBlocks=8;
	
	ABlockSize=sizeof(CFormat)*iBlockSize;


	cfAccumulatorBuffer=malloc(ABlockSize*iNBlocks);

	// seed the random number thingy
	srand(time(NULL));

	printf("low latency convolution processor\n");

	pAccQueue = pipe_new(sizeof(int),8);

	int rc=0;
	
	pipe_consumer_t* pPipeConsumer=pipe_consumer_new(pAccQueue);
	rc=pthread_create (&ptAccumulator, NULL, fAccumulator,(void*)pPipeConsumer);

	if (rc) {
		printf("ERROR creating accumulator %d\n",rc);
		exit(1);
	}

	/*
	 * instantiate a ConvoBuilder structure
	 * and build a convolution thread.
	 */
	sConvoBuilder sCB[iNBlocks];

	for (int i=0;i<iNBlocks;i++) {
		sCB[i].nBlocks=1; // 1 block
		sCB[i].iBlockTarget=i; // where will it put the data?
		sCB[i].ppHandle=pipe_producer_new(pAccQueue);

		rc=pthread_create (&ptProducer, NULL, fConvolution,(void*)&sCB[i]);
	}
	

	if (rc) {
		printf("ERROR creating accumulator %d\n",rc);
		exit(1);
	}


	pipe_free (pAccQueue);

	printf("threads created and running..... \n");
	pthread_join(ptAccumulator,NULL);
	free (cfAccumulatorBuffer);
	exit(0);
}
