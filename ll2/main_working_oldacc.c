#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <stdio.h>
#include <portaudio.h>
#include <string.h>
#include <pthread.h>
#include <getopt.h>
#include <sndfile.h>
#include <fftw3.h>
#include <assert.h>

#include "pqueue.h"
#include "direct.h"
#include "pa_ringbuffer.h"
#include "pa_util.h"

#define SAMPLE_RATE (44100)
#define FRAMESPERBLOCK (256)
#define BUFBLOCKS (16)

pthread_mutex_t pmt;
pthread_mutex_t acc;

int usl=1024; // 512 FPB

SF_INFO		sfinfo_fir1 ;

static int 	verbose_flag;
const char* version="0.1alpha";
char 		szFilter1[1024];
int 		iFilter1=0;
pthread_t 	ptSaw;
void* 		pFir1;
pqueue_t* 	pqTaskQueue;

static void* svTmpBuffer1;
static void* svTmpBuffer2;
volatile int terminate;

PaError err;
PaStream *stream;

/*
 * Convolution Task Structure
 *
 * TODO: Move all buffers out of the kernel structure into the task structure?
 */
typedef struct 	convoTaskStruct {
	int 			blockStart;
	int 			blockCount;
	//int			blockDelay;
	//int 			bufCount;
	int				blockDest;
	int 			firBlock;

	CFormat* 		signal;
	//fftwf_complex*  signalFFT;
	//fftwf_complex*  convoFFT;
	//CFormat* 		convo;
	
	// note this is for the benefit of the priority queue library
	pqueue_pri_t 	priority;
	size_t 			position;
} convoTask;

/*
 * fix up these random array pointers and manage the memory more securely.
 */
typedef struct strPpKernelStruct {

	CFormat* 			ppKernel;
	CFormat* 			ppRingBufferData;
	CFormat*			ppDirectBuffer; // 2 blocks for managing the data from the direct convolution

	CFormat*			accOut;

    PaUtilRingBuffer    convoRB;

	int 				nBlocks;
	int 				nBlocksRoot;


	CFormat*			SignalRBData[4096];
    PaUtilRingBuffer    SignalRBAry [4096];

	CFormat*			ResultRBData[4096];
    PaUtilRingBuffer    ResultRBAry [4096];

	CFormat*			AccumulateRBData[4096];
    PaUtilRingBuffer    AccumulateRBAry [4096];

	CFormat*			pFilterBlockSrc[100]; 	// pointer to ppKernel for each block starting point.
	CFormat*			pFilterBlock[100]; 		// pointer to ppKernel for each block starting point.
	CFormat*			pSignalBlock[100]; 		// pointer to ppKernel for each block starting point.
	CFormat*			pConvoBlock[100]; 		// pointer to ppKernel for each block starting point.

	size_t				blocksPerSegment[100]; // number of frames within for each filter block.

	fftwf_complex*		pFilterFFTW[100];
	fftwf_complex*		pSignalFFTW[100];
	fftwf_complex*		pConvoFFTW[100];

	fftwf_plan			filterFFTBlockPlans[100];
	fftwf_plan			signalFFTBlockPlans[100];
	fftwf_plan			convoFFTBlockPlans[100];

	convoTask*			cTasks[100];

} strPpKernelTd;

static strPpKernelTd *strPpKernel;
	
static __attribute__((noinline)) unsigned next_pow2(unsigned x)
{
	if (x <= 2) return x;
	
	return (1ULL << 32) >> __builtin_clz(x - 1);
}

typedef struct paTestDataStruct
{
	float left_phase;
	float right_phase;
	int framesPerBuffer;

    /* Ring buffer (FIFO) for "communicating" towards audio callback */
    PaUtilRingBuffer    rBufToRT;
    void*               rBufToRTData;

    /* Ring buffer (FIFO) for "communicating" from audio callback */
    PaUtilRingBuffer    rBufFromRT;
    void*               rBufFromRTData;

} paTestData;  


/*
 * Comparison functions for the Priority Queue
 */
int
fPQPriCmp(pqueue_pri_t next, pqueue_pri_t curr) {
	return (next > curr);
}

void 
fPQPriSet(void * element, pqueue_pri_t priority) {
	convoTask *ct=(convoTask*) element;
	ct->priority=priority;
}

pqueue_pri_t
fPQPriGet(void * element) {
	convoTask *ct=(convoTask*) element;
	return ct->priority;
}

size_t
fPQPosGet(void * element) {
	convoTask *ct=(convoTask*) element;
	return ct->position;
}

void
fPQPosSet(void * element, size_t position) {
	convoTask *ct=(convoTask*) element;
	ct->position=position;
}

/*
 * Schedule a convolution task
 *
 * @start: 	 Time at which the task should be scheduled
 * @hblock:  Block of the fir filter to be convovolved
 * @nblocks: Number of blocks to be convovolved (always +1 output due overlap-add)
 * @delay:   Insertion delay
 * @pri:     Task priority
void
scheduleTask(convoTask* ct, int start, int bufcounter, int hblock, int nblocks, int delay, int priority) {
	//printf("schedule task start=%u,bufcounter=%u,firBlock=%u,nblocks=%u,delay=%u, priority=%u\n",start,bufcounter,hblock,nblocks,delay,priority);
	//convoTask* ct=malloc(sizeof(convoTask));
	ct->priority=priority;
	ct->blockStart=start;
	ct->blockCount=nblocks;
	ct->blockDelay=delay;
	ct->bufCount=bufcounter;
	ct->firBlock=hblock;
	//ct->signal=(CFormat*)malloc(sizeof(CFormat)*FRAMESPERBLOCK*ct->blockCount);
	ring_buffer_size_t a = PaUtil_ReadRingBuffer(&strPpKernel->SignalRBAry[ct->firBlock],ct->signal,ct->blockCount);

	if (a!=ct->blockCount) 
		printf("schedule task ringbuffer empty. tried to read %u blocks\n",ct->blockCount);

	//pthread_mutex_lock(&pmt);
	pqueue_insert(pqTaskQueue,(void*)ct);
	//pthread_mutex_unlock(&pmt);
}
 */
/*
 */ 
void
scheduleTaskPair(convoTask* ct[], int h, int bufcounter, int j, int i) {

	if (verbose_flag) printf("STP: h=%d, bufcounter=%d, j=%d, i=%d \n",h,bufcounter,j,i);	
	bufcounter=(bufcounter+1) % strPpKernel->nBlocks;
	
	convoTask* ct1=ct[h+0];
	convoTask* ct2=ct[h+1];

	ct1->priority=i*(h+0);
	ct1->blockStart=((bufcounter + strPpKernel->nBlocks) - j -1 ) % strPpKernel->nBlocks;
	ct1->blockCount=j;
	ct1->blockDest=(bufcounter+(j*1) -1) % strPpKernel->nBlocks;
	ct1->firBlock=h;

	ct2->priority=i*(h+1);
	ct2->blockStart=((bufcounter + strPpKernel->nBlocks) - j -1 ) % strPpKernel->nBlocks;
	ct2->blockCount=j;
	ct2->blockDest=(bufcounter+(j*2) -1) % strPpKernel->nBlocks;
	ct2->firBlock=h+1;

	ring_buffer_size_t a = PaUtil_ReadRingBuffer(&strPpKernel->SignalRBAry[ct1->firBlock],ct1->signal,ct1->blockCount);
	ct2->signal=NULL;

	//if (a!=ct1->blockCount) 
	//	printf("schedule task ringbuffer empty. tried to read %u blocks\n",ct1->blockCount);

	//pthread_mutex_lock(&pmt);
	pqueue_insert(pqTaskQueue,(void*)ct1);
	pqueue_insert(pqTaskQueue,(void*)ct2);
	//pthread_mutex_unlock(&pmt);
}

/*
 * this is the thread function that manages the convolution tasks.
 */
int
fConvolutionTaskProcessor(void* data) {

	void* pqd;

		/*
		 * 1. get a task from the queue
		 *
		 * 2. might need to mask these with a mutex/critical section
		 */
		pqd=pqueue_pop(pqTaskQueue);

		if (pqd==NULL) {

			return 0;

		} else {

			convoTask* ct=(convoTask*)pqd;

			int firBlock=ct->firBlock;

			int blockSize=(2 * FRAMESPERBLOCK * strPpKernel->blocksPerSegment[ct->firBlock]); // two blocks per buffer
			int blockSizeCpx=2 * ((blockSize/2)+1); // note the FFTW buffers are smaller

			/*
			 * log task pulled for this sequence 
			 */
			if (verbose_flag) 
			printf(
				"dqed task : priority=%llu,blockStart=%u,blockCount=%u,blockDest=%u,firBlock=%u\n",
					ct->priority, ct->blockStart, ct->blockCount, ct->blockDest, ct->firBlock);

			if (ct->signal==NULL) {

				// if ct->signal is NULL then we are reusing the complex FFT data from a previous calculation
				firBlock--;

			} else {

				/*
				 * grab the signal block from the task and put it in the temp memory block
				 *
				 * the FFT input block is 2 times normal size so we clear twice the volume. memset deals with bytes.
				 * we only fill half of it when we put the signal in.
				 */
				memset(strPpKernel->pSignalBlock[ct->firBlock],0,			sizeof(CFormat) * 2 * FRAMESPERBLOCK * strPpKernel->blocksPerSegment[ct->firBlock] );
				memcpy(strPpKernel->pSignalBlock[ct->firBlock],ct->signal,	sizeof(CFormat) * 	  FRAMESPERBLOCK * strPpKernel->blocksPerSegment[ct->firBlock] );
	
				/*
				 * if the ct->firblock is odd, then run the fft over it R2C
				 * if the ct->firblock is even then we can just grab the last odd one for this range.
				 */
				fftwf_execute(strPpKernel->signalFFTBlockPlans[ct->firBlock]);
			
			}

			/*
			 * multiply the FFT results to make the convolution
			 */
			fftwf_complex *a,*b,*c;

			for (int i=0 ; i < blockSizeCpx ; i+=2) {
		
				for (int j=0;j<2;j++) {
		
					a=strPpKernel->pSignalFFTW[firBlock]; // <- this will only ever be an ODD number because we don't want to repeat work.
					b=strPpKernel->pFilterFFTW[ct->firBlock];
					c=strPpKernel->pConvoFFTW [ct->firBlock];

					a+=i+j; b+=i+j; c+=i+j;
			
					(*c)[0]=(
						 ((*a)[0] * (*b)[0])
						-((*a)[1] * (*b)[1])
					);
			
					(*c)[1]=(
						 ((*a)[1] * (*b)[0])
						+((*a)[0] * (*b)[1])
					);
				}
			}

			/*
			 * run the inverse fft over it to get back to time domain.
			 */
			fftwf_execute(strPpKernel->convoFFTBlockPlans[ct->firBlock]);

			/*
			 * scale the result back
			 */

			float scale=1.0/(FRAMESPERBLOCK * strPpKernel->blocksPerSegment[ct->firBlock]);

			for (int i=0;i<2 * FRAMESPERBLOCK * strPpKernel->blocksPerSegment[ct->firBlock];i++) {
				(*(strPpKernel->pConvoBlock[ct->firBlock]+i)).left*=scale;
				(*(strPpKernel->pConvoBlock[ct->firBlock]+i)).right*=scale;
			}

			/*
			 * push the results into the accumulator ring buffer
			 *
			 * doing this in a loop block by block. becauase of the accumulator ringbuffer being vertical for each time slot. only one block per slot please.
			 *
			 */
			 int j=0; int i=ct->blockDest;	int bb=strPpKernel->blocksPerSegment[ct->firBlock];
			 while(j<2*ct->blockCount) {
		//	 	if (verbose_flag) printf("writing into acc buffer position %u size %u result block %u - ct->blockCount=%u\n",i,bb,j,ct->blockCount);
				PaUtil_WriteRingBuffer(&strPpKernel->AccumulateRBAry[i],strPpKernel->pConvoBlock[ct->firBlock]+(j*FRAMESPERBLOCK),1);
				i=((i+1+strPpKernel->nBlocks)) % strPpKernel->nBlocks; j++;
			 }
		}
	return 1;
}

/*
 * this is the direct convolution routine
 *
 * this goes right back into the output ringbuffer in the interrupt routine.
 */
void convodirect(CFormat* out, CFormat* in) {
	size_t n;

   	for (n = 0; n < 2*FRAMESPERBLOCK- 1; n++) {
        size_t kmin, kmax, k;

        out[n].left 	= 0;
        out[n].right 	= 0;

       	kmin = (n >= FRAMESPERBLOCK - 1) ? n - (FRAMESPERBLOCK - 1) : 0;
       	kmax = (n <	 FRAMESPERBLOCK - 1) ? n :  FRAMESPERBLOCK - 1;

        for (k = kmin; k <= kmax; k++) {
        	out[n].left 	+= in[k].left 	* strPpKernel->ppKernel[n - k].left;
        	out[n].right 	+= in[k].right 	* strPpKernel->ppKernel[n - k].right;
   		}
	}
}

static paTestData data;
int bufcounter=0;
int bufFirstRun=1;

/* 
 * This routine will be called by the PortAudio engine when audio is needed.
 * It may called at interrupt level on some machines so don't do anything
 * that could mess up the system like calling malloc() or free().
 */ 
static int 
paTestCallback( 
	const void *inputBuffer, 
	void *outputBuffer,
	unsigned long framesPerBuffer,
	const PaStreamCallbackTimeInfo* timeInfo,
	PaStreamCallbackFlags statusFlags,
	void *userData ) {

	CFormat* rbBuf;

	assert(framesPerBuffer == FRAMESPERBLOCK);
	if (verbose_flag) {
	if (statusFlags & paInputUnderflow) {
		printf("callback underflow\n");
	}
	if (statusFlags & paInputOverflow) {
		printf("callback overflow\n");
	}
	if (statusFlags & paOutputUnderflow) {
		printf("callback output underflow\n");
	}
	if (statusFlags & paOutputOverflow) {
		printf("callback output overflow\n");
	}
	}
					  
	/* Cast data passed through stream to our structure. */
	paTestData *data = (paTestData*)userData;
	CFormat *out = (CFormat*)outputBuffer;
	CFormat *in = (CFormat*)inputBuffer;
	
	/* clear the output buffer */
	memset(out,0,sizeof(CFormat) * FRAMESPERBLOCK); //num channels

	/*
	 * Task Scheduler
	 * Very Simple
	 */
	if (verbose_flag) printf("current bufcounter=%u---------------------------------------------------------------------------------------------------------------------------------------------------------------------\n",bufcounter);

	/*
	 * perform the direct convolution and then push the result onto the accumulator buffers
	 *
	 * two blocks of CFormat
	 */
	memset(strPpKernel->ppDirectBuffer,0,sizeof(CFormat) * FRAMESPERBLOCK * 2);
	convodirect(strPpKernel->ppDirectBuffer,in); 
	PaUtil_WriteRingBuffer(&strPpKernel->AccumulateRBAry[bufcounter],strPpKernel->ppDirectBuffer,1);
	PaUtil_WriteRingBuffer(&strPpKernel->AccumulateRBAry[(bufcounter+1) % strPpKernel->nBlocks],strPpKernel->ppDirectBuffer+FRAMESPERBLOCK,1);

	/*
	for (int kk=0;kk<20;kk++) {

		printf("[%d]leftDir[%d]=%f\n",-1,kk,strPpKernel->ppDirectBuffer[kk].left);
		printf("[%d]rightDir[%d]=%f\n",-1,kk,strPpKernel->ppDirectBuffer[kk].right);

	}
	*/

	/*
	 * mix in some of the dry signal before we start mixing in the convolution.
	 */
	float scale=1.0/1;
	for (int i=0;i<FRAMESPERBLOCK;i++) {
		(*(out+i)).left=(*(in+i)).left * scale;
		(*(out+i)).right=(*(in+i)).right * scale;
	}

	/*
	 * nBlocksRoot identifies the number of levels we need to go to for
	 * the convolution layers and lengths (p129)
	 *
	 * nBlocksRoot includes all levels of the filter blocks. 
	 * 
	 * This segment ignores the direct convo and only looks at level pairs, so we lose one and divide by two. nBlocksRoot should be odd.
	 */
	for (int i=0;i<(strPpKernel->nBlocksRoot-1)/2;i++) {

		int h=(i*2)+1; int j=(int)pow(2,i);

		/*
		 * if we hit a zero here we are going to schedule a task
		 */
		if ( bufcounter % j ==0 ) { // this is for debugging - only scheduling tasks for h=0,1,2
			if (!bufFirstRun) {
				scheduleTaskPair(strPpKernel->cTasks, h, bufcounter, j, i); 
			}
		}
	}

	int tCnt=0;
	for (int tsk=0;tsk<4;tsk++)
	{
		if (fConvolutionTaskProcessor(NULL)) {
			tCnt++;
		}
	}
	if (verbose_flag) printf("did %d tasks\n",tCnt);

	for (int i=0;i<(strPpKernel->nBlocksRoot-1)/2;i++) {

		int h=(i*2)+1; 

		/*
		 * Push the blocks onto the ringbuffers after we have done the tasks
		 * Only push odd numbers on the ring buffer as they are duplicates
		 */
		ring_buffer_size_t a = PaUtil_WriteRingBuffer(&strPpKernel->SignalRBAry[h+0],in,1);
		ring_buffer_size_t b = PaUtil_WriteRingBuffer(&strPpKernel->SignalRBAry[h+1],in,1);
	}

	/*
	 * bufcounter contains the index of the accumulator ringbuffer that we need to
	 * flush and output.
	 */
	 memset(strPpKernel->accOut,0,sizeof(CFormat) * FRAMESPERBLOCK);
	 int accCounter=0;
	 while (1) {
	 	ring_buffer_size_t rbst=PaUtil_ReadRingBuffer(&strPpKernel->AccumulateRBAry[bufcounter],strPpKernel->accOut,1); //FRAMESPERBLOCK
		if (rbst == 0) break;
		assert(rbst == 1);
		for (int kk=0;kk<FRAMESPERBLOCK;kk++) {  // or FRAMESPERBLOCK
			out[kk].left+=strPpKernel->accOut[kk].left;
			out[kk].right+=strPpKernel->accOut[kk].right;
			/*
			if (kk<10) {
				printf("[%d]left[%d]=%f\n",accCounter,kk,strPpKernel->accOut[kk].left);
				printf("[%d]right[%d]=%f\n",accCounter,kk,strPpKernel->accOut[kk].right);
			}
			*/
		}
		accCounter++;
	 }
	 if ((accCounter<14) && verbose_flag)
	 	printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>accCounter=%u\n",accCounter);

	 /* 
	  * increment and loop the bufcounter.
	  */
	 bufcounter=(bufcounter+1) % strPpKernel->nBlocks;
	 bufFirstRun=0;

	 return 0;
}

/*
 * Here's where we analyze and prepare the kernel for convolution
 * 
 * @ *pFir1 is the pointer to the kernel buffer
 * @ fir1_buffer_len is the number of kernel frames available
 * @ FramesPerBlock is the number of frames per processing block
 * @ nBufferBlocks is the number of buffer blocks we are working with
 * @ CFormat *** ppKernel is the output variable to store the array of kernel pointers
 */
strPpKernelTd*
prepareKernel(
	CFormat * pFir1, 
	size_t fir1_buffer_len) 
	{

	/*
	 * allocate memory for the kernel structure
	 */
	strPpKernelTd* ppKernel=malloc(sizeof(strPpKernelTd));

	// one temp block for output thread.
	ppKernel->accOut=PaUtil_AllocateMemory(sizeof(CFormat) * FRAMESPERBLOCK);

	/*
	 * work out the nearest power of two that will cover the
	 * entire length of the kernel at the current block settings
	 */ 
	int nKernelBlocks=fir1_buffer_len/FRAMESPERBLOCK;
	if (verbose_flag) printf("FramesPerBlock=%d\n",FRAMESPERBLOCK);
	if (verbose_flag) printf("nKernelBlocks=%d\n",nKernelBlocks); //actual blocks in the filter from the file
	if (verbose_flag) printf("fir1_buffer_len=%zu\n",fir1_buffer_len); //frames in the filter from the file

	/*
	 * calculate size of the filter block - next best power of two.
	 */

	int i=0; long int k=2; int fBlock=0;

	while (1) {
		k+=2*(pow(2,i));
		if (k>nKernelBlocks) break; 
		else i++;
	}
	fBlock=(2*i)+2;

	if (verbose_flag) printf("Kernel Prep: i=%u,k=%ld\n",fBlock,k);

	/*
	 * Allocate a memory block that can be then divided up
	 * poke this back into the ppKernel location
	 */
	ppKernel->ppKernel=(CFormat*)malloc(sizeof(CFormat)*k*FRAMESPERBLOCK);	// allocate the memory block for the kernel
	memset(ppKernel->ppKernel,0,		sizeof(CFormat)*k*FRAMESPERBLOCK);				// clear the block
	memcpy(ppKernel->ppKernel,pFir1,	sizeof(CFormat)*fir1_buffer_len);				// copy in the kernel to the start of the block

	ppKernel->nBlocks=k;
	ppKernel->nBlocksRoot=fBlock+1; // are we counting h[0]? if not, lose the +1

	if (verbose_flag) printf("................ fblock=%u,i=%u\n",fBlock,i);

	/*
	 * Calculate the pointers for the filter blocks
	 */
	i=0; k=2; fBlock=0; int fbLen=0;

	if (verbose_flag) printf("Kernel Pointers: h=%u,start=%u,length=%u\n",i,0,2); 

	while (1) {

		fbLen=pow(2,i);

		fBlock=(2*i)+1;

		ppKernel->pFilterBlockSrc[fBlock]=ppKernel->ppKernel+(FRAMESPERBLOCK*k);
		//ppKernel->filterBlockFrames[fBlock]=fbLen*FRAMESPERBLOCK;
		ppKernel->blocksPerSegment[fBlock]=fbLen;

		if (verbose_flag) printf("Kernel Pointers: h=%u,start=%lu,length=%u\n",fBlock,k,fbLen); 

		k+=fbLen;

		fBlock=(2*i)+2;

		ppKernel->pFilterBlockSrc[fBlock]=ppKernel->ppKernel+(FRAMESPERBLOCK*k);
		//ppKernel->filterBlockFrames[fBlock]=fbLen*FRAMESPERBLOCK;
		ppKernel->blocksPerSegment[fBlock]=fbLen;

		 if (verbose_flag) printf("Kernel Pointers: h=%u,start=%lu,length=%u\n",fBlock,k,fbLen); 

		k+=fbLen;

		if (k>nKernelBlocks) break; 
		else i++;
	}
	
	/*
	 * prepare an array of ringbuffers.
	 * we are not creating a ringbuffer for h[0] so we might need to shift everything back in the array.
	 */
	for (int ri=1;ri<=fBlock;ri++) {

		if (verbose_flag) printf("Building tasks, rings and buffers for ri=%u\n",ri);

		/*
		 * build the tasks for each level/block of the impulse response filter.
		 */ 
		convoTask* ct=(convoTask*) malloc(sizeof(convoTask));
		ct->signal=(CFormat*)malloc(sizeof(CFormat) * FRAMESPERBLOCK * ppKernel->blocksPerSegment[ri]);
		memset(ct->signal,0,sizeof(CFormat) * FRAMESPERBLOCK * ppKernel->blocksPerSegment[ri]);

		ppKernel->cTasks[ri]=ct;

		/*
		 * Initialize the Signal ringbuffer. This is the ringbuffer that contains the signal information that is pushed in
		 * for each level of the filter i/r during the sample callbackroutine.
		 */
		ppKernel->SignalRBData[ri]=PaUtil_AllocateMemory(sizeof(CFormat) * FRAMESPERBLOCK * next_pow2(ppKernel->blocksPerSegment[ri]));
		memset(ppKernel->SignalRBData[ri],0,sizeof(CFormat) * FRAMESPERBLOCK * ppKernel->blocksPerSegment[ri]);
		PaUtil_InitializeRingBuffer(&ppKernel->SignalRBAry[ri],sizeof(CFormat) * FRAMESPERBLOCK,next_pow2(ppKernel->blocksPerSegment[ri]),ppKernel->SignalRBData[ri]);

		/*
		 * Initialize the Result ringbffers. This is where the result gets pushed once the convolution is done, and prior to 
		 * pushing onto the accumulator ringbuffers.
		 */
		ppKernel->ResultRBData[ri]=PaUtil_AllocateMemory(sizeof(CFormat) * FRAMESPERBLOCK * next_pow2(ppKernel->blocksPerSegment[ri]));
		memset(ppKernel->ResultRBData[ri],0,sizeof(CFormat) * FRAMESPERBLOCK * next_pow2(ppKernel->blocksPerSegment[ri]));
		PaUtil_InitializeRingBuffer(&ppKernel->ResultRBAry[ri],sizeof(CFormat) * FRAMESPERBLOCK,next_pow2(ppKernel->blocksPerSegment[ri]),ppKernel->ResultRBData[ri]);

		/*
		 * calculate the blocksize in frames for each of the calculation buffers
		 */
		//int blockSize=(2*ppKernel->filterBlockFrames[ri])-1; // not -1 just for convenience here
		int blockSize=(2 * FRAMESPERBLOCK * ppKernel->blocksPerSegment[ri]); // two blocks per buffer  - number in frames

		/*
		 * initialize the FFT buffers for each filter block
		 */
		ppKernel->pFilterFFTW[ri]=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * 2 * ((blockSize/2)+1)); // note the FFTW buffers are smaller
		ppKernel->pSignalFFTW[ri]=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * 2 * ((blockSize/2)+1));
		ppKernel->pConvoFFTW[ri] =(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * 2 * ((blockSize/2)+1));

		memset(ppKernel->pFilterFFTW[ri],0,(sizeof(fftwf_complex) * 2 * ((blockSize/2)+1)));
		memset(ppKernel->pSignalFFTW[ri],0,(sizeof(fftwf_complex) * 2 * ((blockSize/2)+1)));
		memset(ppKernel->pConvoFFTW[ri], 0,(sizeof(fftwf_complex) * 2 * ((blockSize/2)+1)));

		/*
		 * initialize the input buffers (REAL) prior to the FFT
		 */
		ppKernel->pFilterBlock[ri] =(CFormat*) fftwf_malloc(sizeof(CFormat) * blockSize);
		ppKernel->pSignalBlock[ri] =(CFormat*) fftwf_malloc(sizeof(CFormat) * blockSize);
		ppKernel->pConvoBlock [ri] =(CFormat*) fftwf_malloc(sizeof(CFormat) * blockSize);

		memset(ppKernel->pFilterBlock[ri],0,(sizeof(CFormat) * blockSize));
		memset(ppKernel->pSignalBlock[ri],0,(sizeof(CFormat) * blockSize));
		memset(ppKernel->pConvoBlock [ri],0,(sizeof(CFormat) * blockSize));

		// copy the block of the filter into the filter blocks
		memcpy(ppKernel->pFilterBlock[ri], ppKernel->pFilterBlockSrc[ri], sizeof(CFormat) * FRAMESPERBLOCK * ppKernel->blocksPerSegment[ri]);
	}

	/*
	 * create FFTW Plans for each block in the filter.
	 */
	for (int ri=1;ri<=fBlock;ri++) {

		if (verbose_flag) printf("planning ri=%u - blocks per segment = %zu\n",ri,ppKernel->blocksPerSegment[ri]);

		int n[]={(2 * FRAMESPERBLOCK * ppKernel->blocksPerSegment[ri])-1}; // blocksize TODO Check this is considered number of frames not number of samples?
		
		ppKernel->filterFFTBlockPlans[ri] =
			fftwf_plan_many_dft_r2c(1,n,2,
           	                    	(float *)ppKernel->pFilterBlock[ri],NULL,
									2,1,
									(fftwf_complex*)ppKernel->pFilterFFTW[ri],NULL,
									2,1,
           	                    	FFTW_ESTIMATE);
									
		if (verbose_flag) printf("done ri=%u - filter\n",ri);

		ppKernel->signalFFTBlockPlans[ri] =
			fftwf_plan_many_dft_r2c(1,n,2,
           	                    	(float *)ppKernel->pSignalBlock[ri],NULL,
									2,1,
									(fftwf_complex*)ppKernel->pSignalFFTW[ri],NULL,
									2,1,
           	                    	FFTW_ESTIMATE);

		if (verbose_flag) printf("done ri=%u - signal\n",ri);

		ppKernel->convoFFTBlockPlans[ri] =
			fftwf_plan_many_dft_c2r(1,n,2,
									(fftwf_complex*)ppKernel->pConvoFFTW[ri],NULL,
									2,1,
           	                    	(float *)ppKernel->pConvoBlock[ri],NULL,
									2,1,
           	                    	FFTW_ESTIMATE);

		if (verbose_flag) printf("done ri=%u - convo\n",ri);

		memset((void*)ppKernel->pFilterFFTW[ri],0,n[0]);
		if (verbose_flag) printf("running filter plan ri=%u \n",ri);
		fftwf_execute(ppKernel->filterFFTBlockPlans[ri]);
		if (verbose_flag) printf("done running filter plan ri=%u \n",ri);

	}

	for (int ii=0;ii<ppKernel->nBlocks;ii++) {
		//if (verbose_flag) printf("acc ringbuf ii=%u,nblox=%u\n",ii,next_pow2(ppKernel->nBlocksRoot));

		ppKernel->AccumulateRBData[ii]=PaUtil_AllocateMemory(	sizeof(CFormat) * FRAMESPERBLOCK * next_pow2(ppKernel->nBlocksRoot));
		memset(ppKernel->AccumulateRBData[ii],0,				sizeof(CFormat) * FRAMESPERBLOCK * next_pow2(ppKernel->nBlocksRoot));
		PaUtil_InitializeRingBuffer(&ppKernel->AccumulateRBAry[ii],sizeof(CFormat) * FRAMESPERBLOCK, next_pow2(ppKernel->nBlocksRoot),ppKernel->AccumulateRBData[ii]);

	}

	/*
	 * Create the buffer that manages the direct convolution routine.
	 * 
	 * We need a buffer that is 2 blocks long.
	 */
	 ppKernel->ppDirectBuffer = PaUtil_AllocateMemory(sizeof(CFormat) * FRAMESPERBLOCK * 2);

	 return ppKernel;
}
	
void
memCleanup() {
}

void
intHandler(int dummy) {
	memCleanup();
	Pa_AbortStream(stream);
	Pa_Terminate();
}



int
main(int argc,char** argv) {

	signal(SIGINT,intHandler);


	printf("convolution reverb low latency\n");

	int c,index,aflag=0,bflag=0;

	memset(szFilter1,0,1024);

	while (1) {
		static struct option long_options[] =
		{
			{"verbose",	no_argument, 		&verbose_flag, 	1},
			{"filter",	required_argument, 	0, 				'f'},
			//{"filter2",	required_argument, 	0, 				'g'},
			//{"input",	required_argument,	0,				'i'},
			//{"output",	required_argument,	0,				'o'},
			{"help",	no_argument,		0,				'h'},
			{0,0,0,0}
		};
		int option_index=0;
		c=getopt_long(argc,argv,"f:g:h",long_options,&option_index);
		if (c==-1)
			break;
		switch(c) {
			case 0:
				if (long_options[option_index].flag !=0)
					break;
				printf("option %s", long_options[option_index].name);
				if (optarg)
					printf(" with arg %s", optarg);
				printf("\n");
				break;
			case 'f':
				strncpy(szFilter1,optarg,1024);
				break;
	/*		case 'g':
				strncpy(szFilter2,optarg,1024);
				break;
			case 'i':
				strncpy(szInput,optarg,1024);
				break;
			case 'o':
				strncpy(szOutput,optarg,1024);
				break; */
			case 'h':
				//usage();
				exit(0);
				break;
			default:
				break;
		}
	};

	if (verbose_flag)
		;//entitled();

	// check incoming options

	int bail=0;

	if (strlen(szFilter1)==0) {
		printf("no filter supplied (-f)\n");
		bail=1;
	}

	if (bail)
		exit(1);

	/* A SNDFILE is very much like a FILE in the Standard C library. The
	** sf_open function return an SNDFILE* pointer when they sucessfully
	** open the specified file.
	*/
	SNDFILE	*infile_fir1;

	/* A pointer to an SF_INFO struct is passed to sf_open.
	** On read, the library fills this struct with information about the file.
	** On write, the struct must be filled in before calling sf_open.
	*/

	int		readcount ;
	int		rc_fir1 ;

	const char	*fir1 = szFilter1;

	/* The SF_INFO struct must be initialized before using it.
	*/
	memset (&sfinfo_fir1, 0, sizeof (sfinfo_fir1)) ;

	/*
	 * load in the interleaved FIR impulse response filters
	 * prior to convolution.
	 */

	if (! iFilter1) {
		if (! (infile_fir1 = sf_open (fir1, SFM_READ, &sfinfo_fir1)))
		{	/* Open failed so print an error message. */
			printf ("Not able to open fir1 file %s.\n", fir1) ;
			/* Print the error message from libsndfile. */
			puts (sf_strerror (NULL)) ;
			return 1 ;
		} ;
		if (verbose_flag)
			;//sfinfo_print(&sfinfo_fir1,"FIR 1");
		if (sfinfo_fir1.channels != 2)
		{	printf ("FIR1 filter needs 2 channels\n") ;
			return 1 ;
		} ;
	}

	/*
	 * Identify the length of the FIR filter buffers
	 */

	int fir1_buffer_len=0;
	
	int dat_buffer_len=16384;


	if (!iFilter1) {
		fir1_buffer_len=sfinfo_fir1.frames;
	} else {
		fir1_buffer_len=1024;
	}

	/*
	 * allocate a buffer and read in the entire kernel 
	 */
	pFir1=malloc(fir1_buffer_len*sizeof(CFormat));
	sf_readf_float(infile_fir1,pFir1,fir1_buffer_len);

	/* 
	 * prepare the kernel and chop it up into blocks
	 */
	strPpKernel=prepareKernel(pFir1,fir1_buffer_len);

	free (pFir1);																// free up the original buffer
	if (verbose_flag) {
		printf("strPpKernel->nBlocks=%u\n",strPpKernel->nBlocks);
		printf("strPpKernel->nBlocksRoot=%u\n",strPpKernel->nBlocksRoot);
		printf("FRAMESPERBLOCK=%d\n",FRAMESPERBLOCK);
	}
	
	/* 
	 * Build a temporary buffer for holding output data
	 */
    //svTmpBuffer1 = PaUtil_AllocateMemory(sizeof(CFormat) * FRAMESPERBLOCK * BUFBLOCKS);
    //svTmpBuffer2 = PaUtil_AllocateMemory(sizeof(CFormat) * FRAMESPERBLOCK * BUFBLOCKS);

	//strPpKernel->ppRingBufferData=PaUtil_AllocateMemory(sizeof(CFormat) * strPpKernel->nBlocks * FRAMESPERBLOCK);

    //PaUtil_InitializeRingBuffer(&strPpKernel->convoRB, sizeof(CFormat), strPpKernel->nBlocks * FRAMESPERBLOCK, strPpKernel->ppRingBufferData);

	/*
	 * set up the Priority task queue
	 */
	pqTaskQueue=pqueue_init(
		32,
		fPQPriCmp,
		fPQPriGet,
		fPQPriSet,
		fPQPosGet,
		fPQPosSet);

	/*
	 * Reset our termination flag. This gets set to '1' if we want to bail on the
	 * processing thread
	 */
	terminate=0;

	pthread_mutex_init(&pmt,NULL);
	pthread_mutex_init(&acc,NULL);

	err=Pa_Initialize();
	if (err != paNoError)  goto error;

	err=Pa_OpenDefaultStream( &stream,
							  2, //input channels
							  2, //output channels
							  paFloat32,
							  SAMPLE_RATE,
							  FRAMESPERBLOCK,
							  paTestCallback,
							  &data );

	/*
	 * 2. split kernel into multiple buffer blocks
	 * 3. build convolution engine
	 */

	err=Pa_StartStream(stream);
	if (err != paNoError) goto error;

	//pthread_create(&ptSaw,NULL,fConvolutionTaskProcessor,NULL);
	Pa_Sleep(60000);
	terminate=1;
	printf("before join\n");

	err=Pa_StopStream(stream);
	if (err != paNoError) goto error;
	printf("after stop\n");

	pthread_join(ptSaw,NULL);
	printf("after join\n");

	// flush pq queue
	// destory queue

	while (1) {
		//pthread_mutex_lock(&pmt);
		void* pqdt=pqueue_pop(pqTaskQueue);
		//pthread_mutex_unlock(&pmt);
		if (pqdt==NULL) break;
		convoTask* ct=(convoTask*) pqdt;
		free(ct->signal);
		free(ct);
	}
	pqueue_free(pqTaskQueue);

	//free(data.rBufToRTData);
	//free(data.rBufFromRTData);
	printf("after free1\n");

	free(svTmpBuffer1);
	free(svTmpBuffer2);
	printf("after free2\n");


	err=Pa_Terminate();
	if (err != paNoError) goto error;
	pthread_mutex_destroy(&pmt);
	pthread_mutex_destroy(&acc);

	exit(0);

	error:
	printf("PortAudio Error: %s \n",Pa_GetErrorText(err));
	exit(1);

}
