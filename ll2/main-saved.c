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

#include "gspacc.h"

#define SAMPLE_RATE (44100)
#define FRAMESPERBLOCK (64)

pthread_mutex_t pmt;
pthread_mutex_t acc;

int usl=1024; // 512 FPB

SF_INFO		sfinfo_fir1 ;

static int 	verbose_flag;
const char* version="0.1alpha";
char 		szFilter1[1024];
int 		iFilter1=0;
pthread_t 	ptConvoTask;
void* 		pFir1;
pqueue_t* 	pqTaskQueue;

volatile int terminate;
volatile int bufcounter;

PaError err;
PaStream *stream;

gspacc_t * pGspAcc;

float lrMax=4.0f;

typedef struct convoTaskStruct convoTask;

typedef struct taskSplitStruct {
	int 			firBlock;
	int 			blockDest;
	pthread_t		pt;

	fftwf_complex*	filterFFTW; //TODO malloc
	fftwf_complex*	convoFFTW;
	CFormat*		convo;

	fftwf_plan		convoPlan;

	convoTask* 		cTask; // pointer to parent
} taskSplit;
/*
 * Convolution Task Structure
 *
 * TODO: Move all buffers out of the kernel structure into the task structure?
 */
typedef struct 	convoTaskStruct {
	int 			firBlock;
	int 			blockStart;
	int 			blockCount;

	CFormat* 		signal;

	fftwf_complex*	signalFFTW;

	fftwf_plan		signalPlan;

	taskSplit		task1;
	taskSplit		task2;
	
	// note this is for the benefit of the priority queue library
	pqueue_pri_t 	priority;
	size_t 			position;

	clock_t			cStart;
	clock_t			cEnd;
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

	CFormat*			SignalRBData[16384];
    PaUtilRingBuffer    SignalRBAry [16384];

	CFormat*			pFilterBlockSrc[512]; 	// pointer to ppKernel for each block starting point.
	CFormat*			pFilterBlock[512]; 		// pointer to ppKernel for each block starting point.
	//CFormat*			pSignalBlock[512]; 		// pointer to ppKernel for each block starting point.
	//CFormat*			pConvoBlock[512]; 		// pointer to ppKernel for each block starting point.

	size_t				blocksPerSegment[512]; // number of frames within for each filter block.

	fftwf_complex*		pFilterFFTW[512];
	//fftwf_complex*		pSignalFFTW[512];
	//fftwf_complex*		pConvoFFTW[512];

	fftwf_plan			filterFFTBlockPlans[512];
	//fftwf_plan			signalFFTBlockPlans[512];
	//fftwf_plan			convoFFTBlockPlans[512];

	convoTask*			cTasks[512];

} strPpKernelTd;

static strPpKernelTd *strPpKernel;

typedef struct paTestDataStruct
{
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

int firBlockMap(int firBlock) {
	return floor(firBlock/2);
	
}
void
scheduleTaskPair(convoTask* ct[], int h, int bufcounter, int j, int i) {

	if (verbose_flag) printf("STP: h=%d, bufcounter=%d, j=%d, i=%d \n",h,bufcounter,j,i);	

	bufcounter=(bufcounter+1) % strPpKernel->nBlocks;
	
	convoTask* ctp=ct[h+0];
	
	ctp->priority=i*h; //TODO move to initialization
	ctp->blockStart=((bufcounter + strPpKernel->nBlocks) - j - 1 ) % strPpKernel->nBlocks;
	ctp->blockCount=j; //TODO move to initialization
	ctp->firBlock=h; 	// TODO move to initialization
	ctp->cStart=clock();

	ctp->task1.blockDest=(bufcounter+(j*1) - 1) % strPpKernel->nBlocks;
	ctp->task2.blockDest=(bufcounter+(j*2) - 1) % strPpKernel->nBlocks;
	ctp->task1.firBlock=h; 	// TODO - move this to initialization
	ctp->task2.firBlock=h+1;// TODO - move this to initialization

	ctp->task1.cTask=ctp;
	ctp->task2.cTask=ctp;

	memset(ctp->signal,0,sizeof(CFormat) * FRAMESPERBLOCK * strPpKernel->blocksPerSegment[h] * 2);
	ring_buffer_size_t c = PaUtil_ReadRingBuffer(&strPpKernel->SignalRBAry[h],ctp->signal,ctp->blockCount);

	pqueue_insert(pqTaskQueue,(void*)ctp);

	/*
	convoTask* ct1=ct[h+0];
	convoTask* ct2=ct[h+1];

	ct1->priority=i*(h+0);
	ct1->blockStart=((bufcounter + strPpKernel->nBlocks) - j -1 ) % strPpKernel->nBlocks;
	ct1->blockCount=j;
	ct1->blockDest=(bufcounter+(j*1) -1) % strPpKernel->nBlocks;
	ct1->firBlock=h;
	ct1->cStart=clock();

	ct2->priority=i*(h+1);
	ct2->blockStart=((bufcounter + strPpKernel->nBlocks) - j -1 ) % strPpKernel->nBlocks;
	ct2->blockCount=j;
	ct2->blockDest=(bufcounter+(j*2) -1) % strPpKernel->nBlocks;
	ct2->firBlock=h+1;
	ct2->cStart=clock();

	ring_buffer_size_t a = PaUtil_ReadRingBuffer(&strPpKernel->SignalRBAry[ct1->firBlock],ct1->signal,ct1->blockCount);
	ring_buffer_size_t b = PaUtil_ReadRingBuffer(&strPpKernel->SignalRBAry[ct2->firBlock],ct2->signal,ct2->blockCount);

	pqueue_insert(pqTaskQueue,(void*)ct1);
	pqueue_insert(pqTaskQueue,(void*)ct2);
	*/
}

void  *fConvolutionTaskProcessor (void *);

int
fConvolutionTasker(void* data) {
	void * pqd;
	pqd=pqueue_pop(pqTaskQueue);

	if (pqd==NULL) {
		//printf("fConvolutionTasker pqd=NULL\n");
		return 0;
	} else {
		convoTask * ct=(convoTask*)pqd;
		if (ct->blockCount>=128) {
			pthread_create(&ptConvoTask,NULL,fConvolutionTaskProcessor,(void*)ct);
			return 1;
		} else {
			fConvolutionTaskProcessor((void*)ct);
			return 1;
		}
	}
}

void*
fConvolutionTaskMultiplyJoin(void * tSplit) {

	taskSplit * pTask=(taskSplit*)tSplit;
	int blockSize=(2 * FRAMESPERBLOCK * strPpKernel->blocksPerSegment[pTask->firBlock]); // two blocks per buffer
	int blockSizeCpx=2 * ((blockSize/2)+1); // note the FFTW buffers are smaller

	/*
	 * multiply the FFT results to make the convolution
	 */
	fftwf_complex *a,*b,*c;

	for (int i=0 ; i < blockSizeCpx ; i+=2) {

		for (int j=0;j<2;j++) {

			a=pTask->cTask->signalFFTW;
			b=pTask->filterFFTW;
			c=pTask->convoFFTW;

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
	fftwf_execute(pTask->convoPlan);

	/*
	 * scale the result back
	 */

	float scale=1.0/(FRAMESPERBLOCK * strPpKernel->blocksPerSegment[pTask->firBlock]);
	int g=0;
	for (;g<2 * FRAMESPERBLOCK * strPpKernel->blocksPerSegment[pTask->firBlock];g++) {
		(*(pTask->convo+g)).left*=scale;
		(*(pTask->convo+g)).right*=scale;
	}
		g--; // this handles the -1 sample case with the second block of a circular convolution. Set it to zero and it won't affect the accumulation
		(*(pTask->convo+g)).left*=0.0f;
		(*(pTask->convo+g)).right*=0.0f;

	 gspacc_write(pGspAcc,(void *) pTask->convo,pTask->blockDest,2 * pTask->cTask->blockCount); 
	 
	 pTask->cTask->cEnd=clock();

	 if ((pTask->blockDest+(2*pTask->cTask->blockCount) % strPpKernel->nBlocks) > (bufcounter + strPpKernel->nBlocks)) printf("**TOO LATE! at bufcounter%u destination=%u, blockCount=%u\n",bufcounter,pTask->blockDest,pTask->cTask->blockCount);

	 //printf("task duration = %lu\n",ct->cEnd-ct->cStart);

	return 0;
}

/*
 * this is the thread function that manages the convolution tasks.
 */
void *
fConvolutionTaskProcessor(void* pct) {

			convoTask * ct=(convoTask*) pct;

			//int firBlock=ct->firBlock;
			int blockSize=(2 * FRAMESPERBLOCK * strPpKernel->blocksPerSegment[ct->firBlock]); // two blocks per buffer
			int blockSizeCpx=2 * ((blockSize/2)+1); // note the FFTW buffers are smaller

			/*
			 * log task pulled for this sequence 
			 */
			if (verbose_flag) 
			printf(
				"dqed task : priority=%llu,blockStart=%u,blockCount=%u,blockDest1=%u,blockDest2=%u,firBlock=%u\n",
					ct->priority, ct->blockStart, ct->blockCount, ct->task1.blockDest, ct->task2.blockDest, ct->firBlock);

			fftwf_execute(ct->signalPlan);

			if (ct->blockCount>=128) {
				pthread_create(&ct->task1.pt,NULL,fConvolutionTaskMultiplyJoin,&ct->task1);
				pthread_create(&ct->task2.pt,NULL,fConvolutionTaskMultiplyJoin,&ct->task2);
			} else {
				fConvolutionTaskMultiplyJoin(&ct->task1);
				fConvolutionTaskMultiplyJoin(&ct->task2);
			}


	return NULL;
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
int bufFirstRun=1;

void accCallback(void * data, size_t elementSize, void * dest) {
	CFormat * cfData=(CFormat*)data;
	CFormat * cfDest=(CFormat*)dest;
	int i=0;
	while (1) {
		cfDest[i].left+=cfData[i].left;	
		cfDest[i].right+=cfData[i].right;	
		/*
		if (fabsf(cfDest[i].left)>lrMax) {
			lrMax=fabsf(cfDest[i].left);

		}
		if (fabsf(cfDest[i].right)>lrMax) {
			lrMax=fabsf(cfDest[i].right);
		}
		*/
		if (++i * sizeof(CFormat) >= elementSize) break;
	}
}

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
	if (verbose_flag) printf("current bufcounter=%u---------------------------------------------------------------------------------------------\n",bufcounter);

	/*
	 * perform the direct convolution and then push the result onto the accumulator buffers
	 *
	 * two blocks of CFormat
	 */
	memset(strPpKernel->ppDirectBuffer,0,sizeof(CFormat) * FRAMESPERBLOCK * 2);

	convodirect(strPpKernel->ppDirectBuffer,in); 
	
	gspacc_write(pGspAcc,(void*)strPpKernel->ppDirectBuffer,bufcounter,2); // 2

	/*
	 * mix in some of the dry signal before we start mixing in the convolution.
	 */
	float scale=1.0/4;
	for (int i=0;i<FRAMESPERBLOCK;i++) {
		(*(out+i)).left =(*(in+i)).left  * scale;
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
		if ( bufcounter % j ==0 & h<99) { // this is for debugging - only scheduling tasks for h=0,1,2
			if (!bufFirstRun) {
				scheduleTaskPair(strPpKernel->cTasks, h, bufcounter, j, i); 
			}
		}
	}

	int tCnt=0;
	for (int tsk=0;tsk<5;tsk++)
	{
		if (fConvolutionTasker(NULL)) {
			tCnt++;
		}
	}
	if (verbose_flag) printf("bufcount %d did %d tasks\n",bufcounter,tCnt);

	// save the input blocks into horizontal ringbuffers for all the filter blocks.
	for (int i=1;i<=(strPpKernel->nBlocksRoot-1)/2;i++) {

		int h=(i*2)+1; 

		/*
		 * Push the blocks onto the ringbuffers after we have done the tasks
		 * Only push odd numbers on the ring buffer as they are duplicates
		 *
		 * TODO: sort out what we are going to do with even ringbuffers. they are duplicates.
		 */
		ring_buffer_size_t a = PaUtil_WriteRingBuffer(&strPpKernel->SignalRBAry[h],in,1);
		//ring_buffer_size_t b = PaUtil_WriteRingBuffer(&strPpKernel->SignalRBAry[h+1],in,1);

	}

	// grab the accumulator for this block and process it to the output buffer.

	int accCounter=gspacc_read(pGspAcc, bufcounter, &accCallback, outputBuffer);

	if ((accCounter<14) && verbose_flag)
		printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>accCounter=%u\n",accCounter);

	/* 
	 * increment and loop the bufcounter.
	 */
	bufcounter=(bufcounter+1) % strPpKernel->nBlocks;
	bufFirstRun=0;

	return paContinue;
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
		ppKernel->blocksPerSegment[fBlock]=fbLen;

		if (verbose_flag) printf("Kernel Pointers: h=%u,start=%lu,length=%u\n",fBlock,k,fbLen); 

		k+=fbLen;

		fBlock=(2*i)+2;

		ppKernel->pFilterBlockSrc[fBlock]=ppKernel->ppKernel+(FRAMESPERBLOCK*k);
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

		if (verbose_flag) printf("Building rings and buffers for ri=%u\n",ri);

		int blockSize=(2 * FRAMESPERBLOCK * ppKernel->blocksPerSegment[ri]); // two blocks per buffer  - number in frames

		/*
		 * Initialize the Signal ringbuffer. This is the ringbuffer that contains the signal information that is pushed in
		 * for each level of the filter i/r during the sample callbackroutine.
		 */
		ppKernel->SignalRBData[ri]=PaUtil_AllocateMemory(sizeof(CFormat) * FRAMESPERBLOCK * next_pow2(ppKernel->blocksPerSegment[ri]));
		memset(ppKernel->SignalRBData[ri],0,sizeof(CFormat) * FRAMESPERBLOCK * ppKernel->blocksPerSegment[ri]);
		PaUtil_InitializeRingBuffer(&ppKernel->SignalRBAry[ri],sizeof(CFormat) * FRAMESPERBLOCK,next_pow2(ppKernel->blocksPerSegment[ri]),ppKernel->SignalRBData[ri]);

		/*
		 * initialize the input buffers (REAL) prior to the FFT
		 */
		ppKernel->pFilterBlock[ri] =(CFormat*) fftwf_malloc(sizeof(CFormat) * blockSize);
		memset(ppKernel->pFilterBlock[ri],0,(sizeof(CFormat) * blockSize));

		/*
		 * initialize the FFT buffers for each filter block
		 */
		ppKernel->pFilterFFTW[ri]=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * 2 * ((blockSize/2)+1)); // note the FFTW buffers are smaller
		memset(ppKernel->pFilterFFTW[ri],0,(sizeof(fftwf_complex) * 2 * ((blockSize/2)+1)));

		// copy the block of the filter into the filter blocks
		memcpy(ppKernel->pFilterBlock[ri], ppKernel->pFilterBlockSrc[ri], sizeof(CFormat) * FRAMESPERBLOCK * ppKernel->blocksPerSegment[ri]);
	}

	for (int ri=1;ri<ppKernel->nBlocksRoot;ri+=2) {

		printf("building tasks: task %u\n",ri);
		
		/*
		 * calculate the blocksize in frames for each of the calculation buffers
		 */
		//int blockSize=(2*ppKernel->filterBlockFrames[ri])-1; // not -1 just for convenience here
		int blockSize=(2 * FRAMESPERBLOCK * ppKernel->blocksPerSegment[ri]); // two blocks per buffer  - number in frames
		int n[]={(2 * FRAMESPERBLOCK * ppKernel->blocksPerSegment[ri])-1}; // blocksize TODO Check this is considered number of frames not number of samples?

		/*
		 * build the task pairs for each level/block of the impulse response filter.
		 */ 
		convoTask* ct=(convoTask*) malloc(sizeof(convoTask));

		ct->signal=		     (CFormat*)	fftwf_malloc(sizeof(CFormat) * blockSize);
		ct->signalFFTW=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * 2 * ((blockSize/2)+1));

		memset(ct->signal,0,sizeof(CFormat) * blockSize);
		memset(ct->signalFFTW,0,sizeof(fftwf_complex) * 2 * ((blockSize/2)+1));

		ct->task1.convoFFTW=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * 2 * ((blockSize/2)+1));
		ct->task1.convo=(CFormat*)malloc(sizeof(CFormat) * blockSize);

		ct->task2.convoFFTW=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * 2 * ((blockSize/2)+1));
		ct->task2.convo=(CFormat*)malloc(sizeof(CFormat) * blockSize);

		ct->signalPlan =
			fftwf_plan_many_dft_r2c(1,n,2,
           	                    	(float *)ct->signal,NULL,
									2,1,
									(fftwf_complex*)ct->signalFFTW,NULL,
									2,1,
           	                    	FFTW_MEASURE);

		ct->task1.convoPlan =
			fftwf_plan_many_dft_c2r(1,n,2,
									(fftwf_complex*)ct->task1.convoFFTW,NULL,
									2,1,
           	                    	(float *)ct->task1.convo,NULL,
									2,1,
           	                    	FFTW_MEASURE);

		ct->task2.convoPlan =
			fftwf_plan_many_dft_c2r(1,n,2,
									(fftwf_complex*)ct->task2.convoFFTW,NULL,
									2,1,
           	                    	(float *)ct->task2.convo,NULL,
									2,1,
           	                    	FFTW_MEASURE);

		ct->task1.filterFFTW=ppKernel->pFilterFFTW[ri];
		ct->task2.filterFFTW=ppKernel->pFilterFFTW[ri+1];

		ppKernel->cTasks[ri]=ct;
	}
	
	/*
	 * create and execute filter FFTW Plans for each block in the filter.
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
           	                    	FFTW_MEASURE);
									
		if (verbose_flag) printf("done ri=%u - filter\n",ri);

		memset((void*)ppKernel->pFilterFFTW[ri],0,n[0]);
		if (verbose_flag) printf("running filter plan ri=%u \n",ri);
		fftwf_execute(ppKernel->filterFFTBlockPlans[ri]);
		if (verbose_flag) printf("done running filter plan ri=%u \n",ri);

	}

	//
	// initialize accumulator
	// sizeof(CFormat) * FRAMESPERBLOCK = elementSize
	// ppKernel->nBlocks - Number of blocks in the accumulator
	// ppKernel->nBlocksRoot - depth
	//

	pGspAcc=gspacc_initialize((size_t)(sizeof(CFormat) * FRAMESPERBLOCK),ppKernel->nBlocks,ppKernel->nBlocksRoot+2);

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
	fftwf_init_threads();
	fftwf_plan_with_nthreads(4);
	fftwf_import_wisdom_from_filename(".fftwf");
	strPpKernel=prepareKernel(pFir1,fir1_buffer_len);
	fftwf_export_wisdom_to_filename(".fftwf");

	//abort();

	free (pFir1);																// free up the original buffer
	if (verbose_flag) {
		printf("strPpKernel->nBlocks=%u\n",strPpKernel->nBlocks);
		printf("strPpKernel->nBlocksRoot=%u\n",strPpKernel->nBlocksRoot);
		printf("FRAMESPERBLOCK=%d\n",FRAMESPERBLOCK);
	}
	

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

	//abort();

	//pthread_mutex_init(&pmt,NULL);
	//pthread_mutex_init(&acc,NULL);

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
	Pa_Sleep(120000);
	terminate=1;
	printf("before join\n");

	err=Pa_StopStream(stream);
	if (err != paNoError) goto error;
	printf("after stop\n");

	//pthread_join(ptSaw,NULL);
	printf("after join\n");

	// flush pq queue
	// destory queue

	while (1) {
		void* pqdt=pqueue_pop(pqTaskQueue);
		if (pqdt==NULL) break;
		convoTask* ct=(convoTask*) pqdt;
		free(ct->signal);
		free(ct->signalFFTW);
		free(ct);
	}

	pqueue_free(pqTaskQueue);

	gspacc_cleanup(pGspAcc);

	fftwf_cleanup_threads();
	printf("after gspacc_cleanup\n");

	err=Pa_Terminate();
	if (err != paNoError) goto error;
	pthread_mutex_destroy(&pmt);
	pthread_mutex_destroy(&acc);

	exit(0);

	error:
	printf("PortAudio Error: %s \n",Pa_GetErrorText(err));
	exit(1);

}
