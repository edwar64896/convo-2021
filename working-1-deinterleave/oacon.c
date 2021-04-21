/*
 *
 * works on structs of 4 double
 *
 * [ FIR = M ][L-1 Zeros---------] M+L-1=N
 * [ Data = L--------][M-1 Zeros ] L+M-1=N
 *
 *
 * TODO:
 * Notes - the foward and reverese FFT's work fine on their own, avoiding the convolution multiplcation.
 *
 *
 *
 *
 *
 *
 *
 *
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <fftw3.h>
#include <sndfile.h>
#include <math.h>

#include <oacon.h>


/*
 * initialize structure. we are dealing with set frame sizes (CFormat)
 *
 * blocksize = L+M-1
 */
Convo*
convo_init(int filterSize, int dataSize, const char* title) {

	printf("initializing convo with Filtersize=%u and Datasize=%u with Title=%s\n",filterSize,dataSize,title);

	Convo *handle;
	handle=(Convo*)malloc(sizeof(Convo));

	handle->dbgBlock=0;

#ifdef _DBG_CONVO
	handle->snd_debugfile=NULL;
	memset((void*)&handle->sfi_debugfile,0,sizeof(SF_INFO));
#endif

	handle->title=title;

	handle->filterSize=filterSize; //size in samples of the filter 	=M
	handle->dataSize=dataSize; //size in samples of a data block 	=L

	// this is the block size we work in during FFT/convolution. M+L-1
	handle->blockSize=(filterSize+dataSize)-1; // counted in Samples

	// real to complex transforms have a real buffer length of 'n' and a complex buffer size of floor(n/2)+1
	handle->blockSizeCpx=floor(handle->blockSize/2)+1;

	// this is the size of one FFT stream/channel in bytes
	handle->WbufSize=(size_t)(sizeof(fftw_complex)		* handle->blockSizeCpx); //bufsize in bytes

	// this is the size of one sample stream/channel in bytes
	handle->DbufSize=(size_t)(sizeof(double) 			* handle->blockSize); //bufsize in bytes

	// this is the size of one interleaved sample stream (4 channels) in bytes
	handle->CbufSize=(size_t)(sizeof(CFormat) 			* handle->blockSize); //bufsize in bytes

	// this is the size of one interleaved FFT stream (4 channels) in bytes
	handle->FbufSize=(size_t)(sizeof(fftw_complex) 	* 4	* handle->blockSizeCpx);

	// this is the size of one Overlap interleaved sample buffer
	handle->ObufSize=(size_t)(sizeof(CFormat) 			* (handle->filterSize-1));

	// allocation of FFT buffers (interleaved)
	handle->fftw_fir=	(fftw_complex*)	fftw_malloc(handle->FbufSize);
	handle->fftw_input=	(fftw_complex*)	fftw_malloc(handle->FbufSize);
	handle->fftw_output=(fftw_complex*)	fftw_malloc(handle->FbufSize);

	// clearing of FFT buffers we just allocated.
	memset((void*)handle->fftw_input,	0,handle->FbufSize);
	memset((void*)handle->fftw_output,	0,handle->FbufSize);

	// allocation of sample buffers (interleaved)
	handle->fir=		(CFormat*)		fftw_malloc(handle->CbufSize);
	handle->input=		(CFormat*)		fftw_malloc(handle->CbufSize);
	handle->output=		(CFormat*)		fftw_malloc(handle->CbufSize);
 
 	// allocation of Overlap buffer (interleaved)
	handle->olap=		(CFormat*)		fftw_malloc(handle->ObufSize);

	// allocation/init of non-interleaved channel buffers
	convo_initBufs(handle,&handle->firPtr);
	convo_initBufs(handle,&handle->inputPtr);
	convo_initBufs(handle,&handle->outputPtr);
	convo_initBufs(handle,&handle->olapPtr);

	// allocation/init of non-interleaved FFT buffers
	convo_initFftwBufs(handle,&handle->fftw_firPtr);
	convo_initFftwBufs(handle,&handle->fftw_inputPtr);
	convo_initFftwBufs(handle,&handle->fftw_outputPtr);

	// initialize FFTW threading engine
	//fftw_init_threads();
	//fftw_plan_with_nthreads(4);

	// build the FFT plans. 
	convo_initPlansR2C(handle,&handle->fftwBufPlan,&handle->inputPtr,&handle->fftw_inputPtr);
	convo_initPlansR2C(handle,&handle->fftwFilterPlan,&handle->firPtr,&handle->fftw_firPtr);
	//convo_initPlansC2R(handle,&handle->fftwConvPlan,&handle->fftw_inputPtr,&handle->outputPtr);
	convo_initPlansC2R(handle,&handle->fftwConvPlan,&handle->fftw_outputPtr,&handle->outputPtr);

	//fftw_import_wisdom_from_filename(".fftw");
	/*
	handle->fftwBufPlan = fftw_plan_dft_r2c_2d(
								handle->blockSize, 4, 
                               	(double *)handle->input, 
								(fftw_complex*)handle->fftw_input,
                               	FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
	handle->fftwFilterPlan = fftw_plan_dft_r2c_2d(
								handle->blockSize, 4,
                               	(double *)handle->fir, 
								(fftw_complex*)handle->fftw_fir,
                               	FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
	handle->fftwConvPlan = fftw_plan_dft_c2r_2d(
								handle->blockSize, 4,
                               	(fftw_complex *)handle->fftw_output, 
								(double*)handle->output,
                               	FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
	*/

	// int n[]={handle->blockSize};

/*
fftw_plan fftw_plan_many_dft_r2c(
								int rank, 
								const int *n, 
								int howmany,

                                double *in, 
								const int *inembed,
								int istride, 
								int idist,

								fftw_complex *out, 
								const int *onembed,
								int ostride, 
								int odist,

								unsigned flags
								);
	handle->fftwBufPlan =		
		fftw_plan_many_dft_r2c(	1,n,4,
                               	(double *)handle->input,NULL,
								4,1,
								(fftw_complex*)handle->fftw_input,NULL,
								4,1,
                               	FFTW_ESTIMATE);

	handle->fftwFilterPlan = 
		fftw_plan_many_dft_r2c(	1,n,4,
                               	(double *)handle->fir,NULL,
								4,1,
								(fftw_complex*)handle->fftw_fir,NULL,
								4,1,
                               	FFTW_ESTIMATE);

	handle->fftwConvPlan = 
		fftw_plan_many_dft_c2r(	1,n,4,
                               	(fftw_complex *)handle->fftw_output,NULL,
								4,1,
								(double*)handle->output,NULL,
								4,1,
                               	FFTW_ESTIMATE);

	*/
	//fftw_export_wisdom_to_filename(".fftw");
	
	// calculate scaling factor
	handle->scale=1.0f/(handle->blockSize);
	//handle->scale=1.0f;

	return handle;
}

void
convo_destroyPlans(CFormatPlan* plans) {
	fftw_destroy_plan(plans->A);
	fftw_destroy_plan(plans->B);
	fftw_destroy_plan(plans->C);
	fftw_destroy_plan(plans->D);
}	

void
convo_executePlans(CFormatPlan* plans) {
	fftw_execute(plans->A);
	fftw_execute(plans->B);
	fftw_execute(plans->C);
	fftw_execute(plans->D);
}	


void
convo_initPlansC2R(Convo*handle, CFormatPlan* plans, fftwPtr *in, CFormatPtr* out) {
	plans->A= fftw_plan_dft_c2r_1d(
				handle->blockSize,
				(fftw_complex *)in->A,
				(double *)out->A,
				FFTW_ESTIMATE|FFTW_PRESERVE_INPUT);

	plans->B= fftw_plan_dft_c2r_1d(
				handle->blockSize,
				(fftw_complex *)in->B,
				(double *)out->B,
				FFTW_ESTIMATE|FFTW_PRESERVE_INPUT);

	plans->C= fftw_plan_dft_c2r_1d(
				handle->blockSize,
				(fftw_complex *)in->C,
				(double *)out->C,
				FFTW_ESTIMATE|FFTW_PRESERVE_INPUT);

	plans->D= fftw_plan_dft_c2r_1d(
				handle->blockSize,
				(fftw_complex *)in->D,
				(double *)out->D,
				FFTW_ESTIMATE|FFTW_PRESERVE_INPUT);
}
void
convo_initPlansR2C(Convo*handle,CFormatPlan* plans, CFormatPtr* in, fftwPtr* out) {
	plans->A= fftw_plan_dft_r2c_1d(
				handle->blockSize,
				(double *)in->A,
				(fftw_complex *)out->A,
				FFTW_ESTIMATE|FFTW_PRESERVE_INPUT);

	plans->B= fftw_plan_dft_r2c_1d(
				handle->blockSize,
				(double *)in->B,
				(fftw_complex *)out->B,
				FFTW_ESTIMATE|FFTW_PRESERVE_INPUT);

	plans->C= fftw_plan_dft_r2c_1d(
				handle->blockSize,
				(double *)in->C,
				(fftw_complex *)out->C,
				FFTW_ESTIMATE|FFTW_PRESERVE_INPUT);

	plans->D= fftw_plan_dft_r2c_1d(
				handle->blockSize,
				(double *)in->D,
				(fftw_complex *)out->D,
				FFTW_ESTIMATE|FFTW_PRESERVE_INPUT);
}

void
convo_zeroBuf(Convo * handle, CFormatPtr* in) {
	memset((void*)in->A,0,handle->DbufSize);
	memset((void*)in->B,0,handle->DbufSize);
	memset((void*)in->C,0,handle->DbufSize);
	memset((void*)in->D,0,handle->DbufSize);
}

void
convo_zeroFftwBuf(Convo * handle, fftwPtr* in) {
	memset((void*)in->A,0,handle->WbufSize);
	memset((void*)in->B,0,handle->WbufSize);
	memset((void*)in->C,0,handle->WbufSize);
	memset((void*)in->D,0,handle->WbufSize);
}

void
convo_initFftwBufs(Convo * handle, fftwPtr* in) {
	in->A=(fftw_complex*)fftw_malloc(handle->WbufSize);
	in->B=(fftw_complex*)fftw_malloc(handle->WbufSize);
	in->C=(fftw_complex*)fftw_malloc(handle->WbufSize);
	in->D=(fftw_complex*)fftw_malloc(handle->WbufSize);
}

void
convo_initBufs(Convo * handle, CFormatPtr* in) {
	in->A=(double*)fftw_malloc(handle->DbufSize);
	in->B=(double*)fftw_malloc(handle->DbufSize);
	in->C=(double*)fftw_malloc(handle->DbufSize);
	in->D=(double*)fftw_malloc(handle->DbufSize);
}

void
convo_freeFftwBufs(fftwPtr* in) {
	fftw_free(in->A);
	fftw_free(in->B);
	fftw_free(in->C);
	fftw_free(in->D);
}
void
convo_freeBufs(CFormatPtr* in) {
	fftw_free(in->A);
	fftw_free(in->B);
	fftw_free(in->C);
	fftw_free(in->D);
}

void
convo_transformIn(CFormat* in, CFormatPtr out,int size) {
	for (int i=0;i<size;i++) {
		out.A[i]=in[i].A;	
		out.B[i]=in[i].B;	
		out.C[i]=in[i].C;	
		out.D[i]=in[i].D;	
	}
	
}

void
convo_transformOutScale(Convo*handle,CFormatPtr in,CFormat* out,int size) {
	for (int i=0;i<size;i++) {
		out[i].A=in.A[i]*handle->scale;
		out[i].B=in.B[i]*handle->scale;
		out[i].C=in.C[i]*handle->scale;	
		out[i].D=in.D[i]*handle->scale;	
	}
}

/*
 * Initialize filter block 
 */
void
convo_initFilter(Convo * handle) {

	assert(handle != NULL);
	memset(handle->fir,0,handle->CbufSize);

}


/*
 * if we are not running convo_testFilter, we are loading the filter into the handle->fir buffer in
 * our processing routine.
 */
void 
convo_testFilter(Convo * handle) {
	assert(handle != NULL);
	handle->fir[0].A=(double)1.0;
	handle->fir[0].B=(double)1.0;
	handle->fir[0].C=(double)1.0;
	handle->fir[0].D=(double)1.0;
	//convo_dumpRealBuffer(handle->fir,10);
}

/*
 * process filter block and build FFT
 */
void
convo_buildFilter(Convo * handle) {

	assert(handle != NULL);
	//fftw_execute(handle->fftwFilterPlan);
	convo_transformIn(handle->fir,handle->firPtr,handle->blockSize);
	convo_executePlans(&handle->fftwFilterPlan);
	//convo_dumpImagBuffer1(handle, handle->fftw_firPtr.A,handle->blockSize,"A Filter FFTW");
	//convo_dumpImagBuffer1(handle, handle->fftw_firPtr.C,handle->blockSize,"C Filter FFTW");
}

void
convo_dumpRealBuffer1(Convo* handle, double* buf, int size,const char* title) {
	for (int i=0;i<size;i++) {
		if (i>100) break;
		printf("%s:%s:%u:   %f\n",handle->title,title,i,buf[i]);
	}
}
void
convo_dumpImagBuffer1(Convo * handle, fftw_complex* buf, int size,const char* title) {
	fftw_complex *a;
	for (int i=0;i<size;i+=1) {
		if (i>100) break;
		a=buf+i;
		printf("%s:%s: %u:   %f,%f\n",handle->title,title,i,*a[0],*a[1]);
	}
}
void
convo_dumpRealBuffer(CFormat* buf, int size) {
	for (int i=0;i<size;i++) {
		printf("%u:   %f,%f,%f,%f\n",i,buf[i].A,buf[i].B,buf[i].C,buf[i].D);
	}
}
void
convo_dumpImagBuffer(fftw_complex* buf, int size) {
	fftw_complex *a,*b,*c,*d;
	for (int i=0;i<size*4;i+=4) {
		a=&buf[i];
		b=&buf[i];
		c=&buf[i];
		d=&buf[i];
		printf("%u:   %f,%f,%f,%f,%f,%f,%f,%f\n",i,*a[0],*a[1],*b[0],*b[1],*c[0],*c[1],*d[0],*d[1]);
	}
}

/*
 * this performs the convolution multiplication of spectra
 */
void
convo_multiply4FFTs(Convo * handle) {

	convo_multiplyFFTs(handle,handle->fftw_inputPtr.A,handle->fftw_firPtr.A,handle->fftw_outputPtr.A);
	convo_multiplyFFTs(handle,handle->fftw_inputPtr.B,handle->fftw_firPtr.B,handle->fftw_outputPtr.B);
	convo_multiplyFFTs(handle,handle->fftw_inputPtr.C,handle->fftw_firPtr.C,handle->fftw_outputPtr.C);
	convo_multiplyFFTs(handle,handle->fftw_inputPtr.D,handle->fftw_firPtr.D,handle->fftw_outputPtr.D);

}

/*
 * this performs the convolution multiplication of spectra
 */
void
convo_multiplyFFTs(Convo * handle, fftw_complex * bufIn, fftw_complex * bufFir, fftw_complex * bufOut) {

	assert(handle != NULL);
	assert(bufIn  != NULL);
	assert(bufOut != NULL);
	assert(bufFir != NULL);

	fftw_complex *a,*b,*c;

	for (int i=0 ; i < handle->blockSizeCpx ; i++) {

		a=bufIn+i;
		b=bufFir+i;
		c=bufOut+i;

		(*c)[0]=(
			 ((*a)[0]*(*b)[0])
			-((*a)[1]*(*b)[1])
		);

		(*c)[1]=(
			 ((*a)[1]*(*b)[0])
			+((*a)[0]*(*b)[1])
		);

	}
}

/* 
 * final overlap add stage after the IFFT process
 */
void
convo_overlapAdd(Convo * handle, CFormat * d_out) {

	assert(handle != NULL);
	assert(d_out != NULL);

	// we now add the first half of the convolution solution to the saved overlap block from the previous call
	// to this function. we also normalize here too.
	for (int i=0;i<handle->dataSize;i++) {
		if (i<handle->filterSize-1) {
			d_out[i].A=(handle->output[i].A+handle->olap[i].A);
			d_out[i].B=(handle->output[i].B+handle->olap[i].B);
			d_out[i].C=(handle->output[i].C+handle->olap[i].C);
			d_out[i].D=(handle->output[i].D+handle->olap[i].D);
//			d_out[i].A=(handle->output[i].A+handle->olap[i].A)*handle->scale;
//			d_out[i].B=(handle->output[i].B+handle->olap[i].B)*handle->scale;
//			d_out[i].C=(handle->output[i].C+handle->olap[i].C)*handle->scale;
//			d_out[i].D=(handle->output[i].D+handle->olap[i].D)*handle->scale;
		} else {
			d_out[i].A=(handle->output[i].A);
			d_out[i].B=(handle->output[i].B);
			d_out[i].C=(handle->output[i].C);
			d_out[i].D=(handle->output[i].D);
//			d_out[i].A=(handle->output[i].A)*handle->scale;
//			d_out[i].B=(handle->output[i].B)*handle->scale;
//			d_out[i].C=(handle->output[i].C)*handle->scale;
//			d_out[i].D=(handle->output[i].D)*handle->scale;
		}
	}
}

/*
 * cleanup and free memory structures 
 */
void
convo_cleanup(Convo * handle) {

	assert(handle != NULL);

#ifdef _DBG_CONVO
	sf_close(handle->snd_debugfile);
#endif

	convo_destroyPlans(&handle->fftwBufPlan);
	convo_destroyPlans(&handle->fftwFilterPlan);
	convo_destroyPlans(&handle->fftwConvPlan);

	convo_freeBufs(&handle->firPtr);
	convo_freeBufs(&handle->inputPtr);
	convo_freeBufs(&handle->outputPtr);

	fftw_free(handle->fir);
	fftw_free(handle->input);
	fftw_free(handle->output);

	convo_freeFftwBufs(&handle->fftw_firPtr);
	convo_freeFftwBufs(&handle->fftw_inputPtr);
	convo_freeFftwBufs(&handle->fftw_outputPtr);

	fftw_free(handle->fftw_fir);
	fftw_free(handle->fftw_input);
	fftw_free(handle->fftw_output);

	free(handle);
}

/*
 * perform the convolution process sequence
 */
void
convo_volve(Convo * handle, CFormat *d_out, CFormat *d_in) {

	assert(handle != NULL);
	assert(d_out != NULL);
	assert(d_in != NULL);

	/*
	 * N=M+L-1
	 *
	 * initialize input and output structures
	 *
	 * TODO: CHECK SEE IF WE NEED THESE STILL
	 */
	//memset((void*)handle->fftw_input,0,handle->CbufSize); // this is to reset the contents of the Input buffer
	//memset((void*)handle->fftw_output,0,handle->CbufSize); // this resets the contents of the output buffer.

	// fill half the input buffers with input data
	// note this is an interleaved buffer
	memcpy((void*)handle->input,d_in,sizeof(CFormat)*handle->dataSize); // copy the contents of the transfer buffer into the input buffer.

	// convert from interleaved format to plain format
	convo_zeroBuf(handle,&handle->inputPtr);
	convo_transformIn(handle->input,handle->inputPtr,handle->dataSize);

	// execute the FFT on the input data
	convo_executePlans (&handle->fftwBufPlan);

	//Multiply the data FFT with the existing filter FFT
	convo_multiply4FFTs(handle);

	// inverse FFT the convolution to get the time domain data
	convo_executePlans (&handle->fftwConvPlan);

	// back to interleaved format
	convo_transformOutScale(handle,handle->outputPtr,handle->output,handle->blockSize);

#ifdef _DBG_CONVO
	convo_debugBlock(handle,handle->outputPtr.A,handle->dataSize,"A_buf_aftr_transformOut");
#endif

	// fix up the circular convolution wrap
	convo_overlapAdd(handle, (CFormat *)d_out);

	// We need to copy out the last block and save it for next time.
	memcpy(
		handle->olap,
		handle->output+handle->dataSize,
		handle->ObufSize
	);

}


void 
convo_debugBlock(Convo * handle, double* buffer, size_t size, char* title) {
	char sbuf[50];

#ifdef _DBG_CONVO
	//printf("debugging block \n");
	if (handle->snd_debugfile==NULL) {
		printf("allocating debug SFI debug block\n");

		handle->sfi_debugfile.channels=1;
		handle->sfi_debugfile.format=SF_FORMAT_WAV|SF_FORMAT_PCM_24;
		handle->sfi_debugfile.samplerate=48000;

		sprintf(sbuf,"dbg-%s-%s.wav",handle->title,title);
	
		if (! (handle->snd_debugfile = sf_open (sbuf, SFM_WRITE, &handle->sfi_debugfile)))
		{	printf ("Not able to open output file %s.\n", sbuf) ;
			puts (sf_strerror (NULL)) ;
		} ;

	}

	sf_writef_double(handle->snd_debugfile,(double*)buffer,size);
	//printf("writing debug block of size %zu\n",size);


#endif
}
