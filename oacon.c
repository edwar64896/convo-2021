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
convo_init(int filterSize, int dataSize, const char* title, int verbose) {

#ifdef _DBG_CONVO
	if (verbose)
		printf("initializing convo with Filtersize=%u and Datasize=%u with Title=%s\n",filterSize,dataSize,title);
#endif

	Convo *handle;
	handle=(Convo*)malloc(sizeof(Convo));

#ifdef _DBG_CONVO
	handle->dbgBlock=0;
	handle->snd_debugfile=NULL;
	memset((void*)&handle->sfi_debugfile,0,sizeof(SF_INFO));
#endif

	handle->title=title;
	handle->verbose=verbose;

	handle->filterSize=filterSize; //size in samples of the filter 	=M
	handle->dataSize=dataSize; //size in samples of a data block 	=L

	// this is the block size we work in during FFT/convolution. M+L-1
	handle->blockSize=(filterSize+dataSize)-1; // counted in Samples

	// real to complex transforms have a real buffer length of 'n' and a complex buffer size of floor(n/2)+1
	handle->blockSizeCpx=floor(handle->blockSize/2)+1;

	// this is the size of one FFT stream/channel in bytes
//	handle->WbufSize=(size_t)(sizeof(fftw_complex)		* handle->blockSizeCpx); //bufsize in bytes

	// this is the size of one sample stream/channel in bytes
//	handle->DbufSize=(size_t)(sizeof(double) 			* handle->blockSize); //bufsize in bytes

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

	// initialize FFTW threading engine
	fftw_init_threads();
	fftw_plan_with_nthreads(4);

	int n[]={handle->blockSize};

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

	
	// calculate scaling factor
	handle->scale=1.0f/(handle->blockSize);

	return handle;
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
}

/*
 * process filter block and build FFT
 */
void
convo_buildFilter(Convo * handle) {

	assert(handle != NULL);
	fftw_execute(handle->fftwFilterPlan);
}

#ifdef _DBG_CONVO
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
#endif

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

	for (int i=0 ; i < handle->blockSizeCpx * 4 ; i+=4) {

		for (int j=0;j<4;j++) {

			a=bufIn+i+j;
			b=bufFir+i+j;
			c=bufOut+i+j;
	
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
			d_out[i].A=(handle->output[i].A+handle->olap[i].A)*handle->scale;
			d_out[i].B=(handle->output[i].B+handle->olap[i].B)*handle->scale;
			d_out[i].C=(handle->output[i].C+handle->olap[i].C)*handle->scale;
			d_out[i].D=(handle->output[i].D+handle->olap[i].D)*handle->scale;
		} else {
			d_out[i].A=(handle->output[i].A)*handle->scale;
			d_out[i].B=(handle->output[i].B)*handle->scale;
			d_out[i].C=(handle->output[i].C)*handle->scale;
			d_out[i].D=(handle->output[i].D)*handle->scale;
		}
	}
}

/*
 * cleanup and free memory structures 
 */
void
convo_cleanup(Convo * handle) {

	assert(handle != NULL);

	//cleanup fftw

	fftw_destroy_plan(handle->fftwBufPlan);
	fftw_destroy_plan(handle->fftwFilterPlan);
	fftw_destroy_plan(handle->fftwConvPlan);
	fftw_cleanup_threads();

#ifdef _DBG_CONVO
	sf_close(handle->snd_debugfile);
#endif

	fftw_free(handle->fir);
	fftw_free(handle->input);
	fftw_free(handle->output);

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
	memset((void*)handle->fftw_input,0,handle->CbufSize); // this is to reset the contents of the Input buffer
	memset((void*)handle->fftw_output,0,handle->CbufSize); // this resets the contents of the output buffer.

	// fill half the input buffers with input data
	// note this is an interleaved buffer
	memcpy((void*)handle->input,d_in,sizeof(CFormat)*handle->dataSize); // copy the contents of the transfer buffer into the input buffer.

	// execute the FFT on the input data
	fftw_execute (handle->fftwBufPlan);

	//Multiply the data FFT with the existing filter FFT
	convo_multiplyFFTs(handle,handle->fftw_input,handle->fftw_fir,handle->fftw_output);

	// inverse FFT the convolution to get the time domain data
	fftw_execute (handle->fftwConvPlan);

#ifdef _DBG_CONVO
	//convo_debugBlock(handle,handle->outputPtr.A,handle->dataSize,"A_buf_aftr_transformOut");
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


#ifdef _DBG_CONVO
void 
convo_debugBlock(Convo * handle, double* buffer, size_t size, char* title) {
	char sbuf[50];

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
}
#endif
