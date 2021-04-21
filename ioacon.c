/*
 *
 * works on structs of 4 double
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

#include "ioacon.h"


/*
 * initialize structure. we are dealing with set frame sizes (CFormat)
 * additionally, we are not going to use overlap add with reverse convolution
 * this will be a one-off to find the impulse response
 *
 * All buffers must be equal - input buffer size, output buffer size and the resultant filter size.
 *
 */
iConvo*
iconvo_init(int blockSize) {

	iConvo *handle;

	handle=(iConvo*)malloc(sizeof(iConvo));

	handle->blockSize=blockSize;
	handle->bufferSize=(2 * blockSize) -1; //real buffer size
	handle->bufferSizeCpx=floor(handle->bufferSize/2)+1; //complex buffer sizes

	handle->BbufSize=(size_t)(sizeof(CFormat) 			* handle->blockSize );
	handle->CbufSize=(size_t)(sizeof(CFormat) 			* handle->bufferSize );
	handle->FbufSize=(size_t)(sizeof(fftw_complex) 	* 4	* handle->bufferSizeCpx );

	handle->fftw_input=	(fftw_complex*)	fftw_malloc(handle->FbufSize);
	handle->fftw_output=(fftw_complex*)	fftw_malloc(handle->FbufSize);
	handle->fftw_filter=(fftw_complex*)	fftw_malloc(handle->FbufSize);

	handle->filter=		(CFormat*)		fftw_malloc(handle->CbufSize);
	handle->input=		(CFormat*)		fftw_malloc(handle->CbufSize);
	handle->output=		(CFormat*)		fftw_malloc(handle->CbufSize);
	handle->filterOut=	(CFormat*)		fftw_malloc(handle->BbufSize);

	memset((void*)handle->fftw_filter,	0,handle->FbufSize);
	memset((void*)handle->fftw_input,	0,handle->FbufSize);
	memset((void*)handle->fftw_output,	0,handle->FbufSize);

	memset((void*)handle->filter,		0,handle->CbufSize);
	memset((void*)handle->input,		0,handle->CbufSize);
	memset((void*)handle->output,		0,handle->CbufSize);

	memset((void*)&handle->scale,0,sizeof(CFormat));
	memset((void*)&handle->max,0,sizeof(CFormat));


	//printf("Checking values of buffer %f\n",handle->filter[0].A);

	fftw_init_threads();
	fftw_plan_with_nthreads(4);

	int n[]={handle->bufferSize};

	handle->fftwInputPlan =		
				fftw_plan_many_dft_r2c(
								1,n,4,
                               	(double *)handle->input, NULL,4,1,
								(fftw_complex*)handle->fftw_input, NULL,4,1,
                               	FFTW_ESTIMATE|FFTW_DESTROY_INPUT);

	handle->fftwOutputPlan = 
				fftw_plan_many_dft_r2c(
								1,n,4,
                               	(double *)handle->output, NULL,4,1,
								(fftw_complex*)handle->fftw_output, NULL,4,1,
                               	FFTW_ESTIMATE|FFTW_DESTROY_INPUT);

	handle->fftwFilterPlan = 
				fftw_plan_many_dft_c2r(
								1,n,4,
                               	(fftw_complex *)handle->fftw_filter, NULL,4,1,
								(double*)handle->filter, NULL,4,1,
                               	FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
	return handle;
}

/*
 * Initialize input block 
 */
void
iconvo_initInput(iConvo * handle) {

	assert(handle != NULL);
	memset(handle->input,0,handle->CbufSize);

}
/*
 * Initialize output block 
 */
void
iconvo_initOutput(iConvo * handle) {

	assert(handle != NULL);
	memset(handle->output,0,handle->CbufSize);
}

/*
 * this performs the iconvolution division (inverse) of spectra
 *
 * using formula from http://www.mathportal.org/formulas/algebra/complex.php
 */
void
iconvo_divideFFTs(iConvo * handle) {

	assert(handle != NULL);

	fftw_complex *a,*b,*c;

	for (int i=0 ; i < handle->bufferSizeCpx * 4 ; i++) {

		a=handle->fftw_input+i;
		b=handle->fftw_output+i;
		c=handle->fftw_filter+i;

		(*c)[0]=((((*a)[0]*(*b)[0])+((*b)[1]*(*a)[0]))/(((*a)[0]*(*a)[0])+((*a)[1]*(*a)[1])));
		(*c)[1]=((((*b)[0]*(*a)[0])-((*b)[0]*(*a)[1]))/(((*a)[0]*(*a)[0])+((*a)[1]*(*a)[1])));
	}
}

/*
 * cleanup and free memory structures 
 */
void
iconvo_cleanup(iConvo * handle) {

	assert(handle != NULL);

	fftw_destroy_plan(handle->fftwInputPlan);
	fftw_destroy_plan(handle->fftwFilterPlan);
	fftw_destroy_plan(handle->fftwOutputPlan);

	fftw_free(handle->filterOut);
	fftw_free(handle->filter);
	fftw_free(handle->input);
	fftw_free(handle->output);

	fftw_free(handle->fftw_filter);
	fftw_free(handle->fftw_input);
	fftw_free(handle->fftw_output);

	free(handle);
}

void
iconvo_findmax(iConvo * handle) {

	handle->max.A=0.0f;
	handle->max.B=0.0f;
	handle->max.C=0.0f;
	handle->max.D=0.0f;
	handle->overmax=0.0f;

	CFormat* d_out=handle->filter;

//		if (fabs(d_out[i].A) > (handle->max.A)) { (handle->max.A)=fabs(d_out[i].A); }
//		if (fabs(d_out[i].B) > (handle->max.B)) { (handle->max.B)=fabs(d_out[i].B); }
//		if (fabs(d_out[i].C) > (handle->max.C)) { (handle->max.C)=fabs(d_out[i].C); }
//		if (fabs(d_out[i].D) > (handle->max.D)) { (handle->max.D)=fabs(d_out[i].D); }
	for (int i=0;i<handle->bufferSize;i++) { if (fabs(d_out[i].A) > (handle->overmax)) { (handle->overmax)=fabs(d_out[i].A); } }
	for (int i=0;i<handle->bufferSize;i++) { if (fabs(d_out[i].B) > (handle->overmax)) { (handle->overmax)=fabs(d_out[i].B); } }
	for (int i=0;i<handle->bufferSize;i++) { if (fabs(d_out[i].C) > (handle->overmax)) { (handle->overmax)=fabs(d_out[i].C); } }
	for (int i=0;i<handle->bufferSize;i++) { if (fabs(d_out[i].D) > (handle->overmax)) { (handle->overmax)=fabs(d_out[i].D); } } 

	printf("Max Values : %f,%f,%f,%f \n", handle->max.A,handle->max.B,handle->max.C,handle->max.D);
}

void
iconvo_scale(iConvo * handle) {
	CFormat* d_out=handle->filter;

	iconvo_findmax(handle);

	handle->scale.A=1.0f/handle->overmax;
	handle->scale.B=1.0f/handle->overmax;
	handle->scale.C=1.0f/handle->overmax;
	handle->scale.D=1.0f/handle->overmax;


	for (int i=0;i<handle->bufferSize;i++) {
		d_out[i].A=d_out[i].A*handle->scale.A;
		d_out[i].B=d_out[i].B*handle->scale.B;
		d_out[i].C=d_out[i].C*handle->scale.C;
		d_out[i].D=d_out[i].D*handle->scale.D;
	}
}

/*
 * perform the iconvolution process sequence
 */
void
iconvo_devolve(iConvo * handle) {

	assert(handle != NULL);
	
	printf("inside devolve\n");

	// execute the FFT on the input and output data
	fftw_execute (handle->fftwInputPlan);
	printf("done Input fft\n");
	fftw_execute (handle->fftwOutputPlan);
	printf("done Output fft\n");

	//divide the output FFt by the input FFt to get the Impulse Response Filter
	iconvo_divideFFTs(handle);
	printf("done fft division\n");

	// inverse FFT the iconvolution to get the time domain data
	fftw_execute (handle->fftwFilterPlan);
	printf("done Filter FFT\n");

	iconvo_scale(handle);

	iconvo_split(handle);
}

void
iconvo_split(iConvo * handle) {
	size_t stBuf1=floor(handle->blockSize/2);
	size_t stBuf2=handle->blockSize-floor(handle->blockSize/2);
	memcpy(handle->filterOut,		handle->filter+stBuf2,	stBuf1*sizeof(CFormat));
	memcpy(handle->filterOut+stBuf1,handle->filter,			stBuf2*sizeof(CFormat));

}

void 
iconvo_debugBlock(iConvo * handle, double* buffer, char* title, SF_INFO* sfi) {
	SNDFILE * snd_debugfile;
	char sbuf[50];

	handle->dbgBlock++;

	sprintf(sbuf,"dbg-%s-%03u.wav",title,handle->dbgBlock);
	
	
	if (! (snd_debugfile = sf_open (sbuf, SFM_WRITE, sfi)))
	{	printf ("Not able to open output file %s.\n", sbuf) ;
		puts (sf_strerror (NULL)) ;
	} ;

	sf_writef_double(snd_debugfile,(double*)buffer,handle->bufferSize);

	sf_close(snd_debugfile);

}
