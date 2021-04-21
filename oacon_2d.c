/*
 *
 * works on structs of 4 double
 *
 * [ FIR = M ][L-1 Zeros---------] M+L-1=N
 * [ Data = L--------][M-1 Zeros ] L+M-1=N
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <fftw3.h>
#include <sndfile.h>

#include "oacon.h"


/*
 * initialize structure. we are dealing with set frame sizes (CFormat)
 *
 * blocksize = L+M-1
 */
Convo*
convo_init(int filterSize, int dataSize) {

	Convo *handle;
	handle=(Convo*)malloc(sizeof(Convo));

	handle->dbgBlock=0;

	handle->filterSize=filterSize; //size in samples of the filter 	=M
	handle->dataSize=dataSize; //size in samples of a data block 	=L

	handle->blockSize=(filterSize+dataSize)-1; // counted in Samples

	handle->CbufSize=(size_t)(sizeof(CFormat) 			* handle->blockSize);
	handle->FbufSize=(size_t)(sizeof(fftw_complex) 	* 4	* handle->blockSize);
	handle->ObufSize=(size_t)(sizeof(CFormat) 			* (handle->filterSize-1));

	handle->fftw_fir=	(fftw_complex*)	fftw_malloc(handle->FbufSize);
	handle->fftw_input=	(fftw_complex*)	fftw_malloc(handle->FbufSize);
	handle->fftw_output=(fftw_complex*)	fftw_malloc(handle->FbufSize);

	handle->fir=		(CFormat*)		fftw_malloc(handle->CbufSize);
	handle->input=		(CFormat*)		fftw_malloc(handle->CbufSize);
	handle->output=		(CFormat*)		fftw_malloc(handle->CbufSize);

	handle->olap=		(CFormat*)		fftw_malloc(handle->ObufSize);

	memset((void*)handle->fftw_fir,		0,handle->FbufSize);
	memset((void*)handle->fftw_input,	0,handle->FbufSize);
	memset((void*)handle->fftw_output,	0,handle->FbufSize);

	fftw_init_threads();
	fftw_plan_with_nthreads(4);

	fftw_import_wisdom_from_filename(".fftw");
	handle->fftwBufPlan = fftw_plan_dft_r2c_2d(
								handle->blockSize, 4, 
                               	(double *)handle->input, 
								(fftw_complex*)handle->fftw_input,
                               	FFTW_PATIENT|FFTW_DESTROY_INPUT);
	handle->fftwFilterPlan = fftw_plan_dft_r2c_2d(
								handle->blockSize, 4,
                               	(double *)handle->fir, 
								(fftw_complex*)handle->fftw_fir,
                               	FFTW_PATIENT|FFTW_DESTROY_INPUT);
	handle->fftwConvPlan = fftw_plan_dft_c2r_2d(
								handle->blockSize, 4,
                               	(fftw_complex *)handle->fftw_output, 
								(double*)handle->output,
                               	FFTW_PATIENT|FFTW_DESTROY_INPUT);
	fftw_export_wisdom_to_filename(".fftw");
	
	// calculate scaling factor
	handle->scale=1.0f/(handle->blockSize*8);
	//handle->scale=1.0f;

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

/*
 * this performs the convolution division (inverse) of spectra
 *
 * using formula from http://www.mathportal.org/formulas/algebra/complex.php
 */
void
convo_divideFFTs(Convo * handle) {

	assert(handle != NULL);

	fftw_complex *a,*b,*c;

	for (int i=0 ; i < handle->blockSize * 4 ; i++) {
		a=&handle->fftw_input[i];
		b=&handle->fftw_output[i];
		c=&handle->fftw_fir[i];

		*c[0]=( (((*a[0])*(*b[0]))+((*b[1])*(*a[0])))/(((*a[0])*(*a[0]))+((*a[1])*(*a[1]))) );
		*c[1]=( (((*b[0])*(*a[0]))-((*b[0])*(*a[1])))/(((*a[0])*(*a[0]))+((*a[1])*(*a[1]))) );

	}
}
/*
 * this performs the convolution multiplication of spectra
 */
void
convo_multiplyFFTs(Convo * handle) {

	assert(handle != NULL);

	fftw_complex *a,*b,*c;

	for (int i=0 ; i < handle->blockSize * 4 ; i++) {

		a=&handle->fftw_input[i];
		b=&handle->fftw_output[i];
		c=&handle->fftw_fir[i];

		*b[0]=(
			 ((*a[0])*(*c[0]))
			-((*a[1])*(*c[1]))
		);

		*b[1]=(
			 ((*a[1])*(*c[0]))
			+((*a[0])*(*c[1]))
		);

	}
}

/*
 * this performs the convolution multiplication of spectra
 */
void
convo_multiplyFFTs_old(Convo * handle) {

	assert(handle != NULL);

	for (int i=0 ; i < handle->blockSize * 4 ; i++) {

		handle->fftw_output[i][0]=(
			 ((handle->fftw_input[i][0])*(handle->fftw_fir[i][0]))
			-((handle->fftw_input[i][1])*(handle->fftw_fir[i][1]))
		);

		handle->fftw_output[i][1]=(
			 ((handle->fftw_input[i][1])*(handle->fftw_fir[i][0]))
			+((handle->fftw_input[i][0])*(handle->fftw_fir[i][1]))
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

	fftw_destroy_plan(handle->fftwBufPlan);
	fftw_destroy_plan(handle->fftwFilterPlan);
	fftw_destroy_plan(handle->fftwConvPlan);

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
	 */
	memset((void*)handle->fftw_input,0,handle->CbufSize); // this is to reset the contents of the Input buffer
	memset((void*)handle->fftw_output,0,handle->CbufSize); // this resets the contents of the output buffer.

	// fill half the input buffers with input data
	memcpy((void*)handle->input,d_in,sizeof(BFormat)*handle->dataSize); // copy the contents of the transfer buffer into the input buffer.

	// execute the FFT on the input data
	fftw_execute (handle->fftwBufPlan);

	//Multiply the data FFT with the existing filter FFT
	convo_multiplyFFTs(handle);

	// inverse FFT the convolution to get the time domain data
	fftw_execute (handle->fftwConvPlan);

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
convo_debugBlock(Convo * handle, double* buffer, char* title, SF_INFO* sfi) {
	SNDFILE * snd_debugfile;
	char sbuf[50];

	handle->dbgBlock++;

	sprintf(sbuf,"dbg-%s-%03u.wav",title,handle->dbgBlock);
	
	
	if (! (snd_debugfile = sf_open (sbuf, SFM_WRITE, sfi)))
	{	printf ("Not able to open output file %s.\n", sbuf) ;
		puts (sf_strerror (NULL)) ;
	} ;

	sf_writef_double(snd_debugfile,(double*)buffer,handle->blockSize);

	sf_close(snd_debugfile);

}
