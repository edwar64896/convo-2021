
#ifndef _OACON_H
#define _OACON_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <fftw3.h>

#include <sndfile.h>

struct AFormatStruct {
	double FLU;
	double FRD;
	double BLD;
	double BRU;
};

struct BFormatStruct {
	double W;
	double X;
	double Y;
	double Z;
};

struct CFormatStruct {
	double A;
	double B;
	double C;
	double D;
};

struct CFormatPtrStruct {
	double *A;
	double *B;
	double *C;
	double *D;
};


struct fftwPtrStruct {
	fftw_complex *A;
	fftw_complex *B;
	fftw_complex *C;
	fftw_complex *D;
};


typedef struct AFormatStruct AFormat;
typedef struct BFormatStruct BFormat;
typedef struct CFormatStruct CFormat;
typedef struct CFormatPtrStruct CFormatPtr;
typedef struct fftwPtrStruct fftwPtr;

struct CFormatPlanStruct{
	fftw_plan A;
	fftw_plan B;
	fftw_plan C;
	fftw_plan D;
};

typedef struct CFormatPlanStruct CFormatPlan;

struct ConvolutionStruct {

	const char* title;
		
	fftw_complex *fftw_fir;
	fftw_complex *fftw_input;
	fftw_complex *fftw_output;

	CFormat *fir;
	CFormat *input;
	CFormat *output;
	CFormat *olap;

	CFormatPtr firPtr;
	CFormatPtr inputPtr;
	CFormatPtr outputPtr;
	CFormatPtr olapPtr;
	
	fftwPtr fftw_firPtr;
	fftwPtr fftw_inputPtr;
	fftwPtr fftw_outputPtr;
	fftwPtr fftw_olapPtr;

	size_t CbufSize; // size of the buffer that will hold the FFT Input
	size_t FbufSize;
	size_t ObufSize;
	size_t WbufSize;
	size_t DbufSize;

	size_t filterSize; // =M
	size_t dataSize;   // =L

	double scale;

	CFormatPlan fftwBufPlan;
	CFormatPlan fftwFilterPlan;
	CFormatPlan fftwConvPlan;

	int blockSize;
	int blockSizeCpx;
	int frameSize;

	int dbgBlock;

	SNDFILE * snd_debugfile;
	SF_INFO sfi_debugfile;

};

typedef struct ConvolutionStruct Convo;

Convo* convo_init(int filterSize, int dataSize, const char* title) ;
void convo_initFilter(Convo * handle) ;
void convo_testFilter(Convo * handle) ;
void convo_buildFilter(Convo * handle) ;
void convo_multiplyFFTs(Convo * handle,fftw_complex*,fftw_complex*,fftw_complex*) ;
void convo_overlapAdd(Convo * handle, CFormat * d_out) ;
void convo_cleanup(Convo * handle) ;
void convo_volve(Convo * handle, CFormat *d_out, CFormat *d_in) ;
void convo_debugBlock(Convo * handle, double* buffer, size_t size, char* title) ;
void convo_dumpRealBuffer1(Convo*, double* buf, int size,const char* ) ;
void convo_dumpImagBuffer1(Convo*, fftw_complex* buf, int size,const char*) ;
void convo_dumpRealBuffer(CFormat* buf, int size) ;
void convo_dumpImagBuffer(fftw_complex* buf, int size) ;
void convo_initBufs(Convo *,CFormatPtr* in) ;
void convo_initFftwBufs(Convo *,fftwPtr* in) ;
void convo_freeBufs(CFormatPtr* in) ;
void convo_initPlansR2C(Convo *,CFormatPlan*,CFormatPtr*,fftwPtr*);
void convo_initPlansC2R(Convo *,CFormatPlan*,fftwPtr*,CFormatPtr*);
void convo_executePlans(CFormatPlan*);

#endif //_OACON_H
