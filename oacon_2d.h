
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

typedef struct AFormatStruct AFormat;
typedef struct BFormatStruct BFormat;
typedef struct CFormatStruct CFormat;

struct ConvolutionStruct {
		
	fftw_complex *fftw_fir;
	fftw_complex *fftw_input;
	fftw_complex *fftw_output;

	CFormat *fir;
	CFormat *input;
	CFormat *output;
	CFormat *olap;

	size_t CbufSize; // size of the buffer that will hold the FFT Input
	size_t FbufSize;
	size_t ObufSize;

	size_t filterSize; // =M
	size_t dataSize;   // =L

	double scale;

	fftw_plan fftwBufPlan;
	fftw_plan fftwFilterPlan;
	fftw_plan fftwConvPlan;

	int blockSize;
	int frameSize;

	int dbgBlock;
	
};

typedef struct ConvolutionStruct Convo;

Convo* convo_init(int filterSize, int dataSize) ;
void convo_initFilter(Convo * handle) ;
void convo_testFilter(Convo * handle) ;
void convo_buildFilter(Convo * handle) ;
void convo_multiplyFFTs(Convo * handle) ;
void convo_overlapAdd(Convo * handle, CFormat * d_out) ;
void convo_cleanup(Convo * handle) ;
void convo_volve(Convo * handle, CFormat *d_out, CFormat *d_in) ;
void convo_debugBlock(Convo * handle, double* buffer, char* title, SF_INFO* sfi) ;

#endif //_OACON_H
