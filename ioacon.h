
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

struct iConvolutionStruct {
		
	fftw_complex *fftw_input;
	fftw_complex *fftw_output;
	fftw_complex *fftw_filter;

	CFormat *filterOut;
	CFormat *filter;
	CFormat *input;
	CFormat *output;

	CFormat max;
	CFormat scale;
	double overmax;

	size_t blockSize;
	size_t bufferSize;
	size_t bufferSizeCpx;

	size_t BbufSize; // size of the buffer that will hold the FFT Input
	size_t CbufSize; // size of the buffer that will hold the FFT Input
	size_t FbufSize;
	size_t ObufSize;

	size_t inputSize; // =M
	size_t outputSize;   // =L

//	double scale;

	fftw_plan fftwInputPlan;
	fftw_plan fftwOutputPlan;
	fftw_plan fftwFilterPlan;

	int dbgBlock;
};

typedef struct iConvolutionStruct iConvo;

iConvo* iconvo_init(int bufferSize);
void iconvo_initInput(iConvo * handle) ;
void iconvo_initOutput(iConvo * handle) ;
void iconvo_devolve(iConvo * handle);
void iconvo_divideFFTs(iConvo * handle) ;
void iconvo_split(iConvo * handle) ;
void iconvo_cleanup(iConvo * handle) ;

#endif //_OACON_H
