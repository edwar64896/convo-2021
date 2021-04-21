
#include	<stdio.h>
#include	<string.h>
#include	<stdlib.h>
#include 	<fftw3.h>

/* Include this header file to use functions from libsndfile. */
#include	<sndfile.h>

/*	This will be the length of the buffer used to hold.frames while
**	we process them.
*/
//#define		BUFFER_LEN	1024

/* libsndfile can handle more than 6 channels but we'll restrict it to 6. */
//#define		MAX_CHANNELS	6

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

typedef struct AFormatStruct AFormat;
typedef struct BFormatStruct BFormat;

/* Function prototype. */
static void convovolveA(AFormat *d_out, AFormat *d_in, AFormat* fir, int count) ;
static void convovolveB(BFormat *d_out, BFormat *d_in, BFormat* fir, int count) ;
static void matrix(BFormat *d_out, AFormat *d_in, int count) ;
static void sfinfo_print(SF_INFO* sfinfo,char*);
static void debugBlock(BFormat*,int ,char* ,int ,SF_INFO* ) ;

BFormat *fftw_d_in_buf;
BFormat *fftw_fir_buf;
BFormat *fftw_conv_buf; // this will hold the final convoluted data (inverse FFT after product calculation)
BFormat *fftw_olap_add; // this will hold the overlap amount from the previous invocation of the convolution function.

fftw_complex *fftw_d_in_buf_fft;
fftw_complex *fftw_fir_buf_fft;
fftw_complex *fftw_conv_buf_fft; // this will be used to hold the product of both the data buffer and the FIR FFT's

fftw_plan fftwPlanBuf;
fftw_plan fftwPlanFIR;
fftw_plan fftwPlanCONV;

SF_INFO		sfinfo ;
SF_INFO		sfinfo_out ;
SF_INFO		sfinfo_dbg ;
SF_INFO		sfinfo_fir1 ;
SF_INFO		sfinfo_fir2 ;

int
main (void) {

	/* A SNDFILE is very much like a FILE in the Standard C library. The
	** sf_open function return an SNDFILE* pointer when they sucessfully
	** open the specified file.
	*/
	SNDFILE	*infile, *outfile ;
	SNDFILE	*infile_fir1, *infile_fir2 ;

	/* A pointer to an SF_INFO struct is passed to sf_open.
	** On read, the library fills this struct with information about the file.
	** On write, the struct must be filled in before calling sf_open.
	*/

	int		readcount ;
	int		rc_fir1 ;
	int		rc_fir2 ;

	const char	*fir1 = "fir1.wav" ;
	const char	*fir2 = "fir2.wav" ;

	const char	*infilename = "input.wav" ;
	const char	*outfilename = "output.wav" ;

	/* The SF_INFO struct must be initialized before using it.
	*/
	memset (&sfinfo, 0, sizeof (sfinfo)) ;
	memset (&sfinfo_out, 0, sizeof (sfinfo_out)) ;
	memset (&sfinfo_dbg, 0, sizeof (sfinfo_dbg)) ;
	memset (&sfinfo_fir1, 0, sizeof (sfinfo_fir1)) ;
	memset (&sfinfo_fir2, 0, sizeof (sfinfo_fir2)) ;

	/*
	 * load in the interleaved FIR impulse response filters
	 * prior to convolution.
	 */

	if (! (infile_fir1 = sf_open (fir1, SFM_READ, &sfinfo_fir1)))
	{	/* Open failed so print an error message. */
		printf ("Not able to open fir1 file %s.\n", fir1) ;
		/* Print the error message from libsndfile. */
		puts (sf_strerror (NULL)) ;
		return 1 ;
	} ;
	sfinfo_print(&sfinfo_fir1,"FIR 1");
	if (sfinfo_fir1.channels != 4)
	{	printf ("FIR1 filter needs 4 channels\n") ;
		return 1 ;
	} ;
	if (! (infile_fir2 = sf_open (fir2, SFM_READ, &sfinfo_fir2)))
	{	/* Open failed so print an error message. */
		printf ("Not able to open fir2 file %s.\n", fir2) ;
		/* Print the error message from libsndfile. */
		puts (sf_strerror (NULL)) ;
		return 1 ;
	};
	sfinfo_print(&sfinfo_fir2,"FIR 2");
	if (sfinfo_fir2.channels != 4)
	{	printf ("FIR2 filter needs 4 channels\n") ;
		return 1 ;
	} ;


	/*
	 * Check the FIR buffers match
	 */
	if (sfinfo_fir1.frames != sfinfo_fir2.frames) {
		printf("FIR lengths do not match\n");
		return 1;
	}
	if (sfinfo_fir1.samplerate != sfinfo_fir2.samplerate) {
		printf("FIR samplerates do not match\n");
		return 1;
	}
	if (sfinfo_fir1.channels != sfinfo_fir2.channels) {
		printf("FIR channel counts do not match\n");
		return 1;
	}
	if (sfinfo_fir1.format != sfinfo_fir2.format) {
		printf("FIR formats do not match\n");
		return 1;
	}

	/*
	 * Identify the length of the FIR filter buffers
	 */
	int fir_buffer_len=sfinfo_fir1.frames;


	printf("creating fft plans\n");

	// allocate convolution buffers
	fftw_d_in_buf=(BFormat*)fftw_malloc(sizeof(BFormat)*fir_buffer_len*2);
	fftw_fir_buf=(BFormat*)fftw_malloc(sizeof(BFormat)*fir_buffer_len*2);
	fftw_conv_buf=(BFormat*)fftw_malloc(sizeof(BFormat)*fir_buffer_len*2);

	fftw_olap_add=(BFormat*)fftw_malloc(sizeof(BFormat)*fir_buffer_len);
	memset((void*)fftw_olap_add,0,sizeof(BFormat)*fir_buffer_len);

	/*
	 * N=M+L-1
	 */

	fftw_d_in_buf_fft=(fftw_complex*)fftw_malloc((sizeof(fftw_complex)*4*fir_buffer_len*2)-4);
	fftw_fir_buf_fft=(fftw_complex*)fftw_malloc((sizeof(fftw_complex)*4*fir_buffer_len*2)-4);
	fftw_conv_buf_fft=(fftw_complex*)fftw_malloc((sizeof(fftw_complex)*4*fir_buffer_len*2)-4);

	printf("sizeof fftw_complex=%lu\n",sizeof(fftw_complex));

	/*
	 * Create the FFT Plan for our FFTW Convolution
	 * FFT SIZED at 
	 * M=LEN(FIR)
	 * L=M
	 * N=M+L-1=2M-1
	 */
	fftwPlanBuf=fftw_plan_dft_r2c_2d((fir_buffer_len*2)-1,4, // reversing fir_buffer_len and 4 makes no difference - still noise output.====
                               (double *)fftw_d_in_buf, fftw_d_in_buf_fft,
                               FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
	fftwPlanFIR=fftw_plan_dft_r2c_2d((fir_buffer_len*2)-1,4,
                               (double *)fftw_fir_buf, fftw_fir_buf_fft,
                               FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
	fftwPlanCONV=fftw_plan_dft_c2r_2d((fir_buffer_len*2)-1,4,
                               (fftw_complex *)fftw_conv_buf_fft, (double*)fftw_conv_buf,
                               FFTW_ESTIMATE|FFTW_DESTROY_INPUT);

	printf("done with creating fft plans\n");

	/*
	 * create buffers for the FIR filters and for the processee
	 */

	printf("fir_buffer_len = %u\n\n",fir_buffer_len);
	printf("Sizeof AFormat = %lu\n\n",sizeof(AFormat));

	AFormat* data_fir1 = (AFormat*) malloc(sizeof(AFormat)*fir_buffer_len);
	memset(data_fir1,0,sizeof(AFormat)*fir_buffer_len);

	BFormat* data_fir2 = (BFormat*) malloc(sizeof(BFormat)*fir_buffer_len);
	memset(data_fir2,0,sizeof(BFormat)*fir_buffer_len);

	/*
	 * 2 blocks, all of the size of the FIR Buffer
	 */
	AFormat* data_in = (AFormat*) malloc(sizeof(AFormat)*fir_buffer_len*2);
	memset(data_in,0,sizeof(AFormat)*fir_buffer_len*2);

	/*
	 *  one output block
	 */

	AFormat* data_out1 = (AFormat*) malloc(sizeof(AFormat)*fir_buffer_len*2);
	memset(data_out1,0,sizeof(AFormat)*fir_buffer_len);

	BFormat* data_out2 = (BFormat*) malloc(sizeof(BFormat)*2*fir_buffer_len);
	memset(data_out2,0,sizeof(BFormat)*fir_buffer_len*2);

	BFormat* data_out3 = (BFormat*) malloc(sizeof(BFormat)*fir_buffer_len*2);
	memset(data_out3,0,sizeof(BFormat)*fir_buffer_len*2);

	/*
	 * load the FIR data into the buffers
	 */
	rc_fir1 = sf_readf_double (infile_fir1, (double*)data_fir1, fir_buffer_len);
	rc_fir2 = sf_readf_double (infile_fir2, (double*)data_fir2, fir_buffer_len);

	/* TEMP DEBUG 
	memset(data_fir1,(double)0,sizeof(BFormat)*fir_buffer_len);
	memset(data_fir2,(double)0,sizeof(BFormat)*fir_buffer_len);
	data_fir2[0].W=(double)1.0;
	data_fir2[0].X=(double)1.0;
	data_fir2[0].Y=(double)1.0;
	data_fir2[0].Z=(double)1.0;
	 TEMP DEBUG */

	sf_close(infile_fir1);
	sf_close(infile_fir2);

	/* PREP the Impulse Response FFT */

	memset((void*)fftw_fir_buf,0,(sizeof(BFormat)*fir_buffer_len*2)-1);
	memcpy(fftw_fir_buf,(double*)data_fir2,sizeof(BFormat)*fir_buffer_len);
	fftw_execute (fftwPlanFIR);

	/* OPEN the input file */

	if (! (infile = sf_open (infilename, SFM_READ, &sfinfo)))
	{	/* Open failed so print an error message. */
		printf ("Not able to open input file %s.\n", infilename) ;
		/* Print the error message from libsndfile. */
		puts (sf_strerror (NULL)) ;
		return 1 ;
	} ;
	sfinfo_print(&sfinfo,"INPUT");

	sfinfo_out.samplerate=sfinfo.samplerate;
	sfinfo_out.channels=sfinfo.channels;
	sfinfo_out.format=sfinfo.format;
	sfinfo_dbg=sfinfo_out;
	sfinfo_out.samplerate=48000;
	sfinfo_dbg.samplerate=48000;

	if (sfinfo.channels != 4)
	{	printf ("Not able to process more than 4 channels\n") ;
		return 1 ;
	} ;

	/* Open the output file. */

	if (! (outfile = sf_open (outfilename, SFM_WRITE, &sfinfo_out)))
	{	printf ("Not able to open output file %s.\n", outfilename) ;
		puts (sf_strerror (NULL)) ;
		return 1 ;
	} ;

	int src_buffer_ptr=0;

	/*
	 * convolution routine
	 *
	 * data_in* triple buffer for incoming datafile
	 * data_fir1* buffer for incoming fir1
	 * data_fir2* buffer for incoming fir2
	 */

	/*
	 * read in 2 blocks of the datafile
	 */

	//sf_readf_double(infile,(double*)data_in,fir_buffer_len*2);

	/*
         * keep a track of the blocks in the buffers;
	 */
	int iOut1=0,iOut2=0,iOut3=0,iIn1=3;
	int finishing=0,finishing1=0;
	int iCount=0;
	int iBlockCount=0;

	while (1) {
		printf("Working block %u\n\n",++iCount);
		
		// convovolve one block (len of fir1 buffer)

		/*
		* process PRE filter - FIR1
		*
		*/

		/*
		* 1.    [A][B] - convolution on block A with FIR [F]. This will output two Blocks [C][D]
		*       [C][D]-  
		* 2. [ ][E][F] - shift blocks E and F left...
		*    [E][F][G] - add block C to F and copy block D down to G
		*
		*/

		/*
		* [2][3] 	- data_in (convolution on first block, outputs two blocks)
		* [1][ ]	- data_out1	
		* [ ][1]	- data_out2
		* [ ]		- data_out3
		*/


		convovolveA(data_out1,data_in,data_fir1,fir_buffer_len);
		//debugBlock((BFormat *)data_out1,fir_buffer_len,"data_out1",iCount,&sfinfo_dbg); 

		printf("First Convo block %u\n",iCount);

		/*
		* [2][3][4] 	- data_in
		* [2]		- data_out1	
		* [ ][1]	- data_out2
		* [ ]		- data_out3
		*/

		/* move a block from the end of data_out2 to the start*/

		memcpy(data_out2,data_out2+fir_buffer_len,(sizeof(AFormat)*fir_buffer_len));
		//debugBlock((BFormat *)data_out2,fir_buffer_len,"data_out2",iCount,&sfinfo_dbg); 

		printf("First memcpy block %u\n",iCount);
		/*
		* [2][3][4] 	- data_in
		* [2]		- data_out1	
		* [1][1]	- data_out2
		* [ ]		- data_out3
		*/

		/* process Matrix filter
		* this is where we convert from A-Format to B-Format.
		*
		* samples in data_out2 and data_out3 will be in B-Format
		*
		* we put this block at the end of data_out2
		*/

		matrix(	data_out2 + fir_buffer_len, 
			data_out1, 
			fir_buffer_len);

		debugBlock((BFormat *)(data_out2+fir_buffer_len),fir_buffer_len,"matrix",iCount,&sfinfo_dbg); 

		printf("matrix block %u\n",iCount);

		/*
		* [2][3][4] 	- data_in
		* [2]		- data_out1	
		* [1][2]	- data_out2
		* [ ]		- data_out3
		*/
		/*
		* process POST filter - FIR2
		*/

		convovolveB(data_out3,data_out2,data_fir2,fir_buffer_len);

		debugBlock((BFormat *)(data_out3),fir_buffer_len,"data_out3",iCount,&sfinfo_dbg); 

		printf("second convo block %u\n",iCount);

		/*
		* [2][3][4] 	- data_in
		* [2]		- data_out1	
		* [1][2]	- data_out2
		* [1]		- data_out3
		*/

		// write first processed block to the output file

		if (iCount>3 ) {
			iBlockCount++;
			printf("writingF at 0x%p for %u bytes\n",(void*)data_out3,fir_buffer_len);
			sf_writef_double(outfile,(double*)data_out3,fir_buffer_len);
			sf_write_sync (outfile) ;
			printf("write block %u\n",iBlockCount);
		}


		// shift entire remaining buffer one block to the left

		memcpy(data_in,(data_in+fir_buffer_len),(sizeof(AFormat)*fir_buffer_len));

		printf("second memcpy %u\n",iCount);

		/*
		* [3][4][4] 	- data_in
		* [2]		- data_out1	
		* [1][2]	- data_out2
		* [1]		- data_out3
		*/

		// load in one more new block into the last block of the buffer
			
		//debugBlock((BFormat *)data_in,fir_buffer_len,"data_in",iCount,&sfinfo_dbg); 

		int rc=0;


		if (!finishing) {

			rc=sf_readf_double(infile,(double*)(data_in+fir_buffer_len),fir_buffer_len);
			printf("rc=%u\n\n",rc);

		} else {

			printf("we are finishing soon = %u\n\n",finishing);
			if (finishing==3) {
				printf("finishing on block %u\n",iCount);
				break;
			}

		}

		/*
		* [3][4][5] 	- data_in
		* [2]		- data_out1	
		* [1][2]	- data_out2
		* [1]		- data_out3
		*/

		//    if the block is incomplete or empty, fill it with zeros and flag future finish
		//
		// fir_buffer_len * 4 becuase 4 channels
		//

		if (rc<fir_buffer_len) finishing++;

		for (;rc<(fir_buffer_len);rc++) {
			(data_in+fir_buffer_len+rc)->FLU=(double)0;
			(data_in+fir_buffer_len+rc)->FRD=(double)0;
			(data_in+fir_buffer_len+rc)->BLD=(double)0;
			(data_in+fir_buffer_len+rc)->BRU=(double)0;
		}
	};

	/* Close input and output files. */
	
	printf("closing files \n\n");
	
	sf_close (outfile) ;
	sf_close (infile) ;

	fftw_destroy_plan( fftwPlanBuf);
	fftw_destroy_plan( fftwPlanFIR);
	fftw_destroy_plan( fftwPlanCONV);

	fftw_free(fftw_d_in_buf);
	fftw_free(fftw_fir_buf);
	fftw_free(fftw_conv_buf);
	fftw_free(fftw_olap_add);
	fftw_free(fftw_d_in_buf_fft);
	fftw_free(fftw_fir_buf_fft);
	fftw_free(fftw_conv_buf_fft);

	free( data_fir1 );
	free( data_fir2 );
	free( data_in );
	free( data_out1 );
	free( data_out2 );
	free( data_out3 );

	return 0 ;
} /* main */

static void
convovolveA(AFormat *d_out, AFormat *d_in, AFormat *fir, int count) {
	memcpy(d_out,d_in,(sizeof(AFormat)*count));

}

int fftcnt=0;

// count is the size of the FIR buffer.


static void 
prepFIR(BFormat *fir,int count) {
}




static void
convovolveB(BFormat *d_out, BFormat *d_in, BFormat *fir, int count) {
//	memcpy(d_out,d_in,(sizeof(AFormat)*count));
	//clear the buffer

	//debugBlock((BFormat *)(data_out2+fir_buffer_len),fir_buffer_len,"matrix",iCount,&sfinfo_dbg); 

	/*
	 * N=M+L-1
	 */
	memset((void*)fftw_d_in_buf,0,(sizeof(BFormat)*count*2)-1);
	//memset((void*)fftw_fir_buf,0,(sizeof(BFormat)*count*2)-1);
	memset((void*)fftw_conv_buf,0,(sizeof(BFormat)*count*2)-1);

	printf("memset_done\n");

	// fill half the input buffers with input data
	memcpy(fftw_d_in_buf,d_in,sizeof(BFormat)*count);
	//memcpy(fftw_fir_buf,fir,sizeof(BFormat)*count);

	debugBlock((BFormat *)(fftw_d_in_buf),(count*2)-1,"befrIFFT",fftcnt,&sfinfo_dbg); 

	printf("memcpy_done\n");

	// execute the FFT's
	fftw_execute (fftwPlanBuf);
	//fftw_execute (fftwPlanFIR);

	printf("fft_forward_done\n");
	
	// calculate scaling factor
	double scale=1.0f/(count*64);

	//Multiply the FFT's

	printf("count=%u\n",count);

	// half the block size due this is complex number field
	// rectangular multiplcation due this is complex numbers
	for (int i=0;i<(4*count*2)-4;i++) {
		fftw_conv_buf_fft[i][0]=(
			 ((fftw_d_in_buf_fft[i][0])*(fftw_fir_buf_fft[i][0]))
			-((fftw_d_in_buf_fft[i][1])*(fftw_fir_buf_fft[i][1]))
		);
		fftw_conv_buf_fft[i][1]=(
			 ((fftw_d_in_buf_fft[i][1])*(fftw_fir_buf_fft[i][0]))
			+((fftw_d_in_buf_fft[i][0])*(fftw_fir_buf_fft[i][1]))
		);
	}

	printf("multiply_done\n");

	// inverse FFT to get the convolution
	fftw_execute (fftwPlanCONV);

	// now fftw_conv_buf should have two blocks worth of data. 
	printf("fft_reverse done\n");

	debugBlock((BFormat *)(fftw_conv_buf),(count*2)-1,"aftrIFFT",fftcnt++,&sfinfo_dbg); 

	// we now add the first half of the convolution solution to the saved overlap block from the previous call
	// to this function. we also normalize here too.
	for (int i=0;i<count;i++) {
		if (i<count-1) {
			d_out[i].W=(fftw_conv_buf[i].W+fftw_olap_add[i].W)*scale;
			d_out[i].X=(fftw_conv_buf[i].X+fftw_olap_add[i].X)*scale;
			d_out[i].Y=(fftw_conv_buf[i].Y+fftw_olap_add[i].Y)*scale;
			d_out[i].Z=(fftw_conv_buf[i].Z+fftw_olap_add[i].Z)*scale;
		} else {
			d_out[i].W=(fftw_conv_buf[i].W)*scale;
			d_out[i].X=(fftw_conv_buf[i].X)*scale;
			d_out[i].Y=(fftw_conv_buf[i].Y)*scale;
			d_out[i].Z=(fftw_conv_buf[i].Z)*scale;
		}
	}

	// We need to copy out the last block and save it for next time.

	memcpy(fftw_olap_add,(void*)(fftw_conv_buf+count),(sizeof(BFormat)*count)-4);
}
	
/*
 * OK THIS WORKS!
 *
 W' = FLU+FRD+BLD+BRU
 X' = FLU+FRD-BLD-BRU
 Y' = FLU-FRD+BLD-BRU
 Z' = FLU-FRD-BLD+BRU
 *
 FLU:0
 FRD:1
 BLD:2
 BRU:3
 W':0
 X':1
 Y':2
 Z':3
 *
 *
 */
static void
matrix(BFormat *d_out, AFormat *d_in, int count) {
//	memcpy(d_out,d_in,(sizeof(AFormat)*count));
	//memset(d_out,(double)50,(sizeof(double)*channels*count));

	int ptr;

	for (ptr = 0 ; ptr < count ; ptr += 1) {
		d_out->W=(d_in->FLU+d_in->FRD+d_in->BLD+d_in->BRU)/2.0;
		d_out->X=(d_in->FLU+d_in->FRD-d_in->BLD-d_in->BRU)/2.0;
		d_out->Y=(d_in->FLU-d_in->FRD+d_in->BLD-d_in->BRU)/2.0;
		d_out->Z=(d_in->FLU-d_in->FRD-d_in->BLD+d_in->BRU)/2.0;
		d_in++;d_out++;
	}
}

static void
sfinfo_print(SF_INFO* sfinfo, char* title) {
	printf("%s loaded\n",title);
	printf("Number of channels: %u\n",sfinfo->channels);
	printf("samplerate: %u\n",sfinfo->samplerate);
	printf("format: %u\n",sfinfo->format);

	double secs=(double)sfinfo->frames/(double)sfinfo->samplerate;

	printf("length: %f\n\n",secs);

}


int dbgCnt=0;
static void 
debugBlock(BFormat *buffer,int bufsize,char* title,int block,SF_INFO* sfi) {
	SNDFILE * snd_debugfile;
	char sbuf[50];

	dbgCnt++;

	sprintf(sbuf,"dbg-%s-%u-%03u.wav",title,block,dbgCnt);
	
	
	if (! (snd_debugfile = sf_open (sbuf, SFM_WRITE, sfi)))
	{	printf ("Not able to open output file %s.\n", sbuf) ;
		puts (sf_strerror (NULL)) ;
	} ;

	sf_writef_double(snd_debugfile,(double*)buffer,bufsize);

	sf_close(snd_debugfile);

}
