
#include	<stdio.h>
#include	<string.h>
#include	<stdlib.h>
#include 	<fftw3.h>
#include "oacon.h"

/* Include this header file to use functions from libsndfile. */
#include	<sndfile.h>

static void matrix(BFormat *d_out, AFormat *d_in, int count) ;
static void sfinfo_print(SF_INFO* sfinfo,char*);
static void debugBlock(BFormat*,int ,char* ,int ,SF_INFO* ) ;

SF_INFO		sfinfo ;
SF_INFO		sfinfo_out ;
SF_INFO		sfinfo_dbg ;
SF_INFO		sfinfo_fir1 ;
SF_INFO		sfinfo_fir2 ;

Convo* firHandle1;
Convo* firHandle2;

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
	int dat_buffer_len=fir_buffer_len*40;
	
	firHandle1=convo_init(fir_buffer_len,dat_buffer_len,"FIR1");
	firHandle2=convo_init(fir_buffer_len,dat_buffer_len,"FIR2");

	/* * 2 blocks, all of the size of the FIR Buffer */

	AFormat* data_in = (AFormat*) malloc(sizeof(AFormat)*dat_buffer_len*2);
	memset(data_in,0,sizeof(AFormat)*dat_buffer_len*2);

	/*
	 *  one output block
	 */

	AFormat* data_out1 = (AFormat*) malloc(sizeof(AFormat)*dat_buffer_len*2);
	memset(data_out1,0,sizeof(AFormat)*dat_buffer_len);

	BFormat* data_out2 = (BFormat*) malloc(sizeof(BFormat)*2*dat_buffer_len);
	memset(data_out2,0,sizeof(BFormat)*dat_buffer_len*2);

	BFormat* data_out3 = (BFormat*) malloc(sizeof(BFormat)*dat_buffer_len*2);
	memset(data_out3,0,sizeof(BFormat)*dat_buffer_len*2);

	/*
	 * load the FIR data into the buffers
	 */

	convo_initFilter(firHandle1);
	convo_initFilter(firHandle2);

	rc_fir1 = sf_readf_double (infile_fir1, (double*)firHandle1->fir, firHandle1->filterSize);
	rc_fir2 = sf_readf_double (infile_fir2, (double*)firHandle2->fir, firHandle2->filterSize);

	//convo_testFilter(firHandle1);
	//convo_testFilter(firHandle2);

	/* PREP the Impulse Response FFT */

	convo_buildFilter(firHandle1);
	convo_buildFilter(firHandle2);

	sf_close(infile_fir1);
	sf_close(infile_fir2);

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
	{	printf ("Input file should have 4 channels. This one doesn't. :-(\n") ;
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

	setbuf(stdout,NULL);
	while (1) {
		++iCount;
		printf("#");
		
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

		//debugBlock((BFormat *)data_in,dat_buffer_len,"data_in",iCount,&sfinfo_dbg); 
		convo_volve(firHandle1,(CFormat*)data_out1,(CFormat*)data_in);
		debugBlock((BFormat *)data_out1,dat_buffer_len,"data_out1",iCount,&sfinfo_dbg); 

		/*
		* [2][3][4] 	- data_in
		* [2]		- data_out1	
		* [ ][1]	- data_out2
		* [ ]		- data_out3
		*/

		/* move a block from the end of data_out2 to the start*/

		memcpy(data_out2,data_out2+dat_buffer_len,(sizeof(AFormat)*dat_buffer_len));
		//debugBlock((BFormat *)data_out2,fir_buffer_len,"data_out2",iCount,&sfinfo_dbg); 

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

		matrix(	data_out2 + dat_buffer_len, 
			data_out1, 
			dat_buffer_len);

		debugBlock((BFormat *)(data_out2+dat_buffer_len),dat_buffer_len,"matrix",iCount,&sfinfo_dbg); 


		/*
		* [2][3][4] 	- data_in
		* [2]		- data_out1	
		* [1][2]	- data_out2
		* [ ]		- data_out3
		*/

		/*
		* process POST filter - FIR2
		*/

		convo_volve(firHandle2,(CFormat*)data_out3,(CFormat*)data_out2);

		//debugBlock((BFormat *)(data_out3),dat_buffer_len,"data_out3",iCount,&sfinfo_dbg); 


		/*
		* [2][3][4] 	- data_in
		* [2]		- data_out1	
		* [1][2]	- data_out2
		* [1]		- data_out3
		*/

		// write first processed block to the output file

		if (iCount>3 ) {
			iBlockCount++;
			sf_writef_double(outfile,(double*)data_out3,dat_buffer_len);
			sf_write_sync (outfile) ;
		}


		// shift entire remaining buffer one block to the left

		memcpy(data_in,(data_in+dat_buffer_len),(sizeof(AFormat)*dat_buffer_len));


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

			rc=sf_readf_double(infile,(double*)(data_in+dat_buffer_len),dat_buffer_len);

		} else {

			if (finishing==3) {
				printf("\n\n");
				printf("finishing on block %u\n",iCount);
				break;
			}

		}

		/*
		/ [3][4][5] 	- data_in
		* [2]		- data_out1	
		* [1][2]	- data_out2
		* [1]		- data_out3
		*/

		//    if the block is incomplete or empty, fill it with zeros and flag future finish
		//
		// fir_buffer_len * 4 becuase 4 channels
		//

		if (rc<dat_buffer_len) finishing++;

		for (;rc<(fir_buffer_len);rc++) {
			(data_in+dat_buffer_len+rc)->FLU=(double)0;
			(data_in+dat_buffer_len+rc)->FRD=(double)0;
			(data_in+dat_buffer_len+rc)->BLD=(double)0;
			(data_in+dat_buffer_len+rc)->BRU=(double)0;
		}
	};

	/* Close input and output files. */
	
	printf("closing files \n\n");
	
	sf_close (outfile) ;
	sf_close (infile) ;

	convo_cleanup(firHandle1);
	convo_cleanup(firHandle2);
	
	free( data_in );
	free( data_out1 );
	free( data_out2 );
	free( data_out3 );

	return 0 ;
} /* main */

int fftcnt=0;

/*
 * Convert filtered A Format data to unfiltered B Format data
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
 */
static void
matrix(BFormat *d_out, AFormat *d_in, int count) {

	int ptr;

	for (ptr = 0 ; ptr < count ; ptr += 1) {
		//d_out->W=(d_in->FLU+d_in->FRD+d_in->BLD+d_in->BRU)/2.0;
		d_out->W=(d_in->FLU+d_in->FRD+d_in->BLD+d_in->BRU)/2.0;
		d_out->X=(d_in->FLU+d_in->FRD-d_in->BLD-d_in->BRU)/2.0;
		d_out->Y=(d_in->FLU-d_in->FRD+d_in->BLD-d_in->BRU)/2.0;
		d_out->Z=(d_in->FLU-d_in->FRD-d_in->BLD+d_in->BRU)/2.0;
		d_in++;d_out++;
	}
}

/*
 * debug the SF Info structure
 */
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
