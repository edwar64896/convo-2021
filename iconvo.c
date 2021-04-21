
#include	<stdio.h>
#include	<string.h>
#include	<stdlib.h>
#include 	<fftw3.h>
#include 	"ioacon.h"

/* Include this header file to use functions from libsndfile. */
#include	<sndfile.h>

static void matrix(BFormat *d_out, AFormat *d_in, int count) ;
static void sfinfo_print(SF_INFO* sfinfo,char*);
static void debugBlock(BFormat*,int ,char* ,int ,SF_INFO* ) ;
static void bufPrint(CFormat* buf);

SF_INFO		sfinfo_input ;
SF_INFO		sfinfo_output ;
SF_INFO		sfinfo_filter ;

iConvo* hReverse;

int
main (void) {

	/* A SNDFILE is very much like a FILE in the Standard C library. The
	** sf_open function return an SNDFILE* pointer when they sucessfully
	** open the specified file.
	*/
	SNDFILE	*inputfile, *outputfile, *filterfile;

	/* A pointer to an SF_INFO struct is passed to sf_open.
	** On read, the library fills this struct with information about the file.
	** On write, the struct must be filled in before calling sf_open.
	*/

	int		rc_input ;
	int		rc_output ;


	const char	*inputfilename = "input.wav" ;
	const char	*outputfilename = "output.wav" ;
	const char	*filterfilename = "impulseresponse.wav" ;

	/* The SF_INFO struct must be initialized before using it.
	*/
	memset (&sfinfo_input, 0, sizeof (SF_INFO)) ;
	memset (&sfinfo_output, 0, sizeof (SF_INFO)) ;
	memset (&sfinfo_filter, 0, sizeof (SF_INFO)) ;

	/*
	 * load in the interleaved FIR impulse response filters
	 * prior to convolution.
	 */

	if (! (inputfile = sf_open (inputfilename, SFM_READ, &sfinfo_input)))
	{	/* Open failed so print an error message. */
		printf ("Not able to open input file %s.\n", inputfilename) ;
		/* Print the error message from libsndfile. */
		puts (sf_strerror (NULL)) ;
		return 1 ;
	} ;
	sfinfo_print(&sfinfo_input,"INPUT");
	if (sfinfo_input.channels != 4)
	{	printf ("INPUT FILE needs 4 channels\n") ;
		return 1 ;
	} ;
	if (! (outputfile = sf_open (outputfilename, SFM_READ, &sfinfo_output)))
	{	/* Open failed so print an error message. */
		printf ("Not able to open output file %s.\n", outputfilename) ;
		/* Print the error message from libsndfile. */
		puts (sf_strerror (NULL)) ;
		return 1 ;
	};
	sfinfo_print(&sfinfo_output,"OUTPUT");
	if (sfinfo_output.channels != 4)
	{	printf ("OUTPUT FILE needs 4 channels\n") ;
		return 1 ;
	} ;

	/*
	 * Check the FIR buffers match
	 */
	if (sfinfo_input.frames != sfinfo_output.frames) {
		printf("lengths do not match\n");
	//	return 1;
	}
	if (sfinfo_input.samplerate != sfinfo_output.samplerate) {
		printf("samplerates do not match\n");
		return 1;
	}
	if (sfinfo_input.channels != sfinfo_output.channels) {
		printf("channel counts do not match\n");
		return 1;
	}
	if (sfinfo_input.format != sfinfo_output.format) {
		printf("formats do not match\n");
		return 1;
	}

	/*
	 * Identify the length of the FIR filter buffers
	 */
	int buffer_len=sfinfo_input.frames;
	
	hReverse=iconvo_init(buffer_len);

	/*
	 * load the FIR data into the buffers
	 */

	iconvo_initInput(hReverse);
	iconvo_initOutput(hReverse);

	int blkSize2=hReverse->blockSize/2;
	
	rc_input = sf_readf_double (inputfile, (double*)(hReverse->input+blkSize2), hReverse->blockSize);
	rc_output = sf_readf_double (outputfile, (double*)(hReverse->output+blkSize2), hReverse->blockSize);

	/* PREP the Impulse Response FFT */

	sf_close(inputfile);
	sf_close(outputfile);


	sfinfo_filter.samplerate=sfinfo_input.samplerate;
	sfinfo_filter.channels=sfinfo_input.channels;
	sfinfo_filter.format=sfinfo_input.format;

	if (! (filterfile = sf_open (filterfilename, SFM_WRITE, &sfinfo_filter)))
	{	/* Open failed so print an error message. */
		printf ("Not able to open filter output file %s.\n", filterfilename) ;
		/* Print the error message from libsndfile. */
		puts (sf_strerror (NULL)) ;
		return 1 ;
	} ;


	iconvo_devolve(hReverse);

	sf_writef_double(filterfile,(double*)hReverse->filterOut,hReverse->blockSize);
	sf_write_sync (filterfile) ;

	bufPrint(hReverse->filter);

	printf("closing files \n\n");
	
	sf_close (filterfile) ;

	iconvo_cleanup(hReverse);
	
	return 0 ;
} /* main */

int fftcnt=0;

static void
bufPrint(CFormat* buf) {
	for (int i =0 ;i<10;i++) {
		printf("sample=%u - values: %f,%f,%f,%f\n",i,buf[i].A,buf[i].B,buf[i].C,buf[i].D);
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
