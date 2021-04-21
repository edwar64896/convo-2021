#include	<stdio.h>
#include	<string.h>
#include	<stdlib.h>
#include	<getopt.h>

/* Include this header file to use functions from libsndfile. */
#include	<sndfile.h>


static int verbose_flag;
const char* version="0.2alpha";

char szInput[1024];
SF_INFO		sfinfo ;

static void
sfinfo_print(SF_INFO* sfinfo, char* title) {
	printf("%s loaded\n",title);
	printf("Number of channels: %u\n",sfinfo->channels);
	printf("samplerate: %u\n",sfinfo->samplerate);
	switch (sfinfo->format & SF_FORMAT_TYPEMASK) {
		case SF_FORMAT_WAV:
			printf("SF_FORMAT_WAV\n");	
			break;
		case SF_FORMAT_AIFF:
			printf("SF_FORMAT_AIFF\n");	
			break;
		case SF_FORMAT_AU:
			printf("SF_FORMAT_AU\n");	
			break;
		case SF_FORMAT_RAW:
			printf("SF_FORMAT_RAW\n");	
			break;
		case SF_FORMAT_PAF:
			printf("SF_FORMAT_PAF\n");	
			break;
		case SF_FORMAT_SVX:
			printf("SF_FORMAT_SVX\n");	
			break;
		case SF_FORMAT_NIST:
			printf("SF_FORMAT_NIST\n");	
			break;
		case SF_FORMAT_VOC:
			printf("SF_FORMAT_VOC\n");	
			break;
		case SF_FORMAT_IRCAM:
			printf("SF_FORMAT_IRCAM\n");	
			break;
		case SF_FORMAT_W64:
			printf("SF_FORMAT_W64\n");	
			break;
		case SF_FORMAT_MAT4:
			printf("SF_FORMAT_MAT4\n");	
			break;
		case SF_FORMAT_MAT5:
			printf("SF_FORMAT_MAT5\n");	
			break;
		case SF_FORMAT_PVF:
			printf("SF_FORMAT_PVF\n");	
			break;
		case SF_FORMAT_XI:
			printf("SF_FORMAT_XI\n");	
			break;
		case SF_FORMAT_HTK:
			printf("SF_FORMAT_HTK\n");	
			break;
		case SF_FORMAT_SDS:
			printf("SF_FORMAT_SDS\n");	
			break;
		case SF_FORMAT_AVR:
			printf("SF_FORMAT_AVR\n");	
			break;
		case SF_FORMAT_WAVEX:
			printf("SF_FORMAT_WAVEX\n");	
			break;
		case SF_FORMAT_SD2:
			printf("SF_FORMAT_SD2\n");	
			break;
		case SF_FORMAT_FLAC:
			printf("SF_FORMAT_FLAC\n");	
			break;
		case SF_FORMAT_CAF:
			printf("SF_FORMAT_CAF\n");	
			break;
		case SF_FORMAT_WVE:
			printf("SF_FORMAT_WVE\n");	
			break;
		case SF_FORMAT_OGG:
			printf("SF_FORMAT_OGG\n");	
			break;
		case SF_FORMAT_MPC2K:
			printf("SF_FORMAT_MPC2K\n");	
			break;
		case SF_FORMAT_RF64:
			printf("SF_FORMAT_RF64\n");	
			break;
	}
	switch (sfinfo->format & SF_FORMAT_SUBMASK) {
		case SF_FORMAT_PCM_S8:
			printf("SF_FORMAT_PCM_S8\n");	
			break;
		case SF_FORMAT_PCM_16:
			printf("SF_FORMAT_PCM_16\n");	
			break;
		case SF_FORMAT_PCM_24:
			printf("SF_FORMAT_PCM_24\n");	
			break;
		case SF_FORMAT_PCM_32:
			printf("SF_FORMAT_PCM_32\n");	
			break;
		case SF_FORMAT_PCM_U8:
			printf("SF_FORMAT_PCM_U8\n");	
			break;
		case SF_FORMAT_FLOAT:
			printf("SF_FORMAT_FLOAT\n");	
			break;
		case SF_FORMAT_DOUBLE:
			printf("SF_FORMAT_DOUBLE\n");	
			break;
		case SF_FORMAT_ULAW:
			printf("SF_FORMAT_ULAW\n");	
			break;
		case SF_FORMAT_ALAW:
			printf("SF_FORMAT_ALAW\n");	
			break;
		case SF_FORMAT_IMA_ADPCM:
			printf("SF_FORMAT_IMA_ADPCM\n");	
			break;
		case SF_FORMAT_MS_ADPCM:
			printf("SF_FORMAT_MS_ADPCM\n");	
			break;
		case SF_FORMAT_GSM610:
			printf("SF_FORMAT_GSM610\n");	
			break;
		case SF_FORMAT_VOX_ADPCM:
			printf("SF_FORMAT_VOX_ADPCM\n");	
			break;
		case SF_FORMAT_G721_32:
			printf("SF_FORMAT_G721_32\n");	
			break;
		case SF_FORMAT_G723_24:
			printf("SF_FORMAT_G723_24\n");	
			break;
		case SF_FORMAT_G723_40:
			printf("SF_FORMAT_G723_40\n");	
			break;
		case SF_FORMAT_DWVW_12:
			printf("SF_FORMAT_DWVW_12\n");	
			break;
		case SF_FORMAT_DWVW_24:
			printf("SF_FORMAT_DWVW_24\n");	
			break;
		case SF_FORMAT_DWVW_16:
			printf("SF_FORMAT_DWVW_16\n");	
			break;
		case SF_FORMAT_DWVW_N:
			printf("SF_FORMAT_DWVW_N\n");	
			break;
		case SF_FORMAT_DPCM_8:
			printf("SF_FORMAT_DPCM_8\n");	
			break;
		case SF_FORMAT_DPCM_16:
			printf("SF_FORMAT_DPCM_16\n");	
			break;
		case SF_FORMAT_VORBIS:
			printf("SF_FORMAT_VORBIS\n");	
			break;
		}

	double secs=(double)sfinfo->frames/(double)sfinfo->samplerate;

	printf("frames: %lld\n",sfinfo->frames);
	printf("length: %f\n",secs);

}

void 
usage() {
	printf("Check v0.1  usage...\n\n");
	printf("Will check format of a WAV file for you \n\n");
	printf("-i input file name\n");
	printf("--verbose = SFPrint output (working?)\n");
	printf("-h help = this\n");
}

int
main (int argc, char **argv) {

	int c,index;
	SNDFILE *infile;

	while (1) {
		static struct option long_options[] =
		{
			{"verbose",	no_argument, 		&verbose_flag, 	1},
			{"input",	required_argument,	0,				'i'},
			{"help",	no_argument,		0,				'h'},
			{0,0,0,0}
		};
		int option_index=0;
		c=getopt_long(argc,argv,":hi:",long_options,&option_index);
		if (c==-1)
			break;
		switch(c) {
			case 0:
				if (long_options[option_index].flag !=0)
					break;
				printf("option %s", long_options[option_index].name);
				if (optarg)
					printf(" with arg %s", optarg);
				printf("\n");
				break;
			case 'i':
				strncpy(szInput,optarg,1024);
				printf(" input file %s\n",szInput);
				break;
			case 'h':
				usage();
				exit(0);
				break;
			default:
				break;
		}
	};

	const char	*infilename = szInput ;

	memset (&sfinfo, 0, sizeof (sfinfo)) ;

	if (! (infile = sf_open (infilename, SFM_READ, &sfinfo)))
	{	/* Open failed so print an error message. */
		printf ("Not able to open input file %s.\n", infilename) ;
		/* Print the error message from libsndfile. */
		puts (sf_strerror (NULL)) ;
		return 1 ;
	} ;
	if (verbose_flag)
		sfinfo_print(&sfinfo,"INPUT");

	sf_count_t cnt=0;
	sf_count_t items=0;

	int *buf;
	buf=malloc(sizeof(int) * sfinfo.channels);

	int framescount=0;

	while (1) {
		cnt=sf_readf_int(infile,buf,1);
		if (cnt) framescount++;
		else break;
	}
	printf("frames read = %u\n",framescount);
	return 0;

}
