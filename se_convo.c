    /*
    */
    
    #include	<stdio.h>
    #include	<string.h>
    #include	<stdlib.h>
    
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
    static void process_data (double *data, int count, int channels) ;
    static void convovolveA(AFormat *d_out, AFormat *d_in, AFormat* fir, int count) ;
    static void convovolveB(BFormat *d_out, BFormat *d_in, BFormat* fir, int count) ;
    static void matrix(BFormat *d_out, AFormat *d_in, int count) ;
    static void sfinfo_print(SF_INFO* sfinfo,char*);
    static void debugBlock(BFormat*,int ,char* ,int ,SF_INFO* ) ;
    
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
    	SF_INFO		sfinfo ;
    	SF_INFO		sfinfo_out ;
    	SF_INFO		sfinfo_dbg ;
    	SF_INFO		sfinfo_fir1 ;
    	SF_INFO		sfinfo_fir2 ;
    
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
    	AFormat* data_in = (AFormat*) malloc(sizeof(AFormat)*fir_buffer_len*3);
    	memset(data_in,0,sizeof(AFormat)*fir_buffer_len*3);
    
    	/*
    	 *  one output block
    	 */
    
    	AFormat* data_out1 = (AFormat*) malloc(sizeof(AFormat)*fir_buffer_len);
    	memset(data_out1,0,sizeof(AFormat)*fir_buffer_len);
    
    	BFormat* data_out2 = (BFormat*) malloc(sizeof(BFormat)*2*fir_buffer_len);
    	memset(data_out2,0,sizeof(BFormat)*fir_buffer_len*2);
    
    	BFormat* data_out3 = (BFormat*) malloc(sizeof(BFormat)*fir_buffer_len);
    	memset(data_out3,0,sizeof(BFormat)*fir_buffer_len);
    
    	/*
    	 * load the FIR data into the buffers
    	 */
    	rc_fir1 = sf_readf_double (infile_fir1, (double*)data_fir1, fir_buffer_len);
    	rc_fir2 = sf_readf_double (infile_fir2, (double*)data_fir2, fir_buffer_len);
    
    	sf_close(infile_fir1);
    	sf_close(infile_fir2);
    
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
    		* [2][3][4] 	- data_in
    		* [1]		- data_out1	
    		* [ ][1]	- data_out2
    		* [ ]		- data_out3
    		*/
    
    
    		convovolveA(data_out1,data_in,data_fir1,fir_buffer_len);
    		//debugBlock((BFormat *)data_out1,fir_buffer_len,"data_out1",iCount,&sfinfo); 
    
    		printf("First Convo block %u\n",iCount);
    
    		/*
    		* [2][3][4] 	- data_in
    		* [2]		- data_out1	
    		* [ ][1]	- data_out2
    		* [ ]		- data_out3
    		*/
    
    		/* move a block from the end of data_out2 to the start*/
    
    		memcpy(data_out2,data_out2+fir_buffer_len,(sizeof(AFormat)*fir_buffer_len));
    		//debugBlock((BFormat *)data_out2,fir_buffer_len,"data_out2",iCount,&sfinfo); 
    
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
    
    		//debugBlock((BFormat *)(data_out2+fir_buffer_len),fir_buffer_len,"matrix",iCount,&sfinfo); 
    
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
    			//sf_writef_double(outfile,(double*)data_out3,fir_buffer_len);
    			//sf_write_sync (outfile) ;
    			printf("write block %u\n",iBlockCount);
    			if (iCount==6)
    				sf_close(outfile);
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
    			
    		//debugBlock((BFormat *)data_in,fir_buffer_len,"data_in",iCount,&sfinfo); 
    
    		int rc=0;
    
    
    		if (!finishing) {
    
    			rc=sf_readf_double(infile,(double*)(data_in+fir_buffer_len),fir_buffer_len);
    			sf_writef_double(outfile,(double*)(data_in+fir_buffer_len),fir_buffer_len);
    			sf_write_sync(outfile);
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
    			(data_in+(2*fir_buffer_len)+rc)->FLU=(double)0;
    			(data_in+(2*fir_buffer_len)+rc)->FRD=(double)0;
    			(data_in+(2*fir_buffer_len)+rc)->BLD=(double)0;
    			(data_in+(2*fir_buffer_len)+rc)->BRU=(double)0;
    		}
    	};
    
    	/* Close input and output files. */
    	
    	printf("closing files \n\n");
    	
    	sf_close (outfile) ;
    	sf_close (infile) ;
    
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
    
    /*
    	int chan_in,chan_fir,oset_in,oset_fir;
    
    	for (chan_in = 0 ; chan_in < channels ; chan_in++) {
    		printf("chan_in:%u\n",chan_in);
    		for (oset_in = chan_in ; oset_in < count ; oset_in += channels) {
    			for (oset_fir = chan_fir ; oset_fir < count ; oset_fir += channels) {
    				d_out [oset_in] = d_out[oset_in] + (fir[oset_fir] * d_in[oset_fir] );
    			}
    		}
    	}
    */
    }
    static void
    convovolveB(BFormat *d_out, BFormat *d_in, BFormat *fir, int count) {
    	memcpy(d_out,d_in,(sizeof(BFormat)*count));
    
    /*
    	int chan_in,chan_fir,oset_in,oset_fir;
    
    	for (chan_in = 0 ; chan_in < channels ; chan_in++) {
    		printf("chan_in:%u\n",chan_in);
    		for (oset_in = chan_in ; oset_in < count ; oset_in += channels) {
    			for (oset_fir = chan_fir ; oset_fir < count ; oset_fir += channels) {
    				d_out [oset_in] = d_out[oset_in] + (fir[oset_fir] * d_in[oset_fir] );
    			}
    		}
    	}
    */
    }
    
    /*
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
    	memcpy(d_out,d_in,(sizeof(AFormat)*count));
    	//memset(d_out,(double)50,(sizeof(double)*channels*count));
    
    /*
    	int ptr;
    
    	for (ptr = 0 ; ptr < count ; ptr += 1) {
    		d_out->W=(d_in->FLU+d_in->FRD+d_in->BLD+d_in->BRU)/2.0;
    		d_out->X=(d_in->FLU+d_in->FRD-d_in->BLD-d_in->BRU)/2.0;
    		d_out->Y=(d_in->FLU-d_in->FRD+d_in->BLD-d_in->BRU)/2.0;
    		d_out->Z=(d_in->FLU-d_in->FRD-d_in->BLD+d_in->BRU)/2.0;
    	}
    */
    
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
