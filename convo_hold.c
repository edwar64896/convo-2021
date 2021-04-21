
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
	

