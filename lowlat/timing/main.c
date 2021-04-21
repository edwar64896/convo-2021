#include <fftw3.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>

clock_t tDirect_start,tDirect_diff,tFftw_start,tFftw_diff;

pthread_t ptDirect;
pthread_t ptFftw;

void*
fDirect(void * args) {
	sleep(3);
	/*
	 * stuff goes in here.
	 */
	tDirect_diff=clock()-tDirect_start;
	tDirect_start=clock();
	printf("done fDirect\n");
	pthread_exit(NULL);
	return 0;
}

void*
fFftw(void * args) {

	sleep(2);
	/*
	 * stuff goes in here.
	 */
	tFftw_diff=clock()-tFftw_start;
	tFftw_start=clock();
	printf("done fFftw\n");
	pthread_exit(NULL);
	return 0;
}


int
main(int argc, char ** argv) {
	printf("fftw timing checker\n");

	while (1) {
		int args=1;
		pthread_create (&ptDirect,NULL,fDirect,(void*) &args);
		pthread_create (&ptFftw,NULL,fFftw,(void*) &args);

		/*
		 * wait for both threads to end before doing it all again.
		 */

		pthread_join(ptDirect,NULL);
		pthread_join(ptFftw,NULL);

		printf("fFtw_diff=%lu\n",tFftw_diff);
		printf("fDirect_diff=%lu\n",tDirect_diff);

		int tFftw_msec = tFftw_diff * 1000 / CLOCKS_PER_SEC;
		int tDirect_msec = tDirect_diff * 1000 / CLOCKS_PER_SEC;

		//int tFftw_msec = tFftw_diff ;
		//int tDirect_msec = tDirect_diff;

		//float rt=(float)iBlockSize/48000.0;
		//float ratio=(rt*1000.0)/(float)msec;
		
		//printf("Time taken %d seconds %d milliseconds to process %f (s) audio = %f times real time\n", msec/1000, msec%1000,rt,ratio);

		printf("Fftw Time taken %d seconds %d milliseconds to process\n", tFftw_msec/1000, tFftw_msec%1000);
		printf("Direct Time taken %d seconds %d milliseconds to process\n", tDirect_msec/1000, tDirect_msec%1000);

	}


	exit(0);
}
