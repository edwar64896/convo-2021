#include <stddef.h> 
#include <stdio.h>

struct CFormatStruct {
	float left;
	float right;
};

typedef struct CFormatStruct CFormat;

/*
 *
 * Convolution using direct formula rather than going through FFT.
 *
 * This will be quicker for small buffers than using FFT
 *
 */
void convo_direct(	const 	CFormat Signal[ /* SignalLen */ ], size_t SignalLen,
    								const 	CFormat Kernel[ /* KernelLen */ ], size_t KernelLen,
    												CFormat Result[ /* SignalLen + KernelLen - 1 */ ]) {
	size_t n;

   	for (n = 0; n < SignalLen + KernelLen - 1; n++) {
        size_t kmin, kmax, k;

        Result[n].left = 0;
        Result[n].right = 0;

       	kmin = (n >= KernelLen - 1) ? n - (KernelLen - 1) : 0;
       	kmax = (n < SignalLen - 1) ? n : SignalLen - 1;

        for (k = kmin; k <= kmax; k++) {
        	Result[n].left += Signal[k].left * Kernel[n - k].left;
        	Result[n].right += Signal[k].right * Kernel[n - k].right;
   		}
	}
}
