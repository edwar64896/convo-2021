#include <stddef.h> 
#include <stdio.h>

struct CFormatStruct {
	double A;
	double B;
	double C;
	double D;
};

typedef struct CFormatStruct CFormat;

void convo_direct(	const 	CFormat Signal[ /* SignalLen */ ], size_t SignalLen,
    				const 	CFormat Kernel[ /* KernelLen */ ], size_t KernelLen,
    						CFormat Result[ /* SignalLen + KernelLen - 1 */ ]) {
	size_t n;

   	for (n = 0; n < SignalLen + KernelLen - 1; n++) {
        size_t kmin, kmax, k;

        Result[n].A = 0;
        Result[n].B = 0;
        Result[n].C = 0;
        Result[n].D = 0;

       	kmin = (n >= KernelLen - 1) ? n - (KernelLen - 1) : 0;
       	kmax = (n < SignalLen - 1) ? n : SignalLen - 1;

        for (k = kmin; k <= kmax; k++) {
        	Result[n].A += Signal[k].A * Kernel[n - k].A;
        	Result[n].B += Signal[k].B * Kernel[n - k].B;
        	Result[n].C += Signal[k].C * Kernel[n - k].C;
        	Result[n].D += Signal[k].D * Kernel[n - k].D;
   		}
	}
}
