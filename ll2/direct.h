typedef struct {
	float left;
	float right;
} CFormat;

void convo_direct(	const 	CFormat Signal[ /* SignalLen */ ], size_t SignalLen,
    				const 	CFormat Kernel[ /* KernelLen */ ], size_t KernelLen,
    						CFormat Result[ /* SignalLen + KernelLen - 1 */ ]) ;

