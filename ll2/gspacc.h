#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "pa_util.h"
#include "pa_ringbuffer.h"

#ifndef _GSPACC_H
#define _GSPACC_H

//typedef void (*gspacc_cb) (void *, size_t, void *, float);
// accumulator callback function
typedef void (*gspacc_cb) (void *, size_t, void *);

// accumulator data structure
typedef struct gspacc_t_struct {
	PaUtilRingBuffer	*RingBuffer;
	void **	RingBufferData;
	void*	outBuf;
	size_t elementSize;
	size_t numElements;
	size_t depth;
} gspacc_t;

static __attribute__((noinline)) unsigned next_pow2(unsigned x)
{
	if (x <= 2) return x;
	
	return (1ULL << 32) >> __builtin_clz(x - 1);
}

//initialization function
gspacc_t * gspacc_initialize(size_t elementSize, size_t numElements, size_t depth) ;

// accumulator write
void gspacc_write (gspacc_t * handle, void * data, size_t accBlock, size_t numBlocks) ;

//int gspacc_read (gspacc_t * handle, size_t accBlock, gspacc_cb accCb, void * dest, float) ;

//accumulator read
int gspacc_read (gspacc_t * handle, size_t accBlock, gspacc_cb accCb, void * dest) ;

//accumulator cleanup.
void gspacc_cleanup (gspacc_t * handle) ;

#endif // _GSPACC_H
