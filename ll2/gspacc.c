#include <stdio.h>
#include "gspacc.h"



// ringbuffer handles
// ringbuffer data blocks

// write blocks

// read blocks - accumulate results

// initialize
//	element size (bytes)
//	number of elements (bytes)	

// cleanup


/*
 * Accumulator functions.
 *
 * Accumulator is a function that adds blocks of audio in a convolution reverb scenario.
 * 
 * initialize the accumulator
 */
gspacc_t * gspacc_initialize(size_t elementSize, size_t numElements, size_t depth) {

	// allocate the handle for the data structure
	gspacc_t * handle=malloc(sizeof(gspacc_t));

	// debug log
	printf("elementsize=%zu, numElements=%zu, nep2=%zu, depth=%zu \n",elementSize, numElements, numElementP2, depth);

	// store the parameters
	handle->elementSize=elementSize; //block size in bytes (sizeof (one frame) * number of frames in a block
	handle->numElements=next_pow2(numElements); // number of blocks (ringbuffers) one ringbuffer per block
	handle->depth=next_pow2(depth); // depth of each ringbuffer and add some fat

	handle->outBuf=malloc(elementSize); //output temp buffer

	// allocate some memory to hold the ringbuffer structs - one per block
	handle->RingBuffer=malloc(sizeof(PaUtilRingBuffer) * handle->numElements);

	// clear it
	memset(handle->RingBuffer,0,sizeof(PaUtilRingBuffer) * handle->numElements);

	// allocate some memory to hold the data pointers - one per block
	handle->RingBufferData=malloc(sizeof(void*) * handle->numElements);

	// clear it
	memset(handle->RingBufferData,0,sizeof(void*) * handle->numElements);

	// allocate all the ringbuffers for the blocks
	for (int i=0;i<handle->numElements;i++) {

		// allocate a data block for one ringbuffer
		handle->RingBufferData[i]=PaUtil_AllocateMemory( handle->elementSize * handle->depth ); // depth * block

		// clear it
		memset(handle->RingBufferData[i],0,handle->elementSize * handle->depth);

		// initialize the ringbuffer
		PaUtil_InitializeRingBuffer(&handle->RingBuffer[i], handle->elementSize, handle->depth, handle->RingBufferData[i]);

	}

	// return the handle
	return handle;

}

/*
 * write into the ring buffers.
 *
 *
 * @handle - accumulator structure to write 
 * @data - start of the incoming data block
 * @block - accumulator block to start writing
 * @numBlocks - number of blocks to write
 */
void gspacc_write (gspacc_t * handle, void * data, size_t accBlock, size_t numBlocks) {

	// assert some parameter values for debugging
	assert(handle!=NULL);
	assert(data!=NULL);
	assert(numBlocks>0);

	// loop through each block that has been sent
	int i=0; while (i<numBlocks) {
	//	printf("gspacc_write: handle=%p,data=%p,accBlock=%zu,numBlocks=%zu\n",handle,data,accBlock,numBlocks);

		// write one data block into into a ringbuffer (elementSize already defined)
		ring_buffer_size_t rbst=PaUtil_WriteRingBuffer(&handle->RingBuffer[accBlock],data,1);

		// advance the block pointer and the data pointer
		i++; data+=handle->elementSize;

		// advance the output buffer pointer, wrapping around if necessary
		accBlock = (accBlock + 1) % handle->numElements;
	}
}

/*
 * read through all the levels of each block and process
 *
 * @handle - accumulator handle
 * @accBlock - output buffer block that we are going to fill
 * @accCb - callback for block processing
 * @dest - destination buffer
 */
int gspacc_read (gspacc_t * handle, size_t accBlock, gspacc_cb accCb, void * dest) {

	//printf("gspacc_read: handle=%p,accBlock=%zu,dest=%p\n",handle,accBlock,dest);

	// assert some parameter values for debugging
	assert(handle!=NULL);
	assert(accBlock<handle->numElements);
	assert(accCb!=NULL);
	
	ring_buffer_size_t rbst=0;

	// loop thorugh all the levels in the block
	// reading out all the available blocks and processing them
	int i=0; while (1) {
	
		// clear the temporary output buffer
		memset(handle->outBuf,0,handle->elementSize);

		// read a block from the ringbuffer into the temporary output buffer
		rbst=PaUtil_ReadRingBuffer(&handle->RingBuffer[accBlock],handle->outBuf,1);

		// if we got something - process it through the callback
		if (rbst>0) {
			
			// process the output block thorugh the callback
			(* accCb)(handle->outBuf,handle->elementSize,dest);

			// increment the output counter
			i++;
		} else {
			
			// jump out of the loop - we have run out of data
			break;
		}
	}

	// return the output counter for debugging
	return i;
}

/*
 * cleanup all allocated memory structures
 */ 
void gspacc_cleanup (gspacc_t * handle) {
	
	free(handle->outBuf);

	for (int i=0;i<handle->numElements;i++) {
		free(handle->RingBufferData[i]);
	}

	free(handle->RingBufferData);
	free(handle->RingBuffer);
	free(handle);
}
