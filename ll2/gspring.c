#include<stdlib.h>
#include "gspring.h"

gspring_t *
gspring_initialize(size_t element_size, size_t nElements) {
	gspring_t * handle=malloc(sizeof(gspring_t));
	return handle;
}
