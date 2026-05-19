
#include "matrix.h"

void *amem_1d(  /* Allocate one-dimensional array */
int dim,        /* Upper limit, exclusive */
int el_size )   /* Elements size */
	{
	register void *ptr;
	if( ( ptr = MemCalloc( dim, el_size ) ) == (void*)NULL )
        MemSetError( MERR_1D, "amem_1d" );
	return( ptr );
	}

void fmem_1d(   /* Free one-dimensional array */
void *ptr )     /* Array pointer, allocated by amem_1d */
	{
	if( ptr != (void*)NULL )
		MemFree( ptr );
	}
