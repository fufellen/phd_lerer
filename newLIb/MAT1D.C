
#include "matrix.h"

void *amat_1d(  /* Allocate one-dimensional array */
int lo,         /* Lower limit, inclusive */
int hi,         /* Upper limit, inclusive */
int el_size )   /* Elements size */
	{
	register void *ptr;
	if( ( ptr = MemCalloc( hi-lo+1, el_size ) ) == (void*)NULL )
		{
        MemSetError( MERR_1D, "amat_1d" );
		return( ptr );
		}
	return( (void*)((char*)ptr-lo*el_size) );
	}

void fmat_1d(   /* Free one-dimensional array */
void *ptr,      /* Array pointer, allocated by amat_1d */
int lo,         /* Lower limit, inclusive */
int el_size )   /* Elements size */
	{
	if( ptr != (void*)NULL )
		MemFree( (void*)((char*)ptr+lo*el_size) );
	}
