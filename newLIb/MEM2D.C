
#include "matrix.h"

void **amem_2d( /* Allocate two-dimensional array */
int rows,       /* First upper limit, exclusive */
int cols,       /* Second upper limit, exclusive */
int el_size )   /* Elements size */
	{
	register void **ptr;
	register int i;

	if( ( ptr = (void**)MemCalloc( rows, sizeof( void* ) ) ) == (void**)NULL )
		{
        MemSetError( MERR_2D, "amem_2d" );
		return( ptr );
		}

	for( i = 0; i < rows; ++i )
		if( ( ptr[i] = amem_1d( cols, el_size ) ) == (void*)NULL )
			{
            MemSetError( MERR_2D, "amem_2d" );
			fmem_2d( ptr, i );
            return( (void**)NULL );
			}

	return( ptr );
	}

void fmem_2d(           /* Free two-dimensional array */
register void **ptr,    /* Array pointer, allocated by amem_2d */
int rows )              /* First upper limit, exclusive */
	{
	register int i;

	if( ptr != (void*)NULL )
		{
		for( i = 0; i < rows; ++i )
			fmem_1d( ptr[i] );
		MemFree( ptr );
		}
	}
