
#include "matrix.h"

void **amat_2d( /* Allocate two-dimensional array */
int lo_row,     /* First lower limit, inclusive */
int hi_row,     /* First upper limit, inclusive */
int lo_col,     /* Second lower limit, inclusive */
int hi_col,     /* Second upper limit, inclusive */
int el_size )   /* Elements size */
	{
	register void **ptr;
	register int i;

	if( ( ptr = (void**)amat_1d( lo_row, hi_row, sizeof( void* ) ) ) ==
															(void**)NULL )
		{
        MemSetError( MERR_2D, "amat_2d" );
		return( ptr );
		}

	for( i = lo_row; i <= hi_row; ++i )
		if( ( ptr[i] = amat_1d( lo_col, hi_col, el_size ) ) == (void*)NULL )
			{
			fmat_2d( ptr, lo_row, i - 1, lo_col, el_size );
            MemSetError( MERR_2D, "amat_2d" );
            return( (void**)NULL );
			}

	return( ptr );
	}

void fmat_2d(           /* Free two-dimensional array */
register void **ptr,    /* Array pointer, allocated by amat_2d */
int lo_row,             /* First lower limit, inclusive */
int hi_row,             /* First upper limit, inclusive */
int lo_col,             /* Second lower limit, inclusive */
int el_size )           /* Elements size */
	{
	register int i;

	if( ptr != (void*)NULL )
		{
		for( i = lo_row; i <= hi_row; ++i )
			fmat_1d( ptr[i], lo_col, el_size );
		fmat_1d( ptr, lo_row, sizeof(void*) );
		}
	}
