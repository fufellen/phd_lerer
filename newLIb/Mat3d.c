
#include "matrix.h"

void ***amat_3d(    /* Allocate three-dimensional array */
int lo_depth,       /* First lower limit, inclusive */
int hi_depth,       /* First upper limit, inclusive */
int lo_row,         /* Second lower limit, inclusive */
int hi_row,         /* Second upper limit, inclusive */
int lo_col,         /* Third lower limit, inclusive */
int hi_col,         /* Third upper limit, inclusive */
int el_size )       /* Elements size */
	{
	register void ***ptr;
	register int i;

	if( ( ptr = (void***)amat_1d( lo_depth, hi_depth, sizeof(void*) ) ) ==
															(void***)NULL )
		{
        MemSetError( MERR_3D, "amat_3d" );
		return( ptr );
		}

	for( i = lo_depth; i <= hi_depth; ++i )
		if( ( ptr[i] = amat_2d( lo_row, hi_row, lo_col, hi_col, el_size ) ) ==
																(void*)NULL )
			{
			fmat_3d( ptr, lo_depth, i - 1, lo_row, hi_row, lo_col, el_size );
            MemSetError( MERR_3D, "amat_3d" );
            return( (void*)NULL );
			}
	
	return( ptr );
	}


void fmat_3d(           /* Free three-dimensional array */
register void ***ptr,   /* Array pointer, allocated by amat_3d */
int lo_depth,           /* First lower limit, inclusive */
int hi_depth,           /* First upper limit, inclusive */
int lo_row,             /* Second lower limit, inclusive */
int hi_row,             /* Second upper limit, inclusive */
int lo_col,             /* Third lower limit, inclusive */
int el_size )           /* Elements size */
	{
	register int i;

	if( ptr != (void*)NULL )
		{
		for( i = lo_depth; i <= hi_depth; ++i )
			fmat_2d( ptr[i], lo_row, hi_row, lo_col, el_size );
		fmat_1d( ptr, lo_depth, sizeof(void*) );
		}
	}

