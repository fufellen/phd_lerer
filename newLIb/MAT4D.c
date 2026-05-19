#include "matrix_4.h"

void ****amat_4d( /* Allocate three-dimensional array */
int lo_row,     /* First lower limit, inclusive */
int hi_row,     /* First upper limit, inclusive */
int lo_col,     /* Second lower limit, inclusive */
int hi_col,     /* Second upper limit, inclusive */
int lo_col1,     /* Third lower limit, inclusive */
int hi_col1,     /* Third upper limit, inclusive */
int lo_col2,     /* Fourth lower limit, inclusive */
int hi_col2,     /* Fourth upper limit, inclusive */
int el_size )   /* Elements size */
	{
	register void ****ptr;
	register int i;

	if( ( ptr = (void****)amat_1d( lo_row, hi_row, sizeof( void*** ) ) ) ==
															(void****)NULL )
		{
      //  MemSetError( MERR_2D, "amat_4d" );
		return( ptr );
		}

	for( i = lo_row; i <= hi_row; ++i )
		if( ( ptr[i] = amat_3d( lo_col, hi_col, lo_col1, hi_col1,lo_col2, hi_col2, el_size ) ) == (void***)NULL )
			{
			fmat_4d( ptr, lo_row, i - 1, lo_col, hi_col, lo_col1,hi_col1,lo_col2, el_size );
       //     MemSetError( MERR_2D, "amat_4d" );
            return( (void****)NULL );
			}

	return( ptr );
	}

void fmat_4d(           /* Free four-dimensional array */
register void ****ptr,    /* Array pointer, allocated by amat_34d */
int lo_row,             /* First lower limit, inclusive */
int hi_row,             /* First upper limit, inclusive */
int lo_col,             /* Second lower limit, inclusive */
int hi_col,				/* Second upper limit, inclusive */
int lo_col1,		    /* Third lower limit, inclusive */
int hi_col1,			/* Third upper limit, inclusive */
int lo_col2,			/* Fourth lower limit, inclusive */
int el_size )           /* Elements size */
	{
	register int i;

	if( ptr != (void****)NULL )
		{
		for( i = lo_row; i <= hi_row; ++i )
			fmat_3d( ptr[i],lo_col,hi_col, lo_col1,hi_col1,lo_col2, el_size );
		fmat_1d( ptr, lo_row, sizeof(void****) );
		}
	}
