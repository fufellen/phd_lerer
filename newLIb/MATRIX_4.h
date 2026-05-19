
/*
 * Copyright (c) 1992 Konstantin L.Topanov.
 */

#ifndef _MATRIX_

#define _MATRIX_ 1

#ifndef _MEMCHK_
#include "memchk.h"
#endif

void *amem_1d(          /* Allocate one-dimensional array */
int dim,                /* Upper limit, exclusive */
int el_size );          /* Elements size */

void fmem_1d(           /* Free one-dimensional array */
void *ptr );            /* Array pointer, allocated by amem_1d */

void **amem_2d(         /* Allocate two-dimensional array */
int rows,               /* First upper limit, exclusive */
int cols,               /* Second upper limit, exclusive */
int el_size );          /* Elements size */

void fmem_2d(           /* Free two-dimensional array */
register void **ptr,    /* Array pointer, allocated by amem_2d */
int rows );             /* First upper limit, exclusive */

void ***amem_3d(        /* Allocate three-dimensional array */
int depth,              /* First upper limit, exclusive */
int rows,               /* Second upper limit, exclusive */
int cols,               /* Third upper limit, exclusive */
int el_size );          /* Elements size */

void fmem_3d(           /* Free three-dimensional array */
register void ***ptr,   /* Array pointer, allocated by amem_3d */
int depth,              /* First upper limit, exclusive */
int rows );             /* Second upper limit, exclusive */

void *amat_1d(          /* Allocate one-dimensional array */
int lo,                 /* Lower limit, inclusive */
int hi,                 /* Upper limit, inclusive */
int el_size );          /* Elements size */

void fmat_1d(           /* Free one-dimensional array */
void *ptr,              /* Array pointer, allocated by amat_1d */
int lo,                 /* Lower limit, inclusive */
int el_size );          /* Elements size */

void **amat_2d(         /* Allocate two-dimensional array */
int lo_row,             /* First lower limit, inclusive */
int hi_row,             /* First upper limit, inclusive */
int lo_col,             /* Second lower limit, inclusive */
int hi_col,             /* Second upper limit, inclusive */
int el_size );          /* Elements size */

void fmat_2d(           /* Free two-dimensional array */
register void **ptr,    /* Array pointer, allocated by amat_2d */
int lo_row,             /* First lower limit, inclusive */
int hi_row,             /* First upper limit, inclusive */
int lo_col,             /* Second lower limit, inclusive */
int el_size );          /* Elements size */

void ***amat_3d(        /* Allocate three-dimensional array */
int lo_depth,           /* First lower limit, inclusive */
int hi_depth,           /* First upper limit, inclusive */
int lo_row,             /* Second lower limit, inclusive */
int hi_row,             /* Second upper limit, inclusive */
int lo_col,             /* Third lower limit, inclusive */
int hi_col,             /* Third upper limit, inclusive */
int el_size );          /* Elements size */

void fmat_3d(           /* Free three-dimensional array */
register void ***ptr,   /* Array pointer, allocated by amat_3d */
int lo_depth,           /* First lower limit, inclusive */
int hi_depth,           /* First upper limit, inclusive */
int lo_row,             /* Second lower limit, inclusive */
int hi_row,             /* Second upper limit, inclusive */
int lo_col,             /* Third lower limit, inclusive */
int el_size );          /* Elements size */

void ****amat_4d( /* Allocate three-dimensional array */
int lo_row,     /* First lower limit, inclusive */
int hi_row,     /* First upper limit, inclusive */
int lo_col,     /* Second lower limit, inclusive */
int hi_col,     /* Second upper limit, inclusive */
int lo_col1,     /* Third lower limit, inclusive */
int hi_col1,     /* Third upper limit, inclusive */
int lo_col2,     /* Fourth lower limit, inclusive */
int hi_col2,     /* Fourth upper limit, inclusive */
int el_size );   /* Elements size */

void fmat_4d(           /* Free four-dimensional array */
register void ****ptr,    /* Array pointer, allocated by amat_34d */
int lo_row,             /* First lower limit, inclusive */
int hi_row,             /* First upper limit, inclusive */
int lo_col,             /* Second lower limit, inclusive */
int hi_col,				/* Second upper limit, inclusive */
int lo_col1,		    /* Third lower limit, inclusive */
int hi_col1,			/* Third upper limit, inclusive */
int lo_col2,			/* Fourth lower limit, inclusive */
int el_size );           /* Elements size */

#endif
