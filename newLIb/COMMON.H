
/*
 * Copyright (c) 1992 Konstantin L.Topanov.
 */

#ifndef _COMMON_

#define _COMMON_ 1

#if defined( _QC ) || defined( _MSC_VER )
#if defined(_M_I86SM) || defined(_M_I86MM)
#define _NEARDATA
#endif
#elif defined( __TURBOC__ )
#endif

#if !defined( _MSC_VER ) || ( _MSC_VER < 700 )
#define __inline
#endif

#if defined( __TURBOC__ )
#include "msc2bc.h"
#elif defined( __TSC__ )
#include "msc2tsc.h"
#endif

#ifdef __cplusplus
#define bool int
#define FALSE 0
#define TRUE  1
typedef unsigned char uchar;
#else
typedef enum { FALSE = 0, TRUE = 1 } bool;
typedef unsigned char uchar;
#endif  /* __cplusplus */

typedef unsigned char  byte;
typedef unsigned short ushort;
typedef unsigned int   uint;
typedef unsigned int   cardinal;
typedef unsigned long  ulong;

typedef signed char  schar;
typedef signed short sshort;
typedef signed int   sint;
typedef signed long  slong;

#define Forever for(;;)

typedef enum { LT = -2, LE = -1, EQ = 0, GE = 1, GT = 2 } relation;

#ifndef NULL
#define NULL 0
#endif

#define NULLP 0

extern void (*_ExitErrProc)( void );
extern void ExitError( int ErrCode );
extern void ErrorSig( int n, int line, char *filename );

#endif
