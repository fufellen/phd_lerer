
 /*
 * Copyright (c) 1992 Konstantin L.Topanov.
 */

#ifndef _MEMCHK_

#define _MEMCHK_

#include <malloc.h>     /* typedef size_t */
#include <stdarg.h>     /* typedef va_list */

#ifndef _COMMON_
#include "common.h"
#endif

/*#ifndef _MSC2BC_
#include <msc2bc.h>
#endif*/

#ifndef NOMEMCHK

#ifdef MEMREDEFINE
#define malloc  MemMalloc
#define calloc  MemCalloc
#define realloc MemRealloc
#define free    MemFree
#define strdup  MemStrDup
#endif

#ifdef _INTERNAL_MEMCHK

#include <stdio.h>      /* printf */
#include <stdlib.h>     /* exit */
#include <string.h>     /* str */

#endif  /* _INTERNAL_MEMCHK */

#define MemError_t int
#define MallocError  0x01
#define CallocError  0x02
#define ReallocError 0x04
#define StrDupError  0x08
#define HallocError  0x10
#define MERR_1D      0x20
#define MERR_2D      0x40
#define MERR_3D      0x80

extern bool MemOverflowFlag;    /* No memory but MemReserved not freed yet */
extern bool MemCheckHeapFlag;   /* Check heap flag */
extern bool MemNoExitFlag;      /* Do not call MemFatalError when no memory */

extern void (*MemFatalFunc)( char *str, va_list ap );
extern void (*MemNoMemFunc)( MemError_t Code );
extern void (*MemOverflowFunc)( void );

bool MemCheckPtrList( void );                       /* Check allocated list */
void MemSetError( MemError_t Code, const char *Str );	/* Sets error code */
MemError_t MemError( void );                        /* Returns error code */

void *MemMalloc( size_t size );                     /* malloc */
void *MemCalloc( size_t num, size_t size );         /* calloc */
void *MemRealloc( void *memblock, size_t size );    /* realloc */
char *MemStrDup( const char *Source );              /* strdup */
void MemFree( void *memblock );               /* free */

void *MemHalloc( long ElNum, size_t ElSize );
void MemHfree( void *memblock );

bool MemReserve( size_t Size );
void MemFreeReserve( void );

/* Free the memblock and all blocks allocated after memblock. */
void MemAllFree( void *memblock );

uint MemCount( void );                              /* Total alloc */
uint MemDump( void );                               /* Print alloc */

bool IsHeapOK( void );
void HeapStat( void );
void HeapDump( void );

#else

#define MemMalloc(size) malloc(size)
#define MemCalloc(num,size) calloc(num,size)
#define MemRealloc(memblock,size) realloc(memblock,size)
#define MemStrDup(str) strdup(str)
#define MemFree(memblock) free(memblock)

#define CheckNullPointer     fprintf(stderr,"CheckNullPointer not suppoted.");
#define MemAllFree(memblock) fprintf(stderr,"MemAllFree not suppoted.");
#define MemError()           fprintf(stderr,"MemError not suppoted.");
#define MemCount()           fprintf(stderr,"MemCount not suppoted.");
#define MemDump()            fprintf(stderr,"MemDebug not suppoted.");

#endif /* NOMEMCHK */

#endif  /* avoid redefenition */
