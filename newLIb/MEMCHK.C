
/*
 * Copyright (c) 1992 Konstantin L.Topanov.
 */

/* Skip file if not request */
#ifndef NOMEMCHK

#define _INTERNAL_MEMCHK
#include "memchk.h"

#if defined(_M_I86SM) || defined(_M_I86MM)
#define Far2Void(ptr) ((void*)(int)(long)(ptr))
#else
#define Far2Void(ptr) (ptr)
#endif

#if defined( _WINDOWS ) || defined ( _Windows )
#define __WINDOWS
#endif

/************************* Internal types *********************/

typedef struct _PtrList_t       /* Pointer list */
	{
	struct _PtrList_t *prev;    /* Pointer to prev */
	struct _PtrList_t *next;    /* Pointer to next */
	} PtrList_t;

typedef PtrList_t *pPtrList_t;

typedef struct _HugePtrList_t       /* Pointer list */
	{
	struct _HugePtrList_t *prev;    /* Pointer to prev */
	struct _HugePtrList_t *next;    /* Pointer to next */
	} HugePtrList_t;

typedef HugePtrList_t *pHugePtrList_t;

/************************* Globals ***************************/

bool MemOverflowFlag  = FALSE;  /* No memory but MemReserved not freed yet */
bool MemCheckHeapFlag = FALSE;  /* Check heap flag */
bool MemNoExitFlag    = FALSE;  /* Do not call MemFatalError when no memory */

void (*MemFatalFunc)( char *str, va_list ap )  = NULL;
void (*MemNoMemFunc)( MemError_t Code )        = NULL;
void (*MemOverflowFunc)( void )                = NULL;

/************************* Internal data *********************/

static MemError_t mem_error = FALSE;       /* No memory at all */

static pPtrList_t     PtrListHead       = NULL;
static pHugePtrList_t HugePtrListHead   = NULL;
static void * MemReserved = NULL;

/************************* Local functions *****************/

static void MemFatalError( char *str, ... )
	{
	va_list ap;

	va_start( ap, str );

	if( MemFatalFunc != NULL )
		MemFatalFunc( str, ap );
	else
        {
#ifndef __WINDOWS
		vfprintf( stderr, str, ap );
        HeapStat();
#else
        char Buf[2048];
        vsprintf( Buf, str, ap );
        ::MessageBox( GetActiveWindow(), Buf, NULL, MB_ICONSTOP | MB_OK );
#endif
		}

	va_end( ap );

    exit( 1 );
	}

static void MemHeapCheck( void )
	{
	if( MemCheckHeapFlag )
		{
		if( !MemCheckPtrList() )
			MemFatalError( "\a\nMemHeapCheck. List is not OK." );
		if( !IsHeapOK() )
			MemFatalError( "\a\nMemHeapCheck. Heap is not OK." );
		}
	}

static void * ListInsert( void *ptr, pPtrList_t *PtrListHead )
	{
    if( ptr != NULL )
		{
        ((pPtrList_t)ptr)->prev = NULL;
		((pPtrList_t)ptr)->next = *PtrListHead;
        if( *PtrListHead != NULL )
            (*PtrListHead)->prev = (pPtrList_t)ptr;
        *PtrListHead = (pPtrList_t)ptr;
		ptr = &((pPtrList_t)ptr)[1];
		}
	return( Far2Void( ptr ) );
	}

static void * ListDelete( void *ptr, pPtrList_t *PtrListHead )
	{
    if( ptr != NULL )
		{
		pPtrList_t Prev;
		pPtrList_t Next;

		ptr = &((pPtrList_t)ptr)[-1];
		Prev = ((pPtrList_t)ptr)->prev;
		Next = ((pPtrList_t)ptr)->next;

        if( ( ( Next == NULL ) ||
              ( Next != NULL && Next->prev == ptr ) ) &&
            ( ( Prev == NULL && *PtrListHead == ptr ) ||
              ( Prev != NULL && Prev->next == ptr ) ) )
			{
            if( Prev == NULL )
				{
				*PtrListHead = Next;
                if( Next != NULL )
                    Next->prev = NULL;
				}
			else
				{
				Prev->next = Next;
                if( Next != NULL )
					Next->prev = Prev;
				}
			return( Far2Void( ptr ) );
			}
		}
    return( NULL );
	}

static void * _memMalloc( size_t Size )
	{
	void *Ptr;
	Size += sizeof( PtrList_t );
	MemHeapCheck();
    if( ( Ptr = malloc( Size ) ) == NULL )
		{
		MemFreeReserve();
		Ptr = malloc( Size );
		}
	return( Ptr );
	}

static void * _memCalloc( size_t ElNum, size_t ElSize )
	{
	void *Ptr;
	size_t Size = sizeof( PtrList_t ) + ElSize * ElNum;
	MemHeapCheck();
    if( ( Ptr = calloc( Size, sizeof( char ) ) ) == NULL )
		{
		MemFreeReserve();
		Ptr = calloc( Size, sizeof( char ) );
		}
	return( Ptr );
	}

static void * _memRealloc( void *memblock, size_t Size )
	{
	void *Ptr;
	Size += sizeof( PtrList_t );
	MemHeapCheck();
    if( ( Ptr = realloc( memblock, Size ) ) == NULL )
		{
		MemFreeReserve();
		Ptr = realloc( memblock, Size );
		}
	return( Ptr );
	}

static void * _memHalloc( long ElNum, size_t ElSize )
	{
	const uint Align = sizeof( int );
	long Size = sizeof( HugePtrList_t ) + ElSize * ElNum;
	MemHeapCheck();
#if !defined( __TSC__ )
    return( calloc( Size / Align + ( Size % Align != 0 ), Align ) );
#else
    return( calloc( Size / Align + ( Size % Align != 0 ) * Align ) );
#endif
	}

/************************* Exported functions ***************/

bool MemCheckPtrList( void )
	{
	bool Stat = TRUE;
	pPtrList_t Ptr = PtrListHead;

    if( Ptr != NULL && Ptr->prev != NULL )
		Stat = FALSE;

    while( Stat && Ptr != NULL )
		{
		pPtrList_t Next = Ptr->next;
        if( Next != NULL && Ptr != Next->prev )
			Stat = FALSE;
		else
			Ptr = Next;
		}
	return( Stat );
	}

void MemSetError( MemError_t Code, const char *Str )
	{
	mem_error |= Code;
	if( MemNoMemFunc != NULL )  /* If exist call user defined function */
		MemNoMemFunc( Code );
	else if( !MemNoExitFlag )   /* Check allow exit */
		MemFatalError( "\a\nMemSetError. No memory in '%s'. "
					   "Code %04x. MemError %04x. Exiting...",
					   Str, Code, mem_error );
	}

bool MemReserve( size_t Size )
	{
	MemHeapCheck();
    return( ( MemReserved = _memCalloc( Size, sizeof( char ) ) ) != NULL );
	}

void MemFreeReserve( void )
	{
	MemOverflowFlag = TRUE;
    if( MemReserved != NULL )
		{
		free( MemReserved );
        MemReserved = NULL;
		if( MemOverflowFunc != NULL )
			MemOverflowFunc();
		}
	}

MemError_t MemError( void )
	{
	MemError_t Code = mem_error;
	mem_error = 0;
	MemHeapCheck();
	return( Code );
	}

void *MemMalloc( size_t Size )
	{
	void *ptr;
	if( ( ptr = (void*)ListInsert( _memMalloc( Size ),
                                   &PtrListHead ) ) == NULL )
		MemSetError( MallocError, "MemMalloc" );
	return( ptr );
	}

void *MemCalloc( size_t ElNum, size_t ElSize )
	{
	void *ptr;
	if( ( ptr = (void*)ListInsert( _memCalloc( ElNum, ElSize ),
                                   &PtrListHead ) ) == NULL )
		MemSetError( CallocError, "MemCalloc" );
	return( ptr );
	}

void *MemRealloc( void *memblock, size_t Size )
	{
    void *ptr = NULL;

    if( memblock == NULL )
		{
		if( ( ptr = (void*)ListInsert( _memMalloc( Size ),
                                       &PtrListHead ) ) == NULL )
			MemSetError( ReallocError, "MemRealloc" );
		}
	else if( ( memblock = (void*)ListDelete( memblock,
                                             &PtrListHead ) ) != NULL )
		{
		if( ( ptr = (void*)ListInsert( _memRealloc( memblock, Size ),
                                       &PtrListHead ) ) == NULL )
			MemSetError( ReallocError, "MemRealloc" );
		}
	else
		MemFatalError( "\a\nMemRealloc. Reallocate not allocated pointer." );

	return( ptr );
	}

char *MemStrDup( const char *Source )
	{
    char *ptr = (char*)MemMalloc( strlen( Source ) + 1 );
    if( ptr != NULL )
		strcpy( ptr, Source );
	else
		MemSetError( StrDupError, "MemStrDup" );
	return( ptr );
	}

void MemFree( void *memblock )
    {
    void *ptr = memblock;
	MemHeapCheck();
    if( memblock != NULL )
		{
		if( ( memblock = (void*)ListDelete( memblock,
                                            &PtrListHead ) ) != NULL )
			free( memblock );
		else
            MemFatalError( "\a\nMemFree. Free not allocated pointer %p.",
                           ptr );
		}
	}

void MemAllFree( void *memblock )
	{
	pPtrList_t slot;
	pPtrList_t next;
	void  *Origin = &((pPtrList_t)memblock)[-1];

	MemHeapCheck();

    for( slot = PtrListHead; slot != NULL; slot = next )
		{
		next = slot->next;
		free( Far2Void( slot ) );
		if( slot == Origin )
			{
			PtrListHead = next;
            next->prev = NULL;
			return;
			}
		}

    if( memblock != NULL )
		MemFatalError( "\a\nMemAllFree. All memory has set free. Free not allocated pointer." );
	else
        PtrListHead = NULL;
	}

void *MemHalloc( long ElNum, size_t ElSize )
	{
	void *ptr;
	if( ( ptr = ListInsert( _memHalloc( ElNum, ElSize ),
                            (pPtrList_t*)&HugePtrListHead ) ) == NULL )
		MemSetError( HallocError, "MemHalloc" );
	return( ptr );
	}

void MemHfree( void *memblock )
    {
    void *ptr = memblock;
	MemHeapCheck();
    if( memblock != NULL )
		{
		if( ( memblock = ListDelete( memblock,
                         (pPtrList_t*)&HugePtrListHead ) ) != NULL )
			free( memblock );
		else
            MemFatalError( "\a\nMemHfree. Free not allocated pointer %Fp.",
                           ptr );
		}
	}

uint MemCount( void )
	{
	uint count = 0;                      /* Used slots counter. */
	pPtrList_t slot;

	MemHeapCheck();

    for( slot = PtrListHead; slot != NULL; slot = slot->next )
		++count;

	return( count );
	}

uint MemDump( void )
	{
	uint count = 0;                      /* Used slots counter. */
	pPtrList_t slot;

	MemHeapCheck();

	fprintf( stderr, "\n\nHead %p", PtrListHead );
    for( slot = PtrListHead; slot != NULL; slot = slot->next )
		{
		fprintf( stderr, "\nSlot=%p, Prev=%p, Next=%p, Ptr=%p",
						 slot, slot->prev, slot->next, &slot[1] );
		++count;
		}
	fprintf( stderr, "\nCount = %u\n", count );

	return( count );
	}

#ifndef __WINDOWS

/*
 * Reports on the status returned by _heapchk
 */
void HeapStat( void )
	{
    fprintf( stderr, "\nHeap status: " );
	switch( _heapchk() )
		{
		case _HEAPOK:
			fprintf( stderr, "OK - heap is fine" );
			break;
		case _HEAPEMPTY:
			fprintf( stderr, "OK - empty heap" );
			break;
		case _HEAPEND:
			fprintf( stderr, "OK - end of heap" );
			break;
		case _HEAPBADPTR:
			fprintf( stderr, "ERROR - bad pointer to heap" );
			break;
		case _HEAPBADBEGIN:
			fprintf( stderr, "ERROR - bad start of heap" );
			break;
		case _HEAPBADNODE:
			fprintf( stderr, "ERROR - bad node in heap" );
			break;
		}
	fprintf( stderr, "\n\n" );
	}

bool IsHeapOK( void )
	{
	return( _heapchk() == _HEAPOK );
	}

#if defined( _QC ) || defined( _MSC_VER )
/*
 * Writes each block in the heap.
 */
void HeapDump( void )
    {
	static char *Mess[] = { "Free", "Used" };
	struct _heapinfo hi;
	unsigned heapstatus;

	/* Walk through entries, checking free blocks. */
	hi._pentry = NULL;
	while( (heapstatus = _heapwalk( &hi )) == _HEAPOK )
		{
        fprintf( stderr, "\nPointer=%p, size=%u, stat=%s",
                 hi._pentry, hi._size, Mess[hi._useflag] );
        }
	}
#endif /* _QC _MSC_VER */

#else
bool IsHeapOK( void )
    {
    return TRUE;
    }
#endif /* __WINDOWS */

#ifdef __cplusplus

void *operator new( size_t size )
    {
    return( MemMalloc( size ) );
    }

void operator delete( void *ptr )
    {
    MemFree( ptr );
    }

#endif /* __cplusplus */

#endif /* NOMEMCHK */
