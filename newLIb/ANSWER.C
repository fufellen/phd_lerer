
#include <conio.h>
#include "simple.h"

#ifdef NULL
#undef NULL
#endif
#define NULL 0

int answer( char *code )
	{
	int c;

	while( kbhit() )
		getch();

	for(;;)
		{
		int j;
		int zone;

		c = getch();
		if( code[0] == NULL )
			return( NULL );
		if( code[1] == NULL )
			return( NULL );
		for( j = 0, zone=0; code[j] != NULL; ++j )
			{
			if( code[j] == code[0] )
				++zone;
			else
				if( (int)code[j] == c )
					return( zone );
			}
		}
	}

int yesno()
	{
	switch( answer( " Yy Nn" ) )
		{
		case 1:
			cputs( "Yes." );
			return( 1 );
		case 2:
			cputs( "No." );
			return( 0 );
        }
    return( 0 );
	}
