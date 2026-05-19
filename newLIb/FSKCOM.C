
#include <stdlib.h>
#include <stdio.h>


#define CHAR_COMMENT '/'

void fSkipComment( FILE *stream )
	{
	int c;

	do
		{
		do
			c = fgetc( stream );
		while( c == ' ' || c == '\t' || c == '\n' || c == '\r' );

		if( c == CHAR_COMMENT )
			{
			while( ( c = fgetc( stream ) ) != '\n' && c != CHAR_COMMENT )
			if( c == EOF )
				return;
			}
		else
			ungetc( c, stream );
		}
	while( c == ' ' || c == '\t' || c == '\n' || c == '\r' ||
		   c == CHAR_COMMENT );
	}
