#include <stdio.h>
#include <time.h>
#include <math.h>

#ifndef CLK_TCK
#define CLK_TCK CLOCKS_PER_SEC
#endif


extern FILE *out;

void tim ( clock_t start )
{
	double sec;
	int hour,min;
		sec =  ( clock() - start ) / (double)CLK_TCK;
		hour = (int)( sec / 3600 );
		min = (int)( ( sec - hour * 3600 ) / 60 );
		sec =  sec - hour * 3600 - min * 60 ;
		printf( "\n time is %2d:%2d:%.2lf \n",hour,min,sec);
		fprintf( out,"\n time is %2d:%2d:%.2lf \n",hour,min,sec);

} /* end TIM */
