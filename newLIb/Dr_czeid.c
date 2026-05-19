#include "newlib.h"


void cZeydel( int N, int It, double **RA, double **IA, double *RZ, double *IZ, double eps )
{
	int S, i, j, k;
	double Zn, modA, modX, *RX, *IX;
	
	RX = amat_1d( 1,N, sizeof( double) );
	IX = amat_1d( 1,N, sizeof( double) );
		
	S = 0;
	do
	{
		k = 0;
		for( i=1; i <= N; i++ )
		{
			RX[i] = - RA[i][N+1];
			IX[i] = - IA[i][N+1];
			for( j=1; j <= N; j++ )
			{
				RX[i] += RA[i][j] * RZ[j] - IA[i][j] * IZ[j];
				IX[i] += RA[i][j] * IZ[j] + IA[i][j] * RZ[j];
			}

			modX = sqrt( RX[i] * RX[i] + IX[i] * IX[i] );
			Zn = RA[i][i] * RA[i][i] + IA[i][i] * IA[i][i];
			modA = sqrt( Zn );
			if ( ( modX / modA ) >= eps )
				k = 1;
	
			RX[i] = RZ[i] - ( RX[i] * RA[i][i] + IX[i] * IA[i][i] ) / Zn;
			IX[i] = IZ[i] - ( IX[i] * RA[i][i] - RX[i] * IA[i][i] ) / Zn;

			RZ[i] = RX[i];
			IZ[i] = IX[i];
		}
		S++;
	//	printf( "\n  S = %d", S );
	//	fprintf( out,"\n  S = %d", S );
	}
	while ( k == 1 && S <= It);

	for( i=1; i <= N; i++ )
	{
		RA[i][N+1] = RX[i];
		IA[i][N+1] = IX[i];
	}
	
		if ( S >= It )
		{
			printf ("\n  Иттерации не сходятся");
	
			for( i=1; i <= N; i++ )
				RA[i][N+1] = IA[i][N+1] = 0.0;
		
		}

	fmat_1d( (void*)RX, 1, sizeof( double) );
	fmat_1d( (void*)IX, 1, sizeof( double) );
	
	return;
}

