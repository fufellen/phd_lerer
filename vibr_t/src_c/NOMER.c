
#include "VIBRr.h"


void NOMER( void) //new
{ // maping of 4 indices[(Ex,Ey,Ez),l,m,n] of base function in eq.(2) of presentation into 1 index that becomes index of unknown in 
	int i,mz,mf,mr,m[2],c,P,
		n_hol;
	
	for ( P=0; P <= 1; P++ )
	for ( m[P] = 0,i=1; i <= 3; i++ )
	{ // LOOP over (Ex,Ey,Ez) - coordinates
		c = (i == 1) ? 1 + P : 1 + ( 1 - P);
		
			for ( mf=c; mf <= M_fi; mf++ ) // LOOP over angular-harmonics (mf-1)
			for ( mz=1; mz <= M_z; mz++ )
			{ 
				n_hol = N_r[mz] - 1;
				for ( mr=1; mr <= MR[n_hol]; mr++ )
				{ // LOOP over radial-harmonics
				
					if (P)
						No2[i][mz][mr][mf] = m[P] + MS[0];
					else
						No[i][mz][mr][mf] = m[P];

					if (No2[i][mz][mr][mf] >= MS_N || No[i][mz][mr][mf] >= MS[0])
{
		printf ("\n		           	  %4d         	  %4d  	  %4d  	  %4d  	  %4d  	  %4d", No[i][mr][mf][mz],No2[i][mr][mf][mz],i,mr,mf,mz );
			fprintf( out,"\n		           	  %4d         	  %4d  	  %4d  	  %4d  	  %4d  	  %4d", No[i][mr][mf][mz],No2[i][mr][mf][mz],i,mr,mf,mz );
}

				m[P]++;
			}//mr
		}//mz
	}//i
}//NOMER