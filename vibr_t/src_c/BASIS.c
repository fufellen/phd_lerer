
#include "VIBRr.h"

//double rz0,rz1;
void BASIS( double alf, double bet1) //вычисление Rnp calculate
{ // !!!core function 
	double R,*J_n,*J_n1,x,x1,y,ro,ro2,tet,tet1,U,V,Z,rz,ro2_z,c_tet,s_tet,c,n1,si,co,Rnu,Rnu1; 
	double rz0,rz1;
	int M_r1,nl,n_hol,mz,mf,mr,nm,n,mr1;
	
	J_n = ( double*)amat_1d( -1,M_fi-1, sizeof( double) );
	if (N_La > 1)
		J_n1 = ( double*)amat_1d( -1,M_fi-1, sizeof( double) );
		
	for ( mz=1; mz <= ms; mz++ )
	{
		c = a_[mz] / b_[mz];

		x = a_[mz] * alf;
		y = b_[mz] * bet1;
		ro2 = x * x + y * y;
		ro = sqrt(ro2); // ROpq in presentation p.6
		tet = (ro > 1.0e-10) ? atan2(y,x) : atan2(b_[mz] * sifi,a_[mz] * cofi); // THETApq
				
		n_hol = N_r[mz] - 1;
				
		c_tet = x / ro;
		s_tet = y / ro;

		U = c_tet * c_fi0[n_hol] + c * s_tet * s_fi0[n_hol];
		V = s_tet * c_fi0[n_hol] - c_tet * s_fi0[n_hol] /  c;

		Z = sqrt(U * U + V * V);
		tet= atan2(V,U);
	
		rz = ro * Z;
		ro2_z = rz * rz;
		
		tet1 = tet - FI_r_;  
		
		n = n_hol;
		//for ( n=1;n <= N_hol;n++ )
		for (nl = 1;nl <= N_l[n];nl++ )
		{
			M_r1 = (nl ==1) ? M_r_1 : M_r_;

			Rnu = r_l[n][nl]; 
			Rnu1 = (nl == 1) ? 0.0 :r_l[n][nl-1];
			
			rz0 = rz * Rnu;
			rz1 = rz * Rnu1;
						
			for ( mf=-1; mf < M_fi; mf++ )
			{
				J_n[mf] = jn(mf,rz0); 
				if (nl != 1)
					J_n1[mf] = jn(mf,rz1);
			}
			
			for ( mf=1; mf <= M_fi; mf++ )
			{
				nm = mf - 1;
				n1 = (double)nm;
				
				co = cos(tet1 * n1);
				si = sin(tet1 * n1);
				
				x = rz0 * J_n[nm-1] - J_n[nm] * n1; 
				if (nl != 1)
					x1 = rz1 * J_n1[nm-1] - J_n1[nm] * n1;

				for ( mr1=1; mr1 <= M_r1; mr1++ )
				{
					mr = (nl == 1) ? mr1 : mr1 + (nl - 1) * M_r_  + 1;
					
					y = j[nm][mr][n];

					U = 1.0 / (ro2_z - y * y);
					if (mr <= M_r1)
						R = x * J[nm][mr][n] * U;
					else
						R = (x * J[nm][mr][n] - x1 * Y1[nm][mr][n]) * U;

					R_c2[mz][mf][mr] = R * co;
					R_s2[mz][mf][mr] = R * si;
				}//mr // R_c and R_s stands for Clm [see "Presentation" p.6]
			

			}//mf
		}//n

	}//mz

	fmat_1d( (void*)J_n,  -1, sizeof( double) );
	
	if (N_La > 1)
		fmat_1d( (void*)J_n1,  -1, sizeof( double) );
		
}//Basis