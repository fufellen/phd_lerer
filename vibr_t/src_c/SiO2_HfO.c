void SiO2_(double lam, double *Re_eps, double *Im_eps)
{
	double Lam2;

	Lam2 = lam / 1000.0;
	Lam2 *= Lam2;
	/*
	Dispersionnaya formula dlya SiO2
	n2−1 = 0.6961663λ^2 / (λ^2−0.0684043^2) + 0.4079426λ^2 / (λ^2−0.1162414^2) + 
	0.8974794λ^2 / (λ^2−9.896161^2)
	*/
	
	 *Re_eps = 1.0 + 0.6961663 * Lam2 / (Lam2 - 0.0684043 * 0.0684043) + 
							 0.4079426 * Lam2 / (Lam2 - 0.1162414 * 0.1162414) +
							 0.8974794 * Lam2 / (Lam2 - 9.896161 * 9.896161);
	 *Im_eps = 0.0;
 }//SiO2

void HfO_(double lam, double *Re_eps, double *Im_eps)
{
	double Lam2;

	Lam2 = lam / 1000.0;
	Lam2 *= Lam2;
	/*
		Dispersionnaya formula dlya HfO2
	n2−1 = 1.9558λ^2 / (λ^2−0.154942) + 1.345λ^2 / (λ^2−0.0634^2) + 10.41λ^2 / (λ^2−27.12^2)

	*/
	
	*Re_eps = 1.0 + 1.9558 * Lam2 / (Lam2 - 0.15494 * 0.15494) +
							1.345 * Lam2 / (Lam2 - 0.0634 * 0.0634) +
							10.41 * Lam2 / (Lam2 - 27.12 * 27.12);
	*Im_eps = 0.0;
}//HfO