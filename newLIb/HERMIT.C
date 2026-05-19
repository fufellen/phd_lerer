void Hermit( int N,double x,double *He )
{
	/* Полиномы Эрмита * exp(- x * x) */
	int n;
	double e;
	e = exp(- x * x);
	He[0] = e;
	if (!N)
		return;
	He[1] = 2.0 * x * e;
	for ( n=1; n < N; n++ )
		He[n+1] = 2.0 * ( x * He[n] - He[n-1] * n );
}/* Hermit */

double Arth( double x )
{
	return (0.5 * log( ( 1.0 + x ) / ( 1.0 - x ) ));
}/* Arth */

double th( double x )
{
	double e;
	e = exp( - 2.0 * x );
	return ( ( 1.0 - e ) / ( 1.0 + e ) );
}/* th */

void GEGEN( int N,double x,double lam,double *Ge )
{
	/* Полиномы Гегенбаура */
	int n;
	double en;
	Ge[0] = 1.0;
	if (!N)
		return;
	Ge[1] = 2.0 * x * lam;
	for ( n=0; n < N-1; n++ )
	{
		en = (int)n;
		Ge[n+2] = ( 2.0 * ( lam + en + 1.0 ) * x * Ge[n+1] -
							  ( 2.0 * lam + en ) * Ge[n] ) / (2.0 + en);
	}
}/* GEGEN */
