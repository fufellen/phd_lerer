
void SortInt( int m, int *x )
{
	int t;
	int i, j, k;
	for ( i=2; i <= m; i++ )
	{
		k = i;
		t = x[k];
		for ( j=i-1; j >= 1; j-- )
			if ( x[j] > t )
			{
				k = j;
				x[j+1] = x[j];
			}  /* end j */
		x[k] = t;
	} /* end i */
}

