#include "VIBRr.h"


int CLineEqC1(int n, int nb, scomplex **ca, int nca, scomplex **cb, int ncb)
{
  int info;
  char trans='T';
  int *ipiv;

  ipiv = malloc(nca*sizeof(int));

  zgetrf(&n, &n, ca[0], &nca, ipiv, &info);
  if(info) goto Exit;
  zgetrs(&trans, &n, &nb, ca[0], &nca, ipiv, cb[0], &ncb, &info);
Exit:
  free(ipiv);
  if(info) printf("CLineEqC1: The solution is absent (info=%ld)\n",info);
  return info;
}

/*----------------------------------------------------
   аҐиҐ­ЁҐ бЁбвҐ¬л «Ё­Ґ©­ле  «ЈҐЎа ЁзҐбЄЁе га ў­Ґ­Ё©
   ca - ¬ ваЁж  бЁбвҐ¬л, n - а §¬Ґа ¬ ваЁжл,
   cb - ¬ ббЁў бв®«Ўж®ў бў®Ў®¤­ле з«Ґ­®ў,
   nb - Є®«-ў® бв®«Ўж®ў бў®Ў®¤­ле з«Ґ­®ў
   ђҐ§г«мв в Ї®¬Ґй Ґвбп ў ¬ ббЁў бв®«Ўж®ў Їа ўле з бвҐ© (cb)
------------------------------------------------------*/
int CLineEqC(scomplex **ca, int n, scomplex **cb, int nb)
{
  return CLineEqC1(n, nb, ca, n, cb, n);
}

int DLineEqD1(int n, int nb, double **ca, int nca, double **cb, int ncb)
{
  int info;
  char trans='T';
  int *ipiv;

  ipiv = malloc(nca*sizeof(int));

  dgetrf(&n, &n, ca[0], &nca, ipiv, &info);
  if(info) goto Exit;
  dgetrs(&trans, &n, &nb, ca[0], &nca, ipiv, cb[0], &ncb, &info);
Exit:
  free(ipiv);
  if(info) printf("CLineEqC1: The solution is absent (info=%ld)\n",info);
  return info;
}

int DLineEqFactor(int n,		// matrix size
				  double **ca,	// the matrix
				  int *ipiv)	// int array of size n
// Makes LU factorization of the matrix, array ipiv contains information about interchanged rows
{
  int info;
 
  dgetrf(&n, &n, ca[0], &n, ipiv, &info);
  if(info) printf("DLineEqFactor: The solution is absent (info=%ld)\n",info);
  return info;
}

int DLineEqSolve(int n,			// matrix size
				int nb,			// number of right hand sides
				double **ca,	// the matrix
				double **cb,	// right hand side matrix
				int * ipiv)		// int array of size n from DLineEqFactor
//Solves SLAE with matrix factorized previously by DLineEqFactor
{
  int info;
  char trans='T';
  
  dgetrs(&trans, &n, &nb, ca[0], &n, ipiv, cb[0], &n, &info);

  if(info) printf("DLineEqSolve: The solution is absent (info=%ld)\n",info);
  return info;
}


int CLineEqFactor(int n,		// matrix size
				  scomplex **ca,	// the matrix
				  int *ipiv)	// int array of size n
// Makes LU factorization of the matrix, array ipiv contains information about interchanged rows
{
  int info;
 
  zgetrf(&n, &n, ca[0], &n, ipiv, &info);
  if(info) printf("CLineEqFactor: The solution is absent (info=%ld)\n",info);
  return info;
}

int CLineEqSolve(int n,			// matrix size
				int nb,			// number of right hand sides
				scomplex **ca,	// the matrix
				scomplex **cb,	// right hand side matrix
				int * ipiv)		// int array of size n from CLineEqFactor
//Solves SLAE with matrix factorized previously by CLineEqFactor
{
  int info;
  char trans='T';
  
  zgetrs(&trans, &n, &nb, ca[0], &n, ipiv, cb[0], &n, &info);

  if(info) printf("CLineEqSolve: The solution is absent (info=%ld)\n",info);
  return info;
}



/*----------------------------------------------------
   аҐиҐ­ЁҐ бЁбвҐ¬л «Ё­Ґ©­ле  «ЈҐЎа ЁзҐбЄЁе га ў­Ґ­Ё©
   ca - ¬ ваЁж  бЁбвҐ¬л, n - а §¬Ґа ¬ ваЁжл,
   cb - ¬ ббЁў бв®«Ўж®ў бў®Ў®¤­ле з«Ґ­®ў,
   nb - Є®«-ў® бв®«Ўж®ў бў®Ў®¤­ле з«Ґ­®ў
   ђҐ§г«мв в Ї®¬Ґй Ґвбп ў ¬ ббЁў бв®«Ўж®ў Їа ўле з бвҐ© (cb)
------------------------------------------------------*/
int DLineEqD(double **ca, int n, double **cb, int nb)
{
  return DLineEqD1(n, nb, ca, n, cb, n);
}


scomplex ** AllocC(int arow, int acol)
{
scomplex **p, *p_blas;
int i;

	if ((p_blas = (scomplex *) calloc(arow*acol, sizeof(scomplex) ))) 
	{
		if( (p = (scomplex **) calloc(arow , sizeof(scomplex*) )) == NULL) 
		{
			free(p_blas); 
			printf ("Cannot allocate memory!(AllocC)");
			exit(1);
			return NULL;
		};

		for ( i = 0; i < arow; i++, p_blas += acol) 
			p[i] = p_blas;
		return p;
	}
	else 
	{
	printf ("Cannot allocate memory!(AllocC)");
			exit(1);
	}
		return NULL;
		
}

void DeAllocC(scomplex **p)
{
scomplex  *p_blas;

		p_blas = p[0];
		free(p_blas);
		free(p);
}

void ** AllocV(int arow, int acol, int size)
{
void **p, *p_blas;
int i;

	if ((p_blas = (void *) calloc( arow*acol, size ))) 
	{
		if( (p = (void **) calloc(arow, sizeof(void*) )) == NULL) 
		{
			free(p_blas); 
			printf ("Cannot allocate memory!(AllocV)");
			exit(1);
			return NULL;
		};

		for ( i = 0; i < arow; i++, (char *)p_blas += acol*size) 
			p[i] = p_blas;
		return p;
	}
	else 
	{
	printf ("Cannot allocate memory!(AllocV)");
			exit(1);
	}
		return NULL;
		
}

void DeAllocV(void **p)
{
void  *p_blas;

		p_blas = p[0];
		free(p_blas);
		free(p);
}

void MatrixSumD(int m,int n, double a, double **x, double **y)
{
// матрицы представляют собой массивы m-мерных векторов размера n 
//y=a*x +y
int i, inc = 1;
//	i = m*n;
	// выполняем одной строкой
	for (i = 0; i < m; i++)
		DAXPY(&n,&a,x[i],&inc,y[i],&inc);
}

void MatrixCopyD(int m,int n, double **x, double **y) //копируем x в y
{
	// матрицы представляют собой массивы m-мерных векторов размера n 
	// y = x;
	int i, inc = 1;
//	i = m*n;
	// выполняем одной строкой
	for (i = 0; i < m; i++)
		dcopy(&n,x[i],&inc,y[i],&inc);
}

void MatrixCopyC(int m,int n, scomplex **x, scomplex **y) //копируем x в y
{
	// матрицы представляют собой массивы m-мерных векторов размера n 
	// y = x;
	int i, inc = 1;
//	i = m*n;
	// выполняем одной строкой
	for (i = 0; i < m; i++)
		zcopy(&n,x[i],&inc,y[i],&inc);
}

void MatrixDotD(int m,int n,  double **x, double **y, double *res)
{
// матрицы представляют собой массивы m-мерных векторов размера n 
//процедура осуществляет скалярное произведение m-мерных векторов для всего массива размера n 
int i, j;
double *r;
	
	r = (double *)calloc(n,sizeof(double));

	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
			r[i] += x[j][i]*y[j][i] ;

	j = 1;
	// res = r;
		dcopy(&n,r,&j,res,&j);
	free(r);
}

void MatrixAbsD(int m,int n, double **x,  double *res)
{
// матрицы представляют собой массивы m-мерных векторов размера n 
//процедура находит длины m-мерных векторов для всего массива размера n 
int i, j, bg;

double **y, * r, pow = 2.;
	
	bg = m * n;
	y = (double **)AllocV(m,n,sizeof(double));
	 r = (double *)calloc(n,sizeof(double));

	//возвели в квадрат каждый элемент матрицы
	vdPowx(bg, x[0], pow, y[0]);
		
	// сумма квадратов
	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
			r[i] += y[j][i];

	//извлекли корень
	vdSqrt(n, r , res);

	free(r);
	DeAllocV(y);

/*	
почему то не работает
	for (i = 0; i < n; i++)
	{
		res[i] = dnrm2(&bg,x[0],&n);
	}
*/

}

void MatrixDivD(int n, double *x,  double *y,  double *res)
{
// деление элементов вектора x на элементы вектора y
	
	vdDiv(n,x,y,res);

}

void MatrixScaleD(int m, int n, double **x,  double *y,  double **res)
{
//x -  матрица представляет собой массив m-мерных векторов размера n 
//y - вектор размера n 
//процедура умножает каждый m-мерный вектор из x на элемент вектора y
	
	
int i, j;

		for (j = 0; j < m; j++)
			for (i = 0; i < n; i++)
			res[j][i] = x[j][i]*y[i];
}

void MatrixNormalizeD(int m, int n, double **x,  double *y,  double **res)
{
//x -  матрица представляет собой массив m-мерных векторов размера n 
//y - вектор размера n содержащий длины этих векторов
//процедура умножает каждый m-мерный вектор из x на элемент вектора y
	
	
int i, j;

		for (j = 0; j < m; j++)
			for (i = 0; i < n; i++)
			res[j][i] = x[j][i]/y[i];
}


void VectorScaleD(int m, int n,  double a, double **x)
{
//x -  матрица представляет собой массив m-мерных векторов размера n 

//процедура умножает каждый m-мерный вектор из x на a
	
	
int i, inc = 1  ;
	//s = m*n;
for (i = 0; i < m; i++)
	dscal(&n,&a,x[i],&inc);
}

void VectorScaleC(int m, int n,  scomplex a, scomplex **x)
{
//x -  матрица представляет собой массив m-мерных векторов размера n 

//процедура умножает каждый m-мерный вектор из x на a
	
	
int i, inc = 1  ;
	//s = m*n;
for (i = 0; i < m; i++)
	zscal(&n,&a,x[i],&inc);
}

double VectorMaxAbs(int m, int n, scomplex **x)
{
//x -  матрица представляет собой массив m-мерных векторов размера n 

//процедура ищет максимум каждго m-мерного вектора из x и далее максимальный из них
	
	
int i, inc = 1, j ;
//scomplex *p;
double z_old = 0., z;
	
//	p = x[0];

 //i = m*n;

/*
for (i = 0; i < m; i++)
{
	index = izamax(&n,x[i],&inc);
		z = x[i][index].real * x[i][index].real + x[i][index].imag * x[i][index].imag;
	if (z > z_old) z_old = z;

}
*/
for (j = 0; j < m; j++)
for (i = 0; i < n; i++)
{
		z = x[i][j].real * x[i][j].real + x[i][j].imag * x[i][j].imag;
	if (z > z_old) z_old = z;

}

	return sqrt(z_old);
}

double snorm(scomplex x)
{
	return (x.real * x.real + x.imag * x.imag);
}

int MatrixPack(scomplex ** A, // исходная матрица
			int n,			// размер матрицы
			double tol,		// точность определения нуля
			scomplex * values,// выходной массив со значениями ненулевых элементов матрицы
			int  * columns,		// выходной массив с номерами столбцов ненулевых элементов
			int * rowindex)	// выходной массив с номерами первых ненулевых элементов в массиве values
{
// Converts general sparse matrix to sparse matrix storage format (see MKL 7.1 manual)
	int i,j,k,p;
	
	
	k = 0;
	
	for (i = 0; i < n; i++)
	{
		p = 0;
//		fprintf(out," \n");
		for (j = 0; j < n; j++)
		{
			if ( (snorm(A[i][j]) > tol ) || (snorm(A[j][i]) > tol ) || (i == j) ) 
			{
				// not  zero
				
				values[k] = A[i][j];
					columns[k] = j + 1; // column index of the element in Fortran convention
				k++;// 

				p++; // non-zero counter in the row

				if (p == 1) rowindex[i] = k;

			}
//				fprintf(out,"  %lf ", snorm(A[i][j]));
		} // j
	}// i
	rowindex[n] = k+1;

	values = realloc(values,sizeof(scomplex) * k);
	 columns = realloc(columns,sizeof(int) * k);

/*		fprintf(out,"\n \n");

	for (i = 0; i < k; i++)
	{
		fprintf(out,"  (%lf, %lf) ", values[i].real,values[i].imag);

	}	

		fprintf(out,"\n \n");
for (i = 0; i < k; i++)
{
	fprintf(out,"  %d ", columns[i]);
}	

		fprintf(out,"\n \n \n");

for (i = 0; i <= n; i++)
{
	fprintf(out,"  %d ", rowindex[i]);
}	

fflush(out);
*/
	 fprintf(out,"\n number of non-zeros = %d, sparsity = %lf %", k, (1.-(double)k / n / n) * 100);
	 printf("\n number of non-zeros = %d, sparsity = %lf %", k, (1.-(double)k / n / n) * 100);
	 return k;
}


void ConjGrNR(scomplex **A,	// матрица системы (с точки зрения BLAS она транспонирована At)
			  int n,		// размер системной матрицы
			  scomplex ** B,// матрица правых частей (с точки зрения BLAS все нормально)
			  int m,		// число правых частей
			  scomplex ** G,// матрица начальных значений размер такой же как у B
			  double tol)	// точность контроля сходимости
{
	double *eps,  alpha, *ro1, *ro2, beta, res ;
	scomplex **q, **p,**r,**x, **u, al, bt;

	int size, lda,ldb,ldc, nn, mm, i, incx = 1, counter = 1,w;
	 char transa='N',transb='N';

	 nn = n; mm = m;

	
	eps = (double *)calloc(m, sizeof(double));
	ro1 = (double *)calloc(m, sizeof(double));
	ro2 = (double *)calloc(m, sizeof(double));
/*	alpha = (double *)calloc(m, sizeof(double));
	beta = (double *)calloc(m, sizeof(double));
*/
	q = AllocC(m,n);
	p = AllocC(m,n);
	r = AllocC(m,n);
	x = AllocC(m,n);
	u = AllocC(m,n);

	// r = Ah * B - Ah * A * G
	al.imag = 0.;
	bt.imag = 0.;

		lda = n; ldb = n; ldc = n;

	if (G == NULL)
	{
	// сначала делаем из нашей матрицы At Ah
	size = n * n;
	zlacgv (&size, A[0], &incx);
	
	// r = Ah * B
		bt.real = 0.;
	
		al.real = 1.;
		


		zgemm(&transa, &transb, &nn, &mm, &nn, &al, A[0], &lda, B[0], &ldb, &bt, r[0], &ldc);
	}
	else
	{
	// r = Ah * ( B - A * G )
		bt.real = 1.;
		al.real = -1.;
		


		// r = B
		size = n * m;
		zcopy( &size, B[0], &incx, r[0], &incx );
		// r = r - A * G
		transa = 'T';
		zgemm(&transa, &transb, &nn, &mm, &nn, &al, A[0], &lda, G[0], &ldb, &bt, r[0], &ldc);

		//  делаем из нашей матрицы At Ah
		size = n * n;
		zlacgv(&size, A[0], &incx);
		
		// r = Ah * r;
		bt.real = 0.;
		al.real = 1.;
		
		transa='N';



		zgemm(&transa, &transb, &nn, &mm, &nn, &al, A[0], &lda, r[0], &ldb, &bt, q[0], &ldc);
		size = n * m;
		zcopy( &size, q[0], &incx, r[0], &incx );

		//x = G
		zcopy( &size, G[0], &incx, x[0], &incx );
	}

	for ( i = 0; i < m; i++)
	{
		ro1[i] = 1.;
			eps[i] = 1.;
	}
	
	do
	{
		w = 0;
		for ( i = 0; i < m; i++)// цикл по правым частям		
		{
			
			res = dznrm2( &nn, r[i], &incx );
			ro2[i] = res * res;
			
		

			if (eps[i] < tol) continue;
			{
				w++;
			
				if ( counter == 1)
				{
					// p[i] = r[i]
					zcopy( &nn, r[i], &incx, p[i], &incx );
				}
				else
				{
					// beta = ro2 / ro1;
					beta = ro2[i] / ro1[i];
					/// p[i] = r[i] + beta * p[i]

						// q[i] = r[i]
						zcopy( &nn, r[i], &incx, q[i], &incx );
						
						// q[i] = q[i] + beta * p[i]				
						bt.real = beta;
						zaxpy ( &nn, &bt, p[i], &incx, q[i], &incx );
						// p[i] = q[i]
						zcopy( &nn, q[i], &incx, p[i], &incx );
					
				}
				
				/// q[i] = A * p[i]
					bt.real = 0;	
					al.real = 1.;
					transa='C'; // делаем из Ah A в текущем умножении

					zgemv( &transa, &nn, &nn, &al, A[0], &lda, p[i], &incx, &bt, q[i], &incx );

				///alpha = ro2 / ||q[i]||
					
					res = dznrm2( &nn, q[i], &incx );
					alpha = ro2[i] / res / res;

				/// x[i] = x[i] + alpha * p[i]
					// u[i] = alpha * p[i] делаем так , чтобы в явном виде получить приращение x
						al.real = alpha;
					// u[i] = p[i]
						zcopy( &nn, p[i], &incx, u[i], &incx );
					//u[i] = alpha * u[i]
						zscal ( &nn, &al, u[i], &incx );
						
					/// x[i] = x[i] + u[i]
						al.real = 1;
						zaxpy ( &nn, &al, u[i], &incx, x[i], &incx );

				/// r[i] = r[i] - alpha * Ah * q[i]

					bt.real = 1;	
					al.real = -alpha;
					transa='N'; // 

					zgemv( &transa, &nn, &nn, &al, A[0], &lda, q[i], &incx, &bt, r[i], &incx );
						
/*
					// q[i] = -alpha * q[i]
						al.real = -alpha;
						zscal ( &nn, &al, q[i], &incx );
					
					//r[i] = r[i] + q[i]
						al.real = 1.;
						zaxpy ( &nn, &al, q[i], &incx, r[i], &incx );
	*/			
				//
					beta = dznrm2 ( &nn, x[i], &incx );
					res = dznrm2 ( &nn, u[i], &incx );
					eps[i] = res / beta;

				ro1[i] = ro2[i]; // передаем текущее значение в предыдущее
				
			}
		}

	//	res = dnrm2(&mm,eps,&incx);

	counter ++;

	}while (w);
//while (( res > tol) && w);
	fprintf(out,"\n \n \n");

	fprintf(out,"Number of iterations is  %d ", counter);

		size = n * m;
		zcopy( &size, x[0], &incx, B[0], &incx );
	
	free(eps);
	free(ro1);
	free(ro2);

	DeAllocC(q);
	DeAllocC(p);
	DeAllocC(r);
	DeAllocC(x);
	DeAllocC(u);	
//	lda = 
//	zgemm (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
	

}



void ConjGrNRJacobi(scomplex **A,	// матрица системы (с точки зрения BLAS она транспонирована At)
			  int n,		// размер системной матрицы
			  scomplex ** B,// матрица правых частей (с точки зрения BLAS все нормально)
			  int m,		// число правых частей
			  scomplex ** G,// матрица начальных значений размер такой же как у B
			  double tol)	// точность контроля сходимости
{
	double *eps,   beta,  res ;
	scomplex **q, **p,**r,**x, **u, **z,*d, *ro1, *ro2, al, bt;

	int size, lda,ldb,ldc, nn, mm, i,qq = 0, incx = 1, counter = 1,w ;
	 char transa='N',transb='N';

	 nn = n; mm = m;

	
	eps = (double *)calloc(m, sizeof(double));
	ro1 = (scomplex *)calloc(m, sizeof(scomplex));
	ro2 = (scomplex *)calloc(m, sizeof(scomplex));
	d = (scomplex *)calloc(n, sizeof(scomplex));

/*	alpha = (double *)calloc(m, sizeof(double));
	beta = (double *)calloc(m, sizeof(double));
*/
	q = AllocC(m,n);
	p = AllocC(m,n);
	r = AllocC(m,n);
	x = AllocC(m,n);
	u = AllocC(m,n);
	z = AllocC(m,n); 
//	d = AllocC(m,n);	 

	// Making Jacobi preconditioner
	beta = tol * tol;
	for (i = 0; i < n; i++)
	{
		res = snorm(A[i][i]);
	if (res > beta) d[i].real = 1. / res;
	//		d[i].real = 1.;
	}

	// r = Ah * B - Ah * A * G
	al.imag = 0.;
	bt.imag = 0.;

		lda = n; ldb = n; ldc = n;

	if (G == NULL)
	{
	// сначала делаем из нашей матрицы At Ah
	size = n * n;
	zlacgv (&size, A[0], &incx);
	
	// r = Ah * B
		bt.real = 0.;
	
		al.real = 1.;
		


		zgemm(&transa, &transb, &nn, &mm, &nn, &al, A[0], &lda, B[0], &ldb, &bt, r[0], &ldc);
	}
	else
	{
	// r = Ah * ( B - A * G )
		bt.real = 1.;
		al.real = -1.;
		


		// r = B
		size = n * m;
		zcopy( &size, B[0], &incx, r[0], &incx );
		// r = r - A * G
		transa = 'T';
		zgemm(&transa, &transb, &nn, &mm, &nn, &al, A[0], &lda, G[0], &ldb, &bt, r[0], &ldc);

		//  делаем из нашей матрицы At Ah
		size = n * n;
		zlacgv(&size, A[0], &incx);
		
		// r = Ah * r;
		bt.real = 0.;
		al.real = 1.;
		
		transa='N';



		zgemm(&transa, &transb, &nn, &mm, &nn, &al, A[0], &lda, r[0], &ldb, &bt, q[0], &ldc);
		size = n * m;
		zcopy( &size, q[0], &incx, r[0], &incx );

		//x = G
		zcopy( &size, G[0], &incx, x[0], &incx );
	}

	for ( i = 0; i < m; i++)
	{
	//	ro1[i] = 1.;
			eps[i] = 1.;
	}
	
	do
	{
		w = 0;
		for ( i = 0; i < m; i++)// цикл по правым частям		
		{
			
		
		

			if (eps[i] < tol) continue;

			// Preconditioning

				// z[q] = r[q] * d[q]	
/*			for (qq = 0; qq < n; qq++)
			{
				z[i][qq].real = r[i][qq].real * d[qq]; 
					z[i][qq].imag = r[i][qq].imag * d[qq]; 
			}
*/
			al.imag = bt.imag = 0.;
				bt.real = 0.;
				al.real = 1.;

			zgbmv ( &transb, &nn, &nn, &qq, &qq, &al, d, &incx, r[i], &incx, &bt, z[i], &incx );
	
			zdotc ( &ro2[i], &nn, r[i], &incx, z[i], &incx );
			//res = dznrm2( &nn, r[i], &incx );
		//	ro2[i] = res * res;
			
			{
				w++;
			
				if ( counter == 1)
				{
					// p[i] = z[i]
					zcopy( &nn, z[i], &incx, p[i], &incx );
				}
				else
				{
					// beta = ro2 / ro1;
					zladiv (&bt, &ro2[i] , &ro1[i]);
					/// p[i] = z[i] + beta * p[i]

						// q[i] = z[i]
						zcopy( &nn, z[i], &incx, q[i], &incx );
						
						// q[i] = q[i] + beta * p[i]				
						//bt.real = beta;
						zaxpy ( &nn, &bt, p[i], &incx, q[i], &incx );
						// p[i] = q[i]
						zcopy( &nn, q[i], &incx, p[i], &incx );
					
				}
				
				/// q[i] = A * p[i]
					al.imag = bt.imag = 0.;
					bt.real = 0;	
					al.real = 1.;
					transa='C'; // делаем из Ah A в текущем умножении

					zgemv( &transa, &nn, &nn, &al, A[0], &lda, p[i], &incx, &bt, q[i], &incx );

				///alpha = ro2 / ||q[i]||
					
					res = dznrm2( &nn, q[i], &incx );
					//alpha = ro2[i] / res / res;
						al.real = ro2[i].real / res / res;
						al.imag = ro2[i].imag / res / res;
				/// x[i] = x[i] + alpha * p[i]
					// u[i] = alpha * p[i] делаем так , чтобы в явном виде получить приращение x
						//al.real = alpha;
					// u[i] = p[i]
						zcopy( &nn, p[i], &incx, u[i], &incx );
					//u[i] = alpha * u[i]
						zscal ( &nn, &al, u[i], &incx );
						
					/// x[i] = x[i] + u[i]
							al.imag = bt.imag = 0.;
						bt.real = 1;
						zaxpy ( &nn, &bt, u[i], &incx, x[i], &incx );

				/// r[i] = r[i] - alpha * Ah * q[i]

					bt.real = 1;	
				//	al.real = -alpha;	
					al.real *= -1.;
					transa='N'; // 

					zgemv( &transa, &nn, &nn, &al, A[0], &lda, q[i], &incx, &bt, r[i], &incx );
						

					beta = dznrm2 ( &nn, x[i], &incx );
					res = dznrm2 ( &nn, u[i], &incx );
					eps[i] = res / beta;

				ro1[i].real = ro2[i].real;
				ro1[i].imag = ro2[i].imag;
				// передаем текущее значение в предыдущее
				
			}
		}

	//	res = dnrm2(&mm,eps,&incx);

	counter ++;

	}while (w);
//while (( res > tol) && w);
	fprintf(out,"\n \n \n");

	fprintf(out,"Number of iterations is  %d ", counter);

		size = n * m;
		zcopy( &size, x[0], &incx, B[0], &incx );
	
	free(eps);
	free(ro1);
	free(ro2);
	free(d);
//
	DeAllocC(q);
	DeAllocC(p);
	DeAllocC(r);
	DeAllocC(x);
	DeAllocC(u);
	DeAllocC(z);
//	DeAllocC(d);	
//	lda = 
//	zgemm (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)	

}



void ConjGrM(scomplex **A,	// матрица системы (с точки зрения BLAS она транспонирована At)
			  int n,		// размер системной матрицы
			  scomplex ** B,// матрица правых частей (с точки зрения BLAS все нормально)
			  int m,		// число правых частей
			  scomplex ** G,// матрица начальных значений размер такой же как у B
			  double tol)	// точность контроля сходимости
{
	// обычный метод сопряженных градиентов, рабтает только для эрмитовых(самосопряженных) матриц
	double *eps,    beta,  res ;
	scomplex **q, **p,**r,**x, **u, **z,*d, *ro1, *ro2, al, bt;

	int size, lda,ldb,ldc, nn, mm, i,qq = 0, incx = 1, counter = 1,w ;
	 char transa='N',transb='N';

	 nn = n; mm = m;

	
	eps = (double *)calloc(m, sizeof(double));
	ro1 = (scomplex *)calloc(m, sizeof(scomplex));
	ro2 = (scomplex *)calloc(m, sizeof(scomplex));
	d = (scomplex *)calloc(n, sizeof(scomplex));

/*	alpha = (double *)calloc(m, sizeof(double));
	beta = (double *)calloc(m, sizeof(double));
*/
	q = AllocC(m,n);
	p = AllocC(m,n);
	r = AllocC(m,n);
	x = AllocC(m,n);
	u = AllocC(m,n);
	z = AllocC(m,n); 
//	d = AllocC(m,n);	 

	al.imag = 0.;
	bt.imag = 0.;
	bt.real = 1.;
	// Making Jacobi preconditioner
/**/	beta = tol * tol;

	for (i = 0; i < n; i++)
	{
		// d[i].real = 1. / res;
		 	zladiv (&d[i], &bt , &A[i][i]);
	//		d[i].real = 1.;
	}

	// r =  B -  A * G


		lda = n; ldb = n; ldc = n;

	if (G == NULL)
	{
		// r = B
	
				size = n * m;
		zcopy( &size, B[0], &incx, r[0], &incx );


	
	}
	else
	{
	// r =  ( B - A * G )
		bt.real = 1.;
		al.real = -1.;
		


		// r = B
		size = n * m;
		zcopy( &size, B[0], &incx, r[0], &incx );
		// r = r - A * G
		transa = 'T';
		zgemm(&transa, &transb, &nn, &mm, &nn, &al, A[0], &lda, G[0], &ldb, &bt, r[0], &ldc);


		//x = G
		zcopy( &size, G[0], &incx, x[0], &incx );
	}

	for ( i = 0; i < m; i++)
	{
	//	ro1[i] = 1.;
			eps[i] = 1.;
	}
	
	do
	{
		w = 0;
		for ( i = 0; i < m; i++)// цикл по правым частям		
		{
			
		
			if (eps[i] < tol) continue;

			al.imag = bt.imag = 0.;
				bt.real = 0.;
				al.real = 1.;


			res = dznrm2( &nn, r[i], &incx );
			ro2[i].real = res * res;
			ro2[i].imag = 0.;
			
			{
				w++;
			
				if ( counter == 1)
				{
					// p[i] = z[i]
					zcopy( &nn, r[i], &incx, p[i], &incx );
				}
				else
				{
					// beta = ro2 / ro1;
					zladiv (&bt, &ro2[i] , &ro1[i]);
					/// p[i] = r[i] + beta * p[i]

						// q[i] = r[i]
						zcopy( &nn, r[i], &incx, q[i], &incx );
						
						// q[i] = q[i] + beta * p[i]				
						//bt.real = beta;
						zaxpy ( &nn, &bt, p[i], &incx, q[i], &incx );
						// p[i] = q[i]
						zcopy( &nn, q[i], &incx, p[i], &incx );
					
				}
				
				/// q[i] = A * p[i]
					al.imag = bt.imag = 0.;
					bt.real = 0;	
					al.real = 1.;
					transa='T'; // делаем из At A в текущем умножении

					zgemv( &transa, &nn, &nn, &al, A[0], &lda, p[i], &incx, &bt, q[i], &incx );

				///alpha = ro2 / (p[i] ,q[i])
					// bt = (p[i] ,q[i])
					zdotc ( &bt, &nn, p[i], &incx, q[i], &incx );
					//alpha = ro2 / bt
					zladiv (&al, &ro2[i] , &bt);
					

				/// x[i] = x[i] + alpha * p[i]
				//	// u[i] = alpha * p[i] делаем так , чтобы в явном виде получить приращение x
						//al.real = alpha;
					// u[i] = p[i]
						zcopy( &nn, p[i], &incx, u[i], &incx );
					//u[i] = alpha * u[i]
						zscal ( &nn, &al, u[i], &incx );
						
					/// x[i] = x[i] + u[i]
						//	al.imag = bt.imag = 0.;
						bt.real = 1;
						zaxpy ( &nn, &bt, u[i], &incx, x[i], &incx );

				/// r[i] = r[i] - alpha *  q[i]
					bt.imag = 0.;
					bt.real = 1.;	
				//	al.real = -alpha;	
					al.real *= -1.;
					al.imag *= -1.;
					//transa='N'; // 

					zaxpy ( &nn, &al, q[i], &incx, r[i], &incx );

					//zgemv( &transa, &nn, &nn, &al, A[0], &lda, q[i], &incx, &bt, r[i], &incx );
						

					beta = dznrm2 ( &nn, x[i], &incx );
					res = dznrm2 ( &nn, u[i], &incx );
					eps[i] = res / beta;

				ro1[i].real = ro2[i].real;
				ro1[i].imag = ro2[i].imag;
				// передаем текущее значение в предыдущее
				
			}
		}

	//	res = dnrm2(&mm,eps,&incx);

	counter ++;

	}while (w && (counter < 1000));
//while (( res > tol) && w);
	fprintf(out,"\n \n \n");

	fprintf(out,"Number of iterations is  %d ", counter);

		size = n * m;
		zcopy( &size, x[0], &incx, B[0], &incx );
	
	free(eps);
	free(ro1);
	free(ro2);
	free(d);
//
	DeAllocC(q);
	DeAllocC(p);
	DeAllocC(r);
	DeAllocC(x);
	DeAllocC(u);
	DeAllocC(z);
//	DeAllocC(d);	
//	lda = 
//	zgemm (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)	

}


void BiCGStab(scomplex **A,	// матрица системы (с точки зрения BLAS она транспонирована At)
			  int n,		// размер системной матрицы
			  scomplex ** B,// матрица правых частей (с точки зрения BLAS все нормально)
			  int m,		// число правых частей
			  scomplex ** G,// матрица начальных значений размер такой же как у B
			  double tol)	// точность контроля сходимости
{
	// 
	double *eps,     res, dtol ;
	scomplex **q,*t, **p,**r,**x, **u,*d, **z, **r0,**v, *s,  *ro1, *ro2, al, bt,  *omega, *alpha;

	int size, lda,ldb,ldc, nn, mm, i,qq = 0, incx = 1, counter = 1,w ;
	 char transa='N',transb='N';

	 nn = n; mm = m;

	
	eps = (double *)calloc(m, sizeof(double));
	ro1 = (scomplex *)calloc(m, sizeof(scomplex));
	ro2 = (scomplex *)calloc(m, sizeof(scomplex));
	d = (scomplex *)calloc(n, sizeof(scomplex));
	omega = (scomplex *)calloc(m, sizeof(scomplex));
	alpha = (scomplex *)calloc(m, sizeof(scomplex));
	s = (scomplex *)calloc(n, sizeof(scomplex));
	t = (scomplex *)calloc(n, sizeof(scomplex));

	q = AllocC(m,n);
	p = AllocC(m,n);
	r = AllocC(m,n);
	r0 = AllocC(m,n);
	x = AllocC(m,n);
	u = AllocC(m,n);
	z = AllocC(m,n);
	v = AllocC(m,n); 
	 	 
	 
	al.imag = 0.;
	al.real = 1.;
	bt.imag = 0.;
	bt.real = 1.;


	// r =  B -  A * G


		lda = n; ldb = n; ldc = n;
/*
*/	if (G == NULL)
	{
		// r = B
	
				size = n * m;
		zcopy( &size, B[0], &incx, r[0], &incx );


	
	}
	else
	{
	// r =  ( B - A * G )
		bt.real = 1.;
		al.real = -1.;
		


		// r = B
		size = n * m;
		zcopy( &size, B[0], &incx, r[0], &incx );
		// r = r - A * G
		transa = 'T';
		zgemm(&transa, &transb, &nn, &mm, &nn, &al, A[0], &lda, G[0], &ldb, &bt, r[0], &ldc);


		//x = G
		zcopy( &size, G[0], &incx, x[0], &incx );
	}


		
	for ( i = 0; i < m; i++)
	{

			eps[i] = 1.;

	}

// fix r0
		size = n * m;
			zcopy( &size, r[0], &incx, r0[0], &incx );
/*		fprintf(stdout," \n norm of r0 = %lf \n",dznrm2( &size, r0[0], &incx ));
		fprintf(out," \n norm of r0 = %lf \n",dznrm2( &size, r0[0], &incx ));
*/	//	

		dtol = tol * tol;
	do
	{
		w = 0;

		for ( i = 0; i < m; i++)// цикл по правым частям		
		{
			
		
			if (eps[i] < tol) continue;

			al.imag = bt.imag = 0.;
				bt.real = 0.;
				al.real = 1.;

/*
			res = dznrm2( &nn, r[i], &incx );
			ro2[i].real = res * res;
			ro2[i].imag = 0.;*/
		
				// ro2 = (r0, r)
			zdotc ( &ro2[i], &nn, r0[i], &incx, r[i], &incx );
			if (sqrt(snorm(ro2[i])) < 1.e-40) 
			{
				printf ("\nBiCGStab fails, since ro = (%e, %e) at column i = %d \n",ro2[i].real,ro2[i].imag,i);
				//exit(1);
				eps[i] = 0;
			}
	/**/		
			{
				w++;
			
				if ( counter == 1)
				{
					// p[i] = r[i]
					zcopy( &nn, r[i], &incx, p[i], &incx );
				//	d[i] = dznrm2( &nn, r0[i], &incx );
				}
				else
				{
					// beta = (ro2 / ro1)(alpha / omega);
					zladiv (&bt, &ro2[i] , &ro1[i]);
					zladiv (&al, &alpha[i] , &omega[i]);

					bt.real = bt.real * al.real - bt.imag * al.imag;
					bt.imag = bt.real * al.imag + bt.imag * al.real;

					/// p[i] = r[i] + beta * (p[i] - omega[i] * v[i])
						// p[i] - omega[i] * v[i]
						al.real = -1. * omega[i].real;
						al.imag = -1. * omega[i].imag;

						zaxpy ( &nn, &al, v[i], &incx, p[i], &incx );
						//p[i] *= beta
						zscal ( &nn, &bt, p[i], &incx );
						// p[i] = r[i] + p[i]
						al.real = 1.;
						al.imag = 0.;

						zaxpy ( &nn, &al, r[i], &incx, p[i], &incx );

				}

					

				/// v[i] = A * p[i]
					al.imag = bt.imag = 0.;
					bt.real = 0;	
					al.real = 1.;
					transa='T'; // делаем из At A в текущем умножении

				zgemv( &transa, &nn, &nn, &al, A[0], &lda, p[i], &incx, &bt, v[i], &incx );


				///alpha = ro2 / (r0[i] ,v[i])

					//zdotc ( &ro2[i], &nn, r0[i], &incx, r[i], &incx );
					// bt = (p[i] ,q[i])
					zdotc ( &bt, &nn, r0[i], &incx, v[i], &incx );
					//alpha = ro2 / bt
					zladiv (&alpha[i], &ro2[i] , &bt);
				///s = r[i] - alpha[i] * v[i]
					// s = r[i]
					zcopy( &nn, r[i], &incx, s, &incx );
					// s = s -  alpha[i] * v[i]
						al.real = -1. * alpha[i].real;
						al.imag = -1. * alpha[i].imag;
					zaxpy ( &nn, &al, v[i], &incx, s, &incx );

				/// check norm of s
					res = dznrm2( &nn, s, &incx );
					if (res < tol)
					{
						// x[i] = x[i] + alpha[i] * p[i]
						zaxpy ( &nn, &alpha[i], p[i], &incx, x[i], &incx );	
	
						eps[i] = 0.;
						continue;
					}


				/// t = A * s
					al.imag = bt.imag = 0.;
					bt.real = 0;	
					al.real = 1.;
					transa='T'; // делаем из At A в текущем умножении

					zgemv( &transa, &nn, &nn, &al, A[0], &lda, s, &incx, &bt, t, &incx );

				///omega[i] = (t,s)/(t,t)
					//al = (t,s)
					zdotc ( &al, &nn, t, &incx, s, &incx );
					// bt = (t,t)
					zdotc ( &bt, &nn, t, &incx, t, &incx );
					//omega[i] = al / bt
					zladiv (&omega[i], &al , &bt);
				///x[i] = x[i] + alpha[i] * p[i] + omega[i] * s
					// u[i] = [pi]
					zcopy( &nn,  p[i], &incx,u[i], &incx );
					// u[i] *= alpha[i]
					zscal ( &nn, &alpha[i], u[i], &incx );
					// u[i] = u[i] + omega[i] * s
				
					zaxpy ( &nn, &omega[i], s, &incx, u[i], &incx );	
					// x[i] += u[i]

					al.real = 1.;
					al.imag = 0.;
					zaxpy ( &nn, &al, u[i], &incx, x[i], &incx );	
				///r[i] = s - omega[i] * t
					// r[i] = s
					zcopy( &nn,  s, &incx,r[i], &incx );

					al.real = -1. * omega[i].real;
					al.imag = -1. * omega[i].imag;
					zaxpy ( &nn, &al, t, &incx, r[i], &incx );	

				res = dznrm2( &nn, u[i], &incx );

				eps[i] = res/dznrm2( &nn, x[i], &incx );;// / d[i];

				if (snorm(omega[i]) < dtol) eps[i] = 0.;

				ro1[i].real = ro2[i].real;
				ro1[i].imag = ro2[i].imag;
				// передаем текущее значение в предыдущее
				
			}
		}

	//	res = dnrm2(&mm,eps,&incx);

	counter ++;

	}while (w && (counter < 10000));
//while (( res > tol) && w);
	fprintf(out,"\n \n \n");

	fprintf(out,"Number of iterations is  %d ", counter);
	fflush(out);
		size = n * m;
		zcopy( &size, x[0], &incx, B[0], &incx );
	
	free(eps);
	free(ro1);
	free(ro2);
	free(d);
	free(omega);
	free(alpha);//
	DeAllocC(q);
	DeAllocC(p);
	DeAllocC(r);
	DeAllocC(r0);
	DeAllocC(x);
	DeAllocC(u);
	DeAllocC(z);
	DeAllocC(v);
	free(t);
	free(s);



}

void DILUPivots(int n,			// размер матрицы системы
				scomplex ** A,	// матрица системы
				int w,			// число боковых диагоналей используемых для генерции вектора диагональных элементов
				scomplex * d)	// вектор содержащий обратные величины диагональных элементов DILU preconditioner
{
int incy, size, incx  = 1, i,j;
scomplex al, bt;
	// Making D-ILU preconditioner
	// d[i] = A[i][i]
	incy = n + 1;
	size = n;
		zcopy( &size, A[0], &incy, d, &incx );

		//w = 8;
		bt.real = 1.;
		bt.imag = 0.;
		for (i = 0; i < n; i++)
		{
			al.real = d[i].real;
			al.imag = d[i].imag;
			zladiv(&d[i],&bt,&al); //1/d[i]
			
			for (j = i + 1; j < __min(i + w,n); j++)
			{
				//d[j] = d[j] - A[j][i] * d[i] * A[i][j];
				al.real = A[j][i].real * A[i][j].real - A[j][i].imag * A[i][j].imag ;
				al.imag = A[j][i].real * A[i][j].imag + A[j][i].imag * A[i][j].real ;

				d[j].real -= d[i].real * al.real - d[i].imag * al.imag;
				d[j].imag -= d[i].real * al.imag + d[i].imag * al.real;
			}
		}
}


void DILUSolve(int n,			// размер матрицы системы
				scomplex ** A,	// матрица системы
				int w,			// число боковых диагоналей используемых для генерции вектора диагональных элементов
				scomplex * d,	// вектор содержащий обратные величины диагональных элементов DILU preconditioner
				scomplex * y,	// вектор правой части системы
			   scomplex *x)		// вектор решений
{
int   incx  = 1, i,j;
scomplex al,  *z;
	
z = (scomplex *)calloc(n, sizeof(scomplex));



		for (i = 0; i < n; i++)
		{
			
			al.real = al.imag = 0.; 

			for (j =__max(i - w,0) ; j < i ; j++)
			{
				//al += A[j][i] * s[i] ;
				al.real += A[j][i].real * z[j].real - A[j][i].imag * z[j].imag ;
				al.imag += A[j][i].real * z[j].imag + A[j][i].imag * z[j].real ;

			
			}
			z[i].real = y[i].real - al.real;
			z[i].imag = y[i].imag - al.imag;

				z[i].real = d[i].real * z[i].real - d[i].imag * z[i].imag ;
				z[i].imag = d[i].real * z[i].imag + d[i].imag * z[i].real ;


		}

		for (i = n - 1 ; i >= 0; i--)
		{
			al.real = al.imag = 0.; 
			for (j = i + 1; j < __min(i + w,n); j++)
			{
				//al += A[j][i] * s[i] ;
				al.real += A[j][i].real * z[j].real - A[j][i].imag * z[j].imag ;
				al.imag += A[j][i].real * z[j].imag + A[j][i].imag * z[j].real ;
			}				
			
				al.real = d[i].real * al.real - d[i].imag * al.imag ;
				al.imag = d[i].real * al.imag + d[i].imag * al.real ;

			x[i].real = z[i].real - al.real;
			x[i].imag = z[i].imag - al.imag;

		}

	free(z);
}

void BiCGStabPreCond(scomplex **A,	// матрица системы (с точки зрения BLAS она транспонирована At)
			  int n,		// размер системной матрицы
			  scomplex ** B,// матрица правых частей (с точки зрения BLAS все нормально)
			  int m,		// число правых частей
			  scomplex ** G,// матрица начальных значений размер такой же как у B
			  double tol)	// точность контроля сходимости
{
	// 
	double *eps,     res, dtol;
	scomplex **q,*t, **p,**r,**x, **u,*d, **z, **r0,**v, *s,*st,*pt,  *ro1, *ro2, al, bt,  *omega, *alpha;

	int size, lda,ldb,ldc, nn, mm, i,qq = 0, incx = 1, counter = 1,w, band ;
	 char transa='N',transb='N';

	 nn = n; mm = m;

	
	eps = (double *)calloc(m, sizeof(double));
	ro1 = (scomplex *)calloc(m, sizeof(scomplex));
	ro2 = (scomplex *)calloc(m, sizeof(scomplex));
	d = (scomplex *)calloc(n, sizeof(scomplex));
	omega = (scomplex *)calloc(m, sizeof(scomplex));
	alpha = (scomplex *)calloc(m, sizeof(scomplex));
	s = (scomplex *)calloc(n, sizeof(scomplex));
	t = (scomplex *)calloc(n, sizeof(scomplex));
	st = (scomplex *)calloc(n, sizeof(scomplex));
	pt = (scomplex *)calloc(n, sizeof(scomplex));
	/*	beta = (double *)calloc(m, sizeof(double));
*/
	q = AllocC(m,n);
	p = AllocC(m,n);
	r = AllocC(m,n);
	r0 = AllocC(m,n);
	x = AllocC(m,n);
	u = AllocC(m,n);
	z = AllocC(m,n);
	v = AllocC(m,n); 
	 	 
//	d = AllocC(m,n);	 
	al.imag = 0.;
	al.real = 1.;
	bt.imag = 0.;
	bt.real = 1.;

/*
	// Making D-ILU preconditioner
	// d[i] = A[i][i]
	incy = n + 1;
	size = n;
		zcopy( &size, A[0], &incy, d, &incx );

		w = 8;

		for (i = 0; i < n; i++)
		{
			al.real = d[i].real;
			al.imag = d[i].imag;
			zladiv(&d[i],&bt,&al); //1/d[i]
			
			for (j = i + 1; j < __min(i + w,n); j++)
			{
				//d[j] = d[j] - A[j][i] * d[i] * A[i][j];
				al.real = A[j][i].real * A[i][j].real - A[j][i].imag * A[i][j].imag ;
				al.imag = A[j][i].real * A[i][j].imag + A[j][i].imag * A[i][j].real ;

				d[j].real -= d[i].real * al.real - d[i].imag * al.imag;
				d[j].imag -= d[i].real * al.imag + d[i].imag * al.real;
			}
		}
*/
	band = n;
	DILUPivots(n,A,band,d);
//
	// r =  B -  A * G


		lda = n; ldb = n; ldc = n;
/*
*/	if (G == NULL)
	{
		// r = B
	
				size = n * m;
		zcopy( &size, B[0], &incx, r[0], &incx );


	
	}
	else
	{
	// r =  ( B - A * G )
		bt.real = 1.;
		al.real = -1.;
		


		// r = B
		size = n * m;
		zcopy( &size, B[0], &incx, r[0], &incx );
		// r = r - A * G
		transa = 'T';
		zgemm(&transa, &transb, &nn, &mm, &nn, &al, A[0], &lda, G[0], &ldb, &bt, r[0], &ldc);


		//x = G
		zcopy( &size, G[0], &incx, x[0], &incx );
	}


		
	for ( i = 0; i < m; i++)
	{
	//	ro1[i] = 1.;
			eps[i] = 1.;
		//	zcopy( &nn, d, &incx, r0[i], &incx );
	}

// fix r0
		size = n * m;
			zcopy( &size, r[0], &incx, r0[0], &incx );
		
	//	

		dtol = tol * tol;
	do
	{
		w = 0;

		for ( i = 0; i < m; i++)// цикл по правым частям		
		{
			
		
			if (eps[i] < tol) continue;

			al.imag = bt.imag = 0.;
				bt.real = 0.;
				al.real = 1.;

/*
			res = dznrm2( &nn, r[i], &incx );
			ro2[i].real = res * res;
			ro2[i].imag = 0.;*/
		
				// ro2 = (r0, r)
			zdotc ( &ro2[i], &nn, r0[i], &incx, r[i], &incx );
	/*		if (snorm(ro2[i]) < dtol) 
			{
				printf ("BiCGStab fails, since ro = 0");
				exit(1);
			}
	*/		
			{
				w++;
			
				if ( counter == 1)
				{
					// p[i] = r[i]
					zcopy( &nn, r[i], &incx, p[i], &incx );
				//	d[i] = dznrm2( &nn, r0[i], &incx );
				}
				else
				{
					// beta = (ro2 / ro1)(alpha / omega);
					zladiv (&bt, &ro2[i] , &ro1[i]);
					zladiv (&al, &alpha[i] , &omega[i]);

					bt.real = bt.real * al.real - bt.imag * al.imag;
					bt.imag = bt.real * al.imag + bt.imag * al.real;

					/// p[i] = r[i] + beta * (p[i] - omega[i] * v[i])
						// p[i] - omega[i] * v[i]
						al.real = -1. * omega[i].real;
						al.imag = -1. * omega[i].imag;

						zaxpy ( &nn, &al, v[i], &incx, p[i], &incx );
						//p[i] *= beta
						zscal ( &nn, &bt, p[i], &incx );
						// p[i] = r[i] + p[i]
						al.real = 1.;
						al.imag = 0.;

						zaxpy ( &nn, &al, r[i], &incx, p[i], &incx );

				}
				//1st preconditioning
				// M * pt = p[i]

		//zcopy( &nn, p[i], &incx, pt, &incx );
				DILUSolve(n,A,band,d,p[i],pt);
				/// v[i] = A * p[i]
					al.imag = bt.imag = 0.;
					bt.real = 0;	
					al.real = 1.;
					transa='T'; // делаем из At A в текущем умножении

				//	zgemv( &transa, &nn, &nn, &al, A[0], &lda, p[i], &incx, &bt, v[i], &incx );
						zgemv( &transa, &nn, &nn, &al, A[0], &lda, pt, &incx, &bt, v[i], &incx );

				///alpha = ro2 / (r0[i] ,v[i])

					//zdotc ( &ro2[i], &nn, r0[i], &incx, r[i], &incx );
					// bt = (p[i] ,q[i])
					zdotc ( &bt, &nn, r0[i], &incx, v[i], &incx );
					//alpha = ro2 / bt
					zladiv (&alpha[i], &ro2[i] , &bt);
				///s = r[i] - alpha[i] * v[i]
					// s = r[i]
					zcopy( &nn, r[i], &incx, s, &incx );
					// s = s -  alpha[i] * v[i]
						al.real = -1. * alpha[i].real;
						al.imag = -1. * alpha[i].imag;
					zaxpy ( &nn, &al, v[i], &incx, s, &incx );

				/// check norm of s
					res = dznrm2( &nn, s, &incx );
					if (res < tol)
					{
						// x[i] = x[i] + alpha[i] * p[i]
					//	zaxpy ( &nn, &alpha[i], p[i], &incx, x[i], &incx );	
						zaxpy ( &nn, &alpha[i], pt, &incx, x[i], &incx );	
						eps[i] = 0.;
						continue;
					}

				//2nd preconditioning
				// M * st = s

				DILUSolve(n,A,band,d,s,st);
			//		zcopy( &nn, s, &incx, st, &incx );

				/// t = A * s
					al.imag = bt.imag = 0.;
					bt.real = 0;	
					al.real = 1.;
					transa='T'; // делаем из At A в текущем умножении

					//zgemv( &transa, &nn, &nn, &al, A[0], &lda, s, &incx, &bt, t, &incx );
					zgemv( &transa, &nn, &nn, &al, A[0], &lda, st, &incx, &bt, t, &incx );
				///omega[i] = (t,s)/(t,t)
					//al = (t,s)
					zdotc ( &al, &nn, t, &incx, s, &incx );
					// bt = (t,t)
					zdotc ( &bt, &nn, t, &incx, t, &incx );
					//omega[i] = al / bt
					zladiv (&omega[i], &al , &bt);
				///x[i] = x[i] + alpha[i] * p[i] + omega[i] * s
					//zcopy( &nn,  p[i], &incx,u[i], &incx );
					zcopy( &nn,  pt, &incx,u[i], &incx );
					zscal ( &nn, &alpha[i], u[i], &incx );
					zaxpy ( &nn, &omega[i], st, &incx, u[i], &incx );	
						
					al.imag  = 0.;
					al.real = 1.;

					zaxpy ( &nn, &al, u[i], &incx, x[i], &incx );	

/*
					zaxpy ( &nn, &alpha[i], p[i], &incx, x[i], &incx );	
					zaxpy ( &nn, &omega[i], s, &incx, x[i], &incx );	
*/				///r[i] = s - omega[i] * t
					// r[i] = s
					zcopy( &nn,  s, &incx,r[i], &incx );

					al.real = -1. * omega[i].real;
					al.imag = -1. * omega[i].imag;
					zaxpy ( &nn, &al, t, &incx, r[i], &incx );	

				res = dznrm2( &nn, u[i], &incx );

				eps[i] = res/dznrm2( &nn, x[i], &incx );;// / d[i];
				if (snorm(omega[i]) < dtol) eps[i] = 0.;

				ro1[i].real = ro2[i].real;
				ro1[i].imag = ro2[i].imag;
				// передаем текущее значение в предыдущее
				
			}
		}

	//	res = dnrm2(&mm,eps,&incx);

	counter ++;

	}while (w && (counter < 10000));
//while (( res > tol) && w);
	fprintf(out,"\n \n \n");

	fprintf(out,"Number of iterations is  %d ", counter);
	fflush(out);
		size = n * m;
		zcopy( &size, x[0], &incx, B[0], &incx );
	
	free(eps);
	free(ro1);
	free(ro2);
	free(d);
	free(omega);
	free(alpha);//
	DeAllocC(q);
	DeAllocC(p);
	DeAllocC(r);
	DeAllocC(r0);
	DeAllocC(x);
	DeAllocC(u);
	DeAllocC(z);
	DeAllocC(v);
	free(t);
	free(s);
	free(pt);
	free(st);
	//	DeAllocC(d);	
//	lda = 
//	zgemm (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)	

}



void PrintMatrixC (scomplex ** A, int n,int m, FILE * f)
{
	int i,j;

	fprintf(f," \n\n");
		for (i = 0; i < n ; i++)
			{
	
			fprintf(f," \n");
					
				for (j = 0; j < m; j++)
					fprintf(f," \t(%8.3lf, %8.3lf)", A[i][j].real, A[i][j].imag);

			}
}


void ConjGrNRBlock(scomplex **A,	// матрица системы (с точки зрения BLAS она транспонирована At)
			  int n,		// размер системной матрицы
			  scomplex ** B,// матрица правых частей (с точки зрения BLAS все нормально)
			  int m,		// число правых частей
			  scomplex ** G,// матрица начальных значений размер такой же как у B
			  double tol)	// точность контроля сходимости
{
	double *eps,     res ;
	scomplex **q, **p,**r,**x, **u,alpha, **beta,*ro1, *ro2,al, bt,resc;

	int size, lda,ldb,ldc, nn, mm, i, j,incx = 1, counter = 1,w;
	 char transa='N',transb='N';

	 nn = n; mm = m;

	
	eps = (double *)calloc(m, sizeof(double));
	ro1 = (scomplex *)calloc(m, sizeof(scomplex));
	ro2 = (scomplex *)calloc(m, sizeof(scomplex));
/*	alpha = (double *)calloc(m, sizeof(double));
	beta = (double *)calloc(m, sizeof(double));
*/
	q = AllocC(m,n);
	p = AllocC(m,n);
	r = AllocC(m,n);
	x = AllocC(m,n);
	u = AllocC(m,n);
//	ro2 = AllocC(m,m);
//	ro1 = AllocC(m,m);
	beta = AllocC(m,m);
	// r = Ah * B - Ah * A * G
	al.imag = 0.;
	bt.imag = 0.;

		lda = n; ldb = n; ldc = n;

	if (G == NULL)
	{
	// сначала делаем из нашей матрицы At Ah
	size = n * n;
	zlacgv (&size, A[0], &incx);
	
	// r = Ah * B
		bt.real = 0.;
	
		al.real = 1.;
		


		zgemm(&transa, &transb, &nn, &mm, &nn, &al, A[0], &lda, B[0], &ldb, &bt, r[0], &ldc);
	}
	else
	{
	// r = Ah * ( B - A * G )
		bt.real = 1.;
		al.real = -1.;
		


		// r = B
		size = n * m;
		zcopy( &size, B[0], &incx, r[0], &incx );
		// r = r - A * G
		transa = 'T';
		zgemm(&transa, &transb, &nn, &mm, &nn, &al, A[0], &lda, G[0], &ldb, &bt, r[0], &ldc);

		//  делаем из нашей матрицы At Ah
		size = n * n;
		zlacgv(&size, A[0], &incx);
		
		// r = Ah * r;
		bt.real = 0.;
		al.real = 1.;
		
		transa='N';



		zgemm(&transa, &transb, &nn, &mm, &nn, &al, A[0], &lda, r[0], &ldb, &bt, q[0], &ldc);
		size = n * m;
		zcopy( &size, q[0], &incx, r[0], &incx );

		//x = G
		zcopy( &size, G[0], &incx, x[0], &incx );
	}

	for ( i = 0; i < m; i++)
	{
		//ro1[i] = 1.;
			eps[i] = 1.;
	}
	
	do
	{
		w = 0;
		

	

	size = n * m;
		if ( counter < 2)
				{
					// p[i] = r[i]
				
					zcopy( &size, r[0], &incx, p[0], &incx );
				}
		//		else
				{

		for ( i = 0; i < m; i++)// цикл по правым частям		
		{//
			 	zdotc(&ro2[i], &nn, r[i], &incx, r[i], &incx);
			for ( j = 0; j < m; j++)// цикл по правым частям	
			{
					
				zdotc(&beta[i][j], &nn, p[i], &incx, p[j], &incx);
			}
		}

					// q[i] = p[i]////
						zcopy( &size, p[0], &incx, q[0], &incx );

					for ( i = 0; i < m; i++)// цикл по правым частям
					{
						for ( j = 0; j < i; j++)// цикл по правым частям	
						{
					//			//p[i] -= beta[i][j] * q[j]
							beta[i][j].real *= 1.;
							beta[i][j].imag *= 1.;	
				//	if (i!=j)	
						zaxpy ( &nn, &beta[i][j], q[j], &incx, p[i], &incx );
							
							
							
						}

					


					}
				
					for ( i = 0; i < m; i++)// цикл по правым частям
					{
						// p[i] *= ro2/ro1


						bt.real = dznrm2(&nn,p[i],&incx );
						
							bt.imag = 0.;
							
							zladiv(&al,&ro2[i],&bt);
/**/
						zscal ( &nn, &al, p[i], &incx );
						// p[i] = p[i]+r[i]
						bt.real = 1.; 
						bt.imag = 0.;
						zaxpy ( &nn, &bt, r[i], &incx, p[i], &incx );
					}
				//	zcopy( &size, q[0], &incx, p[0], &incx );
				}

		for ( i = 0; i < m; i++)// цикл по правым частям		
		{
		
			
		

			if (eps[i] < tol) continue;
			{
				w++;
			
				
				
				/// q[i] = A * p[i]
					al.imag = 0.;
					bt.imag = 0.;
					bt.real = 0;	
					al.real = 1.;
					transa='C'; // делаем из Ah A в текущем умножении

					zgemv( &transa, &nn, &nn, &al, A[0], &lda, p[i], &incx, &bt, q[i], &incx );

				///alpha = ro2 / ||q[i]||

					zdotc(&resc, &nn, q[i], &incx, q[i], &incx);

					zladiv(&alpha,&ro2[i],&resc);
				///	alpha = ro2[i] / res / res;
				

				/// x[i] = x[i] + alpha * p[i]
					// u[i] = alpha * p[i] делаем так , чтобы в явном виде получить приращение x
						//al.real = alpha;
					// u[i] = p[i]
						zcopy( &nn, p[i], &incx, u[i], &incx );
					//u[i] = alpha * u[i]
						zscal ( &nn, &alpha, u[i], &incx );
						
					/// x[i] = x[i] + u[i]
						al.real = 1.;
						zaxpy ( &nn, &al, u[i], &incx, x[i], &incx );

				/// r[i] = r[i] - alpha * Ah * q[i]

					bt.real = 1.;	
					al.real = -alpha.real;
					al.imag = -alpha.imag;
					transa='N'; // 

					zgemv( &transa, &nn, &nn, &al, A[0], &lda, q[i], &incx, &bt, r[i], &incx );
						
		
				//
			/*		bt.real = dznrm2 ( &nn, x[i], &incx );
					res = dznrm2 ( &nn, u[i], &incx );
					eps[i] = res / bt.real;

			*/eps[i] = dznrm2 ( &nn, r[i], &incx );;
				
					ro1[i].real = ro2[i].real; // передаем текущее значение в предыдущее
					ro1[i].imag = ro2[i].imag; // передаем текущее значение в предыдущее
			}
		}
//
		res = dnrm2(&mm,eps,&incx);
		

	counter ++;

	}while ( (res > tol) );
	//while (( res > tol) && w);
	//while (w);

	fprintf(out,"\n \n \n");

	fprintf(out,"Number of iterations is  %d ", counter);

		size = n * m;
		zcopy( &size, x[0], &incx, B[0], &incx );
	
	free(eps);
	free(ro1);
	free(ro2);

	DeAllocC(beta);
//	DeAllocC(ro1);
//	DeAllocC(ro2);
	DeAllocC(q);
	DeAllocC(p);
	DeAllocC(r);
	DeAllocC(x);
	DeAllocC(u);	
//	lda = 
//	zgemm (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
	

}

