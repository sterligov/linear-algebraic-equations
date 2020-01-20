#include "lss_28_10.h"
#define DBL_EPSILON 1e-10

char d, e; // глобальные параметры

double check_solution(double *A, double *B, double *X, int n)
{
	int i, j;
	double sum, sumErr;
	double error, cur, temp;
	sumErr = error = 0.0;

	for (i = 0; i < n; i++)
	{
		sum = 0.0;
		for (j = 0; j < n; j++)
		{
			cur = A[i*n + j] * X[j] - error;
			temp = sum + cur;
			error = temp - sum - cur;
			sum = temp;
			//sum += A[i*n + j] * X[j];
		}
		sumErr += fabs(sum - B[i]) * fabs(sum - B[i]);
	}
	return sqrt(sumErr);
}

int get_X(double *A, double *X, int n, int diag)
{
	/* Здесь массив X используется для вычисления вектора,
	вокруг которого будет происходить отражение */
	int reflSz = n - diag;
	int i, j;
	double s = 0.0;
	double norm, error, cur, temp;

	/* т.к. векторы a_i - ||a_i||*e и a_i отличаются только первым элементом
	то можно сначала сложить их общую часть */
	error = 0.0;
	for (i = diag + 1; i < n; i++)
	{
		cur = (A[n*i + diag] * A[n*i + diag]) - error;
		temp = s + cur;
		error = temp - s - cur;
		s = temp;
		//s += (A[n*i+diag] * A[n*i + diag]); 
	}

	if (fabs(s) < DBL_EPSILON)
	{
		return -1; // не нужно отражать
	}

	norm = sqrt(A[n*diag + diag] * A[n*diag + diag] + s); // ||a_i||
	X[0] = norm; // ||a_i||*e
	for (i = diag, j = 0; i < n; i++, j++)
	{
		X[j] = A[n*i + diag] - X[j]; // a_i - ||a_i||*e
	}
	norm = sqrt(X[0] * X[0] + s); // ||a_i - ||a_i||*e||

	for (i = 0; i < reflSz; i++)
	{
		X[i] /= norm;
	}
	return 0;
}

void get_B(double *B, double *X, int diag, int n)
{
	int i;
	int reflSz = n - diag;
	double error, cur, temp;
	double sum = 0.0;
	error = 0.0;
	for (i = 0; i < reflSz; i++) // x_T*B
	{
		cur = X[i] * B[diag + i] - error;
		temp = sum + cur;
		error = temp - sum - cur;
		sum = temp;
	//	sum += X[i] * B[diag + i]; // x_T*B
	}

	for (i = 0; i < reflSz; i++)
	{
		B[diag + i] -= 2 * sum * X[i]; // B - 2*x_T*B*X
	}
}

void get_A(double *A, double *tmp, double *X, int diag, int n)
{
	int i, j;
	int reflSz = n - diag;
	double error, cur, temp;
	for (i = 0; i < reflSz; i++)
	{
		tmp[i] = 0.0;
		error = 0.0;
		for (j = 0; j < reflSz; j++) // x_T*A
		{
			cur = X[j] * A[n*(diag + j) + (diag + i)] - error;
			temp = tmp[i] + cur;
			error = temp - tmp[i] - cur;
			tmp[i] = temp;
		//	tmp[i] += X[j] * A[n*(diag + j) + (diag + i)]; 
		}
	}

	for (i = 0; i < reflSz; i++) // A - 2*x_T*A*x
	{
		for (j = 0; j < reflSz; j++)
		{
			A[n*(diag + i) + (diag + j)] -= 2 * X[i] * tmp[j]; 
		}
	}
} 

int gauss(double *A, double *B, double *X, double *tmp, int n)
{
	int i, j, r, offset;
	double cur, error, temp;
	r = 0; // кол-во ЛЗ строк
	offset = n;
	for (i = 0; i < n; i++)
	{
		if (fabs(A[n*i + i]) < DBL_EPSILON) // нулевой диагональный элемент => строка нулевая
		{
			if (fabs(B[i]) > DBL_EPSILON)
			{
				return 1; // нет решений
			}
			else
			{
				++r;
				X[n - r] = 0.0;
			}
		}
	}
	
	
	if (r != 0) // вырожденная система
	{
		if (d == TRUE)
		{
			printf("WARNING: The system is degenerated.\n");
		}
		for (i = n - 1; i >= 0; i--)
		{
			error = 0.0;
			if ( fabs(A[n*i + i]) > DBL_EPSILON ) // нулевые строки пропускаем
			{
				for (j = i + 1; j < n; j++)
				{
					cur = -(A[i*n + j] * X[j]) - error;
					temp = B[i] + cur;
					error = temp - B[i] - cur;
					B[i] = temp;
				//	B[i] -= A[i*n + j] * X[j];
				}
				X[i] = B[i] / A[i*n + i];
			}
		}
	}
	else // невырожденная система
	{
		for (i = n - 1; i >= 0; i--)
		{
			error = 0.0;
			for (j = i + 1; j < n; j++)
			{
				cur = -(A[i*n + j] * X[j]) - error;
				temp = B[i] + cur;
				error = temp - B[i] - cur;
				B[i] = temp;
				/*B[i] -= A[i*n + j] * X[j];*/
			}
			X[i] = B[i] / A[i*n + i];
		}
	}

	// Переставим элеметы вектора X в соответствии с начальной матрицей
	int *toInt = (int *)(tmp + offset);
	for (i = 0; i < toInt[0]; i++)
	{
		r = toInt[2*i + 1];
		j = toInt[2*i + 2];
		temp = X[r];
		X[r] = X[j];
		X[j] = temp;
	}
	return 0;
}

int lss_28_10(int n, double *A, double *B, double *X, double *tmp)
{
	int i, k;
	int reflSz; // размерность матрицы отражений 
	transform_matrix(n, A, tmp); // на случай переопределенной матрицы
	for (k = 0; k < n - 1; k++)
	{
		/* Массив X используем для вычисления вектора отражений */
		reflSz = n - k;
		for (i = 1; i < reflSz; i++) 
		{
			X[i] = 0;
		}

		if (d == TRUE)
		{
			printf("Iteration №%d:\nComputing of the vector X\n", k);
		}

		if (get_X(A, X, n, k) == -1)
		{
			if (d == TRUE)
			{
				printf("Reflection isn't required\n");
			}
			continue;
		}

		if (d == TRUE)
		{
			printf("Computing of the matrix A = U*A\n");
		}
		get_A(A, tmp, X, k, n);

		if (d == TRUE)
		{
			printf("Computing of the vector B = U*B\n");
		}	
		get_B(B, X, k, n); 
	}
	/* Теперь в массив X будем записывать ответ */
	if (d == TRUE)
	{
		printf("Computing of the vector of solution\n");
	}
	return gauss(A, B, X, tmp, n);
}

size_t lss_memsize_28_10(int n)
{
	size_t sz;
	sz = n*sizeof(double) + (n + 1)*sizeof(int);
	return sz;
}

int read_from_file(int *n, double **A, double **B, const char *fileName)
{
	int i;
	FILE *in;
	
	in = fopen(fileName, "r");
	if (in == NULL)
	{
		return OPEN_ERR;
	}
		
	if (fscanf(in, "%d", n) != 1)
	{
		fclose(in);
		return READ_ERR;
	}
	
	if (*n <= 0)
	{
		fclose(in);
		return N_ERR;
	}
	
	*A = *B = NULL;
	if (alloc_mem(A, (*n)*(*n)*sizeof(double)) || alloc_mem(B, (*n)*sizeof(double)))
	{
		free(*A);
		free(*B);
		fclose(in);
		return ALLOC_ERR;
	}
	
	for (i = 0; i < (*n) * (*n); i++)
	{
		if (fscanf(in, "%lf", *A + i) != 1)
		{
			free(*A);
			free(*B);
			fclose(in);
			return READ_ERR;
		}
	}


	for (i = 0; i < (*n); i++)
	{
		if (fscanf(in, "%lf", *B + i) != 1)
		{
			free(*A);
			free(*B);
			fclose(in);
			return READ_ERR;
		}
	}

	if (fclose(in))
	{
		free(*A);
		free(*B);
		return CLOSE_ERR;
	}

	return 0;
}

int write_in_file(int n, double *X, const char *fileName, char solution)
{
	int i;
	FILE *out;
	out = fopen(fileName, "w");

	if (out == NULL)
	{
		return OPEN_ERR;
	}

	if (solution == EXIST)
	{
		if (fprintf(out, "%d\n", n) < 1)
		{
			fclose(out);
			return PRINT_ERR;
		}

		for (i = 0; i < n; i++)
		{
			if (fprintf(out, "%1.9lf\n", X[i]) < 1)
			{
				fclose(out);
				return PRINT_ERR;
			}
		}
	}
	else
	{
		if (fprintf(out, "0") < 1)
		{
			fclose(out);
			return PRINT_ERR;
		}
	}
	if (fclose(out))
	{
		return CLOSE_ERR;
	}
	return 0;
}

int alloc_mem(double **tmp, size_t n)
{
	*tmp = (double *)malloc(n);
	if (*tmp == NULL)
	{
		return ALLOC_ERR;
	}
	return 0;
}

void print_matrix(double *A, int n, int m)
{
	int i, j;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{
			printf("%1.9lf ", A[i*n + j]);
		}
		printf("\n");
	}
}



void transform_matrix(int n, double *A, double *tmp)
{
	int i, col;
	int offset = n;
	int count = 0; // счетик количества перестановок
	int *toInt = (int *)(tmp + offset);
	for (i = 0; i < n; i++)
	{
		if (is_null_column(n, A, i) == TRUE)
		{
			if ((col = find_column(n, A, i)) != -1)
			{
				++count;
				toInt[count] = i;
				toInt[count + 1] = col;
				swap_column(n, A, i, col);
			}
		}
	}
	toInt[0] = count;
}

int is_null_column(int n, double *A, int col)
{
	int i;
	char isNull = TRUE;
	for (i = 0; i < n; i++)
	{
		if (fabs(A[n*i + col] - 0.0) > DBL_EPSILON)
		{
			isNull = FALSE;
			break;
		}
	}
	return isNull;
}

void swap_column(int n, double *A, int i, int j)
{
	double tmp;
	int k;
	for (k = 0; k < n; k++)
	{
		tmp = A[k*n + i];
		A[k*n + i] = A[k*n + j];
		A[k*n + j] = tmp;
	}
}

int find_column(int n, double *A, int curColumn)
{
	int i;
	for (i = n - 1; i > curColumn; i--)
	{
		if (is_null_column(n, A, i) == FALSE)
		{
			return i;
		}
	}
	return -1; // не нашли
}

int copy_matrix(int n, double *A, double *B, double **initA, double **initB)
{
	int i, err;
	err = alloc_mem(initA, n*n*sizeof(double));
	if (err != 0)
	{
		return err;
	}
	for (i = 0; i < n*n; i++)
	{
		(*initA)[i] = A[i];
	}
	err = alloc_mem(initB, n*sizeof(double));
	if (err != 0)
	{
		return err;
	}
	for (i = 0; i < n; i++)
	{
		(*initB)[i] = B[i];
	}
	return 0;
}

