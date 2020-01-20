#include "lss_28_10.h"

// Обработка ошибок 
void error_processing(int err);

// Разберем командную строку
int cmd_parse(int argc, char **argv, char *p, char *t, char *h, char *c,
	char **inFileName, char **outFileName);

// Печать справки
void print_help();

// Сравнение строк
int strcmp(char *s1, char *s2); 

int main(int argc, char **argv)
{
	double *A, *B, *X, *tmp, *initA, *initB;
	double curTime;
	int n, err, tact;
	char *inFileName, *outFileName;
	char p, t, h, c, solution;
	size_t sAddMem;
	clock_t sTime;
	e = d = p = t = h = c = FALSE;
	outFileName = "lss_28_10_out.txt";
	inFileName = "lss_28_10_in.txt";

	err = cmd_parse(argc, argv, &p, &t, &h, &c, &inFileName, &outFileName);
	if (err != 0)
	{
		if (e == TRUE)
		{
	 		fprintf(stderr, "ERROR: unknown parametr %s\nUse -h for help\n", argv[err]);
		}
		return PARAMETR_ERR;
	}
	if (h == TRUE)
	{
		print_help();
		return 0;
	}
	if (d == TRUE)
	{
		printf("\nReading from file.................\n");
	}
	
	err = read_from_file(&n, &A, &B, inFileName);		
	if (err != 0)
	{
		error_processing(err);
	}
	if (c == TRUE || p == TRUE)
	{
		/* В процессе вычислений матрица A и вектор B изменятся,
		а для сравнения результата нужна исходные матрица и вектор */
		err = copy_matrix(n, A, B, &initA, &initB);	
		if (err != 0)
		{
			if (e == TRUE)
			{
				error_processing(err);
			}
			free(A);
			free(B);
			free(initA);
			free(initB);
		}
	}
	err = alloc_mem(&X, n*sizeof(double));
	if (err != 0)
	{
		error_processing(err);
	}

	sAddMem = lss_memsize_28_10(n);
	err = alloc_mem(&tmp, sAddMem);
	if (err != 0)
	{
		error_processing(err);
	}

	if (d == TRUE)
	{
		printf("Start computing.................\n");
	}
		
	sTime = clock();
	solution = lss_28_10(n, A, B, X, tmp);
	curTime = (clock() - sTime) / CLOCKS_PER_SEC;
	tact = (int)(clock() - sTime);
	
	
	if (d == TRUE)
	{
		printf("End computing.................\n");
	}
	
	if (c == TRUE && solution == EXIST)
	{
		double sumErr;
		sumErr = check_solution(initA, initB, X, n);
		printf("Total error: %1.9lf\n", sumErr);
	}
	if (d == TRUE)
	{
		printf("Writing in file.................\nDone!\n\n");
	}
		
	if (p == TRUE)
	{
		printf("The initial matrix:\n");
		print_matrix(initA, n, n);
		printf("The transformed matrix:\n");
		print_matrix(A, n, n);
	}
		
	if (t == TRUE)
	{
		printf("Execution time:\n %1.9lf sec\n %d time tacts\n", curTime, tact);
	}
		
	err = write_in_file(n, X, outFileName, solution);
	if (err != 0)
	{
		free(A);
		free(B);
		free(X);
		free(tmp);
		error_processing(err);
	}
	free(A);
	free(B);
	free(X);
	free(tmp);
	return ( (solution == EXIST) ? EXIST : NOT_EXIST );
}

void print_help()
{
	printf("\
		Format: ./lss [input_file] [output_file] [options] \n \
		Options:\n \
		-d    print debug messages\n \
		-e    print errors\n \
		-p    print matrix\n \
		-t    print execution time\n \
		-c    print total error\n \
		-h, -?    help\n \
		Standart input filename: lss_28_10_in.txt\n \
		Standart output filename: lss_28_10_out.txt\n");
}


int cmd_parse(int argc, char **argv, char *p, char *t, char *h, char *c,
	char **inFileName, char **outFileName)
{
	int errParametr = 0;
	int i = 1;
	if (argc > 1)
	{
		if (argv[1][0] != '-')
		{
			*inFileName = argv[1];
			++i;
			if (argc > 2)
			{
				if (argv[2][0] != '-')
				{
					*outFileName = argv[2];
					++i;
				}
			}
		}
	}
	for (;i < argc; i++)
	{
		if ( !strcmp(argv[i], "-h") || !strcmp(argv[i], "-?") )
			*h = TRUE;
		else if (!strcmp(argv[i], "-d"))
			d = TRUE;
		else if (!strcmp(argv[i], "-e"))
			e = TRUE;
		else if (!strcmp(argv[i], "-p"))
			*p = TRUE;
		else if (!strcmp(argv[i], "-t"))
			*t = TRUE;
		else if (!strcmp(argv[i], "-c"))
			*c = TRUE;
		else 
		{
			errParametr = i;
		}
	}
	return errParametr;
}


void error_processing(int err)
{
	if (e == TRUE)
	{
		switch (err)
		{
			case N_ERR:
			{
				fprintf(stderr, "ERROR: matrix size N <= 0 or N is too big\n");
			}
			break;

			case PRINT_ERR:
			{
				fprintf(stderr, "ERROR: at printing in file\n");
			}
			break;

			case READ_ERR:
			{
				fprintf(stderr, "ERROR: at reading from file\nCheck the number of coefficients and their value\n");
			}
			break;

			case ALLOC_ERR:
			{
				fprintf(stderr, "ERROR: in allocation memory\n");
			}
			break;

			case CLOSE_ERR:
			{
				fprintf(stderr, "ERROR: at closing file\n");
			}
			break;

			case OPEN_ERR:
			{
				fprintf(stderr, "ERROR: at opening file\n");
			}
			break;
		
			default:
				break;
		}
	}
	exit(err);
}

int strcmp(char *s1, char *s2)
{
	if (s1 == s2)
	{
		return 0;
	}
	else if (!s1 || !s2)
	{
		return 1;
	}
	while (*s1 != '\0' && *s2 != '\0')
	{
		if (*s1 != *s2)
		{
			return 1;
		}
		++s1; 
		++s2;
	}
	return ( (*s1 != *s2) ? 1 : 0);
}


