#ifndef LSS_H
#define LSS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define TRUE 1
#define FALSE 0
#define EXIST 0

enum
{
	NOT_EXIST = -8,
	PARAMETR_ERR,
	N_ERR,
	PRINT_ERR,
	READ_ERR,
	ALLOC_ERR,
	CLOSE_ERR,
	OPEN_ERR
};

extern char d, e;

// Функция решения
int lss_28_10(int n, double *A, double *B, double *X, double *tmp);

// Функция определяющая размер дополнительной памяти
size_t lss_memsize_28_10(int n);

// Перемножение матриц U и A
void get_A(double *A, double *tmp, double *X, int diag, int n);

// Перемножим U и B
void get_B(double *B, double *X, int diag, int n);

// Проверим корректность решения
double check_solution(double *A, double *B, double *X, int n);

// Вектор X = a_i - ||a_i||* / ||a_i - ||a_i||*e||
int get_X(double *A, double *X, int n, int diag);

// Обратный ход Гаусса
int gauss(double *A, double *B, double *X, double *tmp, int n);

// Запись в файл
int write_in_file(int n, double *X, const char *fileName, char solution);

// Считывание с файла
int read_from_file(int *n, double **A, double **B, const char *fileName);

// Выделение памяти
int alloc_mem(double **tmp, size_t n);

// Печать матрицы
void print_matrix(double *A, int n, int m);

// Перестроим матрицу к виду A = (A'|0)
void transform_matrix(int n, double *A, double *tmp);

// Проверка на нулевой столбец
int is_null_column(int n, double *A, int col);

// Меняем столбцы местами
void swap_column(int n, double *A, int i, int j);

// Найдем ненулевой столбец
int find_column(int n, double *A, int curColumn);

// Копирование матриц
int copy_matrix(int n, double *A, double *B, double **initA, double **initB);



#endif
