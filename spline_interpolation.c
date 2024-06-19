#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#if __unix__
#include "gnuplot_i.h"
#endif

// Проверка, является ли прогоночная матрица трехдиагональной
bool isJacobiMatrix(double **a, int n)
{
    if ((fabs(a[0][0]) < fabs(a[0][1])) || (fabs(a[n-1][n-1]) < fabs(a[n-1][n-2])))
	return false;

    for (int i = 1; i < n - 1; i++)
    {
	if (fabs(a[i][i]) < (fabs(a[i][i-1]) + fabs(a[i][i+1])))
	    return false;
    }

    for (int i = 0; i < n; i++)
	if (a[i][i] == 0)
	    return false;

    return true;		      
}

// Реализация прогонки матрицы
double* through_method(double **a, double b[], int n)
{    
    double *x = (double*)calloc(n, sizeof(double));

    if(!isJacobiMatrix(a, n))
    {
	printf("Матрица не является трехдиагональной!");
	return x;
    }
    
    double A[n], B[n];
    for (int i = 0; i < n; i++)
    {
	A[i] = 0;
	B[i] = 0;
	x[i] = 0;
    }

    // Реализуем прямой и обратный ход прогонки

    A[0] = (-a[0][1]) / (a[0][0]);
    B[0] = b[0] / a[0][0]; 

    for (int i = 1; i < n-1; i++)
    {
	double ai = a[i][i-1];
	double bi = a[i][i];
	double ci = a[i][i+1];
	double di = b[i];
	double ei = ai * A[i-1] + bi;
	
	A[i] = -ci / ei;
	B[i] = (di - ai * B[i-1]) / ei;
    }

    x[n-1] = (b[n-1] - a[n-1][n-2] * B[n-2]) / (a[n-1][n-1] + a[n-1][n-2] * A[n-2]);
    
    for (int i = n - 2; i >= 0; i--)
	x[i] = A[i] * x[i+1] + B[i];
    
    return x;
}

// Вывод результатов вычислений прогонки
double * through_matrix_solution(double coeffts[], double b_vals[], int n)
{
    // Решаем СЛАУ методом прогонки
    int k = 0;

    // Задаем матрицу коэффициентов
    double ** matrix = (double**)calloc(n, sizeof(double*));
    
    for (int i = 0; i < n; i++)
	matrix[i] = (double*)calloc(n, sizeof(double));

    // Заполняем её коэффициентами
    for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++)
	    matrix[i][j] = coeffts[k++]; 

    // Получаем результаты прогонки
    double *p = through_method(matrix, b_vals, n);

    // Выводим значения коэффициентов прогонки
    printf("\nРезультаты прогонки:\n");
    for (int i = 0; i < n; i++)
	printf("%f\n", p[i]);

    for (int i = 0; i < n; i++)
	free(matrix[i]);
    
    free(matrix);

    return p;
}

// Задаем параметры вывода графиков
#define NPOINTS 20
#define SECONDS 60

// Выводим полученные графики сплайнов в gnuplot

void plot_graphs(char **tsplines)
{
    #if __unix__
    gnuplot_ctrl *h1, *h2;
    int i;

    h1 = gnuplot_init();
    h2 = gnuplot_init();
    
    gnuplot_set_axislabel(h1, "x", "X");
    gnuplot_set_axislabel(h1, "y", "Y");

    gnuplot_cmd(h1, "set key top right");
    gnuplot_cmd(h1, "set title 'Исходная функция'");

    gnuplot_cmd(h1, "set xrange[3:6]");
    gnuplot_cmd(h1, "set yrange[-1:1]");

    gnuplot_cmd(h1, "set grid xtics ytics mxtics mytics");
    gnuplot_cmd(h1, "y(x) = ((9.0)/(10.0)) * cos(((16.0)/(11.0)) * x)");    
    
    gnuplot_cmd(h1, "plot y(x) with lines title 'Точное решение'");

    gnuplot_set_axislabel(h2, "x", "X");
    gnuplot_set_axislabel(h2, "y", "Y");

    gnuplot_cmd(h2, "set key top right");
    gnuplot_cmd(h2, "set title 'Интерполяция сплайнами'");
    gnuplot_cmd(h2, "set grid xtics ytics mxtics mytics");
    gnuplot_cmd(h2, "set xrange[3:6]");
    gnuplot_cmd(h2, "set yrange[-1:1]");

    for (int i = 0; i < 5; i++)
	gnuplot_cmd(h2, tsplines[i]);	
    
    gnuplot_cmd(h2, "t(x) = ((9.0)/(10.0)) * cos(((16.0)/(11.0)) * x)");
    gnuplot_cmd(h2, "plot t(x) title 'Исходная функция' lw 2 lc rgb 'blue', \
     S1(x) title 'S1' lw 3, \
     S2(x) title 'S2' lw 3, \
     S3(x) title 'S3' lw 3, \
     S4(x) title 'S4' lw 3, \
     S5(x) title 'S5' lw 3");

    sleep(SECONDS);
    gnuplot_close(h1);
    gnuplot_close(h2);
    #endif
}

int main()
{
    // Определим сплайны разбиением функции на отрезки
    double a = 3.0, b = 6.0, n = 5.0;
    double step = (b - a) / n;
    double x[7], f[7], h[7], c_through_coeffts[36], b_through_coeffts[6],
	a_cof[6], d_cof[6], b_cof[6], c[6];
    double *x_values;

    // Создаем массив для хранения аналитически заданных сплайнов
    char **tsplines = (char**)malloc(5 * sizeof(char*));

    // Вычисляем коэффициенты для интерполяции сплайнами
    printf("Таблично заданная функция:\n");
    for (int i = 0; i < 6; i++)
    {
	// Значение на отрезке текущего шага
	x[i] = a + step * i;

	// Значение функции в точке на отрезке
	f[i] = ((9.0)/(10.0)) * cos(((16.0)/(11.0)) * x[i]);

	// Вывод таблицы значений xi и f(xi)
	printf("x[%d]: %f\n", i, x[i]);
	printf("f[%d]: %f\n\n", i, f[i]);
    }

    // Вычисляем значение шага (h0 не существует, поскольку xi-1 не существует)
    for (int i = 1; i < 7; i++)
    {
	h[i] = x[i] - x[i-1];
    }

    // с0 = вторая производная в точке а = 0
    // сn = вторая производная в точке b = 0
    c[0] = 0.0;
    c[5] = 0.0;

    // Естественные граничные условия для первой и последней строк
    c_through_coeffts[0] = 1.0;
    c_through_coeffts[35] = 1.0;
    b_through_coeffts[0] = 0.0;
    b_through_coeffts[5] = 0.0;

    // Заполняем первую и последнюю строки матрицы
    for (int i = 1; i < 7; i++)
	c_through_coeffts[i] = 0.0;

    for (int i = 0; i < 6; i++)
	c_through_coeffts[i+29] = 0.0;

    // Заполняем значащую матрицу 6x6
    int k = 1;
    for (int i = 6; i < 29; i += 7)
    {
	c_through_coeffts[i] = h[k];
	c_through_coeffts[i+1] = 2 * (h[k] + h[k+1]);
	c_through_coeffts[i+2] = h[k+1];
        c_through_coeffts[i+3] = 0.0;
	c_through_coeffts[i+4] = 0.0;
	c_through_coeffts[i+5] = 0.0;
	++k;
    }

    // Выводим прогоночные коэффициенты ci
    printf("Коэффициенты ci:");
    
    for (int i = 0; i < 36; i++)
    {
	if (i % 6 == 0)
	    printf("\n");
	printf("%f ", c_through_coeffts[i]);
    }

    // Считаем прогоночные коэффициенты для свободного столбца
    for (int i = 1; i < 5; i++)
    {
	b_through_coeffts[i] = 6 * (((f[i+1] - f[i]) / h[i+1]) - ((f[i] - f[i-1]) / h[i]));
    }

    printf("\n\n");
    x_values = through_matrix_solution(c_through_coeffts, b_through_coeffts, 6);
    
    for (int i = 0; i < 6; i++)
    {
	c[i] = x_values[i];
    }

    printf("\nКоэффициенты ai, bi, di:\n");

    // Вычисляем коэффициенты ai, bi, di
    for (int i = 1; i < 6; i++)
    {
	a_cof[i] = f[i];
	d_cof[i] = (c[i] - c[i-1]) / h[i];
	b_cof[i] = ((f[i] - f[i-1]) / h[i]) + c[i] * (h[i] / 3) + c[i-1] * (h[i] / 6);
	printf("a[%d]: %f ", i, a_cof[i]);
	printf("d[%d]: %f ", i, d_cof[i]);
	printf("b[%d]: %f ", i, b_cof[i]);
	printf("\n");
    }

    // Выводим полученные сплайны
    printf("\nПолученные сплайны:\n");
    
    for (int i = 1; i < 6; i++)
    {
	char *tspline = (char*)malloc(256 * sizeof(char));
	
	// Создаем буффер для хранения строк со сплайнами
	const int BUFFER_SIZE = 256;
	char buffer[BUFFER_SIZE];
	
	// Заполняем буффер соответствующим сплайном
	snprintf(buffer, BUFFER_SIZE, "S%d(x) = (x > %f && x < %f) ? %f + (%f)*(x - %f) + (%f/2)*(x - %f)**2 + (%f/6)*(x - %f)**3 : NaN\n", 
		 i, a + step * (i-1), a + step * i, a_cof[i], b_cof[i], x[i], c[i], x[i], d_cof[i], x[i]);

	strcpy(tspline, buffer);
	printf("%s\n", tspline);
	tsplines[i-1] = tspline;		
    }

    // Выводим график сплайнов на экран
    plot_graphs(tsplines);

    // Очищаем все сплайны
    for (int i = 0; i < 5; i++)
    {
	free(tsplines[i]);
    }
	
    free(tsplines);
    free(x_values);

    return 0;
}
