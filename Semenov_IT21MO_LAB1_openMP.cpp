
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#define firstV 10
#define secondV firstV *firstV *firstV

void outpMassive(double *d, int n, double t)
{
    printf("%.8f\n", t);
    for (int i = 0; i < n; i++)
    {
        printf("%.7f ", d[i]);
    }
    printf("\n");
}

void programmFunction(double t, double *x, double *dotX, int n)
{

    int i = 0;
    dotX[i] = -6 * x[i] * (x[i + 1] - x[i - 1 + n]) * 0.5 * firstV - (x[i + 2] - 2 * x[i + 1] + 2 * x[i - 1 + n] - x[i - 2 + n]) * 0.5 * secondV;
    i = 1;
    dotX[i] = -6 * x[i] * (x[i + 1] - x[i - 1]) * 0.5 * firstV - (x[i + 2] - 2 * x[i + 1] + 2 * x[i - 1] - x[i - 2 + n]) * 0.5 * secondV;

#pragma omp parallel for
    for (i = 2; i < n - 2; i++)
    {
        dotX[i] = -6 * x[i] * (x[i + 1] - x[i - 1]) * 0.5 * firstV - (x[i + 2] - 2 * x[i + 1] + 2 * x[i - 1] - x[i - 2]) * 0.5 * secondV;
    }
    i = n - 2;
    dotX[i] = -6 * x[i] * (x[i + 1] - x[i - 1]) * 0.5 * firstV - (x[i + 2 - n] - 2 * x[i + 1] + 2 * x[i - 1] - x[i - 2]) * 0.5 * secondV;
    i = n - 1;
    dotX[i] = -6 * x[i] * (x[i + 1 - n] - x[i - 1]) * 0.5 * firstV - (x[i + 2 - n] - 2 * x[i + 1 - n] + 2 * x[i - 1] - x[i - 2]) * 0.5 * secondV;
}

int methodRK(int n, double t, double *x, double h, double result)
{

    if (h <= 0 || (result - t) <= 0)
    {
        return -1;
    }

    double var1[500];
    double var2[500];
    double var3[500];
    double var4[500];
    double result[500];

    while (true)
    {
        if (t > result)
        {
            break;
        }

        printfArray(x, n, t);

        //var1
        programmFunction(t, x, var1, n);
//var2
#pragma omp parallel for
        for (int i = 0; i < n; i++)
        {
            result[i] = x[i] + 0.5 * h * var1[i];
        }
        programmFunction(t + 0.5 * h, result, var2, n);
//var3
#pragma omp parallel for
        for (int i = 0; i < n; i++)
        {
            result[i] = x[i] + 0.5 * h * var2[i];
        }
        programmFunction(t + 0.5 * h, result, var3, n);
//var4
#pragma omp parallel for
        for (int i = 0; i < n; i++)
        {
            result[i] = x[i] + h * var3[i];
        }
        programmFunction(t + h, result, var4, n);
//result
#pragma omp parallel for
        for (int i = 0; i < n; i++)
        {
            x[i] += h * (var1[i] + 2 * var2[i] + 2 * var3[i] + var4[i]) / 6.0;
        }

        t += h;
    }

    free(var1);
    free(var2);
    free(var3);
    free(var4);
    free(result);
    return 0;
}

int main(int argc, char *argv[])
{
    int n = 500;
    double step_time = 0.001;
    double x[500];
    double from = 0.0, to = 10.0;

    double strtPRGM = 0., endPRGM = 0.;
    omp_set_num_threads(4);
    strtPRGM = omp_get_wtime();

    double k = 1.0;
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        double var5 = i * 0.1;
        x[i] = 2.0 * k * k / (cosh(k * (var5 - 25)) * cosh(k * (var5 - 25)));
    }
    methodRK(n, from, x, step_time, to);

    endPRGM = omp_get_wtime();

    printf("\nВремя на выполнение программы: %.16g\n", endPRGM - strtPRGM);

    free(x);
    return 0;
}
