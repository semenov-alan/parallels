#pragma comment(lib, "C:\\Program Files (x86)\\Microsoft SDKs\\MPI\\Lib\\x86\\msmpi.lib")
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "C:\\Program Files (x86)\\Microsoft SDKs\\MPI\\Include\\mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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

void programmFunction(int rank, int begin_p, int esc_p, double t, double *x, double *dotX, int n)
{
    if (rank == 0)
    {
        int i = 0;
        dotX[i] = -6 * x[i] * (x[i + 1] - x[i - 1 + n]) * 0.5 * firstV - (x[i + 2] - 2 * x[i + 1] + 2 * x[i - 1 + n] - x[i - 2 + n]) * 0.5 * secondV;
        i = 1;
        dotX[i] = -6 * x[i] * (x[i + 1] - x[i - 1]) * 0.5 * firstV - (x[i + 2] - 2 * x[i + 1] + 2 * x[i - 1] - x[i - 2 + n]) * 0.5 * secondV;
        for (i = 2; i < esc_p; i++)
        {
            dotX[i] = -6 * x[i] * (x[i + 1] - x[i - 1]) * 0.5 * firstV - (x[i + 2] - 2 * x[i + 1] + 2 * x[i - 1] - x[i - 2]) * 0.5 * secondV;
        }
    }
    if (rank == 1)
    {
        int i = begin_p;
        for (i = begin_p; i < n - 2; i++)
        {
            dotX[i] = -6 * x[i] * (x[i + 1] - x[i - 1]) * 0.5 * firstV - (x[i + 2] - 2 * x[i + 1] + 2 * x[i - 1] - x[i - 2]) * 0.5 * secondV;
        }
        i = n - 2;
        dotX[i] = -6 * x[i] * (x[i + 1] - x[i - 1]) * 0.5 * firstV - (x[i + 2 - n] - 2 * x[i + 1] + 2 * x[i - 1] - x[i - 2]) * 0.5 * secondV;
        i = n - 1;
        dotX[i] = -6 * x[i] * (x[i + 1 - n] - x[i - 1]) * 0.5 * firstV - (x[i + 2 - n] - 2 * x[i + 1 - n] + 2 * x[i - 1] - x[i - 2]) * 0.5 * secondV;
    }
}

int rankRK(int rank, int begin_p, int esc_p, int n, double t, double *x, double h, double finish)
{

    if (h <= 0 || (finish - t) <= 0)
    {
        return -1;
    }

    double *var1, *var2, *var3, *var4;
    double *result;
    var1 = (double *)malloc(n * sizeof(double));
    var2 = (double *)malloc(n * sizeof(double));
    var3 = (double *)malloc(n * sizeof(double));
    var4 = (double *)malloc(n * sizeof(double));
    result = (double *)malloc(n * sizeof(double));

    while (true)
    {
        if (t > finish)
        {
            break;
        }

        if (rank == 0)
        {
            MPI_Status status;
            MPI_Recv(&(x[esc_p]), n - esc_p, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status);
            MPI_Send(&(x[esc_p - 2]), 2, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);
            MPI_Send(&(x[0]), 2, MPI_DOUBLE, 1, 2, MPI_COMM_WORLD);
            printfArray(x, n, t);
        }
        if (rank == 1)
        {
            MPI_Status status, status2;
            MPI_Send(&(x[begin_p]), n - begin_p, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Recv(&(x[begin_p - 2]), 2, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&(x[0]), 2, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status2);
        }

        int edge_n[4] = {n - 1, n - 2, esc_p, esc_p + 1};
        if (rank == 1)
        {
            edge_n[0] = begin_p - 2;
            edge_n[1] = begin_p - 1;
            edge_n[2] = 1;
            edge_n[3] = 0;
        }

        //k1
        programmFunction(rank, begin_p, esc_p, t, x, k1, n);
        //k2
        for (int i = begin_p; i < esc_p; i++)
        {
            result[i] = x[i] + 0.5 * h * k1[i];
        }
        for (int i = 0; i < 4; i++)
        {
            result[edge_n[i]] = x[edge_n[i]] + 0.5 * h * k1[edge_n[i]];
        }
        programmFunction(rank, begin_p, esc_p, t + 0.5 * h, result, k2, n);
        //k3
        for (int i = begin_p; i < esc_p; i++)
        {
            result[i] = x[i] + 0.5 * h * k2[i];
        }
        for (int i = 0; i < 4; i++)
        {
            result[edge_n[i]] = x[edge_n[i]] + 0.5 * h * k2[edge_n[i]];
        }
        programmFunction(rank, begin_p, esc_p, t + 0.5 * h, result, k3, n);
        //k4
        for (int i = begin_p; i < esc_p; i++)
        {
            result[i] = x[i] + h * k3[i];
        }
        for (int i = 0; i < 4; i++)
        {
            result[edge_n[i]] = x[edge_n[i]] + h * k3[edge_n[i]];
        }
        programmFunction(rank, begin_p, esc_p, t + h, result, k4, n);
        //res
        for (int i = begin_p; i < esc_p; i++)
        {
            x[i] += h * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6.0;
        }

        t += h;
    }

    free(var1);
    free(var2);
    free(var3);
    free(var4);
    free(temp);
    return 0;
}

int main(int argc, char *argv[])
{
    system("chcp 65001"); // Кодировка

    double from = 0, to = 5.0;
    double h = 0.001;
    double k = 2;
    double x0 = 50;

    int rank = 0, size = 2;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
        freopen("result.txt", "w", stdout);
    }

    int n = 100 * N;
    int middle_p = int(n / 2);
    int begin_p = 0, esc_p = n;
    switch (rank)
    {
    case 0:
        begin_p = 0;
        esc_p = middle_p;
        break;
    case 1:
        begin_p = middle_p;
        esc_p = n;
        break;
    default:
        MPI_Finalize();
        return 0;
    }

    double *x = (double *)malloc(n * sizeof(double));

    for (int i = begin_p; i < esc_p; i++)
    {
        double xx = i * 0.1;
        x[i] = 2.0 * k * k / (cosh(k * (xx - x0)) * cosh(k * (xx - x0)));
    }

    rankRK(rank, begin_p, esc_p, n, from, x, h, to);

    free(x);
    MPI_Finalize();
}