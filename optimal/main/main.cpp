#include "problem1.h"
#include "problem1k.h"
#include "problem1z.h"
#include "problem1x.h"
#include "problem1x1.h"
#include "problem1kz.h"
#include "problem3.h"
#include "loadedsystems.h"
#include "example1.h"

#include <QtGui/QGuiApplication>
#include <imaging.h>

#include <float.h>

//void qovma(double *a, double *b, double *c, double *d, double *x, unsigned int n, double *e, unsigned int *E, unsigned int L)
//{
//    double *p = (double*)malloc(sizeof(double)*n);
//    double *q = (double*)malloc(sizeof(double)*n);
//    double **k = (double**) malloc(sizeof(double*)*L);
//    for (unsigned int s=0; s<L; s++) k[s] = (double*)malloc(sizeof(double)*n);

//    for (unsigned int i=0; i<n; i++)
//    {
//        if (i == 0)
//        {
//            p[0] = +d[0]/b[0];
//            q[0] = -c[0]/b[0];

//            for (unsigned int s=0; s<L; s++)
//            {
//                k[s][0] = -e[s]/b[0];
//                //if (i%2==0) k[s][i] *= -1.0;
//            }
//        }
//        else if (i == (n-1))
//        {
//            p[i] = +(d[i]-a[i]*p[i-1])/(b[i]+a[i]*q[i-1]);
//            q[i] = 0.0;

//            for (unsigned int s=0; s<L; s++) k[s][i] = 0.0;
//        }
//        else
//        {
//            p[i] = +(d[i]-a[i]*p[i-1])/(b[i]+a[i]*q[i-1]);
//            q[i] = -c[i]/(b[i]+a[i]*q[i-1]);

//            for (unsigned int s=0; s<L; s++)
//            {
//                if (i<(E[s]-1))
//                    k[s][i] = -(a[i]*k[s][i-1])/(b[i]+a[i]*q[i-1]);
//                else
//                    k[s][i] = 0.0;
//            }

//            for (unsigned int s=0; s<L; s++) if (i==E[s]-1) q[i] += -(a[i]*k[s][i-1])/(b[i]+a[i]*q[i-1]);
//        }
//    }

//    const unsigned int j = (unsigned)0-1;
//    for (unsigned int i=n-1; i != j; i--)
//    {
//        if (i==(n-1))
//        {
//            x[i] = p[i];
//        }
//        else
//        {
//            x[i] = p[i] + q[i]*x[i+1];

//            for (unsigned int s=0; s<L; s++)
//            {
//                if (i<=E[s]-1)
//                {
//                    x[i] = x[i] + k[s][i]*x[E[s]];
//                }
//            }
//        }
//    }

//    free(q);
//    free(q);
//    for (unsigned int s=0; s<L; s++) free(k[s]);
//    free(k);


////    printf("%10.6f %10.6f %10.6f %10.6f\n", p[0], q[0], k[0][0], k[1][0]);
////    printf("%10.6f %10.6f %10.6f %10.6f\n", p[1], q[1], k[0][1], k[1][1]);
////    printf("%10.6f %10.6f %10.6f %10.6f\n", p[2], q[2], k[0][2], k[1][2]);
////    printf("%10.6f %10.6f %10.6f %10.6f\n", p[3], q[3], k[0][3], k[1][3]);
////    printf("%10.6f %10.6f %10.6f %10.6f\n", p[4], q[4], k[0][4], k[1][4]);
////    printf("%10.6f %10.6f %10.6f %10.6f\n", p[5], q[5], k[0][5], k[1][5]);
////    printf("%10.6f %10.6f %10.6f %10.6f\n", p[6], q[6], k[0][6], k[1][6]);
////    printf("%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6]);
//}

int main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    //srand(time(NULL));
    Problem1X1::Main(argc, argv);
    //    LoadedSystems ls;

    //    Example1 ex;
    //    ex.calculate();

//    unsigned int n = 7;
//    double a[] = {0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
//    double b[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
//    double c[] = {2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0};
//    double d[] = {10.0, 5.0, 5.0, 5.0, 5.0, 5.0, 3.0};
//    double *x = (double*)malloc(sizeof(double)*n);
//    double e[] = {3.0, 4.0};
//    unsigned int E[] = {3, 4};
//    qovma(a, b, c, d, x, n, e, E, 2);

    return 0;
}
