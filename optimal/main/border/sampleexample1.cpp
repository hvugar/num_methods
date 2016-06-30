#include "sampleexample1.h"

double ****a = NULL;
double ****b = NULL;
double **c  = NULL;

void check_matrix(double** a, double* b, int n)
{
    int i,j,k;
    for (i=0; i<n; i++)
    {
        int c = 0;
        for (j=0; j<n; j++)
        {
            if (a[i][j] != 0) c++;
        }
        if (c==0)
        {
            if (b[i] == 0)
                fprintf(stderr, "Equation has infinity solutions");
            else
                fprintf(stderr, "Equation has not any solution");
        }
    }
}

void gaussian_elimination(double** a, double* b, double* x, int n)
{
    int k;
    for (k=0; k<n-1; k++)
    {
        int i,j;
        for (i=(k+1); i<n; i++)
        {
            if (a[k][k] != 0.0)
            {
                double f = a[i][k] / a[k][k];
                for (j=k; j<n; j++)
                {
                    a[i][j] = a[i][j] - a[k][j] * f;
                }
                b[i] = b[i] - b[k] * f;
            }
            else
            {
               // printf("%d %f\n", k, a[k][k]);
            }
            check_matrix(a, b, n);
        }
    }

    int i;
    for (i=(n-1); i>=0; i--)
    {
        int j;
        for (j=(n-1); j>i; j--) b[i] -= (a[i][j] * x[j]);
        x[i] = b[i] / a[i][i];
    }
}

#define arguments1 \
    double a111 = x[0]; \
    double a112 = x[1]; \
    double a113 = x[2]; \
    double b111 = x[3]; \
    double b112 = x[4]; \
    double b113 = x[5]; \
    double b211 = x[6]; \
    double b212 = x[7]; \
    double b213 = x[8]; \
    double q1   = x[9]; \
    double M    = x[10];

#define arguments2 \
    double a121 = x[0]; \
    double a122 = x[1]; \
    double a123 = x[2]; \
    double b121 = x[3]; \
    double b122 = x[4]; \
    double b123 = x[5]; \
    double b221 = x[6]; \
    double b222 = x[7]; \
    double b223 = x[8]; \
    double q2   = x[9]; \
    double M    = x[10];

#define arguments3 \
    double a131 = x[0]; \
    double a132 = x[1]; \
    double a133 = x[2]; \
    double b131 = x[3]; \
    double b132 = x[4]; \
    double b133 = x[5]; \
    double b231 = x[6]; \
    double b232 = x[7]; \
    double b233 = x[8]; \
    double q3   = x[9]; \
    double M    = x[10];

double _A11(double t) { C_UNUSED(t); return 0.0; }
double _A12(double t) { C_UNUSED(t); return t; }
double _A13(double t) { C_UNUSED(t); return 2.0*t; }
double _A21(double t) { C_UNUSED(t); return 3.0*t; }
double _A22(double t) { C_UNUSED(t); return 0.0; }
double _A23(double t) { C_UNUSED(t); return -1.0; }
double _A31(double t) { C_UNUSED(t); return 1.0; }
double _A32(double t) { C_UNUSED(t); return 2.0; }
double _A33(double t) { C_UNUSED(t); return 0.0; }

double _B111(double t) { C_UNUSED(t); return 0.0; }
double _B112(double t) { C_UNUSED(t); return 1.0; }
double _B113(double t) { C_UNUSED(t); return 0.0; }
double _B121(double t) { C_UNUSED(t); return 0.0; }
double _B122(double t) { C_UNUSED(t); return 0.0; }
double _B123(double t) { C_UNUSED(t); return 0.0; }
double _B131(double t) { C_UNUSED(t); return 0.0; }
double _B132(double t) { C_UNUSED(t); return -1.0; }
double _B133(double t) { C_UNUSED(t); return 0.0; }

double _B211(double t) { C_UNUSED(t); return 0.0; }
double _B212(double t) { C_UNUSED(t); return 0.0; }
double _B213(double t) { C_UNUSED(t); return 0.0; }
double _B221(double t) { C_UNUSED(t); return 0.0; }
double _B222(double t) { C_UNUSED(t); return 0.0; }
double _B223(double t) { C_UNUSED(t); return 1.0; }
double _B231(double t) { C_UNUSED(t); return 0.0; }
double _B232(double t) { C_UNUSED(t); return 0.0; }
double _B233(double t) { C_UNUSED(t); return 0.0; }

double _C1(double t) { return -t*t*t - 6.0*t*t + 2.0*t*(cos(t)+sin(t)-1.0) + cos(t) + 3.87532; }
double _C2(double t) { return -6.0*t*t + t*(5.0-sin(t)) + sin(t) - 1.56836; }
double _C3(double t) { return -2.0*t*t - 2.0*t + 3.0*cos(t) - sin(t) + 1.12468; }

double R01(double t, double *x, unsigned int n)
{
    C_UNUSED(t); C_UNUSED(x); C_UNUSED(n);
    arguments1;
    //return (a111*a111 + a112*a112 + a113*a113) + (b111*b111 + b112*b112 + b113*b113) + (b211*b211 + b212*b212 + b213*b213) + (q1*q1);
    return a111*a111 + a112*a112 + a113*a113;
}

double S01(double t, double *x, unsigned int n)
{
    C_UNUSED(t); C_UNUSED(x); C_UNUSED(n);
    arguments1;
    double k1 = a111*_A11(t) + a112*_A21(t) + a113*_A31(t);
    double k2 = a111*_A12(t) + a112*_A22(t) + a113*_A32(t);
    double k3 = a111*_A13(t) + a112*_A23(t) + a113*_A33(t);

    //    double h1 = b111*_B111(t) + b112*_B112(t) + b113*_B113(t);
    //    double h2 = b111*_B121(t) + b112*_B122(t) + b113*_B123(t);
    //    double h3 = b111*_B131(t) + b112*_B132(t) + b113*_B133(t);

    //    double d1 = b211*_B211(t) + b212*_B212(t) + b213*_B213(t);
    //    double d2 = b211*_B221(t) + b212*_B222(t) + b213*_B223(t);
    //    double d3 = b211*_B231(t) + b212*_B232(t) + b213*_B233(t);

    double a = a111*k1 + a112*k2 + a113*k3;
    //    double a = ((a111*k1 + a112*k2 + a113*k3) + (a111*(h1+d1) + a112*(h2+d2) + a113*(h3+d3)) - q1*(a111*_C1(t) + a112*_C2(t) + a113*_C3(t)));
    double res = a / R01(t, x, n);
    return res;
}

double R02(double t, double *x, unsigned int n)
{
    C_UNUSED(t); C_UNUSED(x); C_UNUSED(n);
    arguments2;
    //return (a121*a121 + a122*a122 + a123*a123) + (b121*b121 + b122*b122 + b123*b123) + (b221*b221 + b222*b222 + b223*b223) + (q2*q2);
    return a121*a121 + a122*a122 + a123*a123;
}

double S02(double t, double *x, unsigned int n)
{
    C_UNUSED(t); C_UNUSED(x); C_UNUSED(n);
    arguments2;
    double k1 = a121*_A11(t) + a122*_A21(t) + a123*_A31(t);
    double k2 = a121*_A12(t) + a122*_A22(t) + a123*_A32(t);
    double k3 = a121*_A13(t) + a122*_A23(t) + a123*_A33(t);

    //    double h1 = b121*_B111(t) + b122*_B112(t) + b123*_B113(t);
    //    double h2 = b121*_B121(t) + b122*_B122(t) + b123*_B123(t);
    //    double h3 = b121*_B131(t) + b122*_B132(t) + b123*_B133(t);

    //    double d1 = b221*_B211(t) + b222*_B212(t) + b223*_B213(t);
    //    double d2 = b221*_B221(t) + b222*_B222(t) + b223*_B223(t);
    //    double d3 = b221*_B231(t) + b222*_B232(t) + b223*_B233(t);

    //    double res = ((a121*k1 + a122*k2 + a123*k3) + (a121*(h1+d1) + a122*(h2+d2) + a123*(h3+d3)) - q2*(a121*_C1(t) + a122*_C2(t) + a123*_C3(t))) / R02(t, x, n);
    double res = (a121*k1 + a122*k2 + a123*k3) / R02(t, x, n);
    return res;
}

double R03(double t, double *x, unsigned int n)
{
    C_UNUSED(t); C_UNUSED(x); C_UNUSED(n);
    arguments3;
    //return (a111*a111 + a112*a112 + a113*a113) + (b111*b111 + b112*b112 + b113*b113) + (b211*b211 + b212*b212 + b213*b213) + (q1*q1);
    return a131*a131 + a132*a132 + a133*a133;
}

double S03(double t, double *x, unsigned int n)
{
    C_UNUSED(t); C_UNUSED(x); C_UNUSED(n);
    arguments3;
    double k1 = a131*_A11(t) + a132*_A21(t) + a133*_A31(t);
    double k2 = a131*_A12(t) + a132*_A22(t) + a133*_A32(t);
    double k3 = a131*_A13(t) + a132*_A23(t) + a133*_A33(t);

    //    double h1 = b131*_B111(t) + b132*_B112(t) + b133*_B113(t);
    //    double h2 = b131*_B121(t) + b132*_B122(t) + b133*_B123(t);
    //    double h3 = b131*_B131(t) + b132*_B132(t) + b133*_B133(t);

    //    double d1 = b231*_B211(t) + b232*_B212(t) + b233*_B213(t);
    //    double d2 = b231*_B221(t) + b232*_B222(t) + b233*_B223(t);
    //    double d3 = b231*_B231(t) + b232*_B232(t) + b233*_B233(t);

    double a = a131*k1 + a132*k2 + a133*k3;
    //    double a = ((a111*k1 + a112*k2 + a113*k3) + (a111*(h1+d1) + a112*(h2+d2) + a113*(h3+d3)) - q1*(a111*_C1(t) + a112*_C2(t) + a113*_C3(t)));
    double res = a / R03(t, x, n);
    return res;
}

double alpha111(double t, double *x, unsigned int n) { arguments1 return S01(t, x, n) * a111 - (a111*_A11(t) + a112*_A21(t) + a113*_A31(t)); }
double alpha112(double t, double *x, unsigned int n) { arguments1 return S01(t, x, n) * a112 - (a111*_A12(t) + a112*_A22(t) + a113*_A32(t)); }
double alpha113(double t, double *x, unsigned int n) { arguments1 return S01(t, x, n) * a113 - (a111*_A13(t) + a112*_A23(t) + a113*_A33(t)); }
double betta111(double t, double *x, unsigned int n) { arguments1 return S01(t, x, n) * b111 - (a111*_B111(t) + a112*_B121(t) + a113*_B131(t)); }
double betta112(double t, double *x, unsigned int n) { arguments1 return S01(t, x, n) * b112 - (a111*_B112(t) + a112*_B122(t) + a113*_B132(t)); }
double betta113(double t, double *x, unsigned int n) { arguments1 return S01(t, x, n) * b113 - (a111*_B113(t) + a112*_B123(t) + a113*_B133(t)); }
double betta211(double t, double *x, unsigned int n) { arguments1 return S01(t, x, n) * b211 - (a111*_B211(t) + a112*_B221(t) + a113*_B231(t)); }
double betta212(double t, double *x, unsigned int n) { arguments1 return S01(t, x, n) * b212 - (a111*_B212(t) + a112*_B222(t) + a113*_B232(t)); }
double betta213(double t, double *x, unsigned int n) { arguments1 return S01(t, x, n) * b213 - (a111*_B213(t) + a112*_B223(t) + a113*_B233(t)); }
double qamma1(double t, double *x, unsigned int n)   { arguments1 return S01(t, x, n) * q1 + (a111*_C1(t) + a112*_C2(t) + a113*_C3(t)); }
double M1(double t, double *x, unsigned int n)       { arguments1 return S01(t, x, n) * M; }

double alpha121(double t, double *x, unsigned int n) { arguments2 return S02(t, x, n) * a121 - (a121*_A11(t) + a122*_A21(t) + a123*_A31(t)); }
double alpha122(double t, double *x, unsigned int n) { arguments2 return S02(t, x, n) * a122 - (a121*_A12(t) + a122*_A22(t) + a123*_A32(t)); }
double alpha123(double t, double *x, unsigned int n) { arguments2 return S02(t, x, n) * a123 - (a121*_A13(t) + a122*_A23(t) + a123*_A33(t)); }
double betta121(double t, double *x, unsigned int n) { arguments2 return S02(t, x, n) * b121 - (a121*_B111(t) + a122*_B121(t) + a123*_B131(t)); }
double betta122(double t, double *x, unsigned int n) { arguments2 return S02(t, x, n) * b122 - (a121*_B112(t) + a122*_B122(t) + a123*_B132(t)); }
double betta123(double t, double *x, unsigned int n) { arguments2 return S02(t, x, n) * b123 - (a121*_B113(t) + a122*_B123(t) + a123*_B133(t)); }
double betta221(double t, double *x, unsigned int n) { arguments2 return S02(t, x, n) * b221 - (a121*_B211(t) + a122*_B221(t) + a123*_B231(t)); }
double betta222(double t, double *x, unsigned int n) { arguments2 return S02(t, x, n) * b222 - (a121*_B212(t) + a122*_B222(t) + a123*_B232(t)); }
double betta223(double t, double *x, unsigned int n) { arguments2 return S02(t, x, n) * b223 - (a121*_B213(t) + a122*_B223(t) + a123*_B233(t)); }
double qamma2(double t, double *x, unsigned int n)   { arguments2 return S02(t, x, n) * q2 + (a121*_C1(t) + a122*_C2(t) + a123*_C3(t)); }
double M2(double t, double *x, unsigned int n)       { arguments2 return S02(t, x, n) * M; }

double alpha131(double t, double *x, unsigned int n) { arguments3 return S03(t, x, n) * a131 - (a131*_A11(t) + a132*_A21(t) + a133*_A31(t)); }
double alpha132(double t, double *x, unsigned int n) { arguments3 return S03(t, x, n) * a132 - (a131*_A12(t) + a132*_A22(t) + a133*_A32(t)); }
double alpha133(double t, double *x, unsigned int n) { arguments3 return S03(t, x, n) * a133 - (a131*_A13(t) + a132*_A23(t) + a133*_A33(t)); }
double betta131(double t, double *x, unsigned int n) { arguments3 return S03(t, x, n) * b131 - (a131*_B111(t) + a132*_B121(t) + a133*_B131(t)); }
double betta132(double t, double *x, unsigned int n) { arguments3 return S03(t, x, n) * b132 - (a131*_B112(t) + a132*_B122(t) + a133*_B132(t)); }
double betta133(double t, double *x, unsigned int n) { arguments3 return S03(t, x, n) * b133 - (a131*_B113(t) + a132*_B123(t) + a133*_B133(t)); }
double betta231(double t, double *x, unsigned int n) { arguments3 return S03(t, x, n) * b231 - (a131*_B211(t) + a132*_B221(t) + a133*_B231(t)); }
double betta232(double t, double *x, unsigned int n) { arguments3 return S03(t, x, n) * b232 - (a131*_B212(t) + a132*_B222(t) + a133*_B232(t)); }
double betta233(double t, double *x, unsigned int n) { arguments3 return S03(t, x, n) * b233 - (a131*_B213(t) + a132*_B223(t) + a133*_B233(t)); }
double qamma3(double t, double *x, unsigned int n)   { arguments3 return S03(t, x, n) * q3 + (a131*_C1(t) + a132*_C2(t) + a133*_C3(t)); }
double M3(double t, double *x, unsigned int n)       { arguments3 return S03(t, x, n) * M; }

void SampleMain()
{
    unsigned int N = 120;
    a = (double****) malloc(sizeof(double***) * 3);
    for (unsigned int k=0; k<3; k++)
    {
        a[k] = (double***) malloc(sizeof(double**) * 3);
        for (unsigned int j=0; j<3; j++)
        {
            a[k][j] = (double**) malloc(sizeof(double*) * 3);
            for (unsigned int i=0; i<3; i++)
            {
                a[k][j][i] = (double*) malloc(sizeof(double) * (N+1));
                for (unsigned int n=0; n<=N; n++) a[k][j][i][n] = 0.0;
            }
        }
    }

    b = (double****) malloc(sizeof(double***) * 2);
    for (unsigned int k=0; k<2; k++)
    {
        b[k] = (double***) malloc(sizeof(double**) * 3);
        for (unsigned int j=0; j<3; j++)
        {
            b[k][j] = (double**) malloc(sizeof(double*) * 3);
            for (unsigned int i=0; i<3; i++)
            {
                b[k][j][i] = (double*) malloc(sizeof(double) * (N+1));
                for (unsigned int n=0; n<=N; n++) b[k][j][i][n] = 0.0;
            }
        }
    }

    c = (double**) malloc(sizeof(double*) * 3);
    for (unsigned int k=0; k<3; k++)
    {
        c[k] = (double*) malloc(sizeof(double) * (N+1));
        for (unsigned int n=0; n<=N; n++) c[k][n] = 0.0;
    }


    SampleMain1();
    puts("---------------------------------------------------------------------------------------------------");
    SampleMain2();
    puts("---------------------------------------------------------------------------------------------------");
    SampleMain3();

    puts("----------------------------");
    FILE *file = stdout;//fopen("data_a.txt", "w");
        for (unsigned int k=0; k<3; k++)
        {
            for (unsigned int j=0; j<3; j++)
            {
                for (unsigned int i=0; i<3; i++)
                {
                    IPrinter::printVector(a[k][j][i], N+1, NULL, 8, 0, 0, file);
                }
            }
            puts("");
        }

    //    for (unsigned int k=0; k<2; k++)
    //    {
    //        for (unsigned int j=0; j<3; j++)
    //        {
    //            for (unsigned int i=0; i<3; i++)
    //            {
    //                IPrinter::printVector(b[k][j][i], N+1, NULL, 8, 0, 0, file);
    //            }
    //        }
    //        puts("");
    //    }
    for (unsigned int k=0; k<3; k++)
    {
        IPrinter::printVector(c[k], N+1, NULL, 8, 0, 0, file);
    }
    puts("");

    //    fclose(file);

    {
        double **A = (double**) malloc(sizeof(double*) * 15);
        for (unsigned int j=0; j<15; j++) A[j] = (double*) malloc(sizeof(double) * 15);
        double *B = (double*) malloc(sizeof(double) * 15);
        double *X = (double*) malloc(sizeof(double) * 15);

        // 0.0
        A[0][0]  = a[0][0][0][0]; A[0][1]  = a[0][0][1][0]; A[0][2]  = a[0][0][2][0];
        A[0][3]  = a[1][0][0][0]; A[0][4]  = a[1][0][1][0]; A[0][5]  = a[1][0][2][0];
        A[0][6]  = a[2][0][0][0]; A[0][7]  = a[2][0][1][0]; A[0][8]  = a[2][0][2][0];
        A[0][9]  = b[0][0][0][0]; A[0][10] = b[0][0][1][0]; A[0][11] = b[0][0][2][0];
        A[0][12] = b[1][0][0][0]; A[0][13] = b[1][0][1][0]; A[0][14] = b[1][0][2][0];

        A[1][0]  = a[0][1][0][0]; A[1][1]  = a[0][1][1][0]; A[1][2]  = a[0][1][2][0];
        A[1][3]  = a[1][1][0][0]; A[1][4]  = a[1][1][1][0]; A[1][5]  = a[1][1][2][0];
        A[1][6]  = a[2][1][0][0]; A[1][7]  = a[2][1][1][0]; A[1][8]  = a[2][1][2][0];
        A[1][9]  = b[0][1][0][0]; A[1][10] = b[0][1][1][0]; A[1][11] = b[0][1][2][0];
        A[1][12] = b[1][1][0][0]; A[1][13] = b[1][1][1][0]; A[1][14] = b[1][1][2][0];

        A[2][0]  = a[0][2][0][0]; A[2][1]  = a[0][2][1][0]; A[2][2]  = a[0][2][2][0];
        A[2][3]  = a[1][2][0][0]; A[2][4]  = a[1][2][1][0]; A[2][5]  = a[1][2][2][0];
        A[2][6]  = a[2][2][0][0]; A[2][7]  = a[2][2][1][0]; A[2][8]  = a[2][2][2][0];
        A[2][9]  = b[0][2][0][0]; A[2][10] = b[0][2][1][0]; A[2][11] = b[0][2][2][0];
        A[2][12] = b[1][2][0][0]; A[2][13] = b[1][2][1][0]; A[2][14] = b[1][2][2][0];

        // 0.25
        A[3][0]  = a[0][0][0][30]; A[3][1]  = a[0][0][1][30]; A[3][2]  = a[0][0][2][30];
        A[3][3]  = a[1][0][0][30]; A[3][4]  = a[1][0][1][30]; A[3][5]  = a[1][0][2][30];
        A[3][6]  = a[2][0][0][30]; A[3][7]  = a[2][0][1][30]; A[3][8]  = a[2][0][2][30];
        A[3][9]  = b[0][0][0][30]; A[3][10] = b[0][0][1][30]; A[3][11] = b[0][0][2][30];
        A[3][12] = b[1][0][0][30]; A[3][13] = b[1][0][1][30]; A[3][14] = b[1][0][2][30];

        A[4][0]  = a[0][1][0][30]; A[4][1]  = a[0][1][1][30]; A[4][2]  = a[0][1][2][30];
        A[4][3]  = a[1][1][0][30]; A[4][4]  = a[1][1][1][30]; A[4][5]  = a[1][1][2][30];
        A[4][6]  = a[2][1][0][30]; A[4][7]  = a[2][1][1][30]; A[4][8]  = a[2][1][2][30];
        A[4][9]  = b[0][1][0][30]; A[4][10] = b[0][1][1][30]; A[4][11] = b[0][1][2][30];
        A[4][12] = b[1][1][0][30]; A[4][13] = b[1][1][1][30]; A[4][14] = b[1][1][2][30];

        A[5][0]  = a[0][2][0][30]; A[5][1]  = a[0][2][1][30]; A[5][2]  = a[0][2][2][30];
        A[5][3]  = a[1][2][0][30]; A[5][4]  = a[1][2][1][30]; A[5][5]  = a[1][2][2][30];
        A[5][6]  = a[2][2][0][30]; A[5][7]  = a[2][2][1][30]; A[5][8]  = a[2][2][2][30];
        A[5][9]  = b[0][2][0][30]; A[5][10] = b[0][2][1][30]; A[5][11] = b[0][2][2][30];
        A[5][12] = b[1][2][0][30]; A[5][13] = b[1][2][1][30]; A[5][14] = b[1][2][2][30];

        // 0.5
        A[6][0]  = 0.0;            A[6][1]  = 0.0;            A[6][2]  = 0.0;
        A[6][3]  = a[1][0][0][60]; A[6][4]  = a[1][0][1][60]; A[6][5]  = a[1][0][2][60];
        A[6][6]  = a[2][0][0][60]; A[6][7]  = a[2][0][1][60]; A[6][8]  = a[2][0][2][60];
        A[6][9]  = b[0][0][0][60]; A[6][10] = b[0][0][1][60]; A[6][11] = b[0][0][2][60];
        A[6][12] = b[1][0][0][60]; A[6][13] = b[1][0][1][60]; A[6][14] = b[1][0][2][60];

        A[7][0]  = 0.0;            A[7][1]  = 0.0;            A[7][2]  = 0.0;
        A[7][3]  = a[1][1][0][60]; A[7][4]  = a[1][1][1][60]; A[7][5]  = a[1][1][2][60];
        A[7][6]  = a[2][1][0][60]; A[7][7]  = a[2][1][1][60]; A[7][8]  = a[2][1][2][60];
        A[7][9]  = b[0][1][0][60]; A[7][10] = b[0][1][1][60]; A[7][11] = b[0][1][2][60];
        A[7][12] = b[1][1][0][60]; A[7][13] = b[1][1][1][60]; A[7][14] = b[1][1][2][60];

        A[8][0]  = 0.0;            A[8][1]  = 0.0;            A[8][2]  = 0.0;
        A[8][3]  = a[1][2][0][60]; A[8][4]  = a[1][2][1][60]; A[8][5]  = a[1][2][2][60];
        A[8][6]  = a[2][2][0][60]; A[8][7]  = a[2][2][1][60]; A[8][8]  = a[2][2][2][60];
        A[8][9]  = b[0][2][0][60]; A[8][10] = b[0][2][1][60]; A[8][11] = b[0][2][2][60];
        A[8][12] = b[1][2][0][60]; A[8][13] = b[1][2][1][60]; A[8][14] = b[1][2][2][60];

        // 0.75
        A[9][0]  = 0.0;             A[9][1]  = 0.0;             A[9][2]   = 0.0;
        A[9][3]  = a[1][0][0][90];  A[9][4]  = a[1][0][1][90];  A[9][5]   = a[1][0][2][90];
        A[9][6]  = a[2][0][0][90];  A[9][7]  = a[2][0][1][90];  A[9][8]   = a[2][0][2][90];
        A[9][9]  = b[0][0][0][90];  A[9][10] = b[0][0][1][90];  A[9][11]  = b[0][0][2][90];
        A[9][12] = b[1][0][0][90];  A[9][13] = b[1][0][1][90];  A[9][14]  = b[1][0][2][90];

        A[10][0]  = 0.0;            A[10][1]  = 0.0;            A[10][2]  = 0.0;
        A[10][3]  = a[1][1][0][90]; A[10][4]  = a[1][1][1][90]; A[10][5]  = a[1][1][2][90];
        A[10][6]  = a[2][1][0][90]; A[10][7]  = a[2][1][1][90]; A[10][8]  = a[2][1][2][90];
        A[10][9]  = b[0][1][0][90]; A[10][10] = b[0][1][1][90]; A[10][11] = b[0][1][2][90];
        A[10][12] = b[1][1][0][90]; A[10][13] = b[1][1][1][90]; A[10][14] = b[1][1][2][90];

        A[11][0]  = 0.0;            A[11][1]  = 0.0;            A[11][2]  = 0.0;
        A[11][3]  = a[1][2][0][90]; A[11][4]  = a[1][2][1][90]; A[11][5]  = a[1][2][2][90];
        A[11][6]  = a[2][2][0][90]; A[11][7]  = a[2][2][1][90]; A[11][8]  = a[2][2][2][90];
        A[11][9]  = b[0][2][0][90]; A[11][10] = b[0][2][1][90]; A[11][11] = b[0][2][2][90];
        A[11][12] = b[1][2][0][90]; A[11][13] = b[1][2][1][90]; A[11][14] = b[1][2][2][90];

        // 1.0
        A[12][0]  = 0.0;             A[12][1]  = 0.0;             A[12][2]  = 0.0;
        A[12][3]  = 0.0;             A[12][4]  = 0.0;             A[12][5]  = 0.0;
        A[12][6]  = a[2][0][0][120]; A[12][7]  = a[2][0][1][120]; A[12][8]  = a[2][0][2][120];
        A[12][9]  = b[0][0][0][120]; A[12][10] = b[0][0][1][120]; A[12][11] = b[0][0][2][120];
        A[12][12] = b[1][0][0][120]; A[12][13] = b[1][0][1][120]; A[12][14] = b[1][0][2][120];

        A[13][0]  = 0.0;             A[13][1]  = 0.0;             A[13][2]  = 0.0;
        A[13][3]  = 0.0;             A[13][4]  = 0.0;             A[13][5]  = 0.0;
        A[13][6]  = a[2][1][0][120]; A[13][7]  = a[2][1][1][120]; A[13][8]  = a[2][1][2][120];
        A[13][9]  = b[0][1][0][120]; A[13][10] = b[0][1][1][120]; A[13][11] = b[0][1][2][120];
        A[13][12] = b[1][1][0][120]; A[13][13] = b[1][1][1][120]; A[13][14] = b[1][1][2][120];

        A[14][0]  = 0.0;             A[14][1] = 0.0;              A[14][2]  = 0.0;
        A[14][3]  = 0.0;             A[14][4] = 0.0;              A[14][5]  = 0.0;
        A[14][6]  = a[2][2][0][120]; A[14][7] = a[2][2][1][120];  A[14][8]  = a[2][2][2][120];
        A[14][9]  = b[0][2][0][120]; A[14][10] = b[0][2][1][120]; A[14][11] = b[0][2][2][120];
        A[14][12] = b[1][2][0][120]; A[14][13] = b[1][2][1][120]; A[14][14] = b[1][2][2][120];

        B[0]  = c[0][0];
        B[1]  = c[1][0];
        B[2]  = c[2][0];
        B[3]  = c[0][30];
        B[4]  = c[1][30];
        B[5]  = c[2][30];
        B[6]  = c[0][60];
        B[7]  = c[1][60];
        B[8]  = c[2][60];
        B[9]  = c[0][90];
        B[10] = c[1][90];
        B[11] = c[1][90];
        B[12] = c[0][120];
        B[13] = c[1][120];
        B[14] = c[2][120];

        for (unsigned int j=0; j<15; j++)
        {
            IPrinter::printVector(A[j], 15, NULL, 15, 0, 0, stdout);
            if (j%3==2) puts("---");
        }


        //        a[0][0] = 2.0;//-1.1342761192;
        //        a[0][1] = 3.0;//+3.1027273713;
        //        a[0][2] = 2.0;//+0.5284992092+1.0723190298;

        //        a[1][0] = 3.0;//-0.6913673895+0.9308570887;
        //        a[1][1] = 4.0;//+0.0967300723;
        //        a[1][2] = -2.0;//+0.7159989023;

        //        a[2][0] = 5.0;//+1.5716460549;
        //        a[2][1] = -3.0;//-2.5552946924-0.5587377192;
        //        a[2][2] = 2.0;//+0.0199430778;

        //        b[0] = 14.0;//+2.4074101105;
        //        b[1] = 5.0;//+1.2087028239;
        //        b[2] = 5.0;//+4.9533129323;

        //        x[0] = 1000.0;
        //        x[1] = 2000.0;
        //        x[2] = 3000.0;

        gaussian_elimination(A, B, X, 15);
        puts("-----------------");
        printf("%f %f %f\n", X[0], X[1], X[2]);
    }
}

void SampleMain1()
{
    {
        double t0 = 0.0;
        double t1 = 0.5;
        unsigned int N = 60;
        double h = (t1 - t0) / N;
        unsigned int n = 11;

        double *x10 = (double*) malloc(sizeof(double) * n);
        // alpha
        x10[0] = 2.0;
        x10[1] = 0.0;
        x10[2] = 0.0;
        // beta1
        x10[3] = 0.0;
        x10[4] = 0.0;
        x10[5] = 0.0;
        // beta2
        x10[6] = 0.0;
        x10[7] = 0.0;
        x10[8] = 0.0;
        // qamma1
        x10[9] = -1.3569663561501328033501918172533;
        // M
        x10[10] = 1.0;

        double **x11 = (double**)malloc(sizeof(double*) * n);
        for (unsigned int i=0; i<n; i++) x11[i] = (double*)malloc(sizeof(double) * N+1);

        ODE1stOrderEquationN eqs[11];

        eqs[0] = &alpha111;
        eqs[1] = &alpha112;
        eqs[2] = &alpha113;

        eqs[3] = &betta111;
        eqs[4] = &betta112;
        eqs[5] = &betta113;

        eqs[6] = &betta211;
        eqs[7] = &betta212;
        eqs[8] = &betta213;

        eqs[9] = &qamma1;
        eqs[10] = &M1;

        runge_kutta_rk4_system(t0, t1, x10, x11, n, N, h, eqs);

        double *a211 = (double*) malloc(sizeof(double) * (N+1));
        double *a212 = (double*) malloc(sizeof(double) * (N+1));
        double *a213 = (double*) malloc(sizeof(double) * (N+1));

        double *a311 = (double*) malloc(sizeof(double) * (N+1));
        double *a312 = (double*) malloc(sizeof(double) * (N+1));
        double *a313 = (double*) malloc(sizeof(double) * (N+1));

        for (unsigned int i=0; i<=N; i++)
        {
            double M = x11[10][i];
            a211[i] = M * 0.0;
            a212[i] = M * 3.0;
            a213[i] = M * 0.0;
            a311[i] = M * 0.0;
            a312[i] = M * 0.0;
            a313[i] = M * 1.0;
        }

        FILE* file = stdout;//fopen("file1.txt", "w");
        unsigned int M = 4;
        //        IPrinter::printVector(x11[0], N+1, "alpha1_11", M, 0, 0, file);
        //        IPrinter::printVector(x11[1], N+1, "alpha1_12", M, 0, 0, file);
        //        IPrinter::printVector(x11[2], N+1, "alpha1_13", M, 0, 0, file);
        //        puts("");
        //        IPrinter::printVector(x11[3], N+1, "betta1_11", M, 0, 0, file);
        //        IPrinter::printVector(x11[4], N+1, "betta1_12", M, 0, 0, file);
        //        IPrinter::printVector(x11[5], N+1, "betta1_13", M, 0, 0, file);
        //        puts("");
        //        IPrinter::printVector(x11[6], N+1, "betta2_11", M, 0, 0, file);
        //        IPrinter::printVector(x11[7], N+1, "betta2_12", M, 0, 0, file);
        //        IPrinter::printVector(x11[8], N+1, "betta2_13", M, 0, 0, file);
        //        puts("");
        IPrinter::printVector(x11[9], N+1, "qamma_1  ", M, 0, 0, file);
        //        IPrinter::printVector(x11[10], N+1, "M        ", M, 0, 0, file);
        //        puts("");
        //        IPrinter::printVector(a211, N+1, "alpha2_11", M, 0, 0, file);
        //        IPrinter::printVector(a212, N+1, "alpha2_12", M, 0, 0, file);
        //        IPrinter::printVector(a213, N+1, "alpha2_13", M, 0, 0, file);
        //        puts("");
        //        IPrinter::printVector(a311, N+1, "alpha3_11", M, 0, 0, file);
        //        IPrinter::printVector(a312, N+1, "alpha3_12", M, 0, 0, file);
        //        IPrinter::printVector(a313, N+1, "alpha3_13", M, 0, 0, file);
        //fclose(file);

        for (unsigned int i=0; i<=N; i++)
        {
            a[0][0][0][i] = x11[0][i];
            a[0][0][1][i] = x11[1][i];
            a[0][0][2][i] = x11[2][i];

            a[1][0][0][i] = a211[i];
            a[1][0][1][i] = a212[i];
            a[1][0][2][i] = a213[i];

            a[2][0][0][i] = a311[i];
            a[2][0][1][i] = a312[i];
            a[2][0][2][i] = a313[i];

            b[0][0][0][i] = x11[3][i];
            b[0][0][1][i] = x11[4][i];
            b[0][0][2][i] = x11[5][i];

            b[1][0][0][i] = x11[6][i];
            b[1][0][1][i] = x11[7][i];
            b[1][0][2][i] = x11[8][i];

            c[0][i] = x11[9][i];
        }
    }

    puts("---------------------------------------------------");

    {
        double t0 = 0.5;
        double t1 = 1.0;
        unsigned int N = 60;
        double h = (t1 - t0) / N;
        unsigned int n = 11;

        double *x10 = (double*) malloc(sizeof(double) * n);
        // alpha
        x10[0] = +1.9380720363+0.0000000000;
        x10[1] = -0.0742330301+2.7572018084;
        x10[2] = -0.4882276514+0.0000000000;
        // beta1
        x10[3] = +0.0000000000;
        x10[4] = -1.0116112226;
        x10[5] = +0.0000000000;
        // beta2
        x10[6] = +0.0000000000;
        x10[7] = +0.0000000000;
        x10[8] = +0.0188066184;
        // qamma1
        x10[9] = +2.6858630270;
        // M
        x10[10] = 1.0;

        double **x11 = (double**)malloc(sizeof(double*) * n);
        for (unsigned int i=0; i<n; i++) x11[i] = (double*)malloc(sizeof(double) * N+1);

        ODE1stOrderEquationN eqs[11];

        eqs[0] = &alpha111;
        eqs[1] = &alpha112;
        eqs[2] = &alpha113;

        eqs[3] = &betta111;
        eqs[4] = &betta112;
        eqs[5] = &betta113;

        eqs[6] = &betta211;
        eqs[7] = &betta212;
        eqs[8] = &betta213;

        eqs[9] = &qamma1;
        eqs[10] = &M1;

        runge_kutta_rk4_system(t0, t1, x10, x11, n, N, h, eqs);

        double *a311 = (double*) malloc(sizeof(double) * (N+1));
        double *a312 = (double*) malloc(sizeof(double) * (N+1));
        double *a313 = (double*) malloc(sizeof(double) * (N+1));

        for (unsigned int i=0; i<=N; i++)
        {
            double M = x11[10][i];
            a311[i] = M * 0.0;
            a312[i] = M * 0.0;
            a313[i] = M * 0.9190672695;
        }

        FILE* file = stdout;//fopen("file1.txt", "w");
        unsigned int M = 4;
        //        IPrinter::printVector(x11[0], N+1, "alpha2_11", M, 0, 0, file);
        //        IPrinter::printVector(x11[1], N+1, "alpha2_12", M, 0, 0, file);
        //        IPrinter::printVector(x11[2], N+1, "alpha2_13", M, 0, 0, file);
        //        puts("");
        //        IPrinter::printVector(x11[3], N+1, "betta1_11", M, 0, 0, file);
        //        IPrinter::printVector(x11[4], N+1, "betta1_12", M, 0, 0, file);
        //        IPrinter::printVector(x11[5], N+1, "betta1_13", M, 0, 0, file);
        //        puts("");
        //        IPrinter::printVector(x11[6], N+1, "betta2_11", M, 0, 0, file);
        //        IPrinter::printVector(x11[7], N+1, "betta2_12", M, 0, 0, file);
        //        IPrinter::printVector(x11[8], N+1, "betta2_13", M, 0, 0, file);
        //        puts("");
        IPrinter::printVector(x11[9], N+1, "qamma_1  ", M, 0, 0, file);
        //        IPrinter::printVector(x11[10], N+1, "M        ", M, 0, 0, file);
        //        puts("");
        //        IPrinter::printVector(a311, N+1, "alpha3_11", M, 0, 0, file);
        //        IPrinter::printVector(a312, N+1, "alpha3_12", M, 0, 0, file);
        //        IPrinter::printVector(a313, N+1, "alpha3_13", M, 0, 0, file);
        //fclose(file);

        for (unsigned int i=0; i<=N; i++)
        {
            a[1][0][0][60+i] = x11[0][i];
            a[1][0][1][60+i] = x11[1][i];
            a[1][0][2][60+i] = x11[2][i];

            a[2][0][0][60+i] = a311[i];
            a[2][0][1][60+i] = a312[i];
            a[2][0][2][60+i] = a313[i];

            b[0][0][0][60+i] = x11[3][i];
            b[0][0][1][60+i] = x11[4][i];
            b[0][0][2][60+i] = x11[5][i];

            b[1][0][0][60+i] = x11[6][i];
            b[1][0][1][60+i] = x11[7][i];
            b[1][0][2][60+i] = x11[8][i];

            c[0][60+i] = x11[9][i];
        }
    }
}

void SampleMain2()
{
    double t0 = 0.0;
    double t1 = 1.0;
    unsigned int N = 120;
    double h = (t1 - t0) / N;
    unsigned int n = 11;

    double *x10 = (double*) malloc(sizeof(double) * n);
    // alpha
    x10[0] = 0.0;
    x10[1] = 1.0;
    x10[2] = 0.0;
    // beta1
    x10[3] = 0.0;
    x10[4] = 0.0;
    x10[5] = 0.0;
    // beta2
    x10[6] = 0.0;
    x10[7] = 0.0;
    x10[8] = 0.0;
    // qamma1
    x10[9] = 3.6829419696157930133050046432606;
    // M
    x10[10] = 1.0;

    double **x11 = (double**)malloc(sizeof(double*) * n);
    for (unsigned int i=0; i<n; i++) x11[i] = (double*)malloc(sizeof(double) * N+1);

    ODE1stOrderEquationN eqs[11];

    eqs[0] = &alpha121;
    eqs[1] = &alpha122;
    eqs[2] = &alpha123;

    eqs[3] = &betta121;
    eqs[4] = &betta122;
    eqs[5] = &betta123;

    eqs[6] = &betta221;
    eqs[7] = &betta222;
    eqs[8] = &betta223;

    eqs[9] = &qamma2;
    eqs[10] = &M2;

    runge_kutta_rk4_system(t0, t1, x10, x11, n, N, h, eqs);

    double *a221 = (double*) malloc(sizeof(double) * (N+1));
    double *a222 = (double*) malloc(sizeof(double) * (N+1));
    double *a223 = (double*) malloc(sizeof(double) * (N+1));

    double *a321 = (double*) malloc(sizeof(double) * (N+1));
    double *a322 = (double*) malloc(sizeof(double) * (N+1));
    double *a323 = (double*) malloc(sizeof(double) * (N+1));

    for (unsigned int i=0; i<=N; i++)
    {
        double M = x11[10][i];
        a221[i] = M * 0.0;
        a222[i] = M * 0.0;
        a223[i] = M * 0.0;
        a321[i] = M * 2.0;
        a322[i] = M * 0.0;
        a323[i] = M * 0.0;
    }

    FILE* file = stdout;//fopen("file1.txt", "w");
    unsigned int K = 8;
    //    IPrinter::printVector(x11[0], N+1, "alpha1_21", K, 0, 0, file);
    //    IPrinter::printVector(x11[1], N+1, "alpha1_22", K, 0, 0, file);
    //    IPrinter::printVector(x11[2], N+1, "alpha1_23", K, 0, 0, file);
    //    puts("");
    //    IPrinter::printVector(x11[3], N+1, "betta1_21", K, 0, 0, file);
    //    IPrinter::printVector(x11[4], N+1, "betta1_22", K, 0, 0, file);
    //    IPrinter::printVector(x11[5], N+1, "betta1_23", K, 0, 0, file);
    //    puts("");
    //    IPrinter::printVector(x11[6], N+1, "betta2_21", K, 0, 0, file);
    //    IPrinter::printVector(x11[7], N+1, "betta2_22", K, 0, 0, file);
    //    IPrinter::printVector(x11[8], N+1, "betta2_23", K, 0, 0, file);
    //    puts("");
    IPrinter::printVector(x11[9], N+1, "qamma_2  ", K, 0, 0, file);
    //    IPrinter::printVector(x11[10], N+1, "M        ", K, 0, 0, file);
    //    puts("");
    //    IPrinter::printVector(a221, N+1, "alpha2_21", K, 0, 0, file);
    //    IPrinter::printVector(a222, N+1, "alpha2_22", K, 0, 0, file);
    //    IPrinter::printVector(a223, N+1, "alpha2_23", K, 0, 0, file);
    //    puts("");
    //    IPrinter::printVector(a321, N+1, "alpha3_21", K, 0, 0, file);
    //    IPrinter::printVector(a322, N+1, "alpha3_22", K, 0, 0, file);
    //    IPrinter::printVector(a323, N+1, "alpha3_23", K, 0, 0, file);
    //fclose(file);

    for (unsigned int i=0; i<=N; i++)
    {
        a[0][1][0][i] = x11[0][i];
        a[0][1][1][i] = x11[1][i];
        a[0][1][2][i] = x11[2][i];

        a[1][1][0][i] = a221[i];
        a[1][1][1][i] = a222[i];
        a[1][1][2][i] = a223[i];

        a[2][1][0][i] = a321[i];
        a[2][1][1][i] = a322[i];
        a[2][1][2][i] = a323[i];

        b[0][1][0][i] = x11[3][i];
        b[0][1][1][i] = x11[4][i];
        b[0][1][2][i] = x11[5][i];

        b[1][1][0][i] = x11[6][i];
        b[1][1][1][i] = x11[7][i];
        b[1][1][2][i] = x11[8][i];

        c[1][i] = x11[9][i];
    }
}

void SampleMain3()
{
    double t0 = 0.0;
    double t1 = 1.0;
    unsigned int N = 120;
    double h = (t1 - t0) / N;
    unsigned int n = 11;

    double *x10 = (double*) malloc(sizeof(double) * n);
    // alpha
    x10[0] = 0.0;
    x10[1] = 0.0;
    x10[2] = 3.0;
    // beta1
    x10[3] = 0.0;
    x10[4] = 0.0;
    x10[5] = 0.0;
    // beta2
    x10[6] = 0.0;
    x10[7] = 0.0;
    x10[8] = 0.0;
    // qamma1
    x10[9] = 3.080604611736279434801873214886;
    // M
    x10[10] = 1.0;

    double **x11 = (double**)malloc(sizeof(double*) * n);
    for (unsigned int i=0; i<n; i++) x11[i] = (double*)malloc(sizeof(double) * N+1);

    ODE1stOrderEquationN eqs[11];

    eqs[0] = &alpha131;
    eqs[1] = &alpha132;
    eqs[2] = &alpha133;

    eqs[3] = &betta131;
    eqs[4] = &betta132;
    eqs[5] = &betta133;

    eqs[6] = &betta231;
    eqs[7] = &betta232;
    eqs[8] = &betta233;

    eqs[9] = &qamma3;
    eqs[10] = &M3;

    runge_kutta_rk4_system(t0, t1, x10, x11, n, N, h, eqs);

    double *a231 = (double*) malloc(sizeof(double) * (N+1));
    double *a232 = (double*) malloc(sizeof(double) * (N+1));
    double *a233 = (double*) malloc(sizeof(double) * (N+1));

    double *a331 = (double*) malloc(sizeof(double) * (N+1));
    double *a332 = (double*) malloc(sizeof(double) * (N+1));
    double *a333 = (double*) malloc(sizeof(double) * (N+1));

    for (unsigned int i=0; i<=N; i++)
    {
        double M = x11[10][i];
        a231[i] = M * 0.0;
        a232[i] = M * 0.0;
        a233[i] = M * 0.0;
        a331[i] = M * 0.0;
        a332[i] = M * -1.0;
        a333[i] = M * 0.0;
    }

    FILE* file = stdout;//fopen("file1.txt", "w");
    unsigned int K = 8;
    //    IPrinter::printVector(x11[0], N+1, "alpha1_31", K, 0, 0, file);
    //    IPrinter::printVector(x11[1], N+1, "alpha1_32", K, 0, 0, file);
    //    IPrinter::printVector(x11[2], N+1, "alpha1_33", K, 0, 0, file);
    //    puts("");
    //    IPrinter::printVector(x11[3], N+1, "betta1_31", K, 0, 0, file);
    //    IPrinter::printVector(x11[4], N+1, "betta1_32", K, 0, 0, file);
    //    IPrinter::printVector(x11[5], N+1, "betta1_33", K, 0, 0, file);
    //    puts("");
    //    IPrinter::printVector(x11[6], N+1, "betta2_31", K, 0, 0, file);
    //    IPrinter::printVector(x11[7], N+1, "betta2_32", K, 0, 0, file);
    //    IPrinter::printVector(x11[8], N+1, "betta2_33", K, 0, 0, file);
    //    puts("");
    IPrinter::printVector(x11[9], N+1, "qamma_3  ", K, 0, 0, file);
    //    IPrinter::printVector(x11[10], N+1, "M        ", K, 0, 0, file);
    //    puts("");
    //    IPrinter::printVector(a231, N+1, "alpha2_31", K, 0, 0, file);
    //    IPrinter::printVector(a232, N+1, "alpha2_32", K, 0, 0, file);
    //    IPrinter::printVector(a233, N+1, "alpha2_33", K, 0, 0, file);
    //    puts("");
    //    IPrinter::printVector(a331, N+1, "alpha3_31", K, 0, 0, file);
    //    IPrinter::printVector(a332, N+1, "alpha3_32", K, 0, 0, file);
    //    IPrinter::printVector(a333, N+1, "alpha3_33", K, 0, 0, file);
    //fclose(file);

    for (unsigned int i=0; i<=N; i++)
    {
        a[0][2][0][i] = x11[0][i];
        a[0][2][1][i] = x11[1][i];
        a[0][2][2][i] = x11[2][i];

        a[1][2][0][i] = a231[i];
        a[1][2][1][i] = a232[i];
        a[1][2][2][i] = a233[i];

        a[2][2][0][i] = a331[i];
        a[2][2][1][i] = a332[i];
        a[2][2][2][i] = a333[i];

        b[0][2][0][i] = x11[3][i];
        b[0][2][1][i] = x11[4][i];
        b[0][2][2][i] = x11[5][i];

        b[1][2][0][i] = x11[6][i];
        b[1][2][1][i] = x11[7][i];
        b[1][2][2][i] = x11[8][i];

        c[2][i] = x11[9][i];
    }
}
