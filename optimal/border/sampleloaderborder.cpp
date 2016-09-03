#include "sampleloaderborder.h"

double A11(double t) { C_UNUSED(t); return 3.0; }
double A12(double t) { C_UNUSED(t); return 1.0; }
double A21(double t) { C_UNUSED(t); return 1.0; }
double A22(double t) { C_UNUSED(t); return 2.0; }

double B11(double t) { C_UNUSED(t); return 1.0; }
double B12(double t) { C_UNUSED(t); return 2.0; }
double B21(double t) { C_UNUSED(t); return 3.0; }
double B22(double t) { C_UNUSED(t); return 1.0; }

double C1(double t) { return t - 3.0*t*t - 0.44; }
double C2(double t) { return -2.0*t - t*t + 0.68; }

double R10(double t, double *x, unsigned int n)
{
    C_UNUSED(t);
    C_UNUSED(n);
    double alpha1_11 = x[0];
    double alpha1_12 = x[1];
    double beta11    = x[2];
    double beta12    = x[3];
    double qamma1    = x[4];
    double M         = x[5];
    return (alpha1_11*alpha1_11 + alpha1_12*alpha1_12) + (beta11*beta11 + beta12*beta12) + qamma1*qamma1;
}

double R20(double t, double *x, unsigned int n)
{
    C_UNUSED(t);
    C_UNUSED(n);
    double alpha1_21 = x[0];
    double alpha1_22 = x[1];
    double beta21    = x[2];
    double beta22    = x[3];
    double qamma2    = x[4];
    double M         = x[5];
    return (alpha1_21*alpha1_21 + alpha1_22*alpha1_22) + (beta21*beta21 + beta22*beta22) + qamma2*qamma2;
}

double S10(double t, double *x, unsigned int n)
{
    double alpha1_11 = x[0];
    double alpha1_12 = x[1];
    double beta11    = x[2];
    double beta12    = x[3];
    double qamma1    = x[4];
    double M         = x[5];
    return (alpha1_11*alpha1_11*A11(t) + alpha1_11*alpha1_12*(A12(t)+A21(t)) + alpha1_12*alpha1_12*A22(t))
            + (alpha1_11*(beta11*B11(t) + beta12*B21(t))+alpha1_12*(beta11*B12(t) + beta12*B22(t)))
            - (qamma1*(alpha1_11*C1(t)+alpha1_12*C2(t))) / R10(t, x, n);
}

double S20(double t, double *x, unsigned int n)
{
    double alpha1_21 = x[0];
    double alpha1_22 = x[1];
    double beta21    = x[2];
    double beta22    = x[3];
    double qamma2    = x[4];
    double M         = x[5];
    return (alpha1_21*alpha1_21*A11(t) + alpha1_21*alpha1_22*(A12(t)+A21(t)) + alpha1_22*alpha1_22*A22(t))
            + (alpha1_21*(beta21*B11(t) + beta22*B21(t))+alpha1_22*(beta21*B12(t) + beta22*B22(t)))
            - (qamma2*(alpha1_21*C1(t)+alpha1_22*C2(t))) / R20(t, x, n);
}

double Alpha1_11(double t, double *x, unsigned int n)
{
    double alpha1_11 = x[0];
    double alpha1_12 = x[1];
    double beta11    = x[2];
    double beta12    = x[3];
    double qamma1    = x[4];
    double M         = x[5];
    return S10(t, x, n) * alpha1_11 - (alpha1_11*A11(t) + alpha1_12*A21(t));
}

double Alpha1_12(double t, double *x, unsigned int n)
{
    double alpha1_11 = x[0];
    double alpha1_12 = x[1];
    double beta11    = x[2];
    double beta12    = x[3];
    double qamma1    = x[4];
    double M         = x[5];
    return S10(t, x, n) * alpha1_12 - (alpha1_11*A12(t) + alpha1_12*A22(t));
}

double Alpha1_21(double t, double *x, unsigned int n)
{
    double alpha1_21 = x[0];
    double alpha1_22 = x[1];
    double beta21    = x[2];
    double beta22    = x[3];
    double qamma2    = x[4];
    double M         = x[5];
    return S20(t, x, n) * alpha1_21 - (alpha1_21*A11(t) + alpha1_22*A21(t));
}

double Alpha1_22(double t, double *x, unsigned int n)
{
    double alpha1_21 = x[0];
    double alpha1_22 = x[1];
    double beta21    = x[2];
    double beta22    = x[3];
    double qamma2    = x[4];
    double M         = x[5];
    return S20(t, x, n) * alpha1_22 - (alpha1_21*A12(t) + alpha1_22*A22(t));
}

double Betta11(double t, double *x, unsigned int n)
{
    double alpha1_11 = x[0];
    double alpha1_12 = x[1];
    double beta11    = x[2];
    double beta12    = x[3];
    double qamma1    = x[4];
    double M         = x[5];
    return S10(t, x, n) * beta11 - (alpha1_11*B11(t) + alpha1_12*B21(t));
}

double Betta12(double t, double *x, unsigned int n)
{
    double alpha1_11 = x[0];
    double alpha1_12 = x[1];
    double beta11    = x[2];
    double beta12    = x[3];
    double qamma1    = x[4];
    double M         = x[5];
    return S10(t, x, n) * beta12 - (alpha1_11*B12(t) + alpha1_12*B22(t));
}

double Betta21(double t, double *x, unsigned int n)
{
    double alpha1_21 = x[0];
    double alpha1_22 = x[1];
    double beta21    = x[2];
    double beta22    = x[3];
    double qamma2    = x[4];
    double M         = x[5];
    return S20(t, x, n) * beta21 - (alpha1_21*B11(t) + alpha1_22*B21(t));
}

double Betta22(double t, double *x, unsigned int n)
{
    double alpha1_21 = x[0];
    double alpha1_22 = x[1];
    double beta21    = x[2];
    double beta22    = x[3];
    double qamma2    = x[4];
    double M         = x[5];
    return S20(t, x, n) * beta22 - (alpha1_21*B12(t) + alpha1_22*B22(t));
}

double Qamma1(double t, double *x, unsigned int n)
{
    double alpha1_11 = x[0];
    double alpha1_12 = x[1];
    double beta11    = x[2];
    double beta12    = x[3];
    double qamma1    = x[4];
    double M         = x[5];
    return S10(t, x, n) * qamma1 + alpha1_11*C1(t) + alpha1_12*C2(t);
}

double Qamma2(double t, double *x, unsigned int n)
{
    double alpha1_21 = x[0];
    double alpha1_22 = x[1];
    double beta21    = x[2];
    double beta22    = x[3];
    double qamma2    = x[4];
    double M         = x[5];
    return S20(t, x, n) * qamma2 + alpha1_21*C1(t) + alpha1_22*C2(t);
}

double M(double t, double *x, unsigned int n)
{
    double alpha1_11 = x[0];
    double alpha1_12 = x[1];
    double beta11    = x[2];
    double beta12    = x[3];
    double qamma1    = x[4];
    double M         = x[5];
    return S10(t, x, n) * M;
}

void SampleLoaderBorder::main()
{
    double h = 0.001;
    unsigned int N = 1000;
    unsigned int n = 6;
    double t0 = 0.0;
    double t1 = 1.0;

    {
        ODE1stOrderEquationN equations1[n];
        equations1[0] = Alpha1_11;
        equations1[1] = Alpha1_12;
        equations1[2] = Betta11;
        equations1[3] = Betta12;
        equations1[4] = Qamma1;
        equations1[5] = M;

        double *x10 = (double*)malloc(sizeof(double) * 6);
        x10[0] = 1.0;
        x10[1] = 0.0;
        x10[2] = 0.0;
        x10[3] = 0.0;
        x10[4] = 2.0;
        x10[5] = 1.0;

        double **x1 = (double **)malloc(sizeof(double*)*n);
        for (unsigned int j=0; j<n; j++) x1[j] = (double*)malloc(sizeof(double)*(N+1));
        runge_kutta_rk4_system(t0, t1, x10, x1, n, N, h, equations1);

        IPrinter::printVector(x1[0], N+1, "a1_11");
        IPrinter::printVector(x1[1], N+1, "a1_12");
        IPrinter::printVector(x1[2], N+1, "b1_11");
        IPrinter::printVector(x1[3], N+1, "b1_12");
        IPrinter::printVector(x1[4], N+1, "q_1  ");
        IPrinter::printVector(x1[5], N+1, "M    ");
    }
    puts("----");
    {
        ODE1stOrderEquationN equations2[n];
        equations2[0] = Alpha1_21;
        equations2[1] = Alpha1_22;
        equations2[2] = Betta21;
        equations2[3] = Betta22;
        equations2[4] = Qamma2;
        equations2[5] = M;

        double *x20 = (double*)malloc(sizeof(double) * 6);
        x20[0] = 1.0;
        x20[1] = 1.0;
        x20[2] = 0.0;
        x20[3] = 0.0;
        x20[4] = 1.0;
        x20[5] = 1.0;

        double **x2 = (double **)malloc(sizeof(double*)*n);
        for (unsigned int j=0; j<n; j++) x2[j] = (double*)malloc(sizeof(double)*(N+1));
        runge_kutta_rk4_system(t0, t1, x20, x2, n, N, h, equations2);

        IPrinter::printVector(x2[0], N+1, "a1_21");
        IPrinter::printVector(x2[1], N+1, "a1_22");
        IPrinter::printVector(x2[2], N+1, "b1_21");
        IPrinter::printVector(x2[3], N+1, "b1_22");
        IPrinter::printVector(x2[4], N+1, "q_2  ");
        IPrinter::printVector(x2[5], N+1, "M    ");
    }

    {
//        ODE1stOrderEquationN equations1[n];
//        equations1[0] = Alpha1_11;
//        equations1[1] = Alpha1_12;
//        equations1[2] = Alpha1_21;
//        equations1[3] = Alpha1_22;
//        equations1[4] = Betta11;
//        equations1[5] = Betta12;
//        equations1[6] = Betta21;
//        equations1[7] = Betta22;
//        equations1[8] = Qamma1;
//        equations1[9] = Qamma2;
//        equations1[10] = M;
    }
}

double RnFunctionR0::fx(const DoubleVector &x)
{
    double t         = x[0];
    double alpha1_11 = x[1];
    double alpha1_12 = x[2];
    //    double alpha1_21 = x[3];
    //    double alpha1_22 = x[4];
    double beta11    = x[3];
    double beta12    = x[4];
    //    double beta21    = x[7];
    //    double beta22    = x[8];
    double qamma1    = x[5];
    //    double qamma2    = x[6];
    double M         = x[6];
    //    double alpha2_11 = x[12];
    //    double alpha2_12 = x[13];
    //    double alpha2_21 = x[14];
    //    double alpha2_22 = x[15];
    //    double alpha3_11 = x[16];
    //    double alpha3_12 = x[17];
    //    double alpha3_21 = x[18];
    //    double alpha3_22 = x[19];
    return (alpha1_11*alpha1_11 + alpha1_12*alpha1_12) + (beta11*beta11 + beta12*beta12) + qamma1*qamma1;
}

double RnFunctionS0::fx(const DoubleVector &x)
{
    double t         = x[0];
    double alpha1_11 = x[1];
    double alpha1_12 = x[2];
    //    double alpha1_21 = x[3];
    //    double alpha1_22 = x[4];
    double beta11    = x[3];
    double beta12    = x[4];
    //    double beta21    = x[7];
    //    double beta22    = x[8];
    double qamma1    = x[5];
    //    double qamma2    = x[6];
    double M         = x[6];
    //    double alpha2_11 = x[12];
    //    double alpha2_12 = x[13];
    //    double alpha2_21 = x[14];
    //    double alpha2_22 = x[15];
    //    double alpha3_11 = x[16];
    //    double alpha3_12 = x[17];
    //    double alpha3_21 = x[18];
    //    double alpha3_22 = x[19];
    RnFunctionR0 R0;
    return (alpha1_11*alpha1_11*A11(t) + alpha1_11*alpha1_12*(A12(t)+A21(t)) + alpha1_12*alpha1_12*A22(t))
            + (alpha1_11*(B11(t)*beta11 + B12(t)*beta12)+alpha1_12*(B21(t)*beta11 + B22(t)*beta12))
            + (qamma1*(alpha1_11*C1(t)+alpha1_12*C2(t))) / R0.fx(x);
}

double RnFunctionAlpha1_11::fx(const DoubleVector &x)
{
    double t         = x[0];
    double alpha1_11 = x[1];
    double alpha1_12 = x[2];
    //    double alpha1_21 = x[3];
    //    double alpha1_22 = x[4];
    double beta11    = x[3];
    double beta12    = x[4];
    //    double beta21    = x[7];
    //    double beta22    = x[8];
    double qamma1    = x[5];
    //    double qamma2    = x[6];
    double M         = x[6];
    //    double alpha2_11 = x[12];
    //    double alpha2_12 = x[13];
    //    double alpha2_21 = x[14];
    //    double alpha2_22 = x[15];
    //    double alpha3_11 = x[16];
    //    double alpha3_12 = x[17];
    //    double alpha3_21 = x[18];
    //    double alpha3_22 = x[19];
    RnFunctionS0 S0;
    return S0.fx(x) * alpha1_11 - alpha1_11*A11(t) - alpha1_12*A12(t);
}

double RnFunctionAlpha1_12::fx(const DoubleVector &x)
{
    double t         = x[0];
    double alpha1_11 = x[1];
    double alpha1_12 = x[2];
    //    double alpha1_21 = x[3];
    //    double alpha1_22 = x[4];
    double beta11    = x[3];
    double beta12    = x[4];
    //    double beta21    = x[7];
    //    double beta22    = x[8];
    double qamma1    = x[5];
    //    double qamma2    = x[6];
    double M         = x[6];
    //    double alpha2_11 = x[12];
    //    double alpha2_12 = x[13];
    //    double alpha2_21 = x[14];
    //    double alpha2_22 = x[15];
    //    double alpha3_11 = x[16];
    //    double alpha3_12 = x[17];
    //    double alpha3_21 = x[18];
    //    double alpha3_22 = x[19];
    RnFunctionS0 S0;
    return S0.fx(x) * alpha1_12 - alpha1_11*A11(t) - alpha1_12*A12(t);
}

//double RnFunctionAlpha1_21::fx(const DoubleVector &x)
//{
//    double t         = x[0];
//    double alpha1_11 = x[1];
//    double alpha1_12 = x[2];
////    double alpha1_21 = x[3];
////    double alpha1_22 = x[4];
//    double beta11    = x[3];
//    double beta12    = x[4];
////    double beta21    = x[7];
////    double beta22    = x[8];
//    double qamma1    = x[5];
////    double qamma2    = x[6];
//    double M         = x[6];
////    double alpha2_11 = x[12];
////    double alpha2_12 = x[13];
////    double alpha2_21 = x[14];
////    double alpha2_22 = x[15];
////    double alpha3_11 = x[16];
////    double alpha3_12 = x[17];
////    double alpha3_21 = x[18];
////    double alpha3_22 = x[19];
//    RnFunctionS0 S0;
//    return S0.fx(x) * alpha1_21 - alpha1_21*A21(t) - alpha1_22*A22(t);
//}

//double RnFunctionAlpha1_22::fx(const DoubleVector &x)
//{
//    double t         = x[0];
//    double alpha1_11 = x[1];
//    double alpha1_12 = x[2];
////    double alpha1_21 = x[3];
////    double alpha1_22 = x[4];
//    double beta11    = x[3];
//    double beta12    = x[4];
////    double beta21    = x[7];
////    double beta22    = x[8];
//    double qamma1    = x[5];
////    double qamma2    = x[6];
//    double M         = x[6];
////    double alpha2_11 = x[12];
////    double alpha2_12 = x[13];
////    double alpha2_21 = x[14];
////    double alpha2_22 = x[15];
////    double alpha3_11 = x[16];
////    double alpha3_12 = x[17];
////    double alpha3_21 = x[18];
////    double alpha3_22 = x[19];
//    RnFunctionS0 S0;
//    return S0.fx(x) * alpha1_22 - alpha1_21*A21(t) - alpha1_22*A22(t);
//}

double RnFunctionBetta11::fx(const DoubleVector &x)
{
    double t         = x[0];
    double alpha1_11 = x[1];
    double alpha1_12 = x[2];
    //    double alpha1_21 = x[3];
    //    double alpha1_22 = x[4];
    double beta11    = x[3];
    double beta12    = x[4];
    //    double beta21    = x[7];
    //    double beta22    = x[8];
    double qamma1    = x[5];
    //    double qamma2    = x[6];
    double M         = x[6];
    //    double alpha2_11 = x[12];
    //    double alpha2_12 = x[13];
    //    double alpha2_21 = x[14];
    //    double alpha2_22 = x[15];
    //    double alpha3_11 = x[16];
    //    double alpha3_12 = x[17];
    //    double alpha3_21 = x[18];
    //    double alpha3_22 = x[19];
    RnFunctionS0 S0;
    return S0.fx(x) * beta11 - alpha1_11*A11(t) - alpha1_12*A12(t);
}

double RnFunctionBetta12::fx(const DoubleVector &x)
{
    double t         = x[0];
    double alpha1_11 = x[1];
    double alpha1_12 = x[2];
    //    double alpha1_21 = x[3];
    //    double alpha1_22 = x[4];
    double beta11    = x[3];
    double beta12    = x[4];
    //    double beta21    = x[7];
    //    double beta22    = x[8];
    double qamma1    = x[5];
    //    double qamma2    = x[6];
    double M         = x[6];
    //    double alpha2_11 = x[12];
    //    double alpha2_12 = x[13];
    //    double alpha2_21 = x[14];
    //    double alpha2_22 = x[15];
    //    double alpha3_11 = x[16];
    //    double alpha3_12 = x[17];
    //    double alpha3_21 = x[18];
    //    double alpha3_22 = x[19];
    RnFunctionS0 S0;
    return S0.fx(x) * beta12 - alpha1_11*A12(t) - alpha1_12*A22(t);
}

//double RnFunctionBetta21::fx(const DoubleVector &x)
//{
//    double t         = x[0];
//    double alpha1_11 = x[1];
//    double alpha1_12 = x[2];
//    //    double alpha1_21 = x[3];
//    //    double alpha1_22 = x[4];
//    double beta11    = x[3];
//    double beta12    = x[4];
//    //    double beta21    = x[7];
//    //    double beta22    = x[8];
//    double qamma1    = x[5];
//    //    double qamma2    = x[6];
//    double M         = x[6];
//    //    double alpha2_11 = x[12];
//    //    double alpha2_12 = x[13];
//    //    double alpha2_21 = x[14];
//    //    double alpha2_22 = x[15];
//    //    double alpha3_11 = x[16];
//    //    double alpha3_12 = x[17];
//    //    double alpha3_21 = x[18];
//    //    double alpha3_22 = x[19];
//    RnFunctionS0 S0;
//    return S0.fx(x) * beta21 - alpha1_21*A21(t) - alpha1_22*A22(t);
//}

//double RnFunctionBetta22::fx(const DoubleVector &x)
//{
//    double t         = x[0];
//    double alpha1_11 = x[1];
//    double alpha1_12 = x[2];
//    //    double alpha1_21 = x[3];
//    //    double alpha1_22 = x[4];
//    double beta11    = x[3];
//    double beta12    = x[4];
//    //    double beta21    = x[7];
//    //    double beta22    = x[8];
//    double qamma1    = x[5];
//    //    double qamma2    = x[6];
//    double M         = x[6];
//    //    double alpha2_11 = x[12];
//    //    double alpha2_12 = x[13];
//    //    double alpha2_21 = x[14];
//    //    double alpha2_22 = x[15];
//    //    double alpha3_11 = x[16];
//    //    double alpha3_12 = x[17];
//    //    double alpha3_21 = x[18];
//    //    double alpha3_22 = x[19];
//    RnFunctionS0 S0;
//    return S0.fx(x) * beta22 - alpha1_21*A22(t) - alpha1_22*A22(t);
//}


double RnFunctionQamma1::fx(const DoubleVector &x)
{
    double t         = x[0];
    double alpha1_11 = x[1];
    double alpha1_12 = x[2];
    //    double alpha1_21 = x[3];
    //    double alpha1_22 = x[4];
    double beta11    = x[3];
    double beta12    = x[4];
    //    double beta21    = x[7];
    //    double beta22    = x[8];
    double qamma1    = x[5];
    //    double qamma2    = x[6];
    double M         = x[6];
    //    double alpha2_11 = x[12];
    //    double alpha2_12 = x[13];
    //    double alpha2_21 = x[14];
    //    double alpha2_22 = x[15];
    //    double alpha3_11 = x[16];
    //    double alpha3_12 = x[17];
    //    double alpha3_21 = x[18];
    //    double alpha3_22 = x[19];
    RnFunctionS0 S0;
    return S0.fx(x) * qamma1 - alpha1_11*C1(t);
}

//double RnFunctionQamma2::fx(const DoubleVector &x)
//{
//    double t         = x[0];
//    double alpha1_11 = x[1];
//    double alpha1_12 = x[2];
//    //    double alpha1_21 = x[3];
//    //    double alpha1_22 = x[4];
//    double beta11    = x[3];
//    double beta12    = x[4];
//    //    double beta21    = x[7];
//    //    double beta22    = x[8];
//    double qamma1    = x[5];
//    //    double qamma2    = x[6];
//    double M         = x[6];
//    //    double alpha2_11 = x[12];
//    //    double alpha2_12 = x[13];
//    //    double alpha2_21 = x[14];
//    //    double alpha2_22 = x[15];
//    //    double alpha3_11 = x[16];
//    //    double alpha3_12 = x[17];
//    //    double alpha3_21 = x[18];
//    //    double alpha3_22 = x[19];
//    RnFunctionS0 S0;
//    return S0.fx(x) * qamma2 - alpha1_12*C2(t);
//}

double RnFunctionM::fx(const DoubleVector &x)
{
    double t         = x[0];
    double alpha1_11 = x[1];
    double alpha1_12 = x[2];
    //    double alpha1_21 = x[3];
    //    double alpha1_22 = x[4];
    double beta11    = x[3];
    double beta12    = x[4];
    //    double beta21    = x[7];
    //    double beta22    = x[8];
    double qamma1    = x[5];
    //    double qamma2    = x[6];
    double M         = x[6];
    //    double alpha2_11 = x[12];
    //    double alpha2_12 = x[13];
    //    double alpha2_21 = x[14];
    //    double alpha2_22 = x[15];
    //    double alpha3_11 = x[16];
    //    double alpha3_12 = x[17];
    //    double alpha3_21 = x[18];
    //    double alpha3_22 = x[19];
    RnFunctionS0 S0;
    return S0.fx(x) * M;
}

//double RnFunctionAlpha2_11::fx(const DoubleVector &x)
//{
//    RnFunctionS0 S0;
//    double t       = x[0];
//    double alpha1_11 = x[1];
//    double alpha1_12 = x[2];
//    double alpha1_21 = x[3];
//    double alpha1_22 = x[4];
//    double beta11  = x[5];
//    double beta12  = x[6];
//    double beta21  = x[7];
//    double beta22  = x[8];
//    double qamma1  = x[9];
//    double qamma2  = x[10];
//    double M       = x[11];
//    double alpha2_11 = x[12];
//    double alpha2_12 = x[13];
//    double alpha2_21 = x[14];
//    double alpha2_22 = x[15];
//    double alpha3_11 = x[16];
//    double alpha3_12 = x[17];
//    double alpha3_21 = x[18];
//    double alpha3_22 = x[19];
//    return M * alpha2_11;
//}

//double RnFunctionAlpha2_12::fx(const DoubleVector &x)
//{
//    double t       = x[0];
//    double alpha1_11 = x[1];
//    double alpha1_12 = x[2];
//    double alpha1_21 = x[3];
//    double alpha1_22 = x[4];
//    double beta11  = x[5];
//    double beta12  = x[6];
//    double beta21  = x[7];
//    double beta22  = x[8];
//    double qamma1  = x[9];
//    double qamma2  = x[10];
//    double M       = x[11];
//    double alpha2_11 = x[12];
//    double alpha2_12 = x[13];
//    double alpha2_21 = x[14];
//    double alpha2_22 = x[15];
//    double alpha3_11 = x[16];
//    double alpha3_12 = x[17];
//    double alpha3_21 = x[18];
//    double alpha3_22 = x[19];
//    return M * alpha2_12;
//}

//double RnFunctionAlpha2_21::fx(const DoubleVector &x)
//{
//    double t       = x[0];
//    double alpha1_11 = x[1];
//    double alpha1_12 = x[2];
//    double alpha1_21 = x[3];
//    double alpha1_22 = x[4];
//    double beta11  = x[5];
//    double beta12  = x[6];
//    double beta21  = x[7];
//    double beta22  = x[8];
//    double qamma1  = x[9];
//    double qamma2  = x[10];
//    double M       = x[11];
//    double alpha2_11 = x[12];
//    double alpha2_12 = x[13];
//    double alpha2_21 = x[14];
//    double alpha2_22 = x[15];
//    double alpha3_11 = x[16];
//    double alpha3_12 = x[17];
//    double alpha3_21 = x[18];
//    double alpha3_22 = x[19];
//    return M * alpha2_21;
//}

//double RnFunctionAlpha2_22::fx(const DoubleVector &x)
//{
//    double t       = x[0];
//    double alpha1_11 = x[1];
//    double alpha1_12 = x[2];
//    double alpha1_21 = x[3];
//    double alpha1_22 = x[4];
//    double beta11  = x[5];
//    double beta12  = x[6];
//    double beta21  = x[7];
//    double beta22  = x[8];
//    double qamma1  = x[9];
//    double qamma2  = x[10];
//    double M       = x[11];
//    double alpha2_11 = x[12];
//    double alpha2_12 = x[13];
//    double alpha2_21 = x[14];
//    double alpha2_22 = x[15];
//    double alpha3_11 = x[16];
//    double alpha3_12 = x[17];
//    double alpha3_21 = x[18];
//    double alpha3_22 = x[19];
//    return M * alpha2_22;
//}

//double RnFunctionAlpha3_11::fx(const DoubleVector &x)
//{
//    double t       = x[0];
//    double alpha1_11 = x[1];
//    double alpha1_12 = x[2];
//    double alpha1_21 = x[3];
//    double alpha1_22 = x[4];
//    double beta11  = x[5];
//    double beta12  = x[6];
//    double beta21  = x[7];
//    double beta22  = x[8];
//    double qamma1  = x[9];
//    double qamma2  = x[10];
//    double M       = x[11];
//    double alpha2_11 = x[12];
//    double alpha2_12 = x[13];
//    double alpha2_21 = x[14];
//    double alpha2_22 = x[15];
//    double alpha3_11 = x[16];
//    double alpha3_12 = x[17];
//    double alpha3_21 = x[18];
//    double alpha3_22 = x[19];
//    return M * alpha2_11;
//}

//double RnFunctionAlpha3_12::fx(const DoubleVector &x)
//{
//    double t       = x[0];
//    double alpha1_11 = x[1];
//    double alpha1_12 = x[2];
//    double alpha1_21 = x[3];
//    double alpha1_22 = x[4];
//    double beta11  = x[5];
//    double beta12  = x[6];
//    double beta21  = x[7];
//    double beta22  = x[8];
//    double qamma1  = x[9];
//    double qamma2  = x[10];
//    double M       = x[11];
//    double alpha2_11 = x[12];
//    double alpha2_12 = x[13];
//    double alpha2_21 = x[14];
//    double alpha2_22 = x[15];
//    double alpha3_11 = x[16];
//    double alpha3_12 = x[17];
//    double alpha3_21 = x[18];
//    double alpha3_22 = x[19];
//    return M * alpha3_12;
//}

//double RnFunctionAlpha3_21::fx(const DoubleVector &x)
//{
//    double t       = x[0];
//    double alpha1_11 = x[1];
//    double alpha1_12 = x[2];
//    double alpha1_21 = x[3];
//    double alpha1_22 = x[4];
//    double beta11  = x[5];
//    double beta12  = x[6];
//    double beta21  = x[7];
//    double beta22  = x[8];
//    double qamma1  = x[9];
//    double qamma2  = x[10];
//    double M       = x[11];
//    double alpha2_11 = x[12];
//    double alpha2_12 = x[13];
//    double alpha2_21 = x[14];
//    double alpha2_22 = x[15];
//    double alpha3_11 = x[16];
//    double alpha3_12 = x[17];
//    double alpha3_21 = x[18];
//    double alpha3_22 = x[19];
//    return M * alpha3_21;
//}

//double RnFunctionAlpha3_22::fx(const DoubleVector &x)
//{
//    double t       = x[0];
//    double alpha1_11 = x[1];
//    double alpha1_12 = x[2];
//    double alpha1_21 = x[3];
//    double alpha1_22 = x[4];
//    double beta11  = x[5];
//    double beta12  = x[6];
//    double beta21  = x[7];
//    double beta22  = x[8];
//    double qamma1  = x[9];
//    double qamma2  = x[10];
//    double M       = x[11];
//    double alpha2_11 = x[12];
//    double alpha2_12 = x[13];
//    double alpha2_21 = x[14];
//    double alpha2_22 = x[15];
//    double alpha3_11 = x[16];
//    double alpha3_12 = x[17];
//    double alpha3_21 = x[18];
//    double alpha3_22 = x[19];
//    return M * alpha3_21;
//}
