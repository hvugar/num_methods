#include "examplem.h"

void ExampleM::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    ExampleM e;

    DoubleMatrix rx(3, e.N+1);
    for (unsigned int k=0; k<=e.N; k++)
    {
        rx[0][k] = e.f(0,k);
        rx[1][k] = e.f(1,k);
        rx[2][k] = e.f(2,k);
    }
    IPrinter::printVector(e.w,e.p,rx.row(0));
    IPrinter::printVector(e.w,e.p,rx.row(1));
    IPrinter::printVector(e.w,e.p,rx.row(2));

    //calculateRX(rx);

//    e.calculate2R2LV1(rx);
//    e.calculate4R2LV1(rx);
//    e.calculate6R2LV1(rx);
}

ExampleM::ExampleM()
{
    h = 0.01;
    N = 100;
    w = 14;
    p = 10;
}

void ExampleM::calculateRX(DoubleMatrix &rx)
{
    rx.at(1,0) = 0.0;

    //rx.clear();
    //rx.resize(3, N+1, 0.0);

    for (unsigned int k=0; k<=N; k++)
    {
        //rx.at(0,0) = 0.0;//f(0,k);
        //rx[1][k] = 0.0;//f(1,k);
        //rx[2][k] = 0.0;//f(2,k);
    }
    IPrinter::printMatrix(w,p,rx);

//    FILE *file1 =fopen("data_rx.txt", "w");
//    IPrinter::printVector(14,10,rx,"rx0",rx.size(),0,0,file1);
//    fclose(file1);
}

//void ExampleM::calculateEta()
//{

//}

double ExampleM::a(unsigned int i, unsigned int j, unsigned int k) const
{
    if (i==0 && j==0) return +2.0;
    if (i==0 && j==1) return +3.0;
    if (i==0 && j==2) return -1.0;

    if (i==1 && j==0) return +4.0;
    if (i==1 && j==1) return +6.0;
    if (i==1 && j==2) return -2.0;

    if (i==2 && j==0) return -1.0;
    if (i==2 && j==1) return +1.0;
    if (i==2 && j==2) return -1.0;

    return NAN;
}

double ExampleM::b(unsigned int i, unsigned int k) const
{
    double t = k*h;

#ifdef SAMPLE_1
    if (i==0) return -(+2.0*sin(20.0*t*t) + 3.0*(cos(10.0*t) - sin(20.0*t)) - 1.0*(t*t*t - sin(8.0*t)*sin(8.0*t))) + (+40.0*t*cos(20.0*t*t));
    if (i==1) return -(+4.0*sin(20.0*t*t) + 6.0*(cos(10.0*t) - sin(20.0*t)) - 2.0*(t*t*t - sin(8.0*t)*sin(8.0*t))) + (-10.0*sin(10.0*t) - 20.0*cos(20.0*t));
    if (i==2) return -(-1.0*sin(20.0*t*t) + 1.0*(cos(10.0*t) - sin(20.0*t)) - 1.0*(t*t*t - sin(8.0*t)*sin(8.0*t))) + (+3.0*t*t - 16.0*cos(8.0*t)*sin(8.0*t));
#endif
#ifdef SAMPLE_2
    if (i==0) return -(+2.0*sin(t*t) + 3.0*(cos(t) - sin(t)) - 1.0*(t*t*t - sin(t)*sin(t))) + (+2.0*t*cos(t*t));
    if (i==1) return -(+4.0*sin(t*t) + 6.0*(cos(t) - sin(t)) - 2.0*(t*t*t - sin(t)*sin(t))) + (-sin(t) - cos(t));
    if (i==2) return -(-1.0*sin(t*t) + 1.0*(cos(t) - sin(t)) - 1.0*(t*t*t - sin(t)*sin(t))) + (+3.0*t*t - 2.0*cos(t)*sin(t));
#endif
    return NAN;
}

double ExampleM::f(unsigned int i, unsigned int k) const
{
    double t = k*h;
#ifdef SAMPLE_1
    if (i==0) return sin(20.0*t*t);
    if (i==1) return cos(10.0*t) - sin(20.0*t);
    if (i==2) return t*t*t - sin(8.0*t)*sin(8.0*t);
#endif
#ifdef SAMPLE_2
    if (i==0) return sin(t*t);
    if (i==1) return cos(t) - sin(t);
    if (i==2) return t*t*t - sin(t)*sin(t);
#endif
    return NAN;
}
