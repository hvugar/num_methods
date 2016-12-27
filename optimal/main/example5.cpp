#include "example5.h"

void Example5::Main(int argc UNUSED_PARAM, char **argv UNUSED_PARAM)
{
    Example5 e;
}

Example5::Example5()
{
    DoubleVector rx(N+1); for (unsigned int i=0; i<=N; i++) rx.at(i) = x(i);
    IPrinter::printVector(14,10,rx);
    IPrinter::printSeperatorLine();

    puts("n=2");
    calculate_n2();
    puts("n=4");
    calculate_n4();
    puts("n=6");
    calculate_n6();
}

void Example5::calculate_n2()
{
    unsigned int K = 2;
    {
        //Right
        DoubleVector nx(N+1); for (unsigned int i=0; i<K; i++) nx.at(i) = x(i);
        for (unsigned int i=K; i<=N; i++)
        {
            double m1 = 1.0 - 1.5*h*a(i) - h*h*b(i);
            double k1 = 2.0 - 2.0*h*a(i);
            double k2 = 0.5*h*a(i) - 1.0;
            double k0 = h*h*c(i);
            nx.at(i) = (k1*nx.at(i-1) + k2*nx.at(i-2) + k0)/m1;
        }
        IPrinter::printVector(14,10,nx);
        nx.clear();
    }
    {
        //Center
        DoubleVector nx(N+1); for (unsigned int i=0; i<K; i++) nx.at(i) = x(i);
        for (unsigned int i=K; i<=N; i++)
        {
            double m1 = 1.0 - 0.5*h*a(i-1);
            double k1 = 2.0 + h*h*b(i-1);
            double k2 = -0.5*h*a(i-1) - 1.0;
            double k0 = h*h*c(i-1);
            nx.at(i) = (k1*nx.at(i-1) + k2*nx.at(i-2) + k0)/m1;
        }
        IPrinter::printVector(14,10,nx);
        nx.clear();
    }
    {
        //Left
        DoubleVector nx(N+1); for (unsigned int i=N; i>N-K; i--) nx.at(i) = x(i);
        for (unsigned int i=N-K; i!=UINT32_MAX; i--)
        {
            double m1 = 1.0 + 1.5*h*a(i) - h*h*b(i);
            double k1 = 2.0 + 2.0*h*a(i);
            double k2 = -0.5*h*a(i) - 1.0;
            double k0 = h*h*c(i);
            nx.at(i) = (k1*nx.at(i+1) + k2*nx.at(i+2) + k0)/m1;
        }
        IPrinter::printVector(14,10,nx);
        nx.clear();
    }
}

void Example5::calculate_n4()
{
    unsigned int K = 4;
    {
        //Right
        DoubleVector nx(N+1); for (unsigned int i=0; i<K; i++) nx.at(i) = x(i);
        for (unsigned int i=K; i<=N; i++)
        {
            double m1 = 70.0 - 50.0*h*a(i) - 24.0*h*h*b(i);
            double k1 = 208.0 - 96.0*h*a(i);
            double k2 = 72.0*h*a(i) - 228.0;
            double k3 = 112.0 - 32.0*h*a(i);
            double k4 = 6.0*h*a(i) - 22.0;
            double k0 = 24.0*h*h*c(i);
            nx.at(i) = (k1*nx.at(i-1) + k2*nx.at(i-2) + k3*nx.at(i-3) + k4*nx.at(i-4) + k0)/m1;
        }
        IPrinter::printVector(14,10,nx);
        nx.clear();
    }
    {
        //Center
        DoubleVector nx(N+1); for (unsigned int i=0; i<K; i++) nx.at(i) = x(i);
        for (unsigned int i=0; i<=N; i++)
        {
//            if (i==0)
//            {
//                double m1 = 70.0 + 50.0*h*a(i) - 24.0*h*h*b(i);
//                double k1 = 96.0*h*a(i) + 208.0;
//                double k2 = -72.0*h*a(i) - 228.0;
//                double k3 = 112.0 + 32.0*h*a(i);
//                double k4 = -22.0 - 6.0*h*a(i);
//                double k0 = 24.0*h*h*c(i);
//                nx.at(i) = (k1*nx.at(i+1) + k2*nx.at(i+2) + k3*nx.at(i+3) + k4*nx.at(i+4) + k0)/m1;
//            }

            //            if (i==K)
            //            {
            //                double m1 = 70.0 - 50.0*h*a(i) - 24.0*h*h*b(i);
            //                double k1 = 208.0 - 96.0*h*a(i);
            //                double k2 = 72.0*h*a(i) - 228.0;
            //                double k3 = 112.0 - 32.0*h*a(i);
            //                double k4 = 6.0*h*a(i) - 22.0;
            //                double k0 = 24.0*h*h*c(i);
            //                nx.at(i) = (k1*nx.at(i-1) + k2*nx.at(i-2) + k3*nx.at(i-3) + k4*nx.at(i-4) + k0)/m1;
            //            }
            //            else if (i==K+1)
            //            {
            //                double m1 = -40.0 + 20.0*h*a(i) - 24.0*h*h*b(i);
            //                double k0 = -6.0*h*a(i)  - 22.0 ;
            //                double k2 = +36.0*h*a(i) - 12.0;
            //                double k3 = -8.0 - 12.0*h*a(i);
            //                double k4 = 2.0*h*a(i) + 2.0;
            //                double k0 = 24.0*h*h*c(i);
            //                nx.at(i) = (k1*nx.at(i-1) + k2*nx.at(i-2) + k3*nx.at(i-3) + k4*nx.at(i-4) + k0)/m1;
            //            }
            //            else
            //            {
            //                double m1 = +2.0*h*a(i-2) - 2.0;
            //                double k1 = +16.0*h*a(i-2) - 32.0;
            //                double k2 = +60.0 + 24.0*h*h*b(i-2);
            //                double k3 = -32.0 - 16.0*h*a(i-2);
            //                double k4 = +2.0*h*a(i-2) + 2.0;
            //                double k0 = +24.0*h*h*c(i-2);
            //                nx.at(i) = (k1*nx.at(i-1) + k2*nx.at(i-2) + k3*nx.at(i-3) + k4*nx.at(i-4) + k0)/m1;
            //            }
        }
        IPrinter::printVector(14,10,nx.mid(0,10));
        nx.clear();

        //        DoubleVector nx(N+1); for (unsigned int i=N; i>N-K; i--) nx.at(i) = x(i);
        //        for (unsigned int i=N-K; i!=UINT32_MAX; i--)
        //        {
        //            double m1 = -2.0*h*a(i+2) - 2.0;
        //            double k1 = -32.0 - 16.0*h*a(i+2);
        //            double k2 = 60.0 + 24.0*h*h*b(i+2);
        //            double k3 = 16.0*h*a(i+2) - 32.0;
        //            double k4 = 2.0 - 2.0*h*a(i+2);
        //            double k0 = 24.0*h*h*c(i+2);
        //            nx.at(i) = (k1*nx.at(i+1) + k2*nx.at(i+2) + k3*nx.at(i+3) + k4*nx.at(i+4) + k0)/m1;
        //        }
        //        IPrinter::printVector(14,10,nx);
        //        nx.clear();
    }
    {
        //Left
        DoubleVector nx(N+1); for (unsigned int i=N; i>N-K; i--) nx.at(i) = x(i);
        for (unsigned int i=N-K; i!=UINT32_MAX; i--)
        {
            double m1 = 70.0 + 50.0*h*a(i) - 24.0*h*h*b(i);
            double k1 = 96.0*h*a(i) + 208.0;
            double k2 = -72.0*h*a(i) - 228.0;
            double k3 = 112.0 + 32.0*h*a(i);
            double k4 = -22.0 - 6.0*h*a(i);
            double k0 = 24.0*h*h*c(i);
            nx.at(i) = (k1*nx.at(i+1) + k2*nx.at(i+2) + k3*nx.at(i+3) + k4*nx.at(i+4) + k0)/m1;
        }
        IPrinter::printVector(14,10,nx);
        nx.clear();
    }
}

void Example5::calculate_n6()
{
    unsigned int K = 6;
    {
        //Right
        DoubleVector nx(N+1); for (unsigned int i=0; i<K; i++) nx.at(i) = x(i);
        for (unsigned int i=K; i<=N; i++)
        {
            double m1 = 812.0 - 441.0*h*a(i) - 180.0*h*h*b(i);
            double k1 = -1080.0*h*a(i) + 3132.0;
            double k2 = +1350.0*h*a(i) - 5265.0;
            double k3 = -1200.0*h*a(i) + 5080.0;
            double k4 = +675.0 *h*a(i) - 2970.0;
            double k5 = -216.0 *h*a(i) + 972.0;
            double k6 = +30.0 * h*a(i) - 137.0;
            double k0 = 180.0*h*h*c(i);
            nx.at(i) = (k1*nx.at(i-1) + k2*nx.at(i-2) + k3*nx.at(i-3) + k4*nx.at(i-4) + k5*nx.at(i-5) + k6*nx.at(i-6) + k0)/m1;
        }
        IPrinter::printVector(14,10,nx);
        nx.clear();
    }
    {
        //Center
        //        DoubleVector nx(N+1); for (unsigned int i=0; i<K; i++) nx.at(i) = x(i);
        //        for (unsigned int i=K; i<=N; i++)
        //        {
        //            double m1 = 2.0*h*a(i-2) - 2.0;
        //            double k1 = 16.0*h*a(i-2) - 32.0;
        //            double k2 = 60.0 + 24.0*h*h*b(i-2);
        //            double k3 = -32.0 - 16.0*h*a(i-2);
        //            double k4 = 2.0*h*a(i-2) + 2.0;
        //            double k0 = 24.0*h*h*c(i-2);
        //            nx.at(i) = (k1*nx.at(i-1) + k2*nx.at(i-2) + k3*nx.at(i-3) + k4*nx.at(i-4) + k0)/m1;
        //        }
        //        IPrinter::printVector(14,10,nx);
        //        nx.clear();

        //        DoubleVector nx(N+1); for (unsigned int i=N; i>N-K; i--) nx.at(i) = x(i);
        //        for (unsigned int i=N-K; i!=UINT32_MAX; i--)
        //        {
        //            double m1 = -2.0*h*a(i+2) - 2.0;
        //            double k1 = -32.0 - 16.0*h*a(i+2);
        //            double k2 = 60.0 + 24.0*h*h*b(i+2);
        //            double k3 = 16.0*h*a(i+2) - 32.0;
        //            double k4 = 2.0 - 2.0*h*a(i+2);
        //            double k0 = 24.0*h*h*c(i+2);
        //            nx.at(i) = (k1*nx.at(i+1) + k2*nx.at(i+2) + k3*nx.at(i+3) + k4*nx.at(i+4) + k0)/m1;
        //        }
        //IPrinter::printVector(14,10,nx);
        //        nx.clear();
    }
    {
        //Left
        DoubleVector nx(N+1); for (unsigned int i=N; i>N-K; i--) nx.at(i) = x(i);
        for (unsigned int i=N-K; i!=UINT32_MAX; i--)
        {
            double m1 = 812.0 + 441.0*h*a(i) - 180.0*h*h*b(i);
            double k1 = +1080.0*h*a(i) + 3132.0;
            double k2 = -1350.0*h*a(i) - 5265.0;
            double k3 = +1200.0*h*a(i) + 5080.0;
            double k4 = -675.0 *h*a(i) - 2970.0;
            double k5 = +216.0 *h*a(i) + 972.0;
            double k6 = -30.0 * h*a(i) - 137.0;
            double k0 = 180.0*h*h*c(i);
            nx.at(i) = (k1*nx.at(i+1) + k2*nx.at(i+2) + k3*nx.at(i+3) + k4*nx.at(i+4) + k5*nx.at(i+5) + k6*nx.at(i+6) + k0)/m1;
        }
        IPrinter::printVector(14,10,nx);
        nx.clear();
    }
}

double Example5::a(unsigned int i) const
{
    double t UNUSED_PARAM = i*h;
#ifdef SAMPLE_1
    return t;
#endif
#ifdef SAMPLE_2
    return 2.0;
#endif
}

double Example5::b(unsigned int i) const
{
    double t UNUSED_PARAM = i*h;
#ifdef SAMPLE_1
    return -3.0;
#endif
#ifdef SAMPLE_2
    return -1.0;
#endif
}

double Example5::c(unsigned int i) const
{
    double t UNUSED_PARAM = i*h;
#ifdef SAMPLE_1
    return t*t + 2.0;
#endif
#ifdef SAMPLE_2
    return exp(2.0*t) - 2.0*cos(t);
#endif
}

double Example5::x(unsigned int i) const
{
    double t UNUSED_PARAM = i*h;
#ifdef SAMPLE_1
    return t*t;
#endif
#ifdef SAMPLE_2
    return exp(2.0*t) + sin(t);
#endif
}
