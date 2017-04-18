#include "example6.h"

void Example6::Main(int argc UNUSED_PARAM, char *arg[] UNUSED_PARAM)
{
    Example6 e;
    e.method1K2();
}

Example6::Example6()
{
}

double Example6::a(unsigned int i, unsigned int j, unsigned int k UNUSED_PARAM) const
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

double Example6::b(unsigned int i, unsigned int k) const
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

double Example6::fx(unsigned int i, unsigned int k) const
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

void Example6::method1K2()
{
    std::vector<unsigned int> s;
    s.push_back(0);
    //s.push_back(25000);
    //s.push_back(N/2);
    s.push_back(N-1);
    s.push_back(N);

    DoubleMatrix rx;
    rx.resize(n,N+1);
    for (unsigned int k=0; k<=N; k++)
    {
        rx.at(0,k) = fx(0,k);
        rx.at(1,k) = fx(1,k);
        rx.at(2,k) = fx(2,k);
    }

    /* initializing betta */
    std::vector<DoubleMatrix> betta;
    betta.resize(N+1);
    for (unsigned int i=0; i<=N; i++)
    {
        betta[i].resize(n,n,0.0);
    }

    DoubleVector eta;
    eta.resize(n,0.0);

    for (unsigned int i=0; i<s.size(); i++)
    {
        betta.at(s[i]).at(0,0) = (rand()%100)*0.1-5.0; betta.at(s[i]).at(0,1) = (rand()%100)*0.1-5.0; betta.at(s[i]).at(0,2) = (rand()%100)*0.1-5.0;
        betta.at(s[i]).at(1,0) = (rand()%100)*0.1-5.0; betta.at(s[i]).at(1,1) = (rand()%100)*0.1-5.0; betta.at(s[i]).at(1,2) = (rand()%100)*0.1-5.0;
        betta.at(s[i]).at(2,0) = (rand()%100)*0.1-5.0; betta.at(s[i]).at(2,1) = (rand()%100)*0.1-5.0; betta.at(s[i]).at(2,2) = (rand()%100)*0.1-5.0;

        IPrinter::print(betta.at(s[i]),n,n);
        IPrinter::printSeperatorLine();

        eta.at(0) = eta.at(0) + betta.at(s[i]).at(0,0)*rx.at(0,s[i]) + betta.at(s[i]).at(0,1)*rx.at(1,s[i]) + betta.at(s[i]).at(0,2)*rx.at(2,s[i]);
        eta.at(1) = eta.at(1) + betta.at(s[i]).at(1,0)*rx.at(0,s[i]) + betta.at(s[i]).at(1,1)*rx.at(1,s[i]) + betta.at(s[i]).at(1,2)*rx.at(2,s[i]);
        eta.at(2) = eta.at(2) + betta.at(s[i]).at(2,0)*rx.at(0,s[i]) + betta.at(s[i]).at(2,1)*rx.at(1,s[i]) + betta.at(s[i]).at(2,2)*rx.at(2,s[i]);
    }

    for (unsigned int k=0; k<=N-2; k++)
    {
        DoubleMatrix d0 = betta[k+0];
        DoubleMatrix d1 = betta[k+1];
        DoubleMatrix d2 = betta[k+2];

        /************************************************** betta k+1 **************************************************/
        /* first row */
        betta[k+1].at(0,0) = -2.0*h*(d0.at(0,0)*a(0,0,k+1) + d0.at(0,1)*a(1,0,k+1) + d0.at(0,2)*a(2,0,k+1)) + d1.at(0,0);
        betta[k+1].at(0,1) = -2.0*h*(d0.at(0,0)*a(0,1,k+1) + d0.at(0,1)*a(1,1,k+1) + d0.at(0,2)*a(2,1,k+1)) + d1.at(0,1);
        betta[k+1].at(0,2) = -2.0*h*(d0.at(0,0)*a(0,2,k+1) + d0.at(0,1)*a(1,2,k+1) + d0.at(0,2)*a(2,2,k+1)) + d1.at(0,2);

        /* second row */
        betta[k+1].at(1,0) = -2.0*h*(d0.at(1,0)*a(0,0,k+1) + d0.at(1,1)*a(1,0,k+1) + d0.at(1,2)*a(2,0,k+1)) + d1.at(1,0);
        betta[k+1].at(1,1) = -2.0*h*(d0.at(1,0)*a(0,1,k+1) + d0.at(1,1)*a(1,1,k+1) + d0.at(1,2)*a(2,1,k+1)) + d1.at(1,1);
        betta[k+1].at(1,2) = -2.0*h*(d0.at(1,0)*a(0,2,k+1) + d0.at(1,1)*a(1,2,k+1) + d0.at(1,2)*a(2,2,k+1)) + d1.at(1,2);

        /* third row */
        betta[k+1].at(2,0) = -2.0*h*(d0.at(2,0)*a(0,0,k+1) + d0.at(2,1)*a(1,0,k+1) + d0.at(2,2)*a(2,0,k+1)) + d1.at(2,0);
        betta[k+1].at(2,1) = -2.0*h*(d0.at(2,0)*a(0,1,k+1) + d0.at(2,1)*a(1,1,k+1) + d0.at(2,2)*a(2,1,k+1)) + d1.at(2,1);
        betta[k+1].at(2,2) = -2.0*h*(d0.at(2,0)*a(0,2,k+1) + d0.at(2,1)*a(1,2,k+1) + d0.at(2,2)*a(2,2,k+1)) + d1.at(2,2);

        /************************************************** betta k+2 **************************************************/
        /* first row */
        betta[k+2].at(0,0) = d0.at(0,0) + d2.at(0,0);
        betta[k+2].at(0,1) = d0.at(0,1) + d2.at(0,1);
        betta[k+2].at(0,2) = d0.at(0,2) + d2.at(0,2);

        /* second row */
        betta[k+2].at(1,0) = d0.at(1,0) + d2.at(1,0);
        betta[k+2].at(1,1) = d0.at(1,1) + d2.at(1,1);
        betta[k+2].at(1,2) = d0.at(1,2) + d2.at(1,2);

        /* third row */
        betta[k+2].at(2,0) = d0.at(2,0) + d2.at(2,0);
        betta[k+2].at(2,1) = d0.at(2,1) + d2.at(2,1);
        betta[k+2].at(2,2) = d0.at(2,2) + d2.at(2,2);

        /***************************************************** eta ****************************************************/
        eta.at(0) = eta.at(0) + 2.0*h*(d0.at(0,0)*b(0,k+1) + d0.at(0,1)*b(1,k+1) + d0.at(0,2)*b(2,k+1));
        eta.at(1) = eta.at(1) + 2.0*h*(d0.at(1,0)*b(0,k+1) + d0.at(1,1)*b(1,k+1) + d0.at(1,2)*b(2,k+1));
        eta.at(2) = eta.at(2) + 2.0*h*(d0.at(2,0)*b(0,k+1) + d0.at(2,1)*b(1,k+1) + d0.at(2,2)*b(2,k+1));
    }

    //    for (unsigned int i=0; i<=N; i++)
    //    {
    //        IPrinter::print(betta.at(i),n,n);
    //        IPrinter::printSeperatorLine();
    //    }

    //    DoubleMatrix M(n*n, n*n);
    //    DoubleVector c(n*n);

    //    M[0][0] = betta.at(N-2).at(0,0); M[0][1] = betta.at(N-2).at(0,1); M[0][2] = betta.at(N-2).at(0,2);
    //    M[0][3] = betta.at(N-1).at(0,0); M[0][4] = betta.at(N-1).at(0,1); M[0][5] = betta.at(N-1).at(0,2);
    //    M[0][6] = betta.at(N-0).at(0,0); M[0][7] = betta.at(N-0).at(0,1); M[0][8] = betta.at(N-0).at(0,2); c[0] = eta[0];

    //    M[1][0] = betta.at(N-2).at(1,0); M[1][1] = betta.at(N-2).at(1,1); M[1][2] = betta.at(N-2).at(1,2);
    //    M[1][3] = betta.at(N-1).at(1,0); M[1][4] = betta.at(N-1).at(1,1); M[1][5] = betta.at(N-1).at(1,2);
    //    M[1][6] = betta.at(N-0).at(1,0); M[1][7] = betta.at(N-0).at(1,1); M[1][8] = betta.at(N-0).at(1,2); c[1] = eta[1];

    //    M[2][0] = betta.at(N-2).at(2,0); M[2][1] = betta.at(N-2).at(2,1); M[2][2] = betta.at(N-2).at(2,2);
    //    M[2][3] = betta.at(N-1).at(2,0); M[2][4] = betta.at(N-1).at(2,1); M[2][5] = betta.at(N-1).at(2,2);
    //    M[2][6] = betta.at(N-0).at(2,0); M[2][7] = betta.at(N-0).at(2,1); M[2][8] = betta.at(N-0).at(2,2); c[2] = eta[2];

    //    M[3][0] = -3.0-2.0*h*a(0,0,N-2); M[3][1] = +0.0-2.0*h*a(0,1,N-2); M[3][2] = +0.0-2.0*h*a(0,2,N-2); M[3][3] = +4.0; M[3][4] = +0.0; M[3][5] = +0.0; M[3][6] = -1.0; M[3][7] = +0.0; M[3][8] = +0.0; c[3] = 2.0*h*b(0,N-2);
    //    M[4][0] = +0.0-2.0*h*a(1,0,N-2); M[4][1] = -3.0-2.0*h*a(1,1,N-2); M[4][2] = +0.0-2.0*h*a(1,2,N-2); M[4][3] = +0.0; M[4][4] = +4.0; M[4][5] = +0.0; M[4][6] = +0.0; M[4][7] = -1.0; M[4][8] = +0.0; c[4] = 2.0*h*b(1,N-2);
    //    M[5][0] = +0.0-2.0*h*a(2,0,N-2); M[5][1] = +0.0-2.0*h*a(2,1,N-2); M[5][2] = -3.0-2.0*h*a(2,2,N-2); M[5][3] = +0.0; M[5][4] = +0.0; M[5][5] = +4.0; M[5][6] = +0.0; M[5][7] = +0.0; M[5][8] = -1.0; c[5] = 2.0*h*b(2,N-2);

    //    M[6][0] = +1.0-2.0*h*a(0,0,N-0); M[6][1] = +0.0-2.0*h*a(0,1,N-0); M[6][2] = +0.0-2.0*h*a(0,2,N-0); M[6][3] = -4.0; M[6][4] = +0.0; M[6][5] = +0.0; M[6][6] = +3.0; M[6][7] = +0.0; M[6][8] = +0.0; c[6] = 2.0*h*b(0,N-0);
    //    M[7][0] = +0.0-2.0*h*a(1,0,N-0); M[7][1] = +1.0-2.0*h*a(1,1,N-0); M[7][2] = +0.0-2.0*h*a(1,2,N-0); M[7][3] = +0.0; M[7][4] = -4.0; M[7][5] = +0.0; M[7][6] = +0.0; M[7][7] = +3.0; M[7][8] = +0.0; c[7] = 2.0*h*b(1,N-0);
    //    M[8][0] = +0.0-2.0*h*a(2,0,N-0); M[8][1] = +0.0-2.0*h*a(2,1,N-0); M[8][2] = +1.0-2.0*h*a(2,2,N-0); M[8][3] = +0.0; M[8][4] = +0.0; M[8][5] = -4.0; M[8][6] = +0.0; M[8][7] = +0.0; M[8][8] = +3.0; c[8] = 2.0*h*b(2,N-0);

    //    for (unsigned int i=0; i<9; i++)
    //    {
    //        printf("%d %14.10f %14.10f\n", i, c[i], M[i][0]*rx.at(0,N-2)+M[i][1]*rx.at(1,N-2)+M[i][2]*rx.at(2,N-2)+
    //                                            M[i][3]*rx.at(0,N-1)+M[i][4]*rx.at(1,N-1)+M[i][5]*rx.at(2,N-1)+
    //                                            M[i][6]*rx.at(0,N-0)+M[i][7]*rx.at(1,N-0)+M[i][8]*rx.at(2,N-0));
    //    }

    DoubleMatrix M(6, 6);
    DoubleVector c(6);

    M[0][0] = betta.at(N-1).at(0,0);
    M[0][1] = betta.at(N-1).at(0,1);
    M[0][2] = betta.at(N-1).at(0,2);
    M[0][3] = betta.at(N-0).at(0,0);
    M[0][4] = betta.at(N-0).at(0,1);
    M[0][5] = betta.at(N-0).at(0,2);
    c[0] = eta[0];

    M[1][0] = betta.at(N-1).at(1,0);
    M[1][1] = betta.at(N-1).at(1,1);
    M[1][2] = betta.at(N-1).at(1,2);
    M[1][3] = betta.at(N-0).at(1,0);
    M[1][4] = betta.at(N-0).at(1,1);
    M[1][5] = betta.at(N-0).at(1,2);
    c[1] = eta[1];

    M[2][0] = betta.at(N-1).at(2,0);
    M[2][1] = betta.at(N-1).at(2,1);
    M[2][2] = betta.at(N-1).at(2,2);
    M[2][3] = betta.at(N-0).at(2,0);
    M[2][4] = betta.at(N-0).at(2,1);
    M[2][5] = betta.at(N-0).at(2,2);
    c[2] = eta[2];


    M[3][0] = 2.0*h*((3.0+2.0*h*a(0,0,N-2))*a(0,0,N-1) + (0.0+2.0*h*a(0,1,N-2))*a(1,0,N-1) + (0.0+2.0*h*a(0,2,N-2))*a(2,0,N-1))+4.0;
    M[3][1] = 2.0*h*((3.0+2.0*h*a(0,0,N-2))*a(0,1,N-1) + (0.0+2.0*h*a(0,1,N-2))*a(1,1,N-1) + (0.0+2.0*h*a(0,2,N-2))*a(2,1,N-1))+0.0;
    M[3][2] = 2.0*h*((3.0+2.0*h*a(0,0,N-2))*a(0,2,N-1) + (0.0+2.0*h*a(0,1,N-2))*a(1,2,N-1) + (0.0+2.0*h*a(0,2,N-2))*a(2,2,N-1))+0.0;
    M[3][3] = -(4.0+2.0*h*a(0,0,N-2));
    M[3][4] = -(0.0+2.0*h*a(0,1,N-2));
    M[3][5] = -(0.0+2.0*h*a(0,2,N-2));
    c[3] = 2.0*h*(b(0,N-2) - ((3.0+2.0*h*a(0,0,N-2))*b(0,N-1)+(0.0+2.0*h*a(0,1,N-2))*b(1,N-1)+(0.0+2.0*h*a(0,2,N-2))*b(2,N-1)));

    M[4][0] = 2.0*h*((0.0+2.0*h*a(1,0,N-2))*a(0,0,N-1) + (3.0+2.0*h*a(1,1,N-2))*a(1,0,N-1) + (0.0+2.0*h*a(1,2,N-2))*a(2,0,N-1))+0.0;
    M[4][1] = 2.0*h*((0.0+2.0*h*a(1,0,N-2))*a(0,1,N-1) + (3.0+2.0*h*a(1,1,N-2))*a(1,1,N-1) + (0.0+2.0*h*a(1,2,N-2)*a(2,1,N-1)))+4.0;
    M[4][2] = 2.0*h*((0.0+2.0*h*a(1,0,N-2))*a(0,2,N-1) + (3.0+2.0*h*a(1,1,N-2))*a(1,2,N-1) + (0.0+2.0*h*a(1,2,N-2)*a(2,2,N-1)))+0.0;
    M[4][3] = -(0.0+2.0*h*a(1,0,N-2));
    M[4][4] = -(4.0+2.0*h*a(1,1,N-2));
    M[4][5] = -(0.0+2.0*h*a(1,2,N-2));
    c[4] = 2.0*h*(b(1,N-2) - ((0.0+2.0*h*a(1,0,N-2))*b(0,N-1)+(3.0+2.0*h*a(1,1,N-2))*b(1,N-1)+(0.0+2.0*h*a(1,2,N-2))*b(2,N-1)));

    M[5][0] = 2.0*h*((0.0+2.0*h*a(2,0,N-2)*a(0,0,N-1)) + (0.0+2.0*h*a(2,1,N-2))*a(1,0,N-1) +(3.0+2.0*h*a(2,2,N-2))*a(2,0,N-1))+0.0;
    M[5][1] = 2.0*h*((0.0+2.0*h*a(2,0,N-2)*a(0,1,N-1)) + (0.0+2.0*h*a(2,1,N-2))*a(1,1,N-1) +(3.0+2.0*h*a(2,2,N-2))*a(2,1,N-1))+0.0;
    M[5][2] = 2.0*h*((0.0+2.0*h*a(2,0,N-2)*a(0,2,N-1)) + (0.0+2.0*h*a(2,1,N-2))*a(1,2,N-1) +(3.0+2.0*h*a(2,2,N-2))*a(2,2,N-1))+4.0;
    M[5][3] = -(0.0+2.0*h*a(2,0,N-2));
    M[5][4] = -(0.0+2.0*h*a(2,1,N-2));
    M[5][5] = -(4.0+2.0*h*a(2,2,N-2));
    c[5] = 2.0*h*(b(2,N-2) - ((0.0+2.0*h*a(2,0,N-2))*b(0,N-1)+(0.0+2.0*h*a(2,1,N-2))*b(1,N-1)+(3.0+2.0*h*a(2,2,N-2))*b(2,N-1)));

//    printf("det %.10f\n", M.determinant());

//    for (unsigned int i=0; i<6; i++)
//    {
//        printf("%d %14.10f %14.10f\n", i, c[i], M[i][0]*rx.at(0,N-1)+M[i][1]*rx.at(1,N-1)+M[i][2]*rx.at(2,N-1)+
//                M[i][3]*rx.at(0,N-0)+M[i][4]*rx.at(1,N-0)+M[i][5]*rx.at(2,N-0));
//    }

    DoubleVector x(6);
    GaussianElimination(M, c, x);

    printf("%14.10f %14.10f %14.10f %14.10f\n", x[0], x[3], fx(0,N-1), fx(0,N-0));
    printf("%14.10f %14.10f %14.10f %14.10f\n", x[1], x[4], fx(1,N-1), fx(1,N-0));
    printf("%14.10f %14.10f %14.10f %14.10f\n", x[2], x[5], fx(2,N-1), fx(2,N-0));

//    for (unsigned int i=0; i<6; i++)
//    {
//        printf("%d %14.10f %14.10f\n", i, c[i], M[i][0]*x[0]+M[i][1]*x[1]+M[i][2]*x[2]+
//                M[i][3]*x[3]+M[i][4]*x[4]+M[i][5]*x[6]);
//        printf("%d %14.10f %14.10f\n", i, c[i], M[i][0]*fx(0,N-1)+M[i][1]*fx(1,N-1)+M[i][2]*fx(2,N-1)+
//                M[i][3]*fx(0,N-0)+M[i][4]*fx(1,N-0)+M[i][5]*fx(2,N-0));
//    }
}
