#include "lbvpode.h"

double N2D1[5][5] =
{
    {-3.0, +4.0, -1.0},
    {+1.0, +0.0, -1.0},
    {+1.0, -4.0, +3.0}
};

double N2D2[5][5] =
{
    {+1.0, -2.0, +1.0},
    {+1.0, -2.0, +1.0},
    {+1.0, -2.0, +1.0}
};

double N4D1[5][5] =
{
    {-25.0, +48.0, -36.0, +16.0, -3.0},
    {-3.0,  -10.0, +18.0, -6.0,  +1.0},
    {+1.0,  -8.0,  +0.0,  +8.0,  -1.0},
    {-1.0,  +6.0,  -18.0, +10.0, +3.0},
    {+3.0, -16.0,  +36.0, -48.0, +25.0}
};

double N4D2[5][5] =
{
    {+70.0, -208.0, +228.0, -112.0, +22.0},
    {+22.0, -40.0,  +12.0,  +8.0,   -2.0},
    {-2.0,  +32.0,  -60.0,  +32.0,  -2.0},
    {-2.0,  +8.0,   +12.0,  -40.0,  +22.0},
    {+22.0, -112.0, +228.0, -208.0, +70.0}
};

double N6D1[7][7] =
{
    {+812.0, -3132.0, +5265.0, -5080.0, +2970.0, -972.0,  +137.0},
    {+137.0, -147.0,  -255.0,  +470.0,  -285.0,  +93.0,   -13.0},
    {-13.0,  +228.0,  -420.0,  +200.0,  +15.0,   -12.0,   +2.0},
    {+2.0,   -27.0,   +270.0,  -490.0,  +270.0,  -27.0,   +2.0},
    {+2.0,   -12.0,   +15.0,   +200.0,  -420.0,  +228.0,  -13.0},
    {-13.0,  +93.0,   -285.0,  +470.0,  -255.0,  -147.0,  +137.0},
    {+137.0, -972.0,  +2970.0, -5080.0, +5265.0, -3132.0, +812.0}
};

double N6D2[7][7] =
{
    {-147.0, +360.0, -450.0,  +400.0, -225.0, +72.0, -10.0},
    {-10.0,  -77.0,  +150.0,  -100.0, +50.0,  -15.0,  +2.0},
    {+2.0,   -24.0,  -35.0,   +80.0,  -30.0,  +8.0,   -1.0},
    {-1.0,   +9.0,   -45.0,   +0.0,   +45.0,  -9.0,   +1.0},
    {+1.0,   -8.0,   +30.0,   -80.0,  +35.0,  +24.0,  -2.0},
    {-2.0,   +15.0,  -50.0,   +100.0, -150.0, +77.0,  +10.0},
    {+10.0,  -72.0,   +225.0, -400.0, +450.0, -360.0, +147.0}
};

void LinearBoundaryValueProblemODE::calculate2N(DoubleVector &x, double h, unsigned int N)
{
    x.clear();
    x.resize(N+1);

    DoubleVector da(N-1);
    DoubleVector db(N-1);
    DoubleVector dc(N-1);
    DoubleVector dd(N-1);
    DoubleVector rx(N-1);

    for (unsigned int i=1; i<=N-1; i++)
    {
        double alpha = r(i)/(h*h) - p(i)/(2.0*h);
        double betta = q(i) - (2.0*r(i))/(h*h);
        double gamma = r(i)/(h*h) + p(i)/(2.0*h);

        da[i-1] = alpha;
        db[i-1] = betta;
        dc[i-1] = gamma;
        dd[i-1] = f(i);
    }

    dd[0]   -= da[0]*boundary(Left);
    dd[N-2] -= dc[N-2]*boundary(Right);

    da[0]   = 0.0;
    dc[N-2] = 0.0;

    tomasAlgorithm(da.data(), db.data(), dc.data(), dd.data(), rx.data(), rx.size());

    x[0] = boundary(Left);
    x[N] = boundary(Right);

    for (unsigned int i=1; i<=N-1; i++)
    {
        x[i] = rx[i-1];
    }

    da.clear();
    db.clear();
    dc.clear();
    dd.clear();
    rx.clear();
}

void LinearBoundaryValueProblemODE::calculateN4L2RD(DoubleVector &x, double h, unsigned int N)
{
    const unsigned int k = 4;
    double m1 = 1.0/(12.0*h);
    double m2 = 1.0/(24.0*h*h);

    x.clear();
    x.resize(N+1, 0.0);

    DoubleMatrix A(k, k, 0.0);
    DoubleVector b(k, 0.0);
    DoubleVector z(k, 0.0);
    DoubleMatrix ems(N-k, k);

    /* border conditions */
    x.at(0) = boundary(Left);
    x.at(N) = boundary(Right);

    double alpha0 = r(1)*m2;
    double betta0 = p(1)*m1;
    A[0][0] = -40.0*alpha0 - 10.0*betta0 + q(1);
    A[0][1] = +12.0*alpha0 + 18.0*betta0;
    A[0][2] = +8.0*alpha0  - 6.0*betta0;
    A[0][3] = -2.0*alpha0  + 1.0*betta0;
    b[0] = f(1) - (22.0*alpha0 - 3.0*betta0)*x.at(0);

    A[0][1] /= A[0][0];
    A[0][2] /= A[0][0];
    A[0][3] /= A[0][0];
    b[0]    /= A[0][0];
    A[0][0] = 1.0;

    ems.at(0,0) = A[0][1];
    ems.at(0,1) = A[0][2];
    ems.at(0,2) = A[0][3];
    ems.at(0,3) = b[0];

    for (unsigned int n=1; n<=N-(k+1); n++)
    {
        // + * * * *
        double alphai = r(n)*m2;
        double bettai = p(n)*m1;
        double g1 = +70.0*alphai  - 25.0*bettai + q(n);
        double g2 = -208.0*alphai + 48.0*bettai;
        double g3 = +228.0*alphai - 36.0*bettai;
        double g4 = -112.0*alphai + 16.0*bettai;
        double g5 = +22.0*alphai  - 3.0*bettai;
        double fi = f(n);

        // * + * * *
        //double alphai = r(n+1)/m2;
        //double bettai = p(n+1)/m1;
        //double g1 = +22.0*alphai - 3.0*bettai;
        //double g2 = -40.0*alphai - 10.0*bettai + q(n+1);
        //double g3 = +12.0*alphai + 18.0*bettai;
        //double g4 = +8.0*alphai  - 6.0*bettai;
        //double g5 = -2.0*alphai  + bettai;
        //double fi = f(n+1);

        // * * + * *
        //double alphai = r(n+2)/m2;
        //double bettai = p(n+2)/m1;
        //double g1 = -2.0*alphai  + 1.0*bettai;
        //double g2 = +32.0*alphai - 8.0*bettai;
        //double g3 = -60.0*alphai + 0.0*bettai + q(n+2);
        //double g4 = +32.0*alphai + 8.0*bettai;
        //double g5 = -2.0*alphai  - 1.0*bettai;
        //double fi = f(n+2);

        g2 /= -g1;
        g3 /= -g1;
        g4 /= -g1;
        g5 /= -g1;
        fi /= +g1;
        g1 = 1.0;

        A[0][0] = A[0][1] + g2;
        A[0][1] = A[0][2] + g3;
        A[0][2] = A[0][3] + g4;
        A[0][3] = g5;
        b[0]    = b[0] - fi;

        A[0][1] /= A[0][0];
        A[0][2] /= A[0][0];
        A[0][3] /= A[0][0];
        b[0]    /= A[0][0];
        A[0][0] = 1.0;

        ems.at(n,0) = A[0][1];
        ems.at(n,1) = A[0][2];
        ems.at(n,2) = A[0][3];
        ems.at(n,3) = b[0];
    }

    double alpha1 = r(N-3)*m2;
    double betta1 = p(N-3)*m1;
    A[1][0] = +22.0*alpha1 - 3.0*betta1;
    A[1][1] = -40.0*alpha1 - 10.0*betta1 + q(N-3);
    A[1][2] = +12.0*alpha1 + 18.0*betta1;
    A[1][3] = +8.0*alpha1  - 6.0*betta1;
    b[1]    = f(N-3) - (-2.0*alpha1 + betta1)*x.at(N);

    double alpha2 = r(N-2)*m2;
    double betta2 = p(N-2)*m1;
    A[2][0] = -2.0*alpha2  + betta2;
    A[2][1] = +32.0*alpha2 - 8.0*betta2;
    A[2][2] = -60.0*alpha2 + q(N-2);
    A[2][3] = +32.0*alpha2 + 8.0*betta2;
    b[2]    = f(N-2) - (-2.0*alpha2 - betta2)*x.at(N);

    double alpha3 = r(N-1)*m2;
    double betta3 = p(N-1)*m1;
    A[3][0] = -2.0*alpha3  - betta3;
    A[3][1] = +8.0*alpha3  + 6.0*betta3;
    A[3][2] = +12.0*alpha3 - 18.0*betta3;
    A[3][3] = -40.0*alpha3 + 10.0*betta3 + q(N-1);
    b[3]    = f(N-1) - (22.0*alpha3 + 3.0*betta3)*x.at(N);

    GaussianElimination(A, b, z);

    x.at(N-1) = z.at(3);
    x.at(N-2) = z.at(2);
    x.at(N-3) = z.at(1);
    x.at(N-4) = z.at(0);

//    for (unsigned int n=N-(k+1); n>=1; n--)
//    {
//        x.at(n) = -ems.at(n-1,0)*x.at(n+1) - ems.at(n-1,1)*x.at(n+2) - ems.at(n-1,2)*x.at(n+3) + ems.at(n-1,3);
//    }

    for (unsigned int n=N-(k+1); n>=1; n--)
    {
        double alphai = r(n)*m2;
        double bettai = p(n)*m1;
        double d0 = +70.0*alphai  - 25.0*bettai + q(n);
        double d1 = -208.0*alphai + 48.0*bettai;
        double d2 = +228.0*alphai - 36.0*bettai;
        double d3 = -112.0*alphai + 16.0*bettai;
        double d4 = +22.0*alphai  - 3.0*bettai;
        double fi = f(n);
        x.at(n) = -d1*x.at(n+1) - d2*x.at(n+2) - d3*x.at(n+3) - d4*x.at(n+4) + fi;
        x.at(n) /= d0;
    }

    ems.clear();
    z.clear();
    b.clear();
    A.clear();
}

void LinearBoundaryValueProblemODE::calculateN4R2LD(DoubleVector &x, double h, unsigned int N)
{
    x.clear();
    x.resize(N+1, 0.0);
    x.at(0) = boundary(Left);
    x.at(N) = boundary(Right);

    DoubleMatrix A(4,4);
    DoubleVector b(4);
    DoubleVector z(4);

    {
        double alpha0 = r(1)/(24.0*h*h);
        double betta0 = p(1)/(12.0*h);
        A[0][0] = -40.0*alpha0 - 10.0*betta0 + q(1);
        A[0][1] = +12.0*alpha0 + 18.0*betta0;
        A[0][2] = +8.0*alpha0  - 6.0*betta0;
        A[0][3] = -2.0*alpha0  + betta0;
        b[0]    = f(1) - (22.0*alpha0 - 3.0*betta0)*x.at(0);
    }
    {
        double alpha1 = r(2)/(24.0*h*h);
        double betta1 = p(2)/(12.0*h);
        A[1][0] = +32.0*alpha1 - 8.0*betta1;
        A[1][1] = -60.0*alpha1 + q(2);
        A[1][2] = +32.0*alpha1 + 8.0*betta1;
        A[1][3] = -2.0*alpha1  - betta1;
        b[1]    = f(2) - (-2.0*alpha1 + betta1)*x.at(0);
    }
    {
        double alpha2 = r(3)/(24.0*h*h);
        double betta2 = p(3)/(12.0*h);
        A[2][0] = +8.0*alpha2  + 6.0*betta2;
        A[2][1] = +12.0*alpha2 - 18.0*betta2;
        A[2][2] = -40.0*alpha2 + 10.0*betta2 + q(3);
        A[2][3] = +22.0*alpha2 + 3.0*betta2;
        b[2]    = f(3) - (-2.0*alpha2 - betta2)*x.at(0);
    }
    {
        double alphaN1 = r(N-1)/(24.0*h*h);
        double bettaN1 = p(N-1)/(12.0*h);
        double c1 = -2.0*alphaN1  - bettaN1;
        double c2 = +8.0*alphaN1  + 6.0*bettaN1;
        double c3 = +12.0*alphaN1 - 18.0*bettaN1;
        double c4 = -40.0*alphaN1 + 10.0*bettaN1 + q(N-1);
        double c = f(N-1) - (22.0*alphaN1 + 3.0*bettaN1)*x.at(N);

        c3 /= c4;
        c2 /= c4;
        c1 /= c4;
        c  /= c4;
        c4 = 1.0;

        for (unsigned int i=N-2; i>=4; i--)
        {
            double alphai = r(i)/(24.0*h*h);
            double bettai = p(i)/(12.0*h);

            double g1 = -2.0*alphai  - bettai;
            double g2 = +8.0*alphai  + 6.0*bettai;
            double g3 = +12.0*alphai - 18.0*bettai;
            double g4 = -40.0*alphai + 10.0*bettai + q(i);
            double g5 = +22.0*alphai + 3.0*bettai;
            double fi = f(i);

            g4 /= -g5;
            g3 /= -g5;
            g2 /= -g5;
            g1 /= -g5;
            fi /= +g5;
            g5 = 1.0;

            c4 = c3 + g4;
            c3 = (c2 + g3)/c4;
            c2 = (c1 + g2)/c4;
            c1 = g1/c4;
            c  = (c - fi)/c4;
            c4 = 1.0;
        }

        A[3][0] = c1;
        A[3][1] = c2;
        A[3][2] = c3;
        A[3][3] = c4;
        b[3]    = c;
    }
    GaussianElimination(A, b, z);

    x.at(1) = z.at(0);
    x.at(2) = z.at(1);
    x.at(3) = z.at(2);
    x.at(4) = z.at(3);
    for (unsigned int i=5; i<=N-1; i++)
    {
        double alpha = r(i)/(24.0*h*h);
        double betta = p(i)/(12.0*h);

        double d0 = +70.0*alpha  + 25.0*betta + q(i);
        double d1 = +208.0*alpha + 48.0*betta;
        double d2 = -228.0*alpha - 36.0*betta;
        double d3 = +112.0*alpha + 16.0*betta;
        double d4 = -22.0*alpha  - 3.0*betta;

        x.at(i) = d1*x.at(i-1) + d2*x.at(i-2) + d3*x.at(i-3) + d4*x.at(i-4) + f(i);
        x.at(i) /= d0;
    }

    A.clear();
    b.clear();
    z.clear();
}

double U(unsigned int i, double h)
{
    double x = i*h;
    return x*x*x;
}

void LinearBoundaryValueProblemODE::calculateN6L2RD(DoubleVector &x, double h, unsigned int N)
{
    const unsigned int k=6;
    double m1 = 1.0/(60.0*h);
    double m2 = 1.0/(180.0*h*h);

    x.clear();
    x.resize(N+1, 0.0);

    DoubleMatrix A(k, k, 0.0);
    DoubleVector b(k, 0.0);
    DoubleVector z(k, 0.0);
    DoubleMatrix ems(N-k, k);

    /* border conditions */
    x.at(0) = boundary(Left);
    x.at(N) = boundary(Right);

    double alpha0 = r(1)*m2;
    double betta0 = p(1)*m1;
    A[0][0] = -147.0*alpha0 - 77.0*betta0 + q(1);
    A[0][1] = -255.0*alpha0 + 150.0*betta0;
    A[0][2] = +470.0*alpha0 - 100.0*betta0;
    A[0][3] = -285.0*alpha0 + 50.0*betta0;
    A[0][4] = +93.0*alpha0  - 15.0*betta0;
    A[0][5] = -13.0*alpha0  + 2.0*betta0;
    b[0]    = f(1) - (+137.0*alpha0 - 10.0*betta0)*x.at(0);

    A[0][1] /= A[0][0];
    A[0][2] /= A[0][0];
    A[0][3] /= A[0][0];
    A[0][4] /= A[0][0];
    A[0][5] /= A[0][0];
    b[0]    /= A[0][0];
    A[0][0] = 1.0;

    ems.at(0,0) = A[0][1];
    ems.at(0,1) = A[0][2];
    ems.at(0,2) = A[0][3];
    ems.at(0,3) = A[0][4];
    ems.at(0,4) = A[0][5];
    ems.at(0,5) = b[0];

    for (unsigned int n=1; n<=N-(k+1); n++)
    {
        // + * * * * * *
        double alphai = r(n)*m2;
        double bettai = p(n)*m1;
        double g1 = +812.0*alphai  - 147.0*bettai + q(n);
        double g2 = -3132.0*alphai + 360.0*bettai;
        double g3 = +5265.0*alphai - 450.0*bettai;
        double g4 = -5080.0*alphai + 400.0*bettai;
        double g5 = +2970.0*alphai - 225.0*bettai;
        double g6 = -972.0*alphai  + 72.0*bettai;
        double g7 = +137.0*alphai  - 10.0*bettai;
        double fi = f(n);

        // * + * * * * *
        //double alphai = r(n+1)*m1;
        //double bettai = p(n+1)*m2;
        //double g1 = +137.0*alphai - 10.0*bettai;
        //double g2 = -147.0*alphai - 77.0*bettai + q(n+1);
        //double g3 = -255.0*alphai + 150.0*bettai;
        //double g4 = +470.0*alphai - 100.0*bettai;
        //double g5 = -285.0*alphai + 50.0*bettai;
        //double g6 = +93.0*alphai  - 15.0*bettai;
        //double g7 = -13.0*alphai  + 2.0*bettai;
        //double fi = f(n+1);

        // * * + * * * *
        //double alphai = r(n+2)*m1;
        //double bettai = p(n+2)*m2;
        //double g1 = -13.0*alphai  + 2.0*bettai;
        //double g2 = +228.0*alphai - 24.0*bettai;
        //double g3 = -420.0*alphai - 35.0*bettai + q(n+2);
        //double g4 = +200.0*alphai + 80.0*bettai;
        //double g5 = +15.0*alphai  - 30.0*bettai;
        //double g6 = -12.0*alphai  + 8.0*bettai;
        //double g7 = +2.0*alphai   - 1.0*bettai;
        //double fi = f(n+2);

        g2 /= -g1;
        g3 /= -g1;
        g4 /= -g1;
        g5 /= -g1;
        g6 /= -g1;
        g7 /= -g1;
        fi /= +g1;
        g1 = 1.0;

        A[0][0] = A[0][1] + g2;
        A[0][1] = A[0][2] + g3;
        A[0][2] = A[0][3] + g4;
        A[0][3] = A[0][4] + g5;
        A[0][4] = A[0][5] + g6;
        A[0][5] = g7;
        b[0]    = b[0] - fi;

        A[0][1] /= A[0][0];
        A[0][2] /= A[0][0];
        A[0][3] /= A[0][0];
        A[0][4] /= A[0][0];
        A[0][5] /= A[0][0];
        b[0]    /= A[0][0];
        A[0][0] = 1.0;

        ems.at(n,0) = A[0][1];
        ems.at(n,1) = A[0][2];
        ems.at(n,2) = A[0][3];
        ems.at(n,3) = A[0][4];
        ems.at(n,4) = A[0][5];
        ems.at(n,5) = b[0];
    }

    double alpha1 = r(N-5)*m2;
    double betta1 = p(N-5)*m1;
    A[1][0] = +137.0*alpha1 - 10.0*betta1;
    A[1][1] = -147.0*alpha1 - 77.0*betta1 + q(N-5);
    A[1][2] = -255.0*alpha1 + 150.0*betta1;
    A[1][3] = +470.0*alpha1 - 100.0*betta1;
    A[1][4] = -285.0*alpha1 + 50.0*betta1;
    A[1][5] = +93.0*alpha1  - 15.0*betta1;
    b[1]    = f(N-5) - (-13.0*alpha1 + 2.0*betta1)*x.at(N);

    double alpha2 = r(N-4)*m2;
    double betta2 = p(N-4)*m1;
    A[2][0] = -13.0*alpha2  + 2.0*betta2;
    A[2][1] = +228.0*alpha2 - 24.0*betta2;
    A[2][2] = -420.0*alpha2 - 35.0*betta2 + q(N-4);
    A[2][3] = +200.0*alpha2 + 80.0*betta2;
    A[2][4] = +15.0*alpha2  - 30.0*betta2;
    A[2][5] = -12.0*alpha2  + 8.0*betta2;
    b[2]    = f(N-4) - (2.0*alpha2 - betta2)*x.at(N);

    double alpha3 = r(N-3)*m2;
    double betta3 = p(N-3)*m1;
    A[3][0] = +2.0*alpha3   - betta3;
    A[3][1] = -27.0*alpha3  + 9.0*betta3;
    A[3][2] = +270.0*alpha3 - 45.0*betta3;
    A[3][3] = -490.0*alpha3 + q(N-3);
    A[3][4] = +270.0*alpha3 + 45.0*betta3;
    A[3][5] = -27.0*alpha3  - 9.0*betta3;
    b[3]    = f(N-3) - (2.0*alpha3 + betta3)*x.at(N);

    double alpha4 = r(N-2)*m2;
    double betta4 = p(N-2)*m1;
    A[4][0] = +2.0*alpha4   + betta4;
    A[4][1] = -12.0*alpha4  - 8.0*betta4;
    A[4][2] = +15.0*alpha4  + 30.0*betta4;
    A[4][3] = +200.0*alpha4 - 80.0*betta4;
    A[4][4] = -420.0*alpha4 + 35.0*betta4 + q(N-2);
    A[4][5] = +228.0*alpha4 + 24.0*betta4;
    b[4]    = f(N-2) - (-13.0*alpha4 - 2.0*betta4)*x.at(N);

    double alpha5 = r(N-1)*m2;
    double betta5 = p(N-1)*m1;
    A[5][0] = -13.0*alpha5  - 2.0*betta5;
    A[5][1] = +93.0*alpha5  + 15.0*betta5;
    A[5][2] = -285.0*alpha5 - 50.0*betta5;
    A[5][3] = +470.0*alpha5 + 100.0*betta5;
    A[5][4] = -255.0*alpha5 - 150.0*betta5;
    A[5][5] = -147.0*alpha5 + 77.0*betta5 + q(N-1);
    b[5]    = f(N-1) - (137.0*alpha5 + 10.0*betta5)*x.at(N);

    GaussianElimination(A, b, z);

    x.at(N-1) = z.at(5);
    x.at(N-2) = z.at(4);
    x.at(N-3) = z.at(3);
    x.at(N-4) = z.at(2);
    x.at(N-5) = z.at(1);
    x.at(N-6) = z.at(0);

//    for (unsigned int n=N-(k+1); n>=1; n--)
//    {
//        x.at(n) = -ems.at(n-1,0)*x.at(n+1)
//                 - ems.at(n-1,1)*x.at(n+2)
//                 - ems.at(n-1,2)*x.at(n+3)
//                 - ems.at(n-1,3)*x.at(n+4)
//                 - ems.at(n-1,4)*x.at(n+5)
//                 + ems.at(n-1,5);
//    }

    for (unsigned int n=N-(k+1); n>=1; n--)
    {
        double alphai = r(n)*m2;
        double bettai = p(n)*m1;
        double d0 = +812.0*alphai  - 147*bettai + q(n);
        double d1 = -3132.0*alphai + 360.0*bettai;
        double d2 = +5265.0*alphai - 450.0*bettai;
        double d3 = -5080.0*alphai + 400.0*bettai;
        double d4 = +2970.0*alphai - 225.0*bettai;
        double d5 = -972.0*alphai  + 72.0*bettai;
        double d6 = +137.0*alphai  - 10.0*bettai;
        x.at(n) = -d1*x.at(n+1) - d2*x.at(n+2) - d3*x.at(n+3) - d4*x.at(n+4) - d5*x.at(n+5) - d6*x.at(n+6) + f(n);
        x.at(n) /= d0;
    }

    ems.clear();
    z.clear();
    b.clear();
    A.clear();
}

void LinearBoundaryValueProblemODE::calculateN6R2LD(DoubleVector &x, double h, unsigned int N)
{
    x.clear();
    x.resize(N+1, 0.0);
    x.at(0) = boundary(Left);
    x.at(N) = boundary(Right);

    unsigned int K=6;
    DoubleMatrix A(K, K, 0.0);
    DoubleVector b(K, 0.0);
    DoubleVector z(K, 0.0);

    double m1 = 1.0/(180.0*h*h);
    double m2 = 1.0/(60.0*h);

    double N6D2[7][7] =
    {
        {+812.0, -3132.0, +5265.0, -5080.0, +2970.0, -972.0,  +137.0},
        {+137.0, -147.0,  -255.0,  +470.0,  -285.0,  +93.0,   -13.0},
        {-13.0,  +228.0,  -420.0,  +200.0,  +15.0,   -12.0,   +2.0},
        {+2.0,   -27.0,   +270.0,  -490.0,  +270.0,  -27.0,   +2.0},
        {+2.0,   -12.0,   +15.0,   +200.0,  -420.0,  +228.0,  -13.0},
        {-13.0,  +93.0,   -285.0,  +470.0,  -255.0,  -147.0,  +137.0},
        {+137.0, -972.0,  +2970.0, -5080.0, +5265.0, -3132.0, +812.0}
    };

    double N6D1[7][7] =
    {
        {-147.0, +360.0, -450.0,  +400.0, -225.0, +72.0, -10.0},
        {-10.0,  -77.0,  +150.0,  -100.0, +50.0,  -15.0,  +2.0},
        {+2.0,   -24.0,  -35.0,   +80.0,  -30.0,  +8.0,   -1.0},
        {-1.0,   +9.0,   -45.0,   +0.0,   +45.0,  -9.0,   +1.0},
        {+1.0,   -8.0,   +30.0,   -80.0,  +35.0,  +24.0,  -2.0},
        {-2.0,   +15.0,  -50.0,   +100.0, -150.0, +77.0,  +10.0},
        {+10.0,  -72.0,   +225.0, -400.0, +450.0, -360.0, +147.0}
    };

    for (unsigned int i=1; i<=5; i++)
    {
        double alpha0 = r(i)*m1;
        double betta0 = p(i)*m2;
        for (unsigned int j=1; j<=6; j++)
        {
            A[i-1][j-1] = N6D2[i][j]*alpha0 + N6D1[i][j]*betta0;
        }
        b[i-1] = f(i) - (N6D2[i][0]*alpha0 + N6D1[i][0]*betta0)*x.at(0);
        A[i-1][i-1] += q(i);
    }

    //        {
    //            double alphaN1 = r(N-1)*m1;
    //            double bettaN1 = p(N-1)*m2;
    //            for (unsigned int j=1; j<=6; j++)
    //            {
    //                A[5][j-1] = N6D2[5][j-1]*alphaN1 + N6D1[5][j-1]*bettaN1;
    //            }
    //            A[5][5] += q(N-1);
    //            b[5] = f(N-1) - (N6D2[5][6]*alphaN1 + N6D1[5][6]*bettaN1)*x.at(N);

    //            A[5][0] = N6D2[5][0]*alphaN1 + N6D1[5][0]*bettaN1;
    //            A[5][1] = N6D2[5][1]*alphaN1 + N6D1[5][1]*bettaN1;
    //            A[5][2] = N6D2[5][2]*alphaN1 + N6D1[5][2]*bettaN1;
    //            A[5][3] = N6D2[5][3]*alphaN1 + N6D1[5][3]*bettaN1;
    //            A[5][4] = N6D2[5][4]*alphaN1 + N6D1[5][4]*bettaN1;
    //            A[5][5] = N6D2[5][5]*alphaN1 + N6D1[5][5]*bettaN1 + q(N-1);
    //            b[5]    = f(N-1) - (N6D2[5][6]*alphaN1 + N6D1[5][6]*bettaN1)*x.at(N);

    //            A[5][4] /= A[5][5];
    //            A[5][3] /= A[5][5];
    //            A[5][2] /= A[5][5];
    //            A[5][1] /= A[5][5];
    //            A[5][0] /= A[5][5];
    //            b[5]    /= A[5][5];
    //            A[5][5] = 1.0;

    //            for (unsigned int i=N-1; i>=7; i--)
    //            {
    //                double alphai = r(i)*m1;
    //                double bettai = p(i)*m2;

    //                double g[8] = {0.0};
    //                for (unsigned int i=1; i<=7; i++)
    //                {
    //                    g[i] = N6D2[6][i-1]*alphai + N6D1[6][i-1]*bettai;
    //                }
    //                g[0]  = f(i);

    //                g[1] = N6D2[6][0]*alphai + N6D1[6][0]*bettai;
    //                g[2] = N6D2[6][1]*alphai + N6D1[6][1]*bettai;
    //                g[3] = N6D2[6][2]*alphai + N6D1[6][2]*bettai;
    //                g[4] = N6D2[6][3]*alphai + N6D1[6][3]*bettai;
    //                g[5] = N6D2[6][4]*alphai + N6D1[6][4]*bettai;
    //                g[6] = N6D2[6][5]*alphai + N6D1[6][5]*bettai;
    //                g[7] = N6D2[6][6]*alphai + N6D1[6][6]*bettai + q(i);
    //                g[0]  = f(i);

    //                for (unsigned int i=6; i>=1; i--)
    //                {
    //                    g[i] /= -g[7];
    //                }
    //                g[0] /= +g[7];
    //                g[7] = 1.0;

    //    //            g[6] /= -g[7];
    //    //            g[5] /= -g[7];
    //    //            g[4] /= -g[7];
    //    //            g[3] /= -g[7];
    //    //            g[2] /= -g[7];
    //    //            g[1] /= -g[7];
    //    //            g[0] /= +g[7];
    //    //            g[7] = 1.0;

    //                A[5][5] =  A[5][4] + g[6];
    //                A[5][4] = (A[5][3] + g[5])/A[5][5];
    //                A[5][3] = (A[5][2] + g[4])/A[5][5];
    //                A[5][2] = (A[5][1] + g[3])/A[5][5];
    //                A[5][1] = (A[5][0] + g[2])/A[5][5];
    //                A[5][0] = g[1]/A[5][5];
    //                b[5]    = (b[5] - g[0])/A[5][5];
    //                A[5][5] = 1.0;
    //            }
    //        }

    //        {
    //            double alpha0 = r(1)*m1;
    //            double betta0 = p(1)*m2;
    //            A[0][0] = N6D2[1][1]*alpha0 + N6D1[1][1]*betta0 + q(1);
    //            A[0][1] = N6D2[1][2]*alpha0 + N6D1[1][2]*betta0;
    //            A[0][2] = N6D2[1][3]*alpha0 + N6D1[1][3]*betta0;
    //            A[0][3] = N6D2[1][4]*alpha0 + N6D1[1][4]*betta0;
    //            A[0][4] = N6D2[1][5]*alpha0 + N6D1[1][5]*betta0;
    //            A[0][5] = N6D2[1][6]*alpha0 + N6D1[1][6]*betta0;
    //            b[0]    = f(1) - (N6D2[1][0]*alpha0 + N6D1[1][0]*betta0)*x.at(0);
    //        }
    //        {
    //            double alpha0 = r(2)*m1;
    //            double betta0 = p(2)*m2;
    //            A[1][0] = N6D2[2][1]*alpha0 + N6D1[2][1]*betta0;
    //            A[1][1] = N6D2[2][2]*alpha0 + N6D1[2][2]*betta0 + q(2);
    //            A[1][2] = N6D2[2][3]*alpha0 + N6D1[2][3]*betta0;
    //            A[1][3] = N6D2[2][4]*alpha0 + N6D1[2][4]*betta0;
    //            A[1][4] = N6D2[2][5]*alpha0 + N6D1[2][5]*betta0;
    //            A[1][5] = N6D2[2][6]*alpha0 + N6D1[2][6]*betta0;
    //            b[1]    = f(2) - (N6D2[2][0]*alpha0 + N6D1[2][0]*betta0)*x.at(0);
    //        }
    //        {
    //            double alpha0 = r(3)*m1;
    //            double betta0 = p(3)*m2;
    //            A[2][0] = N6D2[3][1]*alpha0 + N6D1[3][1]*betta0;
    //            A[2][1] = N6D2[3][2]*alpha0 + N6D1[3][2]*betta0;
    //            A[2][2] = N6D2[3][3]*alpha0 + N6D1[3][3]*betta0 + q(3);
    //            A[2][3] = N6D2[3][4]*alpha0 + N6D1[3][4]*betta0;
    //            A[2][4] = N6D2[3][5]*alpha0 + N6D1[3][5]*betta0;
    //            A[2][5] = N6D2[3][6]*alpha0 + N6D1[3][6]*betta0;
    //            b[2]    = f(3) - (N6D2[3][0]*alpha0 + N6D1[3][0]*betta0)*x.at(0);
    //        }
    //        {
    //            double alpha0 = r(4)*m1;
    //            double betta0 = p(4)*m2;
    //            A[3][0] = N6D2[4][1]*alpha0 + N6D1[4][1]*betta0;
    //            A[3][1] = N6D2[4][2]*alpha0 + N6D1[4][2]*betta0;
    //            A[3][2] = N6D2[4][3]*alpha0 + N6D1[4][3]*betta0;
    //            A[3][3] = N6D2[4][4]*alpha0 + N6D1[4][4]*betta0 + q(4);
    //            A[3][4] = N6D2[4][5]*alpha0 + N6D1[4][5]*betta0;
    //            A[3][5] = N6D2[4][6]*alpha0 + N6D1[4][6]*betta0;
    //            b[3]    = f(4) - (N6D2[4][0]*alpha0 + N6D1[4][0]*betta0)*x.at(0);
    //        }
    //        {
    //            double alpha0 = r(5)*m1;
    //            double betta0 = p(5)*m2;
    //            A[4][0] = N6D2[5][1]*alpha0 + N6D1[5][1]*betta0;
    //            A[4][1] = N6D2[5][2]*alpha0 + N6D1[5][2]*betta0;
    //            A[4][2] = N6D2[5][3]*alpha0 + N6D1[5][3]*betta0;
    //            A[4][3] = N6D2[5][4]*alpha0 + N6D1[5][4]*betta0;
    //            A[4][4] = N6D2[5][5]*alpha0 + N6D1[5][5]*betta0 + q(5);
    //            A[4][5] = N6D2[5][6]*alpha0 + N6D1[5][6]*betta0;
    //            b[4]    = f(5) - (N6D2[5][0]*alpha0 + N6D1[5][0]*betta0)*x.at(0);
    //        }

    //    {
    //        double alphaN1 = r(N-1)*m1;
    //        double bettaN1 = p(N-1)*m2;
    //        A[5][0] = N6D2[5][0]*alphaN1 + N6D1[5][0]*bettaN1;
    //        A[5][1] = N6D2[5][1]*alphaN1 + N6D1[5][1]*bettaN1;
    //        A[5][2] = N6D2[5][2]*alphaN1 + N6D1[5][2]*bettaN1;
    //        A[5][3] = N6D2[5][3]*alphaN1 + N6D1[5][3]*bettaN1;
    //        A[5][4] = N6D2[5][4]*alphaN1 + N6D1[5][4]*bettaN1;
    //        A[5][5] = N6D2[5][5]*alphaN1 + N6D1[5][5]*bettaN1 + q(N-1);
    //        b[5]    = f(N-1) - (N6D2[5][6]*alphaN1 + N6D1[5][6]*bettaN1)*x.at(N);

    //        A[5][4] /= A[5][5];
    //        A[5][3] /= A[5][5];
    //        A[5][2] /= A[5][5];
    //        A[5][1] /= A[5][5];
    //        A[5][0] /= A[5][5];
    //        b[5]    /= A[5][5];
    //        A[5][5] = 1.0;

    //        for (unsigned int i=N-1; i>=7; i--)
    //        {
    //            double alphai = r(i)*m1;
    //            double bettai = p(i)*m2;

    //            double g1 = N6D2[6][0]*alphai + N6D1[6][0]*bettai;
    //            double g2 = N6D2[6][1]*alphai + N6D1[6][1]*bettai;
    //            double g3 = N6D2[6][2]*alphai + N6D1[6][2]*bettai;
    //            double g4 = N6D2[6][3]*alphai + N6D1[6][3]*bettai;
    //            double g5 = N6D2[6][4]*alphai + N6D1[6][4]*bettai;
    //            double g6 = N6D2[6][5]*alphai + N6D1[6][5]*bettai;
    //            double g7 = N6D2[6][6]*alphai + N6D1[6][6]*bettai + q(i);
    //            double fi = f(i);

    //            g6 /= -g7;
    //            g5 /= -g7;
    //            g4 /= -g7;
    //            g3 /= -g7;
    //            g2 /= -g7;
    //            g1 /= -g7;
    //            fi /= +g7;
    //            g7 = 1.0;

    //            A[5][5] =  A[5][4] + g6;
    //            A[5][4] = (A[5][3] + g5)/A[5][5];
    //            A[5][3] = (A[5][2] + g4)/A[5][5];
    //            A[5][2] = (A[5][1] + g3)/A[5][5];
    //            A[5][1] = (A[5][0] + g2)/A[5][5];
    //            A[5][0] = g1/A[5][5];
    //            b[5]    = (b[5] - fi)/A[5][5];
    //            A[5][5] = 1.0;
    //        }
    //    }

    GaussianElimination(A, b, z);

    x.at(1) = z.at(0);
    x.at(2) = z.at(1);
    x.at(3) = z.at(2);
    x.at(4) = z.at(3);
    x.at(5) = z.at(4);
    x.at(6) = z.at(5);
    for (unsigned int i=7; i<=N-1; i++)
    {
        double alphai = r(i)*m1;
        double bettai = p(i)*m2;

        double d0 = N6D2[6][0]*alphai + N6D1[6][0]*bettai;
        double d1 = N6D2[6][1]*alphai + N6D1[6][1]*bettai;
        double d2 = N6D2[6][2]*alphai + N6D1[6][2]*bettai;
        double d3 = N6D2[6][3]*alphai + N6D1[6][3]*bettai;
        double d4 = N6D2[6][4]*alphai + N6D1[6][4]*bettai;
        double d5 = N6D2[6][5]*alphai + N6D1[6][5]*bettai;
        double d6 = N6D2[6][6]*alphai + N6D1[6][6]*bettai + q(i);

        x.at(i) = -d0*x.at(i-6) - d1*x.at(i-5) - d2*x.at(i-4) - d3*x.at(i-3) - d4*x.at(i-2) - d5*x.at(i-1) + f(i);
        x.at(i) /= d6;
    }

    A.clear();
    b.clear();
    z.clear();
}

void LinearBoundaryValueProblemODE::calculate(DoubleVector &x, unsigned int k, double h, unsigned int N, unsigned int direction UNUSED_PARAM)
{
    C_UNUSED(x);
    C_UNUSED(k);
    C_UNUSED(h);
    C_UNUSED(N);
    C_UNUSED(direction);

    //    x.clear();
    //    x.resize(N+1, 0.0);

    //    double p1[7] = {0.0, 0.0, 2.0, 0.0, 12.0, 0.0, 60.0};
    //    double p2[7] = {0.0, 0.0, 1.0, 0.0, 24.0, 0.0, 180.0};

    //    double (*E2)[k+1][k+1] = NULL;
    //    double (*E1)[k+1][k+1] = NULL;

    //    if (k==4)
    //    {
    //        E2 = &N4D2;
    //        E1 = &N4D1;
    //    }
    //    if (k==6)
    //    {
    //        E2 = &N6D2;
    //        E1 = &N6D1;
    //    }

    //    double m1 = p1[k]*h;
    //    double m2 = p2[k]*h*h;

    //    DoubleMatrix A(k, k, 0.0);
    //    DoubleVector b(k, 0.0);
    //    DoubleVector z(k, 0.0);
    //    DoubleMatrix ems(N-k, k);

    //    /* border conditions */
    //    x.at(0) = boundary(Left);
    //    x.at(N) = boundary(Right);

    //    double alpha0 = r(1)/m2;
    //    double betta0 = p(1)/m1;
    //    for (unsigned int j=1; j<=k; j++)
    //    {
    //        A[0][j-1] = (*E2)[1][j]*alpha0 + (*E1)[1][j]*betta0;
    //    }
    //    b[0] = f(1) - ((*E2)[1][0]*alpha0 + (*E1)[1][0]*betta0)*x.at(0);
    //    A[0][0] += q(1);

    //    for (unsigned int j=1; j<=k-1; j++)
    //    {
    //        A[0][j] /= A[0][0];
    //    }
    //    b[0]    /= A[0][0];
    //    A[0][0] = 1.0;

    //    for (unsigned int j=0; j<=k-2; j++)
    //    {
    //        ems.at(0,j) = A[0][j+1];
    //    }
    //    ems.at(0,k-1) = b[0];

    //    for (unsigned int n=1; n<=N-(k+1); n++)
    //    {
    //        // + * * * *
    //        double alphai = r(n)/m2;
    //        double bettai = p(n)/m1;

    //        double g[k+1] = {0.0};
    //        for (unsigned int i=0; i<=k; i++)
    //        {
    //            g[i] = (*E2)[0][i]*alphai + (*E1)[0][i]*bettai;
    //        }
    //        g[0] += q(n);
    //        double fi = f(n);

    //        for (unsigned int i=1; i<=k; i++)
    //        {
    //            g[i] /= -g[0];
    //        }
    //        fi   /= +g[0];
    //        g[0] = 1.0;

    //        for (unsigned int i=0; i<k-1; i++)
    //        {
    //            A[0][i] = A[0][i+1] + g[i+1];
    //        }
    //        A[0][k-1] = g[k];
    //        b[0]    = b[0] - fi;

    //        for (unsigned int j=1; j<=k-1; j++)
    //        {
    //            A[0][j] /= A[0][0];
    //        }
    //        b[0]    /= A[0][0];
    //        A[0][0] = 1.0;

    //        for (unsigned int j=0; j<k-1; j++)
    //        {
    //            ems.at(n,j) = A[0][j+1];
    //        }
    //        ems.at(n,k-1) = b[0];
    //    }

    //    for (unsigned int i=1; i<=k-1; i++)
    //    {
    //        double alpha1 = r(N-(k-i))/m2;
    //        double betta1 = p(N-(k-i))/m1;
    //        for (unsigned int j=0; j<=(k-1); j++)
    //        {
    //            A[i][j] = (*E2)[i][j]*alpha1 + (*E1)[i][j]*betta1;
    //        }
    //        b[i] = f(N-(k-i)) - ((*E2)[i][k]*alpha1 + (*E1)[i][k]*betta1)*x.at(N);
    //        A[i][i] += q(N-(k-i));
    //    }


    //    GaussianElimination(A, b, z);

    //    for (unsigned int i=1; i<=k; i++)
    //    {
    //        x.at(N-i) = z.at(k-i);
    //    }
    //    printf("%14.10f %14.10f %14.10f %14.10f\n", x.at(N-4), x.at(N-3), x.at(N-2), x.at(N-1));

    //    for (unsigned int n=N-(k+1); n>=1; n--)
    //    {
    //        x.at(n) = ems.at(n-1,k-1);
    //        for (unsigned int i=0; i<k-1; i++)
    //        x.at(n) += -ems.at(n-1,i)*x.at(n+i+1);// - ems.at(n-1,1)*x.at(n+2) - ems.at(n-1,2)*x.at(n+3) + ems.at(n-1,3);
    //    }

    //    //    for (unsigned int n=N-(k+1); n>=1; n--)
    //    //    {
    //    //        double alphai = r(n)/m2;
    //    //        double bettai = p(n)/m1;
    //    //        double d0 = +70.0*alphai  - 25.0*bettai + q(n);
    //    //        double d1 = (-208.0*alphai + 48.0*bettai);
    //    //        double d2 = (+228.0*alphai - 36.0*bettai);
    //    //        double d3 = (-112.0*alphai + 16.0*bettai);
    //    //        double d4 = (+22.0*alphai  - 3.0*bettai);
    //    //        double fi = f(n);
    //    //        x.at(n) = -d1*x.at(n+1) - d2*x.at(n+2) - d3*x.at(n+3) - d4*x.at(n+4) + fi;
    //    //        x.at(n) /= d0;
    //    //    }

    //    ems.clear();
    //    z.clear();
    //    b.clear();
    //    A.clear();
}
