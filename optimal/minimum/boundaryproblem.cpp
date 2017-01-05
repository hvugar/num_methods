#include "boundaryproblem.h"
#include "cmethods.h"

void BoundaryProblem::calculate2N(DoubleVector &x, double h, unsigned int N)
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

void BoundaryProblem::calculate4NL2R(DoubleVector &x, double h, unsigned int N)
{
    x.clear();
    x.resize(N+1, 0.0);
    x.at(0) = boundary(Left);
    x.at(N) = boundary(Right);

    DoubleMatrix A(4,4);
    DoubleVector b(4);
    DoubleVector z(4);

    double alpha0 = r(1)/(24.0*h*h);
    double betta0 = p(1)/(12.0*h);
    double c1 = -40.0*alpha0 - 10.0*betta0 + q(1);
    double c2 = +12.0*alpha0 + 18.0*betta0;
    double c3 = +8.0*alpha0  - 6.0*betta0;
    double c4 = -2.0*alpha0  + betta0;
    double c5 = f(1) - (22.0*alpha0 - 3.0*betta0)*x.at(0);

    c2 /= c1;
    c3 /= c1;
    c4 /= c1;
    c5 /= c1;
    c1 = 1.0;

    for (unsigned int i=2; i<=N-4; i++)
    {
        double alphai = r(i)/(24.0*h*h);
        double bettai = p(i)/(12.0*h);
        double q1 = +22.0*alphai - 3.0*bettai;
        double q2 = +40.0*alphai + 10.0*bettai - q(i);
        double q3 = -12.0*alphai - 18.0*bettai;
        double q4 = -8.0*alphai  + 6.0*bettai;
        double q5 = +2.0*alphai - bettai;
        double e  = f(i);

        q2 /= q1;
        q3 /= q1;
        q4 /= q1;
        q5 /= q1;
        e  /= q1;
        q1 = 1.0;

        c1 = c2 + q2;
        c2 = (c3 + q3)/c1;
        c3 = (c4 + q4)/c1;
        c4 = q5/c1;
        c5 = (c5 - e)/c1;
        c1 = 1.0;
    }

    A[0][0] = c1;
    A[0][1] = c2;
    A[0][2] = c3;
    A[0][3] = c4;
    b[0]    = c5;

    double alpha1 = r(N-3)/(24.0*h*h);
    double betta1 = p(N-3)/(12.0*h);
    A[1][0] = +22.0*alpha1 - 3.0*betta1;
    A[1][1] = -40.0*alpha1 - 10.0*betta1 + q(N-3);
    A[1][2] = +12.0*alpha1 + 18.0*betta1;
    A[1][3] = +8.0*alpha1  - 6.0*betta1;
    b[1]    = f(N-3) - (-2.0*alpha1 + betta1)*x.at(N);

    double alpha2 = r(N-2)/(24.0*h*h);
    double betta2 = p(N-2)/(12.0*h);
    A[2][0] = -2.0*alpha2  + betta2;
    A[2][1] = +32.0*alpha2 - 8.0*betta2;
    A[2][2] = -60.0*alpha2 + q(N-2);
    A[2][3] = +32.0*alpha2 + 8.0*betta2;
    b[2]    = f(N-2) - (-2.0*alpha2 - betta2)*x.at(N);

    double alpha3 = r(N-1)/(24.0*h*h);
    double betta3 = p(N-1)/(12.0*h);
    A[3][0] = -2.0*alpha3  - betta3;
    A[3][1] = +8.0*alpha3  + 6.0*betta3;
    A[3][2] = +12.0*alpha3 - 18.0*betta3;
    A[3][3] = -40.0*alpha3 + 10.0*betta3 + q(N-2);
    b[3]    = f(N-1) - (22.0*alpha3 + 3.0*betta3)*x.at(N);

    GaussianElimination(A, b, z);

    printf("%18.15f %18.15f %18.15f %18.15f\n", z[0], z[1], z[2], z[3]);
    //printf("%14.10f %14.10f %14.10f %14.10f %14.10f\n", fx(N-4), fx(N-3), fx(N-2), fx(N-1), fx(N));

    x.at(N-4) = z.at(0);
    x.at(N-3) = z.at(1);
    x.at(N-2) = z.at(2);
    x.at(N-1) = z.at(3);
    for (unsigned int i=N-5; i>=1; i--)
    {
        double d0 = +(70.0* r(i))/(24.0*h*h) - (25.0*p(i))/(12.0*h) + q(i);
        double d1 = +(208.0*r(i))/(24.0*h*h) - (48.0*p(i))/(12.0*h);
        double d2 = -(228.0*r(i))/(24.0*h*h) + (36.0*p(i))/(12.0*h);
        double d3 = +(112.0*r(i))/(24.0*h*h) - (16.0*p(i))/(12.0*h);
        double d4 = -(22.0* r(i))/(24.0*h*h) + (3.0* p(i))/(12.0*h);

        x.at(i) = d1*x.at(i+1) + d2*x.at(i+2) + d3*x.at(i+3) + d4*x.at(i+4) + f(i);
        //x.at(i) = d1*fx(i+1)   + d2*fx(i+2)   + d3*fx(i+3)   + d4*fx(i+4)   + f(i);
        x.at(i) /= d0;
    }

    A.clear();
    b.clear();
    z.clear();
}

void BoundaryProblem::calculate4NR2L(DoubleVector &x, double h, unsigned int N)
{
    x.clear();
    x.resize(N+1, 0.0);
    x.at(0) = boundary(Left);
    x.at(N) = boundary(Right);

    DoubleMatrix A(4,4);
    DoubleVector b(4);
    DoubleVector z(4);

    double alpha0 = r(1)/(24.0*h*h);
    double betta0 = p(1)/(12.0*h);
    A[0][0] = -40.0*alpha0 - 10.0*betta0 + q(1);
    A[0][1] = +12.0*alpha0 + 18.0*betta0;
    A[0][2] = +8.0*alpha0  - 6.0*betta0;
    A[0][3] = -2.0*alpha0  + betta0;
    b[0]    = f(1) - (22.0*alpha0 - 3.0*betta0)*x.at(0);

    double alpha1 = r(2)/(24.0*h*h);
    double betta1 = p(2)/(12.0*h);
    A[1][0] = +32.0*alpha1 - 8.0*betta1;
    A[1][1] = -60.0*alpha1 + q(2);
    A[1][2] = +32.0*alpha1 + 8.0*betta1;
    A[1][3] = -2.0*alpha1  - betta1;
    b[1]    = f(2) - (-2.0*alpha1 + betta1)*x.at(0);

    double alpha2 = r(3)/(24.0*h*h);
    double betta2 = p(3)/(12.0*h);
    A[2][0] = +8.0*alpha2  + 6.0*betta2;
    A[2][1] = +12.0*alpha2 - 18.0*betta2;
    A[2][2] = -40.0*alpha2 + 10.0*betta2 + q(3);
    A[2][3] = +22.0*alpha2 + 3.0*betta2;
    b[2]    = f(3) - (-2.0*alpha2 - betta2)*x.at(0);

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

        double q1 = -2.0*alphai  - bettai;
        double q2 = +8.0*alphai  + 6.0*bettai;
        double q3 = +12.0*alphai - 18.0*bettai;
        double q4 = -40.0*alphai + 10.0*bettai + q(i);
        double q5 = +22.0*alphai + 3.0*bettai;
        double q  = f(i);

        q4 /= -q5;
        q3 /= -q5;
        q2 /= -q5;
        q1 /= -q5;
        q  /= +q5;
        q5 = 1.0;

        c4 = c3 + q4;
        c3 = (c2 + q3)/c4;
        c2 = (c1 + q2)/c4;
        c1 = q1/c4;
        c  = (c - q)/c4;
        c4 = 1.0;
    }

    A[3][0] = c1;
    A[3][1] = c2;
    A[3][2] = c3;
    A[3][3] = c4;
    b[3]    = c;

    GaussianElimination(A, b, z);

    printf("%18.15f %18.15f %18.15f %18.15f\n", z[0], z[1], z[2], z[3]);
    //printf("%14.10f %14.10f %14.10f %14.10f %14.10f\n", fx(1), fx(2), fx(3), fx(4), fx(0));

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
        //x.at(i) = d1*fx(i-1)   + d2*fx(i-2)   + d3*fx(i-3)   + d4*fx(i-4)   + f(i);
        x.at(i) /= d0;
    }

    A.clear();
    b.clear();
    z.clear();
}

void BoundaryProblem::calculate6NL2R(DoubleVector &x, double h, unsigned int N)
{
    x.clear();
    x.resize(N+1, 0.0);
    x.at(0) = boundary(Left);
    x.at(N) = boundary(Right);

    DoubleMatrix A(6,6);
    DoubleVector b(6);
    DoubleVector z(6);

    {
        double alpha0 = r(1)/(180.0*h*h);
        double betta0 = p(1)/(60.0*h);
        double c1 = -147.0*alpha0 - 77.0*betta0 + q(1);
        double c2 = -225.0*alpha0 + 150.0*betta0;
        double c3 = +470.0*alpha0 - 100.0*betta0;
        double c4 = -285.0*alpha0 + 50.0*betta0;
        double c5 = +93.0*alpha0  - 15.0*betta0;
        double c6 = -13.0*alpha0  + 2.0*betta0;
        double d  = f(1) - (137.0*alpha0 - 10.0*betta0)*x.at(0);

        c2 /= c1;
        c3 /= c1;
        c4 /= c1;
        c5 /= c1;
        c6 /= c1;
        d  /= c1;
        c1 = 1.0;

        for (unsigned int i=2; i<=N-6; i++)
        {
            double alphai = r(i)/(180.0*h*h);
            double bettai = p(i)/(60.0*h);

            double q1 = +137.0*alphai - 10.0*bettai;
            double q2 = -147.0*alphai - 77.0*bettai + q(i);
            double q3 = -255.0*alphai + 150.0*bettai;
            double q4 = +470.0*alphai - 100.0*bettai;
            double q5 = -285.0*alphai + 50.0*bettai;
            double q6 = +93.0*alphai  - 15.0*bettai;
            double q7 = -13.0*alphai  + 2.0*bettai;
            double e  = f(i);

            q2 /= -q1;
            q3 /= -q1;
            q4 /= -q1;
            q5 /= -q1;
            q6 /= -q1;
            q7 /= -q1;
            e  /= +q1;
            q1 = 1.0;

            c1 = c2 + q2;
            c2 = (c3 + q3)/c1;
            c3 = (c4 + q4)/c1;
            c4 = (c5 + q5)/c1;
            c5 = (c6 + q6)/c1;
            c6 = q7/c1;
            d  = (d - e)/c1;
            c1 = 1.0;
        }

        A[0][0] = c1;
        A[0][1] = c2;
        A[0][2] = c3;
        A[0][3] = c4;
        A[0][4] = c5;
        A[0][5] = c6;
        b[0]    = d;
    }
    {
        double alpha1 = r(N-5)/(180.0*h*h);
        double betta1 = p(N-5)/(60.0*h);
        A[1][0] = +137.0*alpha1 - 10.0*betta1;
        A[1][1] = -147.0*alpha1 - 77.0*betta1 + q(N-5);
        A[1][2] = -255.0*alpha1 + 150.0*betta1;
        A[1][3] = +470.0*alpha1 - 100.0*betta1;
        A[1][4] = -285.0*alpha1 + 50.0*betta1;
        A[1][5] = +93.0*alpha1  - 15.0*betta1;
        b[1]    = f(N-5) - (-13.0*alpha1 + 2.0*betta1)*x.at(N);
    }
    {
        double alpha2 = r(N-4)/(180.0*h*h);
        double betta2 = p(N-4)/(60.0*h);
        A[2][0] = -13.0*alpha2  + 2.0*betta2;
        A[2][1] = +228.0*alpha2 - 24.0*betta2;
        A[2][2] = -420.0*alpha2 - 35.0*betta2 + q(N-4);
        A[2][3] = +200.0*alpha2 + 80.0*betta2;
        A[2][4] = +15.0*alpha2  - 30.0*betta2;
        A[2][5] = -12.0*alpha2  + 8.0*betta2;
        b[2]    = f(N-4) - (2.0*alpha2 - betta2)*x.at(N);
    }
    {
        double alpha3 = r(N-3)/(180.0*h*h);
        double betta3 = p(N-3)/(60.0*h);
        A[3][0] = +2.0*alpha3   - betta3;
        A[3][1] = -27.0*alpha3  + 9.0*betta3;
        A[3][2] = +270.0*alpha3 - 45.0*betta3;
        A[3][3] = -490.0*alpha3 + q(N-3);
        A[3][4] = +270.0*alpha3 + 45.0*betta3;
        A[3][5] = -27.0*alpha3  - 9.0*betta3;
        b[3]    = f(N-3) - (2.0*alpha3 + betta3)*x.at(N);
    }
    {
        double alpha4 = r(N-2)/(180.0*h*h);
        double betta4 = p(N-2)/(60.0*h);
        A[4][0] = +2.0*alpha4   + betta4;
        A[4][1] = -12.0*alpha4  - 8.0*betta4;
        A[4][2] = +15.0*alpha4  + 30.0*betta4;
        A[4][3] = +200.0*alpha4 - 80.0*betta4;
        A[4][4] = -420.0*alpha4 + 35.0*betta4 + q(N-2);
        A[4][5] = +228.0*alpha4 + 24.0*betta4;
        b[4]    = f(N-2) - (-13.0*alpha4 - 2.0*betta4)*x.at(N);
    }
    {
        double alpha5 = r(N-1)/(180.0*h*h);
        double betta5 = p(N-1)/(60.0*h);
        A[5][0] = -13.0*alpha5  - 2.0*betta5;
        A[5][1] = +93.0*alpha5  + 15.0*betta5;
        A[5][2] = -285.0*alpha5 - 50.0*betta5;
        A[5][3] = +470.0*alpha5 + 100.0*betta5;
        A[5][4] = -225.0*alpha5 - 150.0*betta5;
        A[5][5] = -147.0*alpha5 + 77.0*betta5 + q(N-1);
        b[5]    = f(N-1) - (137.0*alpha5 + 10.0*betta5)*x.at(N);
    }

    IPrinter::printSeperatorLine();
    IPrinter::print(A,A.rows(),A.cols(), 20, 10);
    IPrinter::print(b, b.size(), 20, 10);
    IPrinter::printSeperatorLine();

    GaussianElimination(A, b, z);

    printf("%14.10f %14.10f %14.10f %14.10f %14.10f %14.10f\n", z[0], z[1], z[2], z[3], z[4], z[5]);
    //printf("%14.10f %14.10f %14.10f %14.10f %14.10f\n", fx(N-4), fx(N-3), fx(N-2), fx(N-1), fx(N));

    //    x.at(N-4) = z.at(0);
    //    x.at(N-3) = z.at(1);
    //    x.at(N-2) = z.at(2);
    //    x.at(N-1) = z.at(3);
    //    for (unsigned int i=N-5; i>=1; i--)
    //    {
    //        double d0 = +(70.0* r(i))/(24.0*h*h) - (25.0*p(i))/(12.0*h) + q(i);
    //        double d1 = +(208.0*r(i))/(24.0*h*h) - (48.0*p(i))/(12.0*h);
    //        double d2 = -(228.0*r(i))/(24.0*h*h) + (36.0*p(i))/(12.0*h);
    //        double d3 = +(112.0*r(i))/(24.0*h*h) - (16.0*p(i))/(12.0*h);
    //        double d4 = -(22.0* r(i))/(24.0*h*h) + (3.0* p(i))/(12.0*h);

    //        x.at(i) = d1*x.at(i+1) + d2*x.at(i+2) + d3*x.at(i+3) + d4*x.at(i+4) + f(i);
    //        //x.at(i) = d1*fx(i+1)   + d2*fx(i+2)   + d3*fx(i+3)   + d4*fx(i+4)   + f(i);
    //        x.at(i) /= d0;
    //    }

    A.clear();
    b.clear();
    z.clear();
}
