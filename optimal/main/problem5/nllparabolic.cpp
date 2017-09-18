#include "nllparabolic.h"

void NLLIParabolicIBVP::Main(int argc, char *argv[])
{
    NLLIParabolicIBVP nllp;
}

NLLIParabolicIBVP::NLLIParabolicIBVP()
{
    pu = NULL;
    grid.setTimeDimension(Dimension(0.1, 0, 10));
    grid.addSpaceDimension(Dimension(0.1, 0, 10));

    DoubleMatrix m;
    solveEquationM4(m, 1.0);
    puts("end");
}

NLLIParabolicIBVP::~NLLIParabolicIBVP()
{}

void NLLIParabolicIBVP::solveEquation(DoubleMatrix &u, double a)
{
    double hx = grid.spaceDimension(Dimension::DimensionX).step();
    double ht = grid.timeDimension().step();
    int N = grid.spaceDimension(Dimension::DimensionX).sizeN();
    int M = grid.timeDimension().sizeN();

    u.clear();
    u.resize(M+1, N+1);

    pu = &u;

    for (int n=0; n<=N; n++)
    {
        SpaceNodePDE sn; sn.x = n*hx; sn.i = n;
        u(0,n) = initial(sn);
    }

    DoubleVector u0(N-1);
    DoubleVector ru(N-1, 0.0);
    for (int m=1; m<=1; m++)
    {
        TimeNodePDE tn; tn.t = m*ht; tn.i = m;
        cur_m = m;

        for (int n=0; n<N-1; n++) u0[n] = u(m-1, n+1);

        calculateNewtonMethod(u0, ru, 0.01, 0.0001);

        SpaceNodePDE sn0; sn0.x = 0*hx; sn0.i = 0;
        u(m,0) = boundary(sn0, tn, BoundaryType::Left);

        SpaceNodePDE snN; snN.x = N*hx; sn0.i = N;
        u(m,N) = boundary(snN, tn, BoundaryType::Right);

        for (int n=0; n<N-1; n++) u(m,n+1) = ru[n];

    }
    ru.clear();
    u0.clear();

    //    IPrinter::printMatrix(u);
    IPrinter::printVector(u.row(0));
    IPrinter::printVector(u.row(1));
}

void NLLIParabolicIBVP::solveEquationM1(DoubleMatrix &u, double a)
{
    double hx = grid.spaceDimension(Dimension::DimensionX).step();
    double ht = grid.timeDimension().step();
    int N = grid.spaceDimension(Dimension::DimensionX).sizeN();
    int M = grid.timeDimension().sizeN();

    u.clear();
    u.resize(M+1, N+1);

    for (int n=0; n<=N; n++)
    {
        SpaceNodePDE sn; sn.x = n*hx; sn.i = n;
        u(0,n) = initial(sn);
    }

    SpaceNodePDE node1; node1.x = 2*hx; node1.i = 2;
    SpaceNodePDE node2; node2.x = 5*hx; node2.i = 5;

    SpaceNodePDE sn0; sn0.x = 0*hx; sn0.i = 0;
    SpaceNodePDE snN; snN.x = N*hx; snN.i = N;

    for (int m=1; m<=M; m++)
    {
        TimeNodePDE tn0; tn0.t = (m+0)*ht; tn0.i = m+0;
        TimeNodePDE tn1; tn1.t = (m-1)*ht; tn1.i = m-1;

        double alpha = (a*a)/(hx*hx);
        double betta = 1.0/ht - 2.0*alpha;
        for (int n=0; n<=N; n++)
        {
            SpaceNodePDE sn; sn.x = n*hx; sn.i = n;

            if (n==0)
                u(m,n) = boundary(sn0, tn0, Left);
            else if (n==N)
                u(m,n) = boundary(snN, tn0, Right);
            else
                u(m,n) = (alpha*u(m-1,n-1) + betta*u(m-1,n) + alpha*u(m-1,n+1) + f(sn, tn1)
                          + g(sn,tn1,1)*U(node1, tn1)*U(node1, tn1)
                          + g(sn,tn1,2)*U(node2, tn1)*U(node2, tn1))*ht;

        }
    }

    IPrinter::print(u);
}

void NLLIParabolicIBVP::solveEquationM2(DoubleMatrix &u, double a)
{
    double hx = grid.spaceDimension(Dimension::DimensionX).step();
    double ht = grid.timeDimension().step();
    int N = grid.spaceDimension(Dimension::DimensionX).sizeN();
    int M = grid.timeDimension().sizeN();

    u.clear();
    u.resize(M+1, N+1);

    DoubleVector a1(N-1); for (int n=0; n<N-1; n++) a1[n] = -(a*a)/(hx*hx); a1[0] = 0.0;
    DoubleVector b1(N-1); for (int n=0; n<N-1; n++) b1[n] = +(1.0/ht + 2.0*(a*a)/(hx*hx));
    DoubleVector c1(N-1); for (int n=0; n<N-1; n++) c1[n] = -(a*a)/(hx*hx); c1[N-2] = 0.0;
    DoubleVector d1(N-1, 0.0);

    for (int n=0; n<=N; n++)
    {
        SpaceNodePDE sn; sn.x = n*hx; sn.i = n;
        u(0,n) = initial(sn);
    }

    SpaceNodePDE node1; node1.x = 2*hx; node1.i = 2;
    SpaceNodePDE node2; node2.x = 5*hx; node2.i = 5;

    SpaceNodePDE sn0; sn0.x = 0*hx; sn0.i = 0;
    SpaceNodePDE snN; snN.x = N*hx; snN.i = N;

    for (int m=1; m<=M; m++)
    {
        TimeNodePDE tn0; tn0.t = (m+0)*ht; tn0.i = m+0;
        TimeNodePDE tn1; tn1.t = (m-1)*ht; tn1.i = m-1;

        for (int n=1; n<=N-1; n++)
        {
            SpaceNodePDE sn; sn.x = n*hx; sn.i = n;
            d1[n-1] = (1.0/ht)*u.at(m-1,n) + f(sn,tn0)
                    + g(sn,tn0,1)*U(node1, tn1)*U(node1, tn1)
                    + g(sn,tn0,2)*U(node2, tn1)*U(node2, tn1);

            if (n==1)
            {
                u.at(m,0) = boundary(sn0, tn0, Left);
                d1[n-1] += ((a*a)/(hx*hx))*u.at(m,0);
            }
            if (n==N-1)
            {
                u(m,N) = boundary(snN, tn0, Right);
                d1[n-1] += ((a*a)/(hx*hx))*u.at(m,N);
            }
        }

        DoubleVector rx(N-1);
        tomasAlgorithm(a1.data(), b1.data(), c1.data(), d1.data(), rx.data(), N-1);
        for (int n=1; n<=N-1; n++)
        {
            u.at(m,n) = rx[n-1];
        }
    }

    IPrinter::print(u);
    //IPrinter::print(u.row(0));
    //IPrinter::print(u.row(1));
}

void NLLIParabolicIBVP::solveEquationM4(DoubleMatrix &u, double a)
{
    double hx = grid.spaceDimension(Dimension::DimensionX).step();
    double ht = grid.timeDimension().step();
    int N = grid.spaceDimension(Dimension::DimensionX).sizeN();
    int M = grid.timeDimension().sizeN();
    int L = 2;

    u.clear();
    u.resize(M+1, N+1);

    DoubleVector V(N-1);
    DoubleMatrix W(N-1, L, 0.0);

    // initial condition
    for (int n=0; n<=N; n++)
    {
        SpaceNodePDE sn; sn.x = n*hx; sn.i = n;
        u[0][n] = initial(sn);
    }
    for (int n=1; n<=N-1; n++)
    {
        V[n-1] = u[0][n];

        for (int s=1; s<=L; s++) W.at(n-1,s-1) = 0.0;
    }

    for (int m=1; m<=M; m++)
    {
        TimeNodePDE tn; tn.t = m*ht; tn.i = m;

        /******** Boundary condition ************/
        SpaceNodePDE sn0; sn0.x = 0*hx; sn0.i = 0;
        SpaceNodePDE snN; snN.x = N*hx; snN.i = N;

        u.at(m,0) = boundary(sn0, tn, BoundaryType::Left);
        u.at(m,N) = boundary(snN, tn, BoundaryType::Right);

        /****************************************/

        /////////////////////////////////////////////////////////////
        DoubleVector a1(N-1); for (int n=0; n<N-1; n++) a1[n] = +(a*a)/(hx*hx); a1[0] = 0.0;
        DoubleVector b1(N-1); for (int n=0; n<N-1; n++) b1[n] = -(1.0/ht + 2.0*(a*a)/(hx*hx));
        DoubleVector c1(N-1); for (int n=0; n<N-1; n++) c1[n] = +(a*a)/(hx*hx); c1[N-2] = 0.0;
        DoubleVector d1(N-1); for (int n=0; n<N-1; n++) d1[n] = +0.0;

        /////////////////////////////////////////////////////////
        DoubleMatrix L2(N-1, N-1);
        for (int i=0; i<N-1; i++) for (int j=0; j<N-1; j++) L2.at(i,j) = 0.0;
        for (int n=0; n<N-1; n++)                           L2.at(n,n) = 1.0/ht;

        ///////////////////////////////////////////////////////////
//        DoubleVector F(N-1);

//        for (int n=1; n<=N-1; n++)
//        {
//            SpaceNodePDE sn; sn.x = n*hx; sn.i = n;
//            F.at(n-1) = -f(sn,tn);
//        }

//        DoubleVector d1 = F - (DoubleVector)(L2*V);

//        d1[0]   -= (a*a)/(hx*hx)*u.at(m,0);
//        d1[N-2] -= (a*a)/(hx*hx)*u.at(m,N);

        for (int n=1; n<=N-1; n++)
        {
            SpaceNodePDE sn; sn.x = n*hx; sn.i = n;
            d1[n-1] = -f(sn, tn) - (1.0/ht)*V[n-1];
        }
        d1[0]   -= (a*a)/(hx*hx)*u.at(m,0);
        d1[N-2] -= (a*a)/(hx*hx)*u.at(m,N);

        DoubleVector rx(N-1);
        tomasAlgorithm(a1.data(), b1.data(), c1.data(), d1.data(), rx.data(), N-1);

        //IPrinter::print(rx);
        for (int n=1; n<=N-1; n++) V[n-1] = rx[n-1];

        ////////////////////////////////////////////////////////////////////////

        DoubleMatrix B(N-1, L);
        for (int n=1; n<=N-1; n++)
        {
            SpaceNodePDE sn; sn.x = n*hx; sn.i = n;
            for (int s=1; s<=L; s++) B.at(n-1,s-1) = g(sn, tn, s);
        }

        //IPrinter::printSeperatorLine();
        //IPrinter::print(B);
        //IPrinter::printSeperatorLine();

        DoubleMatrix A = -1.0*B - L2*W;
        //IPrinter::printSeperatorLine();
        //IPrinter::print(A);
        //IPrinter::printSeperatorLine();

        for (int s=1; s<=L; s++)
        {
            for (int n=1; n<=N-1; n++)
                d1.at(n-1) = A.at(n-1,s-1);
            //d1[0]   -= (a*a)/(hx*hx)*u.at(m,0);
            //d1[N-2] -= (a*a)/(hx*hx)*u.at(m,N);

            DoubleVector rx(N-1);
            tomasAlgorithm(a1.data(), b1.data(), c1.data(), d1.data(), rx.data(), N-1);

            for (int n=1; n<=N-1; n++) W.at(n-1,s-1) = rx[n-1];
        }
        A.clear();

        //IPrinter::printSeperatorLine();
        //IPrinter::print(W);
        //IPrinter::printSeperatorLine();

        for (int n=1; n<=N-1; n++)
        {
            SpaceNodePDE node1; node1.x = 2*hx; node1.i = 2;
            SpaceNodePDE node2; node2.x = 5*hx; node2.i = 5;

            u.at(m, n) = V[n-1]
                    + W.at(n-1,0)*U(node1, tn)*U(node1, tn)
                    + W.at(n-1,1)*U(node2, tn)*U(node2, tn);
        }
    }
    IPrinter::print(u);
    //IPrinter::print(u.row(0));
    //IPrinter::print(u.row(1));
}

double NLLIParabolicIBVP::initial(const SpaceNodePDE &sn UNUSED_PARAM) const
{
    return 0.0;
}

double NLLIParabolicIBVP::boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM, BoundaryType boundary) const
{
    double t = tn.t; C_UNUSED(t);

    if (boundary == BoundaryType::Left)  return 0.0;
    if (boundary == BoundaryType::Right) return t;

    return NAN;
}

double NLLIParabolicIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double x = sn.x; C_UNUSED(x);
    double t = tn.t; C_UNUSED(t);

    return x*x - 2.0*t - 0.0081*t*t*g(sn,tn,1) - 0.1296*t*t*g(sn,tn,2);
}

double NLLIParabolicIBVP::g(const SpaceNodePDE &sn, const TimeNodePDE &tn, unsigned int s) const
{
    double x = sn.x; C_UNUSED(x);
    double t = tn.t; C_UNUSED(t);
    if (s==1) return x+t;
    if (s==2) return x*t;

    return NAN;
}

double NLLIParabolicIBVP::U(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double x = sn.x; C_UNUSED(x);
    double t = tn.t; C_UNUSED(t);

    return x*x*t;
}

double NLLIParabolicIBVP::fx(const DoubleVector &x, unsigned int num) const
{
    Dimension spaceDimension = grid.spaceDimension(Dimension::DimensionX);
    Dimension timeDimension = grid.timeDimension();

    int N = spaceDimension.sizeN();
    int M = timeDimension.sizeN();

    double hx = spaceDimension.step();
    double ht = timeDimension.step();
    double a = 1.0;

    double alpha = a*a*ht/(hx*hx);

    int ts[] = { 2, 5 };

    TimeNodePDE tn;
    tn.t = cur_m*ht;
    tn.i = cur_m;

    SpaceNodePDE sn0; sn0.x = 0.0;  sn0.i = 0;
    SpaceNodePDE snN; snN.x = N*hx; snN.i = N;

    if (num == 0)
    {
        SpaceNodePDE sn1; sn1.x = hx;  sn1.i = 1;

        return -(1.0+2.0*alpha)*x[0] + alpha*x[1]
                + ht*g(sn1,tn,1)*x[ts[0]]*x[ts[0]]
                + ht*g(sn1,tn,2)*x[ts[1]]*x[ts[1]]
                + (ht*f(sn1,tn)+ (*pu)(cur_m-1,sn1.i)+alpha*boundary(sn0,tn,BoundaryType::Left));
    }
    else if ((int)num == N-2)
    {
        SpaceNodePDE snN1; snN1.x = (N-1)*hx; snN1.i = N-1;

        return alpha*x[N-3] -(1.0+2.0*alpha)*x[N-2]
                + ht*g(snN1,tn,1)*x[ts[0]]*x[ts[0]]
                + ht*g(snN1,tn,2)*x[ts[1]]*x[ts[1]]
                + (ht*f(snN1,tn)+ (*pu)(cur_m-1,snN1.i)+alpha*boundary(snN,tn,BoundaryType::Right));
    }
    else
    {
        SpaceNodePDE sn; sn.x = (num+1)*hx; sn.i = num+1;

        return  alpha*x[num-1] - (1.0+2.0*alpha)*x[num] + alpha*x[num+1]
                + ht*g(sn,tn,1)*x[ts[0]]*x[ts[0]]
                + ht*g(sn,tn,2)*x[ts[1]]*x[ts[1]]
                + (ht*f(sn,tn)+ (*pu)(cur_m-1,sn.i));
    }

    return NAN;
}


