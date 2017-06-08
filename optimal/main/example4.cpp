#include "example4.h"
#include "rungekutta.h"

void Example4::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Example4 e;
    e.h = 0.01;
    e.N = 100;
    e.F = e.N/10;
    e.n = 3;
    e.K = 4;
    e.w = 12;
    e.p = 8;

    //--------------------------------------------------------------------------
    fputs("Real process solution:\n",stdout);
    IPrinter::printSeperatorLine(NULL,'-', stdout);
    DoubleMatrix rx;
    e.calculateRX(rx);
    IPrinter::printVector(e.w,e.p,rx.row(0),"x1: ");
    IPrinter::printVector(e.w,e.p,rx.row(1),"x2: ");
    IPrinter::printVector(e.w,e.p,rx.row(2),"x3: ");
    IPrinter::printSeperatorLine(NULL,'-', stdout);

    //--------------------------------------------------------------------------

    std::vector<unsigned int> *s = new std::vector<unsigned int>[e.K];
    //    if (e.K==2)
    //    {
    //        s[0].push_back(0); s[0].push_back(2*e.F); s[0].push_back(5*e.F); s[0].push_back(8*e.F); s[0].push_back(9*e.F); s[0].push_back(10*e.F);
    //        s[1].push_back(0); s[1].push_back(1);     s[1].push_back(2);
    //    }
    if (e.K==4)
    {
        s[0].push_back(0); /*s[0].push_back(5*e.F);*/ s[0].push_back(10*e.F);
        s[1].push_back(0); s[1].push_back(1);     s[1].push_back(2);     s[1].push_back(3);     s[1].push_back(4);
        s[2].push_back(0); s[2].push_back(1);     s[2].push_back(2);     s[2].push_back(3);     s[2].push_back(4);
        s[3].push_back(0); s[3].push_back(1);     s[3].push_back(2);     s[3].push_back(3);     s[3].push_back(4);
    }
    //    if (e.K==6)
    //    {
    //        s[0].push_back(0); s[0].push_back(2*e.F); s[0].push_back(5*e.F); s[0].push_back(6*e.F); s[0].push_back(7*e.F); s[0].push_back(8*e.F);s[0].push_back(999); s[0].push_back(10*e.F);
    //        s[1].push_back(0); s[1].push_back(1);     s[1].push_back(2);     s[1].push_back(3);     s[1].push_back(4);     s[1].push_back(5);    s[1].push_back(6);
    //        s[2].push_back(0); s[2].push_back(1);     s[2].push_back(2);     s[2].push_back(3);     s[2].push_back(4);     s[2].push_back(5);    s[2].push_back(6);
    //        s[3].push_back(0); s[3].push_back(1);     s[3].push_back(2);     s[3].push_back(3);     s[3].push_back(4);     s[3].push_back(5);    s[3].push_back(6);
    //    }

    //        IPrinter::printSeperatorLine(NULL,'-',stdout);
    //        DoubleVector x(e.n*e.K);
    //        for (unsigned int i=0; i<e.K; i++)
    //        {
    //            for (unsigned int j=0; j<e.n; j++)
    //            {
    //                x.at(i*e.n+j) = e.fx(j+1,i);
    //            }
    //        }
    //        IPrinter::print(x,x.size(),e.w,e.p,stdout);
    IPrinter::printSeperatorLine(NULL,'-',stdout);

    //    {
    //        class CauchyProblem1 : public CauchyProblem
    //        {
    //        public:
    //            CauchyProblem1() {}

    //            virtual double f(double t, const DoubleVector &x) const
    //            {
    //                return 2.0*x[0] + 3.0*x[1] - x[2]
    //                        - (2.0*sin(20.0*t*t) + 3.0*(cos(10.0*t) - sin(20.0*t)) - (t*t*t-sin(8.0*t)*sin(8.0*t)))
    //                        + 40.0*t*cos(20.0*t*t);
    //            }
    //        };

    //        class CauchyProblem2 : public CauchyProblem
    //        {
    //        public:
    //            CauchyProblem2() {}

    //            virtual double f(double t, const DoubleVector &x) const
    //            {
    //                return 4.0*x[0] + 6.0*x[1] - 2.0*x[2]
    //                        - (4.0*sin(20.0*t*t) + 6.0*(cos(10.0*t) - sin(20.0*t)) - 2.0*(t*t*t-sin(8.0*t)*sin(8.0*t)))
    //                        - 10.0*sin(10.0*t) - 20.0*cos(20.0*t);
    //            }
    //        };

    //        class CauchyProblem3 : public CauchyProblem
    //        {
    //        public:
    //            CauchyProblem3() {}

    //            virtual double f(double t, const DoubleVector &x) const
    //            {
    //                return -1.0*x[0] + 1.0*x[1] - 1.0*x[2]
    //                        - (-sin(20.0*t*t) + (cos(10.0*t) - sin(20.0*t)) -(t*t*t-sin(8.0*t)*sin(8.0*t)))
    //                        + 3.0*t*t - 16.0*cos(8.0*t)*sin(8.0*t);
    //            }
    //        };

    //        std::vector<CauchyProblem*> ps;

    //        CauchyProblem1 *cp1 = new CauchyProblem1;
    //        cp1->y0 = e.fx(1,0);
    //        cp1->x0 = 0.0;
    //        ps.push_back(cp1);

    //        CauchyProblem2 *cp2 = new CauchyProblem2;
    //        cp2->y0 = e.fx(2,0);
    //        cp2->x0 = 0.0;
    //        ps.push_back(cp2);

    //        CauchyProblem3 *cp3 = new CauchyProblem3;
    //        cp3->y0 = e.fx(3,0);
    //        cp3->x0 = 0.0;
    //        ps.push_back(cp3);

    //        DoubleMatrix mx;
    //        CauchyProblem::rungeKutta(ps, 0.0, e.h, e.N, mx);
    //        IPrinter::printVector(e.w,e.p,mx.row(0),"x1: ");
    //        IPrinter::printVector(e.w,e.p,mx.row(1),"x2: ");
    //        IPrinter::printVector(e.w,e.p,mx.row(2),"x3: ");
    //    }


    //    {
    //        //--------------------------------------------------------------------------
    //        e.K = 1;
    //        printf("Initial first 3 points: K: %d \n", e.K);
    //        DoubleMatrix nx1;
    //        e.calculateNX(rx,nx1);
    //        IPrinter::printVector(e.w,e.p,nx1.row(0),"x1: ");
    //        IPrinter::printVector(e.w,e.p,nx1.row(1),"x2: ");
    //        IPrinter::printVector(e.w,e.p,nx1.row(2),"x3: ");
    //        IPrinter::printSeperatorLine(NULL,'-', stdout);
    //    }

    //    {
    //        //--------------------------------------------------------------------------
    //        e.K = 2;
    //        printf("Initial first 3 points: K: %d \n", e.K);
    //        DoubleMatrix nx1;
    //        e.calculateNX(rx,nx1);
    //        IPrinter::printVector(e.w,e.p,nx1.row(0),"x1: ");
    //        IPrinter::printVector(e.w,e.p,nx1.row(1),"x2: ");
    //        IPrinter::printVector(e.w,e.p,nx1.row(2),"x3: ");
    //        IPrinter::printSeperatorLine(NULL,'-', stdout);

    //        //        //--------------------------------------------------------------------------
    //        //        puts("Method #1");
    //        //        IPrinter::printSeperatorLine(NULL,'-', stdout);
    //        //        e.calculateM1(s, rx, nx1);
    //        //        double norm1 = 0.0; for (unsigned int i=0; i<=e.N; i++) norm1 += (rx.row(0).at(i)-nx1.row(0).at(i))*(rx.row(0).at(i)-nx1.row(0).at(i)); printf("norm1: %f\n", sqrt(norm1));
    //        //        double norm2 = 0.0; for (unsigned int i=0; i<=e.N; i++) norm2 += (rx.row(1).at(i)-nx1.row(1).at(i))*(rx.row(1).at(i)-nx1.row(1).at(i)); printf("norm2: %f\n", sqrt(norm2));
    //        //        double norm3 = 0.0; for (unsigned int i=0; i<=e.N; i++) norm3 += (rx.row(2).at(i)-nx1.row(2).at(i))*(rx.row(2).at(i)-nx1.row(2).at(i)); printf("norm3: %f\n", sqrt(norm3));
    //    }

    //    {
    //--------------------------------------------------------------------------
    e.K = 4;
    printf("Initial first 3 points: K: %d \n", e.K);
    DoubleMatrix nx2;
    e.calculateNX(rx,nx2);
    //        IPrinter::printVector(e.w,e.p,nx2.row(0),"x1: ");
    //        IPrinter::printVector(e.w,e.p,nx2.row(1),"x2: ");
    //        IPrinter::printVector(e.w,e.p,nx2.row(2),"x3: ");
    //        IPrinter::printSeperatorLine(NULL,'-', stdout);
    //        printf("+++ %.8f %.10f %.10f %.10f\n", e.b(1,2), e.a(1,1,2), e.a(1,2,2), e.a(1,3,2));

    //            //--------------------------------------------------------------------------
    puts("Method #1");
    IPrinter::printSeperatorLine(NULL,'-', stdout);
    e.calculateM1(s, rx, nx2);
    //            double norm1 = 0.0; for (unsigned int i=0; i<=e.N; i++) norm1 += (rx.row(0).at(i)-nx2.row(0).at(i))*(rx.row(0).at(i)-nx2.row(0).at(i)); printf("norm1: %f\n", sqrt(norm1));
    //            double norm2 = 0.0; for (unsigned int i=0; i<=e.N; i++) norm2 += (rx.row(1).at(i)-nx2.row(1).at(i))*(rx.row(1).at(i)-nx2.row(1).at(i)); printf("norm2: %f\n", sqrt(norm2));
    //            double norm3 = 0.0; for (unsigned int i=0; i<=e.N; i++) norm3 += (rx.row(2).at(i)-nx2.row(2).at(i))*(rx.row(2).at(i)-nx2.row(2).at(i)); printf("norm3: %f\n", sqrt(norm3));
    //    }

    //    {
    //        //--------------------------------------------------------------------------
    //        printf("Initial first 3 points: K: %d \n", e.K);
    //        DoubleMatrix nx3;
    //        e.calculateNX(rx,nx3);
    //        IPrinter::printVector(e.w,e.p,nx3.row(0),"x1: ");
    //        IPrinter::printVector(e.w,e.p,nx3.row(1),"x2: ");
    //        IPrinter::printVector(e.w,e.p,nx3.row(2),"x3: ");
    //        IPrinter::printSeperatorLine(NULL,'-', stdout);

    ////        //--------------------------------------------------------------------------
    ////        puts("Method #1");
    ////        IPrinter::printSeperatorLine(NULL,'-', stdout);
    ////        e.calculateM1(s, rx, nx3);
    ////        double norm1 = 0.0; for (unsigned int i=0; i<=e.N; i++) norm1 += (rx.row(0).at(i)-nx3.row(0).at(i))*(rx.row(0).at(i)-nx3.row(0).at(i)); printf("norm1: %f\n", sqrt(norm1));
    ////        double norm2 = 0.0; for (unsigned int i=0; i<=e.N; i++) norm2 += (rx.row(1).at(i)-nx3.row(1).at(i))*(rx.row(1).at(i)-nx3.row(1).at(i)); printf("norm2: %f\n", sqrt(norm2));
    ////        double norm3 = 0.0; for (unsigned int i=0; i<=e.N; i++) norm3 += (rx.row(2).at(i)-nx3.row(2).at(i))*(rx.row(2).at(i)-nx3.row(2).at(i)); printf("norm3: %f\n", sqrt(norm3));
    //    }
    //--------------------------------------------------------------------------
    //    puts("Method #2");
    //    IPrinter::printSeperatorLine(NULL,'-', stdout);
    //    e.calculateM2(s, rx, nx);
}

void Example4::qovmaM1R2L(const std::vector<DoubleMatrix> &gamma, DoubleVector &eta, const std::vector<unsigned int> &s, std::vector<DoubleMatrix> &BETTA)
{
    unsigned int L = s.size();
    std::vector<DoubleMatrix> betta(N+1);

    unsigned int cur = L-1;
    unsigned int start = s[cur];
    for (unsigned int i=0; i<K; i++)
    {
        if (start-i == s[cur])
        {
            betta[start-i] = gamma[cur--];
        }
        else
        {
            betta[start-i].resize(n,n,0.0);
        }
    }
    stdDoubleMatrixVector A;
    initAMatrices(A);
    for (unsigned int k=start; k>=K; k--)
    {
        updateAMatrices(A,k);
        eta = eta - betta[k]*A[0];
        for (unsigned int i=1; i<K; i++)
        {
            betta[k-i] = betta[k]*A[i] + betta[k-i];
        }
        betta[k-K] = betta[k]*A[K];
        if ( k == (s[cur]+K) )
        {
            betta[k-K] = betta[k-K] + gamma[cur];
            cur--;
        }
    }
    clearAMatrices(A);
    BETTA.clear();
    BETTA.resize(K);
    for (unsigned int i=0; i<K; i++) BETTA[i] = betta[i];

    for (unsigned int i=0; i<=N; i++) betta.at(i).clear();
    betta.clear();
}

void Example4::initAMatrices(std::vector<DoubleMatrix> &A)
{
    A.clear();
    A.resize(K+1);
    A[0].resize(n,1);
    for (unsigned int i=1; i<=K; i++) A[i].resize(n,n);
}

void Example4::updateAMatrices(std::vector<DoubleMatrix> &A, unsigned int k)
{
    if (K==1)
    {
        A[1].at(0,0) = h*a(1,1,k-1)+1.0; A[1].at(0,1) = h*a(1,2,k-1);      A[1].at(0,2) = h*a(1,3,k-1);
        A[1].at(1,0) = h*a(2,1,k-1);     A[1].at(1,1) = h*a(2,2,k-1)+1.0;  A[1].at(1,2) = h*a(2,3,k-1);
        A[1].at(2,0) = h*a(3,1,k-1);     A[1].at(2,1) = h*a(3,2,k-1);      A[1].at(2,2) = h*a(3,3,k-1)+1.0;

        A[0].at(0,0) = h*b(1,k-1);
        A[0].at(1,0) = h*b(2,k-1);
        A[0].at(2,0) = h*b(3,k-1);
    }

    if (K==2)
    {
        A[1].at(0,0) = (2.0*h*a(1,1,k-1)); A[1].at(0,1) = (2.0*h*a(1,2,k-1));  A[1].at(0,2) = (2.0*h*a(1,3,k-1));
        A[1].at(1,0) = (2.0*h*a(2,1,k-1)); A[1].at(1,1) = (2.0*h*a(2,2,k-1));  A[1].at(1,2) = (2.0*h*a(2,3,k-1));
        A[1].at(2,0) = (2.0*h*a(3,1,k-1)); A[1].at(2,1) = (2.0*h*a(3,2,k-1));  A[1].at(2,2) = (2.0*h*a(3,3,k-1));

        A[2].at(0,0) = +1.0;               A[2].at(0,1) = +0.0;                A[2].at(0,2) = +0.0;
        A[2].at(1,0) = +0.0;               A[2].at(1,1) = +1.0;                A[2].at(1,2) = +0.0;
        A[2].at(2,0) = +0.0;               A[2].at(2,1) = +0.0;                A[2].at(2,2) = +1.0;

        A[0].at(0,0) = (2.0*h*b(1,k-1));
        A[0].at(1,0) = (2.0*h*b(2,k-1));
        A[0].at(2,0) = (2.0*h*b(3,k-1));
    }

    //    if (K==4)
    //    {
    //        A[1].at(0,0) = +16.0/3.0;              A[1].at(0,1) = +0.0;                  A[1].at(0,2) = +0.0;
    //        A[1].at(1,0) = +0.0;                   A[1].at(1,1) = +16.0/3.0;             A[1].at(1,2) = +0.0;
    //        A[1].at(2,0) = +0.0;                   A[1].at(2,1) = +0.0;                  A[1].at(2,2) = +16.0/3.0;

    //        A[2].at(0,0) = -12.0;                  A[2].at(0,1) = +0.0;                  A[2].at(0,2) = +0.0;
    //        A[2].at(1,0) = +0.0;                   A[2].at(1,1) = -12.0;                 A[2].at(1,2) = +0.0;
    //        A[2].at(2,0) = +0.0;                   A[2].at(2,1) = +0.0;                  A[2].at(2,2) = -12.0;

    //        A[3].at(0,0) = +16.0;                  A[3].at(0,1) = +0.0;                  A[3].at(0,2) = +0.0;
    //        A[3].at(1,0) = +0.0;                   A[3].at(1,1) = +16.0;                 A[3].at(1,2) = +0.0;
    //        A[3].at(2,0) = +0.0;                   A[3].at(2,1) = +0.0;                  A[3].at(2,2) = +16.0;

    //        A[4].at(0,0) = -(12.0*h*a(1,1,k-4)+25.0)/3.0; A[4].at(0,1) = -4.0*h*a(1,2,k-4);             A[4].at(0,2) = -4.0*h*a(1,3,k-4);
    //        A[4].at(1,0) = -4.0*h*a(2,1,k-4);             A[4].at(1,1) = -(12.0*h*a(2,2,k-4)+25.0)/3.0; A[4].at(1,2) = -4.0*h*a(2,3,k-4);
    //        A[4].at(2,0) = -4.0*h*a(3,1,k-4);             A[4].at(2,1) = -4.0*h*a(3,2,k-4);             A[4].at(2,2) = -(12.0*h*a(3,3,k-4)+25.0)/3.0;

    //        A[0].at(0,0) = -4.0*h*b(1,k-4);
    //        A[0].at(1,0) = -4.0*h*b(2,k-4);
    //        A[0].at(2,0) = -4.0*h*b(3,k-4);
    //    }

    //    if (K==4)
    //    {
    //        A[1].at(0,0) = +6.0;                   A[1].at(0,1) = +0.0;                  A[1].at(0,2) = +0.0;
    //        A[1].at(1,0) = +0.0;                   A[1].at(1,1) = +6.0;                  A[1].at(1,2) = +0.0;
    //        A[1].at(2,0) = +0.0;                   A[1].at(2,1) = +0.0;                  A[1].at(2,2) = +6.0;

    //        A[2].at(0,0) = -18.0;                  A[2].at(0,1) = +0.0;                  A[2].at(0,2) = +0.0;
    //        A[2].at(1,0) = +0.0;                   A[2].at(1,1) = -18.0;                 A[2].at(1,2) = +0.0;
    //        A[2].at(2,0) = +0.0;                   A[2].at(2,1) = +0.0;                  A[2].at(2,2) = -18.0;

    //        A[3].at(0,0) = 12.0*h*a(1,1,k-3)+10.0; A[3].at(0,1) = 12.0*h*a(1,2,k-3);      A[3].at(0,2) = 12.0*h*a(1,3,k-3);
    //        A[3].at(1,0) = 12.0*h*a(2,1,k-3);      A[3].at(1,1) = 12.0*h*a(2,2,k-3)+10.0; A[3].at(1,2) = 12.0*h*a(2,3,k-3);
    //        A[3].at(2,0) = 12.0*h*a(3,1,k-3);      A[3].at(2,1) = 12.0*h*a(3,2,k-3);      A[3].at(2,2) = 12.0*h*a(3,3,k-3)+10.0;

    //        A[4].at(0,0) = +3.0;                   A[4].at(0,1) = +0.0;                  A[4].at(0,2) = +0.0;
    //        A[4].at(1,0) = +0.0;                   A[4].at(1,1) = +3.0;                  A[4].at(1,2) = +0.0;
    //        A[4].at(2,0) = +0.0;                   A[4].at(2,1) = +0.0;                  A[4].at(2,2) = +3.0;

    //        A[0].at(0,0) = 12.0*h*b(1,k-3);
    //        A[0].at(1,0) = 12.0*h*b(2,k-3);
    //        A[0].at(2,0) = 12.0*h*b(3,k-3);
    //    }

    if (K==4)
    {
        A[1].at(0,0) = +8.0;                   A[1].at(0,1) = +0.0;                  A[1].at(0,2) = +0.0;
        A[1].at(1,0) = +0.0;                   A[1].at(1,1) = +8.0;                  A[1].at(1,2) = +0.0;
        A[1].at(2,0) = +0.0;                   A[1].at(2,1) = +0.0;                  A[1].at(2,2) = +8.0;

        A[2].at(0,0) = -12.0*h*a(1,1,k-2);     A[2].at(0,1) = -12.0*h*a(1,2,k-2);    A[1].at(0,2) = -12.0*h*a(1,3,k-2);
        A[2].at(1,0) = -12.0*h*a(2,1,k-2);     A[2].at(1,1) = -12.0*h*a(2,2,k-2);    A[1].at(1,2) = -12.0*h*a(2,3,k-2);
        A[2].at(2,0) = -12.0*h*a(3,1,k-2);     A[2].at(2,1) = -12.0*h*a(3,2,k-2);    A[1].at(2,2) = -12.0*h*a(3,3,k-2);

        A[3].at(0,0) = -8.0;                   A[3].at(0,1) = +0.0;                  A[3].at(0,2) = +0.0;
        A[3].at(1,0) = +0.0;                   A[3].at(1,1) = -8.0;                  A[3].at(1,2) = +0.0;
        A[3].at(2,0) = +0.0;                   A[3].at(2,1) = +0.0;                  A[3].at(2,2) = -8.0;

        A[4].at(0,0) = +1.0;                   A[4].at(0,1) = +0.0;                  A[4].at(0,2) = +0.0;
        A[4].at(1,0) = +0.0;                   A[4].at(1,1) = +1.0;                  A[4].at(1,2) = +0.0;
        A[4].at(2,0) = +0.0;                   A[4].at(2,1) = +0.0;                  A[4].at(2,2) = +1.0;

        A[0].at(0,0) = -12.0*h*b(1,k-2);
        A[0].at(1,0) = -12.0*h*b(2,k-2);
        A[0].at(2,0) = -12.0*h*b(3,k-2);
    }

    //    if (K==4)
    //    {
    //        A[1].at(0,0) = (12.0*h*a(1,1,k-1)-10.0)/3.0; A[1].at(0,1) = (12.0*h*a(1,2,k-1))/3.0;      A[1].at(0,2) = (12.0*h*a(1,3,k-1))/3.0;
    //        A[1].at(1,0) = (12.0*h*a(2,1,k-1))/3.0;      A[1].at(1,1) = (12.0*h*a(2,2,k-1)-10.0)/3.0; A[1].at(1,2) = (12.0*h*a(2,3,k-1))/3.0;
    //        A[1].at(2,0) = (12.0*h*a(3,1,k-1))/3.0;      A[1].at(2,1) = (12.0*h*a(3,2,k-1))/3.0;      A[1].at(2,2) = (12.0*h*a(3,3,k-1)-10.0)/3.0;

    //        A[2].at(0,0) = +6.00;                  A[2].at(0,1) = +0.00;                  A[2].at(0,2) = +0.00;
    //        A[2].at(1,0) = +0.00;                  A[2].at(1,1) = +6.00;                  A[2].at(1,2) = +0.00;
    //        A[2].at(2,0) = +0.00;                  A[2].at(2,1) = +0.00;                  A[2].at(2,2) = +6.00;

    //        A[3].at(0,0) = -2.00;                  A[3].at(0,1) = +0.00;                  A[3].at(0,2) = +0.00;
    //        A[3].at(1,0) = +0.00;                  A[3].at(1,1) = -2.00;                  A[3].at(1,2) = +0.00;
    //        A[3].at(2,0) = +0.00;                  A[3].at(2,1) = +0.00;                  A[3].at(2,2) = -2.00;

    //        A[4].at(0,0) = +1.0/3.0;               A[4].at(0,1) = +0.00;                  A[4].at(0,2) = +0.00;
    //        A[4].at(1,0) = +0.00;                  A[4].at(1,1) = +1.0/3.0;               A[4].at(1,2) = +0.00;
    //        A[4].at(2,0) = +0.00;                  A[4].at(2,1) = +0.00;                  A[4].at(2,2) = +1.0/3.0;

    //        A[0].at(0,0) = +4.0*h*b(1,k-1);
    //        A[0].at(1,0) = +4.0*h*b(2,k-1);
    //        A[0].at(2,0) = +4.0*h*b(3,k-1);
    //    }

    //    if (K==4)
    //    {
    //        A[1].at(0,0) = 0.48*h*a(1,1,k-0)+1.92; A[1].at(0,1) = 0.48*h*a(1,2,k-0);      A[1].at(0,2) = 0.48*h*a(1,3,k-0);
    //        A[1].at(1,0) = 0.48*h*a(2,1,k-0);      A[1].at(1,1) = 0.48*h*a(2,2,k-0)+1.92; A[1].at(1,2) = 0.48*h*a(2,3,k-0);
    //        A[1].at(2,0) = 0.48*h*a(3,1,k-0);      A[1].at(2,1) = 0.48*h*a(3,2,k-0);      A[1].at(2,2) = 0.48*h*a(3,3,k-0)+1.92;

    //        A[2].at(0,0) = -1.44;                  A[2].at(0,1) = +0.00;                  A[2].at(0,2) = +0.00;
    //        A[2].at(1,0) = +0.00;                  A[2].at(1,1) = -1.44;                  A[2].at(1,2) = +0.00;
    //        A[2].at(2,0) = +0.00;                  A[2].at(2,1) = +0.00;                  A[2].at(2,2) = -1.44;

    //        A[3].at(0,0) = +0.64;                  A[3].at(0,1) = +0.00;                  A[3].at(0,2) = +0.00;
    //        A[3].at(1,0) = +0.00;                  A[3].at(1,1) = +0.64;                  A[3].at(1,2) = +0.00;
    //        A[3].at(2,0) = +0.00;                  A[3].at(2,1) = +0.00;                  A[3].at(2,2) = +0.64;

    //        A[4].at(0,0) = -0.12;                  A[4].at(0,1) = +0.00;                  A[4].at(0,2) = +0.00;
    //        A[4].at(1,0) = +0.00;                  A[4].at(1,1) = -0.12;                  A[4].at(1,2) = +0.00;
    //        A[4].at(2,0) = +0.00;                  A[4].at(2,1) = +0.00;                  A[4].at(2,2) = -0.12;

    //        A[0].at(0,0) = 0.48*h*b(1,k-0);
    //        A[0].at(1,0) = 0.48*h*b(2,k-0);
    //        A[0].at(2,0) = 0.48*h*b(3,k-0);
    //    }

    //    if (K==6)
    //    {
    //        A[1].at(0,0) = (60.0*h*a(1,1,k-1)+360.0)/147.0; A[1].at(0,1) = (60.0*h*a(1,2,k-1))/147.0;      A[1].at(0,2) = (60.0*h*a(1,3,k-1))/147.0;
    //        A[1].at(1,0) = (60.0*h*a(2,1,k-1))/147.0;       A[1].at(1,1) = (60.0*h*a(2,2,k-1)+360.0)/147.0; A[1].at(1,2) = (60.0*h*a(2,3,k-1))/147.0;
    //        A[1].at(2,0) = (60.0*h*a(3,1,k-1))/147.0;       A[1].at(2,1) = (60.0*h*a(3,2,k-1))/147.0;      A[1].at(2,2) = (60.0*h*a(3,3,k-1)+360.0)/147.0;

    //        A[2].at(0,0) = -450.0/147.0;                    A[2].at(0,1) = +0.00;                          A[2].at(0,2) = +0.00;
    //        A[2].at(1,0) = +0.00;                           A[2].at(1,1) = -450.0/147.0;                   A[2].at(1,2) = +0.00;
    //        A[2].at(2,0) = +0.00;                           A[2].at(2,1) = +0.00;                          A[2].at(2,2) = -450.0/147.0;

    //        A[3].at(0,0) = +400.0/147.0;                    A[3].at(0,1) = +0.00;                          A[3].at(0,2) = +0.00;
    //        A[3].at(1,0) = +0.00;                           A[3].at(1,1) = +400.0/147.0;                   A[3].at(1,2) = +0.00;
    //        A[3].at(2,0) = +0.00;                           A[3].at(2,1) = +0.00;                          A[3].at(2,2) = +400.0/147.0;

    //        A[4].at(0,0) = -225.0/147.0;                    A[4].at(0,1) = +0.00;                          A[4].at(0,2) = +0.00;
    //        A[4].at(1,0) = +0.00;                           A[4].at(1,1) = -225.0/147.0;                   A[4].at(1,2) = +0.00;
    //        A[4].at(2,0) = +0.00;                           A[4].at(2,1) = +0.00;                          A[4].at(2,2) = -225.0/147.0;

    //        A[5].at(0,0) = +72.0/147.0;                     A[5].at(0,1) = +0.00;                          A[5].at(0,2) = +0.00;
    //        A[5].at(1,0) = +0.00;                           A[5].at(1,1) = +72.0/147.0;                    A[5].at(1,2) = +0.00;
    //        A[5].at(2,0) = +0.00;                           A[5].at(2,1) = +0.00;                          A[5].at(2,2) = +72.0/147.0;

    //        A[6].at(0,0) = -10.0/147.0;                     A[6].at(0,1) = +0.00;                          A[6].at(0,2) = +0.00;
    //        A[6].at(1,0) = +0.00;                           A[6].at(1,1) = -10.0/147.0;                    A[6].at(1,2) = +0.00;
    //        A[6].at(2,0) = +0.00;                           A[6].at(2,1) = +0.00;                          A[6].at(2,2) = -10.0/147.0;

    //        A[0].at(0,0) = (60.0*h*b(1,k-1))/147.0;
    //        A[0].at(1,0) = (60.0*h*b(2,k-1))/147.0;
    //        A[0].at(2,0) = (60.0*h*b(3,k-1))/147.0;
    //    }
}

void Example4::clearAMatrices(std::vector<DoubleMatrix> &A)
{
    for (unsigned int i=1; i<=K; i++) A[i].clear();
    A[0].clear();
    A.clear();
}

double Example4::fx(unsigned int i, unsigned int k) const
{
    double t = k*h;
    C_UNUSED(t);

#ifdef SAMPLE_1
    if (i==1) return sin(2.0*t) + t*t;
    if (i==2) return 3.0*t;
    if (i==3) return cos(2.0*t) - sin(t);
#endif
#ifdef SAMPLE_2
    if (i==1) return t*t+t;
    if (i==2) return 2.0*t;
    if (i==3) return 3.0*t*t;
#endif
#ifdef SAMPLE_3
    if (i==1) return sin(t);
    if (i==2) return cos(t);
    if (i==3) return 2.0*sin(t)+3.0*cos(t);
#endif
#ifdef SAMPLE_4
    if (i==1) return sin(20.0*t*t);
    if (i==2) return cos(10.0*t) - sin(20.0*t);
    if (i==3) return t*t*t - sin(8.0*t)*sin(8.0*t);
#endif
#ifdef SAMPLE_5
    if (i==1) return sin(t);
    if (i==2) return cos(t) - sin(t);
    if (i==3) return t*t*t - sin(t);
#endif
    return 0.0;
}

double Example4::a(unsigned int i, unsigned int j, unsigned int k) const
{
    double t = k*h;
    C_UNUSED(t);

#ifdef SAMPLE_1
    if (i==1 && j==1) return +3.0;
    if (i==1 && j==2) return -t;
    if (i==1 && j==3) return +2.0;

    if (i==2 && j==1) return +1.0;
    if (i==2 && j==2) return +1.0;
    if (i==2 && j==3) return +3.0;

    if (i==3 && j==1) return -2.0;
    if (i==3 && j==2) return +t;
    if (i==3 && j==3) return +1.0;
#endif

#ifdef SAMPLE_2
    if (i==1 && j==1) return -3.0;
    if (i==1 && j==2) return +1.0;
    if (i==1 && j==3) return +1.0;

    if (i==2 && j==1) return +3.0;
    if (i==2 && j==2) return -1.5;
    if (i==2 && j==3) return -1.0;

    if (i==3 && j==1) return +6.0;
    if (i==3 && j==2) return +3.0;
    if (i==3 && j==3) return -2.0;
#endif

#ifdef SAMPLE_3
    if (i==1 && j==1) return +2.0;
    if (i==1 && j==2) return +3.0;
    if (i==1 && j==3) return -1.0;

    if (i==2 && j==1) return +4.0;
    if (i==2 && j==2) return +6.0;
    if (i==2 && j==3) return -2.0;

    if (i==3 && j==1) return -1.0;
    if (i==3 && j==2) return +1.0;
    if (i==3 && j==3) return -1.0;
#endif

#ifdef SAMPLE_4
    if (i==1 && j==1) return +2.0;
    if (i==1 && j==2) return +3.0;
    if (i==1 && j==3) return -1.0;

    if (i==2 && j==1) return +4.0;
    if (i==2 && j==2) return +6.0;
    if (i==2 && j==3) return -2.0;

    if (i==3 && j==1) return -1.0;
    if (i==3 && j==2) return +1.0;
    if (i==3 && j==3) return -1.0;
#endif

#ifdef SAMPLE_5
    if (i==1 && j==1) return +2.0;
    if (i==1 && j==2) return +3.0;
    if (i==1 && j==3) return -1.0;

    if (i==2 && j==1) return +4.0;
    if (i==2 && j==2) return +6.0;
    if (i==2 && j==3) return -2.0;

    if (i==3 && j==1) return -1.0;
    if (i==3 && j==2) return +1.0;
    if (i==3 && j==3) return -1.0;
#endif
    return NAN;
}

double Example4::b(unsigned int i, unsigned int k) const
{
    double t = k*h;
    //    C_UNUSED(t);

#ifdef SAMPLE_1
    if (i==1) return -(+3.0*(sin(2.0*t) + t*t) -   t*(3.0*t) + 2.0*(cos(2.0*t) - sin(t))) + (2.0*cos(2.0*t) + 2.0*t);
    if (i==2) return -(+1.0*(sin(2.0*t) + t*t) + 1.0*(3.0*t) + 3.0*(cos(2.0*t) - sin(t))) + (3.0);
    if (i==3) return -(-2.0*(sin(2.0*t) + t*t) +   t*(3.0*t) + 1.0*(cos(2.0*t) - sin(t))) + (-2.0*sin(2.0*t) - cos(t));
#endif
    //#ifdef SAMPLE_1
    //    if (i==1) return 2.0*t + 2.0*sin(t) - 3.0*sin(2.0*t);
    //    if (i==2) return 3.0*sin(t) - sin(2.0*t) - 3.0*cos(2.0*t) - t*t - 3.0*t + 3.0;
    //    if (i==3) return sin(t) - cos(t) - cos(2.0*t) - t*t;
    //#endif
#ifdef SAMPLE_2
    if (i==1) return 3.0*t + 1.0;
    if (i==2) return 2.0;
    if (i==3) return -6.0*t;
#endif
#ifdef SAMPLE_3
    if (n==1) return cos(t);
    if (n==2) return -sin(t);
    if (n==3) return 4.0*cos(t);
#endif
#ifdef SAMPLE_4
    if (i==1) return - (2.0*sin(20.0*t*t) + 3.0*(cos(10.0*t) - sin(20.0*t)) - (t*t*t-sin(8.0*t)*sin(8.0*t)))
            + 40.0*t*cos(20.0*t*t);
    if (i==2) return - (4.0*sin(20.0*t*t) + 6.0*(cos(10.0*t) - sin(20.0*t)) - 2.0*(t*t*t-sin(8.0*t)*sin(8.0*t)))
            - 10.0*sin(10.0*t) - 20.0*cos(20.0*t);
    if (i==3) return - (-sin(20.0*t*t) + (cos(10.0*t) - sin(20.0*t)) -(t*t*t-sin(8.0*t)*sin(8.0*t)))
            + 3.0*t*t - 16.0*cos(8.0*t)*sin(8.0*t);

    //    if (i==1) return -(+2.0*sin(20.0*t*t) + 3.0*(cos(10.0*t) - sin(20.0*t)) - 1.0*(t*t*t - sin(8.0*t)*sin(8.0*t))) + (+40.0*t*cos(20.0*t*t));
    //    if (i==2) return -(+4.0*sin(20.0*t*t) + 6.0*(cos(10.0*t) - sin(20.0*t)) - 2.0*(t*t*t - sin(8.0*t)*sin(8.0*t))) + (-10.0*sin(10.0*t) - 20.0*cos(20.0*t));
    //    if (i==3) return -(-1.0*sin(20.0*t*t) + 1.0*(cos(10.0*t) - sin(20.0*t)) - 1.0*(t*t*t - sin(8.0*t)*sin(8.0*t))) + (+3.0*t*t - 16.0*cos(8.0*t)*sin(8.0*t));
#endif
#ifdef SAMPLE_5
    if (i==1) return -(+2.0*sin(t) + 3.0*(cos(t) - sin(t)) - 1.0*(t*t*t - sin(t))) + (+cos(t*t));
    if (i==2) return -(+4.0*sin(t) + 6.0*(cos(t) - sin(t)) - 2.0*(t*t*t - sin(t))) + (-sin(t) - cos(t));
    if (i==3) return -(-1.0*sin(t) + 1.0*(cos(t) - sin(t)) - 1.0*(t*t*t - sin(t))) + (+3.0*t*t - sin(t));
#endif
    return NAN;
}

void Example4::calculateNX(const DoubleMatrix &rx, DoubleMatrix &nx)
{
    h = 0.01;
    N = 100;
    double *x1 = new double[N+1]; for (unsigned int i=0; i<=N; i++) x1[i] = 0.0;
    double *x2 = new double[N+1]; for (unsigned int i=0; i<=N; i++) x2[i] = 0.0;
    double *x3 = new double[N+1]; for (unsigned int i=0; i<=N; i++) x3[i] = 0.0;

    x1[0] = fx(1,0); x2[0] = fx(2,0); x3[0] = fx(3,0);
    x1[1] = fx(1,1); x2[1] = fx(2,1); x3[1] = fx(3,1);
    x1[2] = fx(1,2); x2[2] = fx(2,2); x3[2] = fx(3,2);
    x1[3] = fx(1,3); x2[3] = fx(2,3); x3[3] = fx(3,3);
    //x1[4] = fx(1,4); x2[4] = fx(2,4); x3[4] = fx(3,4);

    //printf("%16.12f %16.12f %16.12f %16.12f %16.12f %16.12f\n", x1[0], fx(1,0), x2[0], fx(2,0), x3[0], fx(3,0));
    //printf("%16.12f %16.12f %16.12f %16.12f %16.12f %16.12f\n", x1[1], fx(1,1), x2[1], fx(2,1), x3[1], fx(3,1));
    //printf("%16.12f %16.12f %16.12f %16.12f %16.12f %16.12f\n", x1[2], fx(1,2), x2[2], fx(2,2), x3[2], fx(3,2));
    //printf("%16.12f %16.12f %16.12f %16.12f %16.12f %16.12f\n", x1[3], fx(1,3), x2[3], fx(2,3), x3[3], fx(3,3));
    //printf("%16.12f %16.12f %16.12f %16.12f %16.12f %16.12f\n", x1[4], fx(1,4), x2[4], fx(2,4), x3[4], fx(3,4));
    //puts("");

    //    for (unsigned int k=1; k<=N; k++)
    //    {
    //        x1[k] = x1[k-1] + h*(+2.0*x1[k-1]+3.0*x2[k-1]-1.0*x3[k-1] + b(1,k-1));
    //        x2[k] = x2[k-1] + h*(+4.0*x1[k-1]+6.0*x2[k-1]-2.0*x3[k-1] + b(2,k-1));
    //        x3[k] = x3[k-1] + h*(-1.0*x1[k-1]+1.0*x2[k-1]-1.0*x3[k-1] + b(3,k-1));
    //        if (k<10)
    //            printf("%16.12f %16.12f %16.12f %16.12f %16.12f %16.12f\n", x1[k], fx(1,k), x2[k], fx(2,k), x3[k], fx(3,k));
    //    }

    //    for (unsigned int k=2; k<=N; k++)
    //    {
    //        x1[k] = x1[k-2] + 2.0*h*(+2.0*x1[k-1]+3.0*x2[k-1]-1.0*x3[k-1] + b(1,k-1));
    //        x2[k] = x2[k-2] + 2.0*h*(+4.0*x1[k-1]+6.0*x2[k-1]-2.0*x3[k-1] + b(2,k-1));
    //        x3[k] = x3[k-2] + 2.0*h*(-1.0*x1[k-1]+1.0*x2[k-1]-1.0*x3[k-1] + b(3,k-1));
    //        if (k<12) printf("%16.12f %16.12f %16.12f %16.12f %16.12f %16.12f\n", x1[k], fx(1,k), x2[k], fx(2,k), x3[k], fx(3,k));
    //    }

    for (unsigned int k=4; k<=N; k++)
    {
        x1[k] = 8.0*x1[k-1] - 12.0*h*(+2.0*x1[k-2]+3.0*x2[k-2]-1.0*x3[k-2] + b(1,k-2)) - 8.0*x1[k-3] + x1[k-4];
        x2[k] = 8.0*x2[k-1] - 12.0*h*(+4.0*x1[k-2]+6.0*x2[k-2]-2.0*x3[k-2] + b(2,k-2)) - 8.0*x2[k-3] + x2[k-4];
        x3[k] = 8.0*x3[k-1] - 12.0*h*(-1.0*x1[k-2]+1.0*x2[k-2]-1.0*x3[k-2] + b(3,k-2)) - 8.0*x3[k-3] + x3[k-4];
        //if (k<K+10) printf("%16.12f %16.12f %16.12f %16.12f %16.12f %16.12f\n", x1[k], fx(1,k), x2[k], fx(2,k), x3[k], fx(3,k));
    }
    puts("");
    //printf("%16.12f %16.12f %16.12f %16.12f %16.12f %16.12f\n", x1[N], fx(1,N), x2[N], fx(2,N), x3[N], fx(3,N));

    //exit(0);
    //return;

    nx.clear();

    nx.resize(n,N+1,0.0);

    if (K==1)
    {
        nx.setColumn(0, rx.col(0));
    }

    if (K==2)
    {
        nx.setColumn(0, rx.col(0));
        nx.setColumn(1, rx.col(1));
    }

    if (K==4)
    {
        nx.setColumn(0, rx.col(0));
        nx.setColumn(1, rx.col(1));
        nx.setColumn(2, rx.col(2));
        nx.setColumn(3, rx.col(3));

        //        double aa = 8.0*nx.at(0,3) - 12.0*h*(a(1,1,2)*nx.at(0,2)+a(1,2,2)*nx.at(1,2)+a(1,3,2)*nx.at(2,2))
        //                - 8.0*nx.at(0,1) + nx.at(0,0) - 12.0*h*b(1,2);
        //        printf("--- %.10f %.10f\n", aa, fx(1,4));

        //printf("++++++++ %.12f %.12f %.12f %.12f %.12f %.12f\n", nx.at(0,0), fx(1,0), nx.at(1,0), fx(2,0), nx.at(2,0), fx(3,0));
        //printf("++++++++ %.12f %.12f %.12f %.12f %.12f %.12f\n", nx.at(0,1), fx(1,1), nx.at(1,1), fx(2,1), nx.at(2,1), fx(3,1));
        //printf("++++++++ %.12f %.12f %.12f %.12f %.12f %.12f\n", nx.at(0,2), fx(1,2), nx.at(1,2), fx(2,2), nx.at(2,2), fx(3,2));
        //printf("++++++++ %.12f %.12f %.12f %.12f %.12f %.12f\n", nx.at(0,3), fx(1,3), nx.at(1,3), fx(2,3), nx.at(2,3), fx(3,3));
    }

    if (K==6)
    {
        nx.setColumn(0, rx.col(0));
        nx.setColumn(1, rx.col(1));
        nx.setColumn(2, rx.col(2));
        nx.setColumn(3, rx.col(3));
        nx.setColumn(4, rx.col(4));
        nx.setColumn(5, rx.col(5));
    }

    std::vector<DoubleMatrix> A;
    initAMatrices(A);
    for (unsigned int k=K; k<=N; k++)
    {
        updateAMatrices(A,k);
        if (K==1)
        {
            DoubleVector ck = A[1]*nx.col(k-1) + A[0];
            nx.setColumn(k,ck);
        }
        if (K==2)
        {
            DoubleVector ck = A[1]*nx.col(k-1) + A[2]*nx.col(k-2) + A[0];
            nx.setColumn(k,ck);
        }
        if (K==4)
        {
            double aa1 = 8.0*nx.at(0,k-1) - 12.0*h*(a(1,1,k-2)*nx.at(0,k-2)+a(1,2,k-2)*nx.at(1,k-2)+a(1,3,k-2)*nx.at(2,k-2)) - 8.0*nx.at(0,k-3) + nx.at(0,k-4) - 12.0*h*b(1,k-2);
            double aa2 = 8.0*nx.at(1,k-1) - 12.0*h*(a(2,1,k-2)*nx.at(0,k-2)+a(2,2,k-2)*nx.at(1,k-2)+a(2,3,k-2)*nx.at(2,k-2)) - 8.0*nx.at(1,k-3) + nx.at(1,k-4) - 12.0*h*b(2,k-2);
            double aa3 = 8.0*nx.at(2,k-1) - 12.0*h*(a(3,1,k-2)*nx.at(0,k-2)+a(3,2,k-2)*nx.at(1,k-2)+a(3,3,k-2)*nx.at(2,k-2)) - 8.0*nx.at(2,k-3) + nx.at(2,k-4) - 12.0*h*b(3,k-2);
            nx.at(0,k) = aa1;
            nx.at(1,k) = aa2;
            nx.at(2,k) = aa3;

            //            DoubleVector ck;
            //            ck << aa1 << aa2 << aa3;

            //            DoubleVector ck = A[1]*nx.col(k-1) + A[2]*nx.col(k-2) + A[3]*nx.col(k-3) + A[4]*nx.col(k-4) + A[0];
            //            nx.setColumn(k,ck);
            //if (k<K+10)
            //    printf("++++++++ %.12f %.12f %.12f %.12f %.12f %.12f\n",
            //           nx.at(0,k), fx(1,k), nx.at(1,k), fx(2,k), nx.at(2,k), fx(3,k));
        }
        if (K==6)
        {
            DoubleVector ck = A[1]*nx.col(k-1) + A[2]*nx.col(k-2) + A[3]*nx.col(k-3) + A[4]*nx.col(k-4) + A[5]*nx.col(k-5) + A[6]*nx.col(k-6) + A[0];
            nx.setColumn(k,ck);
        }
    }
    clearAMatrices(A);
}

void Example4::calculateM1(const std::vector<unsigned int> *s, const DoubleMatrix &rx UNUSED_PARAM, const DoubleMatrix &nx UNUSED_PARAM)
{
    DoubleMatrix M(K*n, K*n, 0.0);
    DoubleVector B(K*n,0.0);

    if (K==6)
    {
        calculateM1BE(0,s[0],rx,nx,M,B);
        calculateM1BE(1,s[1],rx,nx,M,B);
        calculateM1BE(2,s[2],rx,nx,M,B);
        calculateM1BE(3,s[3],rx,nx,M,B);
        calculateM1BE(4,s[4],rx,nx,M,B);
        calculateM1BE(5,s[5],rx,nx,M,B);
    }

    if (K==4)
    {
        calculateM1BE(0,s[0],rx,nx,M,B);
        calculateM1BE(1,s[1],rx,nx,M,B);
        calculateM1BE(2,s[2],rx,nx,M,B);
        calculateM1BE(3,s[3],rx,nx,M,B);
    }
    if (K==2)
    {
        calculateM1BE(0,s[0],rx,nx,M,B);
        calculateM1BE(1,s[1],rx,nx,M,B);
    }

    printf("det: %14.10f\n",M.determinant());

    DoubleVector x(n*K);
    GaussianElimination(M,B,x);

    IPrinter::printSeperatorLine(NULL,'-',stdout);
    IPrinter::print(x,x.size(),w,p,stdout);
    IPrinter::printSeperatorLine(NULL,'-',stdout);

    DoubleMatrix cx(n,K);
    if (K==6)
    {
        cx.at(0,0) = x[0]; cx.at(0,1) = x[5]; cx.at(0,2) = x[10]; cx.at(0,3) = x[15];
        cx.at(1,0) = x[1]; cx.at(1,1) = x[6]; cx.at(1,2) = x[11]; cx.at(1,3) = x[16];
        cx.at(2,0) = x[2]; cx.at(2,1) = x[7]; cx.at(2,2) = x[12]; cx.at(2,3) = x[17];
        cx.at(3,0) = x[3]; cx.at(3,1) = x[8]; cx.at(3,2) = x[13]; cx.at(3,3) = x[18];
        cx.at(4,0) = x[4]; cx.at(4,1) = x[9]; cx.at(4,2) = x[14]; cx.at(4,3) = x[19];
    }
    if (K==4)
    {
        cx.at(0,0) = x[0]; cx.at(0,1) = x[3]; cx.at(0,2) = x[6]; cx.at(0,3) = x[9];
        cx.at(1,0) = x[1]; cx.at(1,1) = x[4]; cx.at(1,2) = x[7]; cx.at(1,3) = x[10];
        cx.at(2,0) = x[2]; cx.at(2,1) = x[5]; cx.at(2,2) = x[8]; cx.at(2,3) = x[11];
    }
    if (K==2)
    {
        cx.at(0,0) = x[0]; cx.at(0,1) = x[3];
        cx.at(1,0) = x[1]; cx.at(1,1) = x[4];
        cx.at(2,0) = x[2]; cx.at(2,1) = x[5];
    }

    DoubleMatrix nx1;
    calculateNX(cx,nx1);
    IPrinter::printVector(w,p,nx1.row(0),"x1: ");
    IPrinter::printVector(w,p,nx1.row(1),"x2: ");
    IPrinter::printVector(w,p,nx1.row(2),"x3: ");
    IPrinter::printSeperatorLine(NULL,'-', stdout);
}

void Example4::calculateM1BE(unsigned int c, const std::vector<unsigned int> s,  const DoubleMatrix &rx, const DoubleMatrix &nx UNUSED_PARAM, DoubleMatrix &M, DoubleVector &B)
{
    unsigned int L = s.size();
    std::vector<DoubleMatrix> GAMMA(L);
    DoubleVector ETA(n,0.0);
    fillGamma(GAMMA, ETA, c, K);
    if (c == 0) { for (unsigned int i=0; i<L; i++) ETA = GAMMA[i]*rx.col(s[i]) + ETA; }
    printf("%.10f %.10f %.10f\n", ETA[0], ETA[1], ETA[2]);

    std::vector<DoubleMatrix> betta;
    qovmaM1R2L(GAMMA, ETA, s, betta);

    for (unsigned int i=0; i<n; i++)
    {
        for (unsigned int j=0; j<n; j++)
        {
            if (K==6)
            {
                M[c*n+i][0*n+j] = betta[0][i][j];
                M[c*n+i][1*n+j] = betta[1][i][j];
                M[c*n+i][2*n+j] = betta[2][i][j];
                M[c*n+i][3*n+j] = betta[3][i][j];
                M[c*n+i][4*n+j] = betta[4][i][j];
                M[c*n+i][5*n+j] = betta[5][i][j];
            }
            if (K==4)
            {
                M[c*n+i][0*n+j] = betta[0][i][j];
                M[c*n+i][1*n+j] = betta[1][i][j];
                M[c*n+i][2*n+j] = betta[2][i][j];
                M[c*n+i][3*n+j] = betta[3][i][j];
            }
            if (K==2)
            {
                M[c*n+i][0*n+j] = betta[0][i][j];
                M[c*n+i][1*n+j] = betta[1][i][j];
            }
        }
        B.at(c*n+i) = ETA.at(i);
    }
    ETA.clear();
    for (unsigned int i=0; i<L; i++) GAMMA[i].clear();
    GAMMA.clear();
    //    for (unsigned int i=0; i<=N; i++) betta[i].clear();
    //    betta.clear();
}

void Example4::calculateM2(const std::vector<unsigned int> *s, const DoubleMatrix &rx UNUSED_PARAM, const DoubleMatrix &nx UNUSED_PARAM)
{
    /* real solution vectors */
    std::vector<stdDoubleMatrixVector> P(K);
    for (unsigned int i=0; i<K; i++) P[i].resize(N+1);

    stdDoubleMatrixVector Q(N+1);
    calculatePQ(P, Q);

    /* numerical solution vectors */
    DoubleMatrix nx1;
    calculateNS(nx1,rx,P,Q);
    IPrinter::printVector(w,p,nx1.row(0),"x1: ");
    IPrinter::printVector(w,p,nx1.row(1),"x2: ");
    IPrinter::printVector(w,p,nx1.row(2),"x3: ");
    IPrinter::printSeperatorLine(NULL,'-', stdout);
    nx1.clear();

    /* find x0, x1, x2, x3 */

    DoubleMatrix M(K*n, K*n, 0.0);
    DoubleVector B(K*n,0.0);

    calculateM2BE(0,s[0],nx,M,B,P,Q);
    calculateM2BE(1,s[1],nx,M,B,P,Q);
    calculateM2BE(2,s[2],nx,M,B,P,Q);
    calculateM2BE(3,s[3],nx,M,B,P,Q);

    printf("det: %14.10f\n",M.determinant());

    DoubleVector x(n*K);
    GaussianElimination(M, B, x);

    IPrinter::printSeperatorLine(NULL,'-',stdout);
    IPrinter::print(x,x.size(),w,p,stdout);
    IPrinter::printSeperatorLine(NULL,'-',stdout);

    DoubleMatrix rx1(n,K);
    rx1.at(0,0) = x[0]; rx1.at(0,1) = x[3]; rx1.at(0,2) = x[6]; rx1.at(0,3) = x[9];
    rx1.at(1,0) = x[1]; rx1.at(1,1) = x[4]; rx1.at(1,2) = x[7]; rx1.at(1,3) = x[10];
    rx1.at(2,0) = x[2]; rx1.at(2,1) = x[5]; rx1.at(2,2) = x[8]; rx1.at(2,3) = x[11];

    DoubleMatrix cx;
    calculateNS(cx,rx1,P,Q);
    IPrinter::printVector(w,p,cx.row(0),"x1: ");
    IPrinter::printVector(w,p,cx.row(1),"x2: ");
    IPrinter::printVector(w,p,cx.row(2),"x3: ");
    IPrinter::printSeperatorLine(NULL,'-', stdout);
    cx.clear();

    for (unsigned int i=0; i<=N; i++)
    {
        P[3][i].clear();
        P[2][i].clear();
        P[1][i].clear();
        P[0][i].clear();
        Q[i].clear();
    }
    P[3].clear();
    P[2].clear();
    P[1].clear();
    P[0].clear();
    Q.clear();
}

void Example4::calculateM2BE(unsigned int c, const std::vector<unsigned int> s, const DoubleMatrix &nx, DoubleMatrix &M, DoubleVector &B, const std::vector<stdDoubleMatrixVector> &P, const stdDoubleMatrixVector &Q)
{
    unsigned int L = s.size();
    std::vector<DoubleMatrix> GAMMA(L);
    DoubleVector B1(n,0.0);

    //    for (unsigned int i=0; i<L; i++)
    //    {
    //        GAMMA[i].resize(n,n,0.0);
    //        GAMMA[i].randomData();
    //        B1 = GAMMA[i]*nx.col(s[i]) + B1;
    //    }

    fillGamma(GAMMA, B1, c, K);
    if (c == 0) { for (unsigned int i=0; i<L; i++) B1 = GAMMA[i]*nx.col(s[i]) + B1; }

    DoubleMatrix U3(n,n,0.0);
    DoubleMatrix U2(n,n,0.0);
    DoubleMatrix U1(n,n,0.0);
    DoubleMatrix U0(n,n,0.0);
    DoubleMatrix V0(n,1,0.0);

    for (unsigned int i=0; i<L; i++)
    {
        if (s[i] == 3) U3 = U3 + GAMMA[i];
        else if (s[i] == 2) U2 = U2 + GAMMA[i];
        else if (s[i] == 1) U1 = U1 + GAMMA[i];
        else if (s[i] == 0) U0 = U0 + GAMMA[i];
        else
        {
            U3 = U3 + GAMMA[i]*P[3][s[i]];
            U2 = U2 + GAMMA[i]*P[2][s[i]];
            U1 = U1 + GAMMA[i]*P[1][s[i]];
            U0 = U0 + GAMMA[i]*P[0][s[i]];
            V0 = V0 + GAMMA[i]*Q[s[i]];
        }
    }
    B1 = B1 - V0;

    for (unsigned int i=0; i<n; i++)
    {
        for (unsigned int j=0; j<n; j++)
        {
            M[c*n+i][0*n+j] = U0[i][j];
            M[c*n+i][1*n+j] = U1[i][j];
            M[c*n+i][2*n+j] = U2[i][j];
            M[c*n+i][3*n+j] = U3[i][j];
        }
        B.at(c*n+i) = B1.at(i);
    }
}

void Example4::calculateRX(DoubleMatrix &rx)
{
    rx.clear();
    rx.resize(n,N+1,0.0);
    for (unsigned int k=0; k<=N; k++)
    {
        rx.at(0,k) = fx(1,k);
        rx.at(1,k) = fx(2,k);
        rx.at(2,k) = fx(3,k);
    }
}

void Example4::calculateNS(DoubleMatrix &nx, const DoubleMatrix &rx, const std::vector<stdDoubleMatrixVector> &P, const stdDoubleMatrixVector &Q)
{
    if (nx.empty())
    {
        nx.resize(n,N+1,0.0);

        nx.setColumn(0, rx.col(0));
        nx.setColumn(1, rx.col(1));
        nx.setColumn(2, rx.col(2));
        nx.setColumn(3, rx.col(3));

        for (unsigned int k=K; k<=N; k++)
        {
            DoubleVector ck = P[3][k]*nx.col(3) + P[2][k]*nx.col(2) + P[1][k]*nx.col(1) + P[0][k]*nx.col(0) + Q[k];
            nx.setColumn(k,ck);
        }
    }
}

void Example4::calculatePQ(std::vector<stdDoubleMatrixVector> &P, stdDoubleMatrixVector &Q)
{
    /* calculating P,Q matrices */
    stdDoubleMatrixVector A;
    initAMatrices(A);
    for (unsigned int k=K; k<=N; k++)
    {
        updateAMatrices(A,k);

        if (k==K)
        {
            P[3][k] = A[1];
            P[2][k] = A[2];
            P[1][k] = A[3];
            P[0][k] = A[4];
            Q[k]  = A[0];
        }
        else if (k==K+1)
        {
            P[3][k] = A[1]*P[3][k-1] + A[2];
            P[2][k] = A[1]*P[2][k-1] + A[3];
            P[1][k] = A[1]*P[1][k-1] + A[4];
            P[0][k] = A[1]*P[0][k-1];
            Q[k]  = A[1]*Q[k-1] + A[0];
        }
        else if (k==K+2)
        {
            P[3][k] = A[1]*P[3][k-1] + A[2]*P[3][k-2] + A[3];
            P[2][k] = A[1]*P[2][k-1] + A[2]*P[2][k-2] + A[4];
            P[1][k] = A[1]*P[1][k-1] + A[2]*P[1][k-2];
            P[0][k] = A[1]*P[0][k-1] + A[2]*P[0][k-2];
            Q[k]  = A[1]*Q[k-1] + A[2]*Q[k-2] + A[0];
        }
        else if (k==K+3)
        {
            P[3][k] = A[1]*P[3][k-1] + A[2]*P[3][k-2] + A[3]*P[3][k-3] + A[4];
            P[2][k] = A[1]*P[2][k-1] + A[2]*P[2][k-2] + A[3]*P[2][k-3];
            P[1][k] = A[1]*P[1][k-1] + A[2]*P[1][k-2] + A[3]*P[1][k-3];
            P[0][k] = A[1]*P[0][k-1] + A[2]*P[0][k-2] + A[3]*P[0][k-3];
            Q[k]  = A[1]*Q[k-1] + A[2]*Q[k-2] + A[3]*Q[k-3] + A[0];
        }
        if (k>=2*K)
        {
            P[3][k] = A[1]*P[3][k-1] + A[2]*P[3][k-2] + A[3]*P[3][k-3] + A[4]*P[3][k-4];
            P[2][k] = A[1]*P[2][k-1] + A[2]*P[2][k-2] + A[3]*P[2][k-3] + A[4]*P[2][k-4];
            P[1][k] = A[1]*P[1][k-1] + A[2]*P[1][k-2] + A[3]*P[1][k-3] + A[4]*P[1][k-4];
            P[0][k] = A[1]*P[0][k-1] + A[2]*P[0][k-2] + A[3]*P[0][k-3] + A[4]*P[0][k-4];
            Q[k]  = A[1]*Q[k-1] + A[2]*Q[k-2] + A[3]*Q[k-3] + A[4]*Q[k-4] + A[0];
        }
    }
    A[4].clear();
    A[3].clear();
    A[2].clear();
    A[1].clear();
    A[0].clear();
    A.clear();
    /* calculating P,Q matrices */
}

void Example4::fillGamma(stdDoubleMatrixVector &GAMMA, DoubleVector &ETA, unsigned int s, unsigned int k UNUSED_PARAM)
{
    unsigned int L = GAMMA.size();
    if (K == 4)
    {
        if (s == 0)
        {
            for (unsigned int i=0; i<L; i++)
            {
                GAMMA[i].resize(n,n,1.0);
                //GAMMA[i].randomData();
                for (unsigned int i1=0; i1<n; i1++)
                {
                    for (unsigned int i2=0; i2<n; i2++)
                    {
                        GAMMA[i].at(i1,i2) = (rand()%100)*0.1;//sin(i1+i)+cos(i2+i);
                    }
                }

                IPrinter::print(GAMMA[i], GAMMA[i].rows(), GAMMA[i].cols());
                IPrinter::printSeperatorLine();
                //eta = eta + DoubleVector(GAMMA[i]*nx.col(s[i]));
                //                ETA = ETA + DoubleVector(GAMMA[i]*nx.col(s[i]));
            }

            //            for (unsigned int i=0; i<L; i++)
            //            {
            //                ETA = ETA + DoubleVector(GAMMA[i]*nx.col(s[i]));
            //            }


        }
        if (s == 1)
        {
            GAMMA[0].resize(n,n,0.0);
            GAMMA[0].at(0,0) = GAMMA[0].at(1,1) = GAMMA[0].at(2,2) = -3.0;

            GAMMA[1].resize(n,n,0.0);
            GAMMA[1].at(0,0) = -12.0*h*a(1,1,1)-10.0; GAMMA[1].at(0,1) = -12.0*h*a(1,2,1);      GAMMA[1].at(0,2) = -12.0*h*a(1,3,1);
            GAMMA[1].at(1,0) = -12.0*h*a(2,1,1);      GAMMA[1].at(1,1) = -12.0*h*a(2,2,1)-10.0; GAMMA[1].at(1,2) = -12.0*h*a(2,3,1);
            GAMMA[1].at(2,0) = -12.0*h*a(3,1,1);      GAMMA[1].at(2,1) = -12.0*h*a(3,2,1);      GAMMA[1].at(2,2) = -12.0*h*a(3,3,1)-10.0;

            GAMMA[2].resize(n,n,0.0);
            GAMMA[2].at(0,0) = GAMMA[2].at(1,1) = GAMMA[2].at(2,2) = +18.0;

            GAMMA[3].resize(n,n,0.0);
            GAMMA[3].at(0,0) = GAMMA[3].at(1,1) = GAMMA[3].at(2,2) = -6.0;

            GAMMA[4].resize(n,n,0.0);
            GAMMA[4].at(0,0) = GAMMA[4].at(1,1) = GAMMA[4].at(2,2) = +1.0;

            ETA.at(0) = 12.0*h*b(1,1);
            ETA.at(1) = 12.0*h*b(2,1);
            ETA.at(2) = 12.0*h*b(3,1);
        }
        if (s == 2)
        {
            GAMMA[0].resize(n,n,0.0);
            GAMMA[0].at(0,0) = GAMMA[0].at(1,1) = GAMMA[0].at(2,2) = +1.0;

            GAMMA[1].resize(n,n,0.0);
            GAMMA[1].at(0,0) = GAMMA[1].at(1,1) = GAMMA[1].at(2,2) = -8.0;

            GAMMA[2].resize(n,n,0.0);
            GAMMA[2].at(0,0) = -12.0*h*a(1,1,2); GAMMA[2].at(0,1) = -12.0*h*a(1,2,2); GAMMA[2].at(0,2) = -12.0*h*a(1,3,2);
            GAMMA[2].at(1,0) = -12.0*h*a(2,1,2); GAMMA[2].at(1,1) = -12.0*h*a(2,2,2); GAMMA[2].at(1,2) = -12.0*h*a(2,3,2);
            GAMMA[2].at(2,0) = -12.0*h*a(3,1,2); GAMMA[2].at(2,1) = -12.0*h*a(3,2,2); GAMMA[2].at(2,2) = -12.0*h*a(3,3,2);

            GAMMA[3].resize(n,n,0.0);
            GAMMA[3].at(0,0) = GAMMA[3].at(1,1) = GAMMA[3].at(2,2) = +8.0;

            GAMMA[4].resize(n,n,0.0);
            GAMMA[4].at(0,0) = GAMMA[4].at(1,1) = GAMMA[4].at(2,2) = -1.0;

            ETA.at(0) = 12.0*h*b(1,2);
            ETA.at(1) = 12.0*h*b(2,2);
            ETA.at(2) = 12.0*h*b(3,2);
        }
        if (s == 3)
        {
            GAMMA[0].resize(n,n,0.0);
            GAMMA[0].at(0,0) = GAMMA[0].at(1,1) = GAMMA[0].at(2,2) = -1.0;

            GAMMA[1].resize(n,n,0.0);
            GAMMA[1].at(0,0) = GAMMA[1].at(1,1) = GAMMA[1].at(2,2) = +6.0;

            GAMMA[2].resize(n,n,0.0);
            GAMMA[2].at(0,0) = GAMMA[2].at(1,1) = GAMMA[2].at(2,2) = -18.0;

            GAMMA[3].resize(n,n,0.0);
            GAMMA[3].at(0,0) = -12.0*h*a(1,1,3)+10.0; GAMMA[3].at(0,1) = -12.0*h*a(1,2,3);      GAMMA[3].at(0,2) = -12.0*h*a(1,3,3);
            GAMMA[3].at(1,0) = -12.0*h*a(2,1,3);      GAMMA[3].at(1,1) = -12.0*h*a(2,2,3)+10.0; GAMMA[3].at(1,2) = -12.0*h*a(2,3,3);
            GAMMA[3].at(2,0) = -12.0*h*a(3,1,3);      GAMMA[3].at(2,1) = -12.0*h*a(3,2,3);      GAMMA[3].at(2,2) = -12.0*h*a(3,3,3)+10.0;

            GAMMA[4].resize(n,n,0.0);
            GAMMA[4].at(0,0) = GAMMA[4].at(1,1) = GAMMA[4].at(2,2) = +3.0;

            ETA.at(0) = 12.0*h*b(1,3);
            ETA.at(1) = 12.0*h*b(2,3);
            ETA.at(2) = 12.0*h*b(3,3);
        }
    }

    if (K == 2)
    {
        if (s == 0)
        {
            for (unsigned int i=0; i<L; i++)
            {
                GAMMA[i].resize(n,n,1.0);
                for (unsigned int i1=0; i1<n; i1++)
                {
                    for (unsigned int i2=0; i2<n; i2++)
                    {
                        GAMMA[i].at(i1,i2) = (rand()%100)*0.1;
                    }
                }
            }
        }
        if (s == 1)
        {
            GAMMA[0].resize(n,n,0.0);
            GAMMA[0].at(0,0) = GAMMA[0].at(1,1) = GAMMA[0].at(2,2) = -1.0;

            GAMMA[1].resize(n,n,0.0);
            GAMMA[1].at(0,0) = -2.0*h*a(1,1,2); GAMMA[1].at(0,1) = -2.0*h*a(1,2,2); GAMMA[1].at(0,2) = -2.0*h*a(1,3,2);
            GAMMA[1].at(1,0) = -2.0*h*a(2,1,2); GAMMA[1].at(1,1) = -2.0*h*a(2,2,2); GAMMA[1].at(1,2) = -2.0*h*a(2,3,2);
            GAMMA[1].at(2,0) = -2.0*h*a(3,1,2); GAMMA[1].at(2,1) = -2.0*h*a(3,2,2); GAMMA[1].at(2,2) = -2.0*h*a(3,3,2);

            GAMMA[2].resize(n,n,0.0);
            GAMMA[2].at(0,0) = GAMMA[2].at(1,1) = GAMMA[2].at(2,2) = +1.0;

            ETA.at(0) = 2.0*h*b(1,2);
            ETA.at(1) = 2.0*h*b(2,2);
            ETA.at(2) = 2.0*h*b(3,2);
        }
    }
}

//void Example4::calculateNX(const stdDoubleMatrixVector &rx, DoubleVector &x1, DoubleVector &x2, DoubleVector &x3, stdDoubleVectorVector &nx)
//{
//    nx.clear();
//    nx.resize(N+1);

//    for (unsigned int i=0; i<K; i++)
//    {
//        nx.at(i) = rx.at(i);
//        x1 << nx.at(i).at(0);
//        x2 << nx.at(i).at(1);
//        x3 << nx.at(i).at(2);
//    }

//    std::vector<DoubleMatrix> A;
//    initAMatrices(A);
//    for (unsigned int k=K; k<=N; k++)
//    {
//        updateAMatrices(A,k);
//        if (K==2)
//        {
//            nx.at(k) = A[1]*nx.at(k-1) + A[2]*nx.at(k-2) + A[0];
//        }
//        if (K==4)
//        {
//            nx.at(k) = A[1]*nx.at(k-1) + A[2]*nx.at(k-2) + A[3]*nx.at(k-3) + A[4]*nx.at(k-4) + A[0];
//        }
//        if (K==6)
//        {
//            nx.at(k) = A[1]*nx.at(k-1) + A[2]*nx.at(k-2) + A[3]*nx.at(k-3) + A[4]*nx.at(k-4) + A[5]*nx.at(k-5) + A[6]*nx.at(k-6) + A[0];
//        }

//        x1 << nx.at(k).at(0);
//        x2 << nx.at(k).at(1);
//        x3 << nx.at(k).at(2);
//    }
//    clearAMatrices(A);
//}
