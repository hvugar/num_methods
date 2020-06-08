#include "lode1o.h"
#include "nlode1o.h"
#include <math.h>
#include <float.h>

#include "../matrix2d.h"
#include "../printer.h"
#include "../cmethods.h"
#include "../function.h"
#include "../gradient_cjt.h"
#include "../gradient_sd.h"

NonLocalCondition::NonLocalCondition() {}

NonLocalCondition::NonLocalCondition(unsigned int i, const PointNodeODE &node, const DoubleMatrix &m) : i(i), n(node), m(m)  {}

NonLocalCondition::~NonLocalCondition()
{
    m.clear();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

IFirstOrderLinearODEIVP::IFirstOrderLinearODEIVP() {}

IFirstOrderLinearODEIVP::IFirstOrderLinearODEIVP(const IFirstOrderLinearODEIVP &) {}

IFirstOrderLinearODEIVP& IFirstOrderLinearODEIVP::operator=(const IFirstOrderLinearODEIVP &other)
{
    if (this == &other) { return *this; }
    return *this;
}

IFirstOrderLinearODEIVP::~IFirstOrderLinearODEIVP() {}

void IFirstOrderLinearODEIVP::discritize(const std::vector<NonLocalCondition> &co, std::vector<NonLocalCondition> &cn, unsigned int k) const
{
    const auto cnd_size = static_cast<unsigned int>( co.size() );
    const auto h = dimension().step();
    const auto min = dimension().min();
    const auto max = dimension().max();
    //const auto sze = dimension().size();

    if (k==0)
    {
        for (unsigned int i=0; i<co.size(); i++) cn.push_back(co[i]);
        return;
    }

    unsigned int idx = 0;
    for (unsigned int s=0; s<cnd_size; s++)
    {
        const NonLocalCondition &c = co[s];
        const DoubleMatrix &m = c.m;
        const PointNodeODE &n = c.n;

        double min_x = min*h;
        double max_x = max*h;
        for (auto i=min; i<=max; i++)
        {
            double cx = i*h;
            auto dh = fabs(n.x - cx);

            if (k==2)
            {
                if (dh <= h) cn.push_back(NonLocalCondition(idx++, PointNodeODE(cx, i), (1.0 - dh/h)*m));
            }

            if (k==3) // TO-DO
            {
                double h20 = h*h;
                double h21 = 1.0/h20;
                double h22 = 1.0/(2.0*h20);

                if (min_x < n.x && n.x < min_x+h)
                {
                    if (i==min+0) cn.push_back(NonLocalCondition(idx++, PointNodeODE(cx, i), ((h-dh)*(2.0*h-dh)) * h22 * m));
                    if (i==min+1) cn.push_back(NonLocalCondition(idx++, PointNodeODE(cx, i), ((h-dh)*(h+dh)) * h21 * m));
                    if (i==min+2) cn.push_back(NonLocalCondition(idx++, PointNodeODE(cx, i), ((2.0*h-dh)*(h-dh)) * h22 * m));

                }
                else if (max_x-h < n.x && n.x < max_x)
                {
                    if (i==max-2) cn.push_back(NonLocalCondition(idx++, PointNodeODE(cx, i), ((h-dh)*(2.0*h-dh)) * h22 * m));
                    if (i==max-1) cn.push_back(NonLocalCondition(idx++, PointNodeODE(cx, i), ((h-dh)*(h+dh)) * h21 * m));
                    if (i==max-0) cn.push_back(NonLocalCondition(idx++, PointNodeODE(cx, i), ((2.0*h-dh)*(h-dh)) * h22 * m));
                }
                else
                {
                    if (dh <= h && cx <= n.x)               cn.push_back(NonLocalCondition(idx++, PointNodeODE(cx, i), ((2.0*h-dh)*(h-dh)) * h22 * m));
                    if (dh <= h && cx >= n.x)               cn.push_back(NonLocalCondition(idx++, PointNodeODE(cx, i), ((h-dh)*(h+dh)) * h21 * m));
                    if (dh > h && dh <= 2.0*h && cx >= n.x) cn.push_back(NonLocalCondition(idx++, PointNodeODE(cx, i), ((2.0*h-dh)*(h-dh)) * h22 * m));
                }
            }

            if (k==4)
            {
                double h30 = h*h*h;
                double h32 = (1.0/(2.0*h30));
                double h36 = (1.0/(6.0*h30));

                if (min_x < n.x && n.x < min_x+h)
                {
                    if (i==min+0) cn.push_back(NonLocalCondition(idx++, PointNodeODE(cx, i), ((2.0*h-dh)*(h-dh)*(3.0*h-dh)) * h36 * m));
                    if (i==min+1) cn.push_back(NonLocalCondition(idx++, PointNodeODE(cx, i), ((2.0*h+dh)*(h-dh)*(h+dh)) * h32 * m));
                    if (i==min+2) cn.push_back(NonLocalCondition(idx++, PointNodeODE(cx, i), ((2.0*h-dh)*(h-dh)*(h+dh)) * h32 * m));
                    if (i==min+3) cn.push_back(NonLocalCondition(idx++, PointNodeODE(cx, i), ((2.0*h-dh)*(h-dh)*(3.0*h-dh)) * h36 * m));
                }
                else if (max_x-h < n.x && n.x < max_x)
                {
                    if (i==max-3) cn.push_back(NonLocalCondition(idx++, PointNodeODE(cx, i), ((2.0*h-dh)*(h-dh)*(3.0*h-dh)) * h36 * m));
                    if (i==max-2) cn.push_back(NonLocalCondition(idx++, PointNodeODE(cx, i), ((2.0*h-dh)*(h-dh)*(h+dh)) * h32 * m));
                    if (i==max-1) cn.push_back(NonLocalCondition(idx++, PointNodeODE(cx, i), ((2.0*h+dh)*(h-dh)*(h+dh)) * h32 * m));
                    if (i==max-0) cn.push_back(NonLocalCondition(idx++, PointNodeODE(cx, i), ((2.0*h-dh)*(h-dh)*(3.0*h-dh)) * h36 * m));
                }
                else
                {
                    if (dh <= h)               cn.push_back(NonLocalCondition(idx++, PointNodeODE(cx, i), ((2.0*h-dh)*(h-dh)*(h+dh)) * h32 * m));
                    if (dh > h && dh <= 2.0*h) cn.push_back(NonLocalCondition(idx++, PointNodeODE(cx, i), ((2.0*h-dh)*(h-dh)*(3.0*h-dh)) * h36 * m));
                }
            }
        }
    }
}

void IFirstOrderLinearODEIVP::transferOfCondition(const std::vector<NonLocalCondition> &co, const DoubleVector &d, std::vector<DoubleVector> &x, unsigned int k) const
{
    if (co.size() < 2) throw ExceptionODE(1);

    for (unsigned int i=0; i<co.size(); i++)
    {
        const DoubleMatrix &m = co[i].m;

        if (!m.squareMatrix()) throw ExceptionODE(2);

        if (m.rows() != count()) throw ExceptionODE(3);

        if (d.length() != count()) throw ExceptionODE(4);
    }

    //const unsigned int min = static_cast<unsigned int>( dimension().min() );
    //const unsigned int max = static_cast<unsigned int>( dimension().max() );
    const unsigned int sze = static_cast<unsigned int>( dimension().size() );
    const unsigned int end = static_cast<unsigned int>( sze-k+1 );
    const double h = dimension().step();
    const size_t M = count();

    //const unsigned int L = static_cast<unsigned int>(C.size());

    std::vector<NonLocalCondition> C;
    discritize(co, C);

    DoubleMatrix * const D = new DoubleMatrix[sze+1];
    for (unsigned int i=0; i<=sze; i++) D[i].resize(M,M);

    for (unsigned int i=0; i<C.size(); i++)
    {
        unsigned int index = static_cast<unsigned int>(C[i].n.i);
        D[index] += C[i].m;
    }

    DoubleMatrix **betta = new DoubleMatrix*[end+1];

    for (unsigned int i=0; i<=end; i++)
    {
        betta[i] = new DoubleMatrix[k+1];
        if (i==0)
        {
            betta[i][0] = d;
            for (unsigned int j=1; j<=k; j++) betta[i][j] = D[j-1];
        }
        else
        {

            DoubleMatrix *alpha = new DoubleMatrix[k+1];
            PointNodeODE node((i-1)*h,static_cast<int>(i-1));

            if (k==2)
            {
                DoubleMatrix mx(M, M);
                for (unsigned int r=0; r<M; r++)
                {
                    for (unsigned int c=0; c<M; c++)
                    {
                        mx[r][c] = 2.0*h*A(node, r+1, c+1);
                    }
                    mx[r][r] += 3.0;
                }
                mx.inverse();
                alpha[0].resize(M, 1, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[0])[m][0] = -2.0*h*B(node,m+1); } alpha[0] = mx*alpha[0];
                alpha[1].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[1])[m][m] = +4.0; }               alpha[1] = mx*alpha[1];
                alpha[2].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[2])[m][m] = -1.0; }               alpha[2] = mx*alpha[2];
                mx.clear();
            }

            if (k==4)
            {
                DoubleMatrix mx(M, M);
                for (unsigned int r=0; r<M; r++)
                {
                    for (unsigned int c=0; c<M; c++)
                    {
                        mx[r][c] = 12.0*h*A(node, r+1, c+1);
                    }
                    mx[r][r] += 25.0;
                }
                mx.inverse();
                alpha[0].resize(M, 1, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[0])[m][0] = -12.0*h*B(node,m+1); } alpha[0] = mx*alpha[0];
                alpha[1].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[1])[m][m] = +48.0; }               alpha[1] = mx*alpha[1];
                alpha[2].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[2])[m][m] = -36.0; }               alpha[2] = mx*alpha[2];
                alpha[3].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[3])[m][m] = +16.0; }               alpha[3] = mx*alpha[3];
                alpha[4].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[4])[m][m] = -3.0; }                alpha[4] = mx*alpha[4];
                mx.clear();

            }

            if (k==6)
            {
                DoubleMatrix mx(M, M);
                for (unsigned int r=0; r<M; r++)
                {
                    for (unsigned int c=0; c<M; c++)
                    {
                        mx[r][c] = 60.0*h*A(node, r+1, c+1);
                    }
                    mx[r][r] += 147.0;
                }
                mx.inverse();
                alpha[0].resize(M, 1, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[0])[m][0] = -60.0*h*B(node,m+1); } alpha[0] = mx*alpha[0];
                alpha[1].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[1])[m][m] = +360.0; }              alpha[1] = mx*alpha[1];
                alpha[2].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[2])[m][m] = -450.0; }              alpha[2] = mx*alpha[2];
                alpha[3].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[3])[m][m] = +400.0; }              alpha[3] = mx*alpha[3];
                alpha[4].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[4])[m][m] = -225.0; }              alpha[4] = mx*alpha[4];
                alpha[5].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[5])[m][m] = +72.0;  }              alpha[5] = mx*alpha[5];
                alpha[6].resize(M, M, 0.0); for (unsigned int m=0; m<M; m++) { (alpha[6])[m][m] = -10.0;  }              alpha[6] = mx*alpha[6];
                mx.clear();
            }

            betta[i][0] = betta[i-1][0] - betta[i-1][1]*alpha[0];
            for (unsigned int j=1; j<=k-1; j++)
            {
                betta[i][j] = betta[i-1][j+1] + betta[i-1][1]*alpha[j];
            }
            betta[i][k] = betta[i-1][1]*alpha[k] + D[k+(i-1)];

            for (unsigned int j=0; j<=k; j++) alpha[j].clear();
            delete [] alpha;
        }
    }

    DoubleMatrix F((k+1)*M, (k+1)*M);
    DoubleVector g((k+1)*M);

    if (k==2)
    {
        double t1 = (sze-1)*h; PointNodeODE node1(t1,static_cast<int>((sze-1)));
        double t0 = (sze-0)*h; PointNodeODE node0(t0,static_cast<int>((sze-0)));
        for (unsigned int r=0; r<M; r++)
        {
            for (unsigned int c=0; c<M; c++)
            {
                F[0*M+r][0*M+c] = 0.0;
                F[0*M+r][1*M+c] = betta[end][1][r][c];
                F[0*M+r][2*M+c] = betta[end][2][r][c];

                F[1*M+r][0*M+c] = 0.0;                     if (r==c) F[1*M+r][0*M+c] += -1.0;
                F[1*M+r][1*M+c] = -2.0*h*A(node1,r+1,c+1); if (r==c) F[1*M+r][1*M+c] += +0.0;
                F[1*M+r][2*M+c] = 0.0;                     if (r==c) F[1*M+r][2*M+c] += +1.0;

                F[2*M+r][0*M+c] = 0.0;                     if (r==c) F[2*M+r][0*M+c] += +1.0;
                F[2*M+r][1*M+c] = 0.0;                     if (r==c) F[2*M+r][1*M+c] += -4.0;
                F[2*M+r][2*M+c] = -2.0*h*A(node0,r+1,c+1); if (r==c) F[2*M+r][2*M+c] += +3.0;
            }
            g[0*M+r] = betta[end][0][r][0];
            g[1*M+r] = 2.0*h*B(node1,r+1);
            g[2*M+r] = 2.0*h*B(node0,r+1);
        }
    }

    if (k==4)
    {
        double t3 = (sze-3)*h; PointNodeODE node3(t3,static_cast<int>((sze-3)));
        double t2 = (sze-2)*h; PointNodeODE node2(t2,static_cast<int>((sze-2)));
        double t1 = (sze-1)*h; PointNodeODE node1(t1,static_cast<int>((sze-1)));
        double t0 = (sze-0)*h; PointNodeODE node0(t0,static_cast<int>((sze-0)));

        for (unsigned int r=0; r<M; r++)
        {
            for (unsigned int c=0; c<M; c++)
            {
                F[0*M+r][0*M+c] = 0.0;
                F[0*M+r][1*M+c] = betta[end][1][r][c];
                F[0*M+r][2*M+c] = betta[end][2][r][c];
                F[0*M+r][3*M+c] = betta[end][3][r][c];
                F[0*M+r][4*M+c] = betta[end][4][r][c];

                F[1*M+r][0*M+c] = 0.0;                       if (r==c) F[1*M+r][0*M+c] += -3.0;
                F[1*M+r][1*M+c] = -12.0*h*A(node3, r+1,c+1); if (r==c) F[1*M+r][1*M+c] += -10.0;
                F[1*M+r][2*M+c] = 0.0;                       if (r==c) F[1*M+r][2*M+c] += +18.0;
                F[1*M+r][3*M+c] = 0.0;                       if (r==c) F[1*M+r][3*M+c] += -6.0;
                F[1*M+r][4*M+c] = 0.0;                       if (r==c) F[1*M+r][4*M+c] += +1.0;

                F[2*M+r][0*M+c] = 0.0;                       if (r==c) F[2*M+r][0*M+c] += +1.0;
                F[2*M+r][1*M+c] = 0.0;                       if (r==c) F[2*M+r][1*M+c] += -8.0;
                F[2*M+r][2*M+c] = -12.0*h*A(node2, r+1,c+1); if (r==c) F[2*M+r][2*M+c] += +0.0;
                F[2*M+r][3*M+c] = 0.0;                       if (r==c) F[2*M+r][3*M+c] += +8.0;
                F[2*M+r][4*M+c] = 0.0;                       if (r==c) F[2*M+r][4*M+c] += -1.0;

                F[3*M+r][0*M+c] = 0.0;                       if (r==c) F[3*M+r][0*M+c] += -1.0;
                F[3*M+r][1*M+c] = 0.0;                       if (r==c) F[3*M+r][1*M+c] += +6.0;
                F[3*M+r][2*M+c] = 0.0;                       if (r==c) F[3*M+r][2*M+c] += -18.0;
                F[3*M+r][3*M+c] = -12.0*h*A(node1, r+1,c+1); if (r==c) F[3*M+r][3*M+c] += +10.0;
                F[3*M+r][4*M+c] = 0.0;                       if (r==c) F[3*M+r][4*M+c] += +3.0;

                F[4*M+r][0*M+c] = 0.0;                       if (r==c) F[4*M+r][0*M+c] += +3.0;
                F[4*M+r][1*M+c] = 0.0;                       if (r==c) F[4*M+r][1*M+c] += -16.0;
                F[4*M+r][2*M+c] = 0.0;                       if (r==c) F[4*M+r][2*M+c] += +36.0;
                F[4*M+r][3*M+c] = 0.0;                       if (r==c) F[4*M+r][3*M+c] += -48.0;
                F[4*M+r][4*M+c] = -12.0*h*A(node0, r+1,c+1); if (r==c) F[4*M+r][4*M+c] += +25.0;
            }
            g[0*M+r] = betta[end][0][r][0];
            g[1*M+r] = +12.0*h*B(node3, r+1);
            g[2*M+r] = +12.0*h*B(node2, r+1);
            g[3*M+r] = +12.0*h*B(node1, r+1);
            g[4*M+r] = +12.0*h*B(node0, r+1);
        }

    }

    if (k==6)
    {
        double t5 = (sze-5)*h; PointNodeODE node5(t5,static_cast<int>((sze-5)));
        double t4 = (sze-4)*h; PointNodeODE node4(t4,static_cast<int>((sze-4)));
        double t3 = (sze-3)*h; PointNodeODE node3(t3,static_cast<int>((sze-3)));
        double t2 = (sze-2)*h; PointNodeODE node2(t2,static_cast<int>((sze-2)));
        double t1 = (sze-1)*h; PointNodeODE node1(t1,static_cast<int>((sze-1)));
        double t0 = (sze-0)*h; PointNodeODE node0(t0,static_cast<int>((sze-0)));

        for (unsigned int r=0; r<M; r++)
        {
            for (unsigned int c=0; c<M; c++)
            {
                F[0*M+r][0*M+c] = 0.0;
                F[0*M+r][1*M+c] = betta[end][1][r][c];
                F[0*M+r][2*M+c] = betta[end][2][r][c];
                F[0*M+r][3*M+c] = betta[end][3][r][c];
                F[0*M+r][4*M+c] = betta[end][4][r][c];
                F[0*M+r][5*M+c] = betta[end][5][r][c];
                F[0*M+r][6*M+c] = betta[end][6][r][c];

                F[1*M+r][0*M+c] = 0.0;                      if (r==c) F[1*M+r][0*M+c] += -10.0;
                F[1*M+r][1*M+c] = -60.0*h*A(node5,r+1,c+1); if (r==c) F[1*M+r][1*M+c] += -77.0;
                F[1*M+r][2*M+c] = 0.0;                      if (r==c) F[1*M+r][2*M+c] += +150.0;
                F[1*M+r][3*M+c] = 0.0;                      if (r==c) F[1*M+r][3*M+c] += -100.0;
                F[1*M+r][4*M+c] = 0.0;                      if (r==c) F[1*M+r][4*M+c] += +50.0;
                F[1*M+r][5*M+c] = 0.0;                      if (r==c) F[1*M+r][5*M+c] += -15.0;
                F[1*M+r][6*M+c] = 0.0;                      if (r==c) F[1*M+r][6*M+c] += +2.0;

                F[2*M+r][0*M+c] = 0.0;                      if (r==c) F[2*M+r][0*M+c] += +2.0;
                F[2*M+r][1*M+c] = 0.0;                      if (r==c) F[2*M+r][1*M+c] += -24.0;
                F[2*M+r][2*M+c] = -60.0*h*A(node4,r+1,c+1); if (r==c) F[2*M+r][2*M+c] += -35.0;
                F[2*M+r][3*M+c] = 0.0;                      if (r==c) F[2*M+r][3*M+c] += +80.0;
                F[2*M+r][4*M+c] = 0.0;                      if (r==c) F[2*M+r][4*M+c] += -30.0;
                F[2*M+r][5*M+c] = 0.0;                      if (r==c) F[2*M+r][5*M+c] += +8.0;
                F[2*M+r][6*M+c] = 0.0;                      if (r==c) F[2*M+r][6*M+c] += -1.0;

                F[3*M+r][0*M+c] = 0.0;                      if (r==c) F[3*M+r][0*M+c] += -1.0;
                F[3*M+r][1*M+c] = 0.0;                      if (r==c) F[3*M+r][1*M+c] += +9.0;
                F[3*M+r][2*M+c] = 0.0;                      if (r==c) F[3*M+r][2*M+c] += -45.0;
                F[3*M+r][3*M+c] = -60.0*h*A(node3,r+1,c+1); if (r==c) F[3*M+r][3*M+c] += +0.0;
                F[3*M+r][4*M+c] = 0.0;                      if (r==c) F[3*M+r][4*M+c] += +45.0;
                F[3*M+r][5*M+c] = 0.0;                      if (r==c) F[3*M+r][5*M+c] += -9.0;
                F[3*M+r][6*M+c] = 0.0;                      if (r==c) F[3*M+r][6*M+c] += +1.0;

                F[4*M+r][0*M+c] = 0.0;                      if (r==c) F[4*M+r][0*M+c] += +1.0;
                F[4*M+r][1*M+c] = 0.0;                      if (r==c) F[4*M+r][1*M+c] += -8.0;
                F[4*M+r][2*M+c] = 0.0;                      if (r==c) F[4*M+r][2*M+c] += +30.0;
                F[4*M+r][3*M+c] = 0.0;                      if (r==c) F[4*M+r][3*M+c] += -80.0;
                F[4*M+r][4*M+c] = -60.0*h*A(node2,r+1,c+1); if (r==c) F[4*M+r][4*M+c] += +35.0;
                F[4*M+r][5*M+c] = 0.0;                      if (r==c) F[4*M+r][5*M+c] += +24.0;
                F[4*M+r][6*M+c] = 0.0;                      if (r==c) F[4*M+r][6*M+c] += -2.0;

                F[5*M+r][0*M+c] = 0.0;                      if (r==c) F[5*M+r][0*M+c] += -2.0;
                F[5*M+r][1*M+c] = 0.0;                      if (r==c) F[5*M+r][1*M+c] += +15.0;
                F[5*M+r][2*M+c] = 0.0;                      if (r==c) F[5*M+r][2*M+c] += -50.0;
                F[5*M+r][3*M+c] = 0.0;                      if (r==c) F[5*M+r][3*M+c] += +100.0;
                F[5*M+r][4*M+c] = 0.0;                      if (r==c) F[5*M+r][4*M+c] += -150.0;
                F[5*M+r][5*M+c] = -60.0*h*A(node1,r+1,c+1); if (r==c) F[5*M+r][5*M+c] += +77.0;
                F[5*M+r][6*M+c] = 0.0;                      if (r==c) F[5*M+r][6*M+c] += +10.0;

                F[6*M+r][0*M+c] = 0.0;                      if (r==c) F[6*M+r][0*M+c] += +10.0;
                F[6*M+r][1*M+c] = 0.0;                      if (r==c) F[6*M+r][1*M+c] += -72.0;
                F[6*M+r][2*M+c] = 0.0;                      if (r==c) F[6*M+r][2*M+c] += +225.0;
                F[6*M+r][3*M+c] = 0.0;                      if (r==c) F[6*M+r][3*M+c] += -400.0;
                F[6*M+r][4*M+c] = 0.0;                      if (r==c) F[6*M+r][4*M+c] += +450.0;
                F[6*M+r][5*M+c] = 0.0;                      if (r==c) F[6*M+r][5*M+c] += -360.0;
                F[6*M+r][6*M+c] = -60.0*h*A(node0,r+1,c+1); if (r==c) F[6*M+r][6*M+c] += +147.0;
            }
            g[0*M+r] = betta[end][0][r][0];
            g[1*M+r] = +60.0*h*B(node5,r+1);
            g[2*M+r] = +60.0*h*B(node4,r+1);
            g[3*M+r] = +60.0*h*B(node3,r+1);
            g[4*M+r] = +60.0*h*B(node2,r+1);
            g[5*M+r] = +60.0*h*B(node1,r+1);
            g[6*M+r] = +60.0*h*B(node0,r+1);
        }
    }

    DoubleVector xf((k+1)*M);
    LinearEquation::GaussianElimination(F, g, xf);

    F.clear();
    g.clear();

    x.clear();
    x.resize(sze+1); for (unsigned int n=0; n<=sze; n++) x[n].resize(M);

    unsigned int s = xf.length()-M;
    unsigned int e = xf.length()-1;
    for (unsigned int n=sze; n>=sze-k; n--)
    {
        x[n] = xf.mid(s, e);
        s -= M;
        e -= M;
    }
    xf.clear();

    /*********************** Finding previous values ************************/
    unsigned int stop = static_cast<unsigned int>(0)-1;
    for (unsigned int n=(sze-k)-1; n!=stop; n--)
    {
        x[n] = betta[n][0];
        for (unsigned int j=2; j<=k; j++) x[n] -= betta[n][j]*x[n+j-1];
        for (unsigned int j=n+k; j<=sze; j++) x[n] -= D[j]*x[j];
        betta[n][1].inverse();
        x[n] = betta[n][1]*x[n];
    }
    /************************************************************************/

    /* Cleaning betta vector matrices */
    for (unsigned int n=0; n<=end; n++)
    {
        for (unsigned int i=0; i<=k; i++) betta[n][k].clear();
        delete [] betta[n];
    }
    delete [] betta;
    /************************************************************************/

    //    delete [] D;
    C.clear();
}

void IFirstOrderLinearODEIVP::transferOfConditionN(const std::vector<NonLocalCondition> &co, const DoubleVector &d, std::vector<DoubleVector> &x, unsigned int k, unsigned int schema) const
{
    if (co.size() < 2) throw ExceptionODE(1);

    for (unsigned int i=0; i<co.size(); i++)
    {
        const DoubleMatrix &m = co[i].m;

        if (!m.squareMatrix()) throw ExceptionODE(2);

        if (m.rows() != count()) throw ExceptionODE(3);

        if (d.length() != count()) throw ExceptionODE(4);
    }

    const int min = dimension().min();
    const int max = dimension().max();
    const unsigned int size = static_cast<unsigned int>( dimension().size() );
    const unsigned int end = static_cast<unsigned int>( size-(k+1) );
    const double h = dimension().step();
    const unsigned int M = count();

    std::vector<NonLocalCondition> C = co;
    //discritize(co, C);

    DoubleMatrix *betta = new DoubleMatrix[size];
    for (unsigned int i=0; i<size; i++) betta[i].resize(M, M, 0.0);
    DoubleMatrix gamma = d;

    for (unsigned int i=0; i<C.size(); i++)
    {
        unsigned int j = static_cast<unsigned int>(C[i].n.i);
        betta[j] += C[i].m;
        printf("**** %8d %12.6f \n", j, betta[j][0][0]);
    }

    DoubleMatrix **alpha = new DoubleMatrix*[end];

    if (k==4 /*&& schema == 1*/)
    {
        IPrinter::printSeperatorLine("0");
        for (unsigned int r=0; r<M; r++) { for (unsigned int i=0; i<=k; i++) { for (unsigned int c=0; c<M; c++) { printf("%12.6f ", betta[i][r][c]); } printf(" | "); } printf("%12.6f\n", gamma[r][0]); }
        IPrinter::printSeperatorLine();
    }

    for (unsigned int i=0; i<end; i++)
    {
        alpha[i] = new DoubleMatrix[k+1];
        DoubleMatrix *pAlpha = alpha[i];

        if (k==2)
        {
            const unsigned int i0 = i+0; PointNodeODE node0(i0*h, static_cast<int>(i0) );
            const unsigned int i1 = i+1; PointNodeODE node1(i1*h, static_cast<int>(i1) );
            const unsigned int i2 = i+2; PointNodeODE node2(i2*h, static_cast<int>(i2) );

            pAlpha[0].resize(M, 1, 0.0);
            pAlpha[1].resize(M, M, 0.0);
            pAlpha[2].resize(M, M, 0.0);

            for (unsigned int r=0; r<M; r++)
            {
                (pAlpha[0])[r][0] = +0.4*h*(B(node2, r+1)-4.0*B(node1, r+1));
                (pAlpha[1])[r][r] = +0.8;
                (pAlpha[2])[r][r] = +0.2;
                for (unsigned int c=0; c<M; c++)
                {
                    pAlpha[1][r][c] -= 1.6*h*A(node1, r+1, c+1);
                    pAlpha[2][r][c] += 0.4*h*A(node2, r+1, c+1);
                }
            }
        }

        if (k==4)
        {
            const unsigned int i0 = i+0; PointNodeODE node0(i0*h, static_cast<int>(i0) );
            const unsigned int i1 = i+1; PointNodeODE node1(i1*h, static_cast<int>(i1) );
            const unsigned int i2 = i+2; PointNodeODE node2(i2*h, static_cast<int>(i2) );
            const unsigned int i3 = i+3; PointNodeODE node3(i3*h, static_cast<int>(i3) );
            const unsigned int i4 = i+4; PointNodeODE node4(i4*h, static_cast<int>(i4) );

            if (schema == 0)
            {
                DoubleMatrix mx(M, M);
                for (unsigned int r=0; r<M; r++)
                {
                    for (unsigned int c=0; c<M; c++)
                    {
                        mx[r][c] = 12.0*h*A(node0, r+1, c+1);
                    }
                    mx[r][r] += 25.0;
                }
                mx.inverse();

                pAlpha[1].resize(M, M, 0.0);
                pAlpha[2].resize(M, M, 0.0);
                pAlpha[3].resize(M, M, 0.0);
                pAlpha[4].resize(M, M, 0.0);
                pAlpha[0].resize(M, 1, 0.0);

                for (unsigned int r=0; r<M; r++)
                {
                    (pAlpha[1])[r][r] = +20.0;
                    (pAlpha[2])[r][r] =  +0.0;
                    (pAlpha[3])[r][r] = -20.0;
                    (pAlpha[4])[r][r] = +25.0;
                    (pAlpha[0])[r][0] = -12.0*h*(B(node0, r+1)+B(node1, r+1)+B(node2, r+1)+B(node3, r+1)+B(node4, r+1));

                    for (unsigned int c=0; c<M; c++)
                    {
                        pAlpha[1][r][c] -= 12.0*h*A(node1, r+1, c+1);
                        pAlpha[2][r][c] -= 12.0*h*A(node2, r+1, c+1);
                        pAlpha[3][r][c] -= 12.0*h*A(node3, r+1, c+1);
                        pAlpha[4][r][c] -= 12.0*h*A(node4, r+1, c+1);
                    }
                }

                pAlpha[1] = mx*pAlpha[1];
                pAlpha[2] = mx*pAlpha[2];
                pAlpha[3] = mx*pAlpha[3];
                pAlpha[4] = mx*pAlpha[4];
                pAlpha[0] = mx*pAlpha[0];

                mx.clear();
            }

            if (schema == 1)
            {
                DoubleMatrix mx(M, M);
                for (unsigned int r=0; r<M; r++)
                {
                    mx[r][r] = 5.0;
                }
                mx.inverse();

                pAlpha[0].resize(M, 1, 0.0);
                pAlpha[1].resize(M, M, 0.0);
                pAlpha[2].resize(M, M, 0.0);
                pAlpha[3].resize(M, M, 0.0);
                pAlpha[4].resize(M, M, 0.0);

                for (unsigned int r=0; r<M; r++)
                {
                    (pAlpha[1])[r][r] = +4.0;
                    (pAlpha[2])[r][r] = +0.0;
                    (pAlpha[3])[r][r] = -4.0;
                    (pAlpha[4])[r][r] = +5.0;
                    (pAlpha[0])[r][0] = -12.0*h*(B(node1, r+1)-B(node2, r+1)+B(node3, r+1));

                    for (unsigned int c=0; c<M; c++)
                    {
                        pAlpha[1][r][c] -= 12.0*h*A(node1, r+1, c+1);
                        pAlpha[2][r][c] += 12.0*h*A(node2, r+1, c+1);
                        pAlpha[3][r][c] -= 12.0*h*A(node3, r+1, c+1);
                        pAlpha[4][r][c] +=  0.0;
                    }
                }

                pAlpha[0] = mx*pAlpha[0];
                pAlpha[1] = mx*pAlpha[1];
                pAlpha[2] = mx*pAlpha[2];
                pAlpha[3] = mx*pAlpha[3];
                pAlpha[4] = mx*pAlpha[4];

                mx.clear();
            }

            if (schema == 2)
            {
                DoubleMatrix mx(M, M);
                for (unsigned int r=0; r<M; r++)
                {
                    mx[r][r] = 5.3;
                }
                mx.inverse();

                pAlpha[0].resize(M, 1, 0.0);
                pAlpha[1].resize(M, M, 0.0);
                pAlpha[2].resize(M, M, 0.0);
                pAlpha[3].resize(M, M, 0.0);
                pAlpha[4].resize(M, M, 0.0);

                for (unsigned int r=0; r<M; r++)
                {
                    (pAlpha[1])[r][r] = +3.0;
                    (pAlpha[2])[r][r] = +1.8;
                    (pAlpha[3])[r][r] = -4.6;
                    (pAlpha[4])[r][r] = +5.1;
                    (pAlpha[0])[r][0] = -12.0*h*(1.1*B(node1, r+1)-B(node2, r+1)+B(node3, r+1));

                    for (unsigned int c=0; c<M; c++)
                    {
                        pAlpha[1][r][c] += -13.2*h*A(node1, r+1, c+1);
                        pAlpha[2][r][c] += +12.0*h*A(node2, r+1, c+1);
                        pAlpha[3][r][c] += -12.0*h*A(node3, r+1, c+1);
                        pAlpha[4][r][c] +=   0.0;
                    }
                }

                pAlpha[0] = mx*pAlpha[0];
                pAlpha[1] = mx*pAlpha[1];
                pAlpha[2] = mx*pAlpha[2];
                pAlpha[3] = mx*pAlpha[3];
                pAlpha[4] = mx*pAlpha[4];

                mx.clear();
            }
        }

        if (k==6)
        {
        }

        for (unsigned int j=1; j<=k; j++)
        {
            betta[i+j] = betta[i+j] + betta[i]*pAlpha[j];
        }
        gamma = gamma - betta[i]*pAlpha[0];

        //        DoubleMatrix mx = betta[i+1]; mx.inverse();
        //        for (unsigned int j=i+1; j<=size-1; j++)
        //        {
        //            betta[j] = mx*betta[j];
        //        }
        //        gamma = mx*gamma;


        //        if (k==4 && ( i%100==0 ||  i==end-1 )  /*&& schema == 1*/)
        //        {
        //            IPrinter::printSeperatorLine(std::to_string(i+1).data(), '-');
        //            for (unsigned int r=0; r<M; r++) { for (unsigned int j=1; j<=k+1; j++) { for (unsigned int c=0; c<M; c++) { printf("%12.6f ", betta[i+j][r][c]); } printf(" | ");  }
        //                printf("%12.6f\n", gamma[r][0]); }
        ////            if (i==0)
        ////            {
        ////                IPrinter::printSeperatorLine();
        ////                for (unsigned int r=0; r<M; r++) { for (unsigned int j=1; j<=k; j++) { for (unsigned int c=0; c<M; c++) { printf("%12.6f ", pAlpha[j][r][c]); } printf(" | "); } puts(""); }
        ////            }
        //        }


        //////////////////////////////////////////////////////////////////////////

        if (betta[i+1][0][0]>100000.0)
        {
            double norm1 = 0.0;
            for (unsigned int i=1; i<size; i++)
            {
                norm1 += betta[i][0][0]*betta[i][0][0];
            }
            //norm1 += gamma[0][0]*gamma[0][0];
            norm1 = sqrt(norm1);

            for (unsigned int i=1; i<size; i++)
            {
                betta[i][0][0] /= norm1;
            }
            gamma[0][0] /= norm1;
        }

        ////////////////////////////////////////////////////////////////////////////

        //        if (k==4 /*&& ( i%100==0 || i==end-1 ) && schema == 1*/)
        //        {
        //            IPrinter::printSeperatorLine(std::to_string(i+1).data(), '-');
        //            for (unsigned int r=0; r<M; r++) { for (unsigned int j=1; j<=k+1; j++) { for (unsigned int c=0; c<M; c++) { printf("%12.6f ", betta[i+j][r][c]); } printf(" | ");  }
        //                printf("%12.6f | %20.14f\n", gamma[r][0], betta[size-1][0][0]); }
        //            if (i==0)
        //            {
        //                IPrinter::printSeperatorLine();
        //                for (unsigned int r=0; r<M; r++) { for (unsigned int j=1; j<=k; j++) { for (unsigned int c=0; c<M; c++) { printf("%20.14f ", pAlpha[j][r][c]); } printf("|"); } puts(""); }
        //            }
        //        }
    }

    DoubleMatrix F((k+1)*M, (k+1)*M);
    DoubleVector g((k+1)*M);

    if (k==2)
    {
        double t1 = (size-2)*h; PointNodeODE node1(t1,static_cast<int>((size-2)));
        double t0 = (size-1)*h; PointNodeODE node0(t0,static_cast<int>((size-1)));

        for (unsigned int r=0; r<M; r++)
        {
            for (unsigned int c=0; c<M; c++)
            {
                F[0*M+r][0*M+c] = betta[size-3][r][c];
                F[0*M+r][1*M+c] = betta[size-2][r][c];
                F[0*M+r][2*M+c] = betta[size-1][r][c];

                F[1*M+r][0*M+c] = 0.0;                     if (r==c) F[1*M+r][0*M+c] += -1.0;
                F[1*M+r][1*M+c] = -2.0*h*A(node1,r+1,c+1); if (r==c) F[1*M+r][1*M+c] += +0.0;
                F[1*M+r][2*M+c] = 0.0;                     if (r==c) F[1*M+r][2*M+c] += +1.0;

                F[2*M+r][0*M+c] = 0.0;                     if (r==c) F[2*M+r][0*M+c] += +1.0;
                F[2*M+r][1*M+c] = 0.0;                     if (r==c) F[2*M+r][1*M+c] += -4.0;
                F[2*M+r][2*M+c] = -2.0*h*A(node0,r+1,c+1); if (r==c) F[2*M+r][2*M+c] += +3.0;
            }
            g[0*M+r] = gamma[r][0];
            g[1*M+r] = 2.0*h*B(node1,r+1);
            g[2*M+r] = 2.0*h*B(node0,r+1);
        }
    }

    if (k==4)
    {
        double t4 = (size-5)*h; PointNodeODE node4(t4, static_cast<int>((size-5)));
        double t3 = (size-4)*h; PointNodeODE node3(t3, static_cast<int>((size-4)));
        double t2 = (size-3)*h; PointNodeODE node2(t2, static_cast<int>((size-3)));
        double t1 = (size-2)*h; PointNodeODE node1(t1, static_cast<int>((size-2)));
        double t0 = (size-1)*h; PointNodeODE node0(t0, static_cast<int>((size-1)));

        for (unsigned int r=0; r<M; r++)
        {
            for (unsigned int c=0; c<M; c++)
            {
                F[0*M+r][0*M+c] = betta[size-5][r][c];
                F[0*M+r][1*M+c] = betta[size-4][r][c];
                F[0*M+r][2*M+c] = betta[size-3][r][c];
                F[0*M+r][3*M+c] = betta[size-2][r][c];
                F[0*M+r][4*M+c] = betta[size-1][r][c];

                printf("*************** %8d %20.14f \n", size-5, betta[size-5][0][0]);
                printf("*************** %8d %20.14f \n", size-4, betta[size-4][0][0]);
                printf("*************** %8d %20.14f \n", size-3, betta[size-3][0][0]);
                printf("*************** %8d %20.14f \n", size-2, betta[size-2][0][0]);
                printf("*************** %8d %20.14f \n", size-1, betta[size-1][0][0]);

                //                F[0*M+r][0*M+c] = 0.0;
                //                F[0*M+r][1*M+c] = betta[size-4][r][c];
                //                F[0*M+r][2*M+c] = betta[size-3][r][c];
                //                F[0*M+r][3*M+c] = betta[size-2][r][c];
                //                F[0*M+r][4*M+c] = betta[size-1][r][c];

                //                F[0*M+r][0*M+c] = -12.0*h*A(node4, r+1,c+1); if (r==c) F[0*M+r][0*M+c] += -25.0;
                //                F[0*M+r][1*M+c] = 0.0;                       if (r==c) F[0*M+r][1*M+c] += +48.0;
                //                F[0*M+r][2*M+c] = 0.0;                       if (r==c) F[0*M+r][2*M+c] += -36.0;
                //                F[0*M+r][3*M+c] = 0.0;                       if (r==c) F[0*M+r][3*M+c] += +16.0;
                //                F[0*M+r][4*M+c] = 0.0;                       if (r==c) F[0*M+r][4*M+c] += -3.0;

                //F[1*M+r][0*M+c] = -12.0*h*A(node4, r+1,c+1); if (r==c) F[1*M+r][0*M+c] += -25.0;
                //F[1*M+r][1*M+c] = 0.0;                       if (r==c) F[1*M+r][1*M+c] += +48.0;
                //F[1*M+r][2*M+c] = 0.0;                       if (r==c) F[1*M+r][2*M+c] += -36.0;
                //F[1*M+r][3*M+c] = 0.0;                       if (r==c) F[1*M+r][3*M+c] += +16.0;
                //F[1*M+r][4*M+c] = 0.0;                       if (r==c) F[1*M+r][4*M+c] += -3.0;

                F[1*M+r][0*M+c] = 0.0;                       if (r==c) F[1*M+r][0*M+c] +=  -3.0;
                F[1*M+r][1*M+c] = -12.0*h*A(node3, r+1,c+1); if (r==c) F[1*M+r][1*M+c] += -10.0;
                F[1*M+r][2*M+c] = 0.0;                       if (r==c) F[1*M+r][2*M+c] += +18.0;
                F[1*M+r][3*M+c] = 0.0;                       if (r==c) F[1*M+r][3*M+c] +=  -6.0;
                F[1*M+r][4*M+c] = 0.0;                       if (r==c) F[1*M+r][4*M+c] +=  +1.0;

                F[2*M+r][0*M+c] = 0.0;                       if (r==c) F[2*M+r][0*M+c] += +1.0;
                F[2*M+r][1*M+c] = 0.0;                       if (r==c) F[2*M+r][1*M+c] += -8.0;
                F[2*M+r][2*M+c] = -12.0*h*A(node2, r+1,c+1); if (r==c) F[2*M+r][2*M+c] += +0.0;
                F[2*M+r][3*M+c] = 0.0;                       if (r==c) F[2*M+r][3*M+c] += +8.0;
                F[2*M+r][4*M+c] = 0.0;                       if (r==c) F[2*M+r][4*M+c] += -1.0;

                F[3*M+r][0*M+c] = 0.0;                       if (r==c) F[3*M+r][0*M+c] +=  -1.0;
                F[3*M+r][1*M+c] = 0.0;                       if (r==c) F[3*M+r][1*M+c] +=  +6.0;
                F[3*M+r][2*M+c] = 0.0;                       if (r==c) F[3*M+r][2*M+c] += -18.0;
                F[3*M+r][3*M+c] = -12.0*h*A(node1, r+1,c+1); if (r==c) F[3*M+r][3*M+c] += +10.0;
                F[3*M+r][4*M+c] = 0.0;                       if (r==c) F[3*M+r][4*M+c] +=  +3.0;

                F[4*M+r][0*M+c] = 0.0;                       if (r==c) F[4*M+r][0*M+c] +=  +3.0;
                F[4*M+r][1*M+c] = 0.0;                       if (r==c) F[4*M+r][1*M+c] += -16.0;
                F[4*M+r][2*M+c] = 0.0;                       if (r==c) F[4*M+r][2*M+c] += +36.0;
                F[4*M+r][3*M+c] = 0.0;                       if (r==c) F[4*M+r][3*M+c] += -48.0;
                F[4*M+r][4*M+c] = -12.0*h*A(node0, r+1,c+1); if (r==c) F[4*M+r][4*M+c] += +25.0;
            }
            g[0*M+r] = gamma[r][0];
            //            g[0*M+r] = +12.0*h*B(node4, r+1);
            g[1*M+r] = +12.0*h*B(node3, r+1);
            g[2*M+r] = +12.0*h*B(node2, r+1);
            g[3*M+r] = +12.0*h*B(node1, r+1);
            g[4*M+r] = +12.0*h*B(node0, r+1);
        }
    }

    if (k==6)
    {
        double t5 = (size-6)*h; PointNodeODE node5(t5,static_cast<int>((size-5)));
        double t4 = (size-5)*h; PointNodeODE node4(t4,static_cast<int>((size-4)));
        double t3 = (size-4)*h; PointNodeODE node3(t3,static_cast<int>((size-3)));
        double t2 = (size-3)*h; PointNodeODE node2(t2,static_cast<int>((size-2)));
        double t1 = (size-2)*h; PointNodeODE node1(t1,static_cast<int>((size-1)));
        double t0 = (size-1)*h; PointNodeODE node0(t0,static_cast<int>((size-0)));

        for (unsigned int r=0; r<M; r++)
        {
            for (unsigned int c=0; c<M; c++)
            {
                F[0*M+r][0*M+c] = betta[size-7][r][c];
                F[0*M+r][1*M+c] = betta[size-6][r][c];
                F[0*M+r][2*M+c] = betta[size-5][r][c];
                F[0*M+r][3*M+c] = betta[size-4][r][c];
                F[0*M+r][4*M+c] = betta[size-3][r][c];
                F[0*M+r][5*M+c] = betta[size-2][r][c];
                F[0*M+r][6*M+c] = betta[size-1][r][c];

                F[1*M+r][0*M+c] = 0.0;                      if (r==c) F[1*M+r][0*M+c] += -10.0;
                F[1*M+r][1*M+c] = -60.0*h*A(node5,r+1,c+1); if (r==c) F[1*M+r][1*M+c] += -77.0;
                F[1*M+r][2*M+c] = 0.0;                      if (r==c) F[1*M+r][2*M+c] += +150.0;
                F[1*M+r][3*M+c] = 0.0;                      if (r==c) F[1*M+r][3*M+c] += -100.0;
                F[1*M+r][4*M+c] = 0.0;                      if (r==c) F[1*M+r][4*M+c] += +50.0;
                F[1*M+r][5*M+c] = 0.0;                      if (r==c) F[1*M+r][5*M+c] += -15.0;
                F[1*M+r][6*M+c] = 0.0;                      if (r==c) F[1*M+r][6*M+c] += +2.0;

                F[2*M+r][0*M+c] = 0.0;                      if (r==c) F[2*M+r][0*M+c] += +2.0;
                F[2*M+r][1*M+c] = 0.0;                      if (r==c) F[2*M+r][1*M+c] += -24.0;
                F[2*M+r][2*M+c] = -60.0*h*A(node4,r+1,c+1); if (r==c) F[2*M+r][2*M+c] += -35.0;
                F[2*M+r][3*M+c] = 0.0;                      if (r==c) F[2*M+r][3*M+c] += +80.0;
                F[2*M+r][4*M+c] = 0.0;                      if (r==c) F[2*M+r][4*M+c] += -30.0;
                F[2*M+r][5*M+c] = 0.0;                      if (r==c) F[2*M+r][5*M+c] += +8.0;
                F[2*M+r][6*M+c] = 0.0;                      if (r==c) F[2*M+r][6*M+c] += -1.0;

                F[3*M+r][0*M+c] = 0.0;                      if (r==c) F[3*M+r][0*M+c] += -1.0;
                F[3*M+r][1*M+c] = 0.0;                      if (r==c) F[3*M+r][1*M+c] += +9.0;
                F[3*M+r][2*M+c] = 0.0;                      if (r==c) F[3*M+r][2*M+c] += -45.0;
                F[3*M+r][3*M+c] = -60.0*h*A(node3,r+1,c+1); if (r==c) F[3*M+r][3*M+c] += +0.0;
                F[3*M+r][4*M+c] = 0.0;                      if (r==c) F[3*M+r][4*M+c] += +45.0;
                F[3*M+r][5*M+c] = 0.0;                      if (r==c) F[3*M+r][5*M+c] += -9.0;
                F[3*M+r][6*M+c] = 0.0;                      if (r==c) F[3*M+r][6*M+c] += +1.0;

                F[4*M+r][0*M+c] = 0.0;                      if (r==c) F[4*M+r][0*M+c] += +1.0;
                F[4*M+r][1*M+c] = 0.0;                      if (r==c) F[4*M+r][1*M+c] += -8.0;
                F[4*M+r][2*M+c] = 0.0;                      if (r==c) F[4*M+r][2*M+c] += +30.0;
                F[4*M+r][3*M+c] = 0.0;                      if (r==c) F[4*M+r][3*M+c] += -80.0;
                F[4*M+r][4*M+c] = -60.0*h*A(node2,r+1,c+1); if (r==c) F[4*M+r][4*M+c] += +35.0;
                F[4*M+r][5*M+c] = 0.0;                      if (r==c) F[4*M+r][5*M+c] += +24.0;
                F[4*M+r][6*M+c] = 0.0;                      if (r==c) F[4*M+r][6*M+c] += -2.0;

                F[5*M+r][0*M+c] = 0.0;                      if (r==c) F[5*M+r][0*M+c] += -2.0;
                F[5*M+r][1*M+c] = 0.0;                      if (r==c) F[5*M+r][1*M+c] += +15.0;
                F[5*M+r][2*M+c] = 0.0;                      if (r==c) F[5*M+r][2*M+c] += -50.0;
                F[5*M+r][3*M+c] = 0.0;                      if (r==c) F[5*M+r][3*M+c] += +100.0;
                F[5*M+r][4*M+c] = 0.0;                      if (r==c) F[5*M+r][4*M+c] += -150.0;
                F[5*M+r][5*M+c] = -60.0*h*A(node1,r+1,c+1); if (r==c) F[5*M+r][5*M+c] += +77.0;
                F[5*M+r][6*M+c] = 0.0;                      if (r==c) F[5*M+r][6*M+c] += +10.0;

                F[6*M+r][0*M+c] = 0.0;                      if (r==c) F[6*M+r][0*M+c] += +10.0;
                F[6*M+r][1*M+c] = 0.0;                      if (r==c) F[6*M+r][1*M+c] += -72.0;
                F[6*M+r][2*M+c] = 0.0;                      if (r==c) F[6*M+r][2*M+c] += +225.0;
                F[6*M+r][3*M+c] = 0.0;                      if (r==c) F[6*M+r][3*M+c] += -400.0;
                F[6*M+r][4*M+c] = 0.0;                      if (r==c) F[6*M+r][4*M+c] += +450.0;
                F[6*M+r][5*M+c] = 0.0;                      if (r==c) F[6*M+r][5*M+c] += -360.0;
                F[6*M+r][6*M+c] = -60.0*h*A(node0,r+1,c+1); if (r==c) F[6*M+r][6*M+c] += +147.0;
            }
            g[0*M+r] = gamma[0][r];
            g[1*M+r] = +60.0*h*B(node5,r+1);
            g[2*M+r] = +60.0*h*B(node4,r+1);
            g[3*M+r] = +60.0*h*B(node3,r+1);
            g[4*M+r] = +60.0*h*B(node2,r+1);
            g[5*M+r] = +60.0*h*B(node1,r+1);
            g[6*M+r] = +60.0*h*B(node0,r+1);
        }
    }

    printf("det %14.8f\n", F.determinant());

    IPrinter::printSeperatorLine();
    IPrinter::print(F, F.rows(), F.cols(), 12, 6);
    IPrinter::printSeperatorLine();

    DoubleVector xf((k+1)*M);
    LinearEquation::GaussianElimination(F, g, xf);

    //    double a = 0.53582679*betta[size-5][0][0]+0.72896863*betta[size-4][0][0]+0.87630668*betta[size-3][0][0]+0.96858316*betta[size-2][0][0]+1.00000000*betta[size-1][0][0];
    //    double c =      xf[0]*betta[size-5][0][0]+     xf[1]*betta[size-4][0][0]+     xf[2]*betta[size-3][0][0]+     xf[3]*betta[size-2][0][0]+     xf[4]*betta[size-1][0][0];
    //    double b = gamma[0][0];
    //    printf("----- %14.8f %14.8f %14.8f\n", a, c, b);
    //    printf("***** %14.8f %14.8f %14.8f %14.8f %14.8f\n", xf[4], xf[3], xf[2], xf[1], xf[0]);


    F.clear();
    g.clear();

    if (k==4){
        //    xf[4] = 100.000000;
        //    xf[3] =  99.800100;
        //    xf[2] =  99.600400;
        //    xf[1] =  99.400900;
        //    xf[0] =  99.201600;

        //        xf[4] = 1.000000;
        //        xf[3] =  0.980100;
        //        xf[2] =  0.960400;
        //        xf[1] =  0.940900;
        //        xf[0] =  0.921600;

        //    printf("***** %14.8f %14.8f %14.8f %14.8f %14.8f\n", xf[4], xf[3], xf[2], xf[1], xf[0]);

        for (unsigned int i=0; i<k; i++)
        {
            printf("***** ");
            for (unsigned int r=1; r<=M; r++)
            {
                printf("%14.8f ", xf[i*M+M-r]);
            }
            printf("\n");
        }

        //        printf("***** %14.8f %14.8f %14.8f\n", xf[2], xf[1], xf[0]);
        //        printf("***** %14.8f %14.8f %14.8f\n", xf[5], xf[4], xf[3]);
        //        printf("***** %14.8f %14.8f %14.8f\n", xf[8], xf[7], xf[6]);
        //        printf("***** %14.8f %14.8f %14.8f\n", xf[11], xf[10], xf[9]);
        //        printf("***** %14.8f %14.8f %14.8f\n", xf[14], xf[13], xf[12]);

    }

    //    if (k==2) {
    //    xf[2] = 100.000000;
    //    xf[1] =  99.800100;
    //    xf[0] =  99.600400;
    //}

    x.clear();
    x.resize(size); for (unsigned int n=0; n<size; n++) x[n].resize(M);

    unsigned int s = xf.length()-M;
    unsigned int e = xf.length()-1;

    for (unsigned int n=size-1; n>=size-k-1; n--)
    {
        x[n] = xf.mid(s, e);
        s -= M;
        e -= M;
    }
    xf.clear();

    if (k==2)
    {
        for (unsigned int i=size-(k+1); i!=0; i--)
        {
            x[i-1] = alpha[i-1][1]*x[i] + alpha[i-1][2]*x[i+1] + alpha[i-1][0];
        }
    }

    if (k==4)
    {
        for (unsigned int i=size-(k+1); i!=0; i--)
        {
            x[i-1] = alpha[i-1][1]*x[i] + alpha[i-1][2]*x[i+1] + alpha[i-1][3]*x[i+2] + alpha[i-1][4]*x[i+3] + alpha[i-1][0];
        }
    }

    for (unsigned int i=0; i<=end-1; i++)
    {
        for (unsigned int j=0; j<=k; j++) alpha[i][j].clear();
        delete [] alpha[i];
    }
    delete [] alpha;

    for (unsigned int i=0; i<size; i++) betta[i].clear();
    delete [] betta;
    gamma.clear();
}

void IFirstOrderLinearODEIVP::transferOfConditionM(const std::vector<NonLocalCondition> &C, const DoubleVector &d, std::vector<DoubleVector> &x, unsigned int k) const
{
    const unsigned int size = dimension().size();
    const unsigned int M = count();
    const double h = dimension().step();

    std::vector<NonLocalCondition> Cd = C;
    //    discritize(C, Cd, 4);

    struct MyRnFunction : public RnFunction, public IGradient, public IPrinter
    {
        MyRnFunction(const IFirstOrderLinearODEIVP &ode, const std::vector<NonLocalCondition> &C, const DoubleVector &d, unsigned int size, double h, unsigned int k, unsigned int M)
            : ode(ode), C(C), d(d), size(size), h(h), k(k), M(M) {  }

        const IFirstOrderLinearODEIVP &ode;
        const std::vector<NonLocalCondition> &C;
        const DoubleVector d;
        unsigned int size;
        double h;
        unsigned int k;
        unsigned int M;

        virtual double fx(const DoubleVector &x) const
        {
            double res = 0.0;

            {
                for (unsigned int i=0; i<size-k; i++) // size-k
                {
                    const PointNodeODE node(i*h, static_cast<int>(i));
                    const unsigned int offset = i*M;

                    for (unsigned int r=0; r<M; r++)
                    {
                        double fr = -25.0 * x[offset+0*M+r] + 48.0 * x[offset+1*M+r] - 36.0 * x[offset+2*M+r] + 16.0 * x[offset+3*M+r] - 3.0 * x[offset+4*M+r] - 12.0*h*ode.B(node, r+1);
                        for (unsigned int c=0; c<M; c++) { fr += (-12.0*h*ode.A(node, r+1, c+1))*x[offset+0*M+c]; }
                        res += fr*fr;
                    }
                }

                const unsigned int offset = (size-(k+1))*M;
                const unsigned int i = size-(k+1);

                const PointNodeODE node1((i+1)*h, static_cast<int>(i+1));
                const PointNodeODE node2((i+2)*h, static_cast<int>(i+2));
                const PointNodeODE node3((i+3)*h, static_cast<int>(i+3));
                for (unsigned int r=0; r<M; r++)
                {
                    double fr1 = -3.0 * x[offset+0*M+r] -10.0 * x[offset+1*M+r] + 18.0 * x[offset+2*M+r] -  6.0 * x[offset+3*M+r] +  1.0 * x[offset+4*M + r] -12.0*h*ode.B(node1, r+1);
                    double fr2 = +1.0 * x[offset+0*M+r] - 8.0 * x[offset+1*M+r] +  8.0 * x[offset+3*M+r] -  1.0 * x[offset+4*M+r] - 12.0 *h*ode.B(node2, r+1);
                    double fr3 = -1.0 * x[offset+0*M+r] + 6.0 * x[offset+1*M+r] - 18.0 * x[offset+2*M+r] + 10.0 * x[offset+3*M+r] +  3.0 * x[offset+4*M + r] - 12.0*h*ode.B(node3, r+1);

                    for (unsigned int c=0; c<M; c++)
                    {
                        fr1 += (-12.0*h*ode.A(node1, r+1, c+1))*x[offset+1*M+c];
                        fr2 += (-12.0*h*ode.A(node2, r+1, c+1))*x[offset+2*M+c];
                        fr3 += (-12.0*h*ode.A(node3, r+1, c+1))*x[offset+3*M+c];
                    }

                    res += fr1*fr1;
                    res += fr2*fr2;
                    res += fr3*fr3;

                }

                {
                    const PointNodeODE node((i+4)*h, static_cast<int>(i+4));

                    for (unsigned int r=0; r<M; r++)
                    {
                        double fr = 0.0;

                        for (unsigned int s=0; s<C.size(); s++)
                        {
                            const NonLocalCondition &Cs = C[s];

                            for (unsigned int c=0; c<M; c++) { fr += (Cs.m[r][c])*x[Cs.n.i*M+c]; }
                        }

                        fr += -d[r];
                        res += fr*fr;
                    }
                }
            }

            return res;
        }

        virtual void gradient(const DoubleVector &x, DoubleVector &g) const
        {
            const DoubleMatrix E = DoubleMatrix::IdentityMatrix(M);
            const unsigned int N = x.length();

            g.resize(N, 0.0);
            for (unsigned int i=0; i<g.length(); i++) g[i] = 0.0;

            if (k==2)
            {
            }

            if (k==4)
            {

                const PointNodeODE nodeP0(0.0*h, static_cast<int>(0));
                const PointNodeODE nodeP1(1.0*h, static_cast<int>(1));
                const PointNodeODE nodeP2(2.0*h, static_cast<int>(2));
                const PointNodeODE nodeP3(3.0*h, static_cast<int>(3));
                const PointNodeODE nodeP4(4.0*h, static_cast<int>(4));

                /********************************************************* loop from 0 to 4 *********************************************************/
                for (unsigned int r=0; r<M; r++)
                {
                    const unsigned int r1 = r+1;

                    double fr0 = -25.0 * x[0*M+r] + 48.0 * x[1*M+r] - 36.0 * x[2*M+r] + 16.0 * x[3*M+r] - 3.0 * x[4*M+r] - 12.0 * h*ode.B(nodeP0, r1);
                    double fr1 = -25.0 * x[1*M+r] + 48.0 * x[2*M+r] - 36.0 * x[3*M+r] + 16.0 * x[4*M+r] - 3.0 * x[5*M+r] - 12.0 * h*ode.B(nodeP1, r1);
                    double fr2 = -25.0 * x[2*M+r] + 48.0 * x[3*M+r] - 36.0 * x[4*M+r] + 16.0 * x[5*M+r] - 3.0 * x[6*M+r] - 12.0 * h*ode.B(nodeP2, r1);
                    double fr3 = -25.0 * x[3*M+r] + 48.0 * x[4*M+r] - 36.0 * x[5*M+r] + 16.0 * x[6*M+r] - 3.0 * x[7*M+r] - 12.0 * h*ode.B(nodeP3, r1);
                    double fr4 = -25.0 * x[4*M+r] + 48.0 * x[5*M+r] - 36.0 * x[6*M+r] + 16.0 * x[7*M+r] - 3.0 * x[8*M+r] - 12.0 * h*ode.B(nodeP4, r1);

                    for (unsigned int c=0; c<M; c++)
                    {
                        const unsigned int c1 = c+1;

                        fr0 += (-12.0*h*ode.A(nodeP0, r1, c1)) * x[0*M+c];
                        fr1 += (-12.0*h*ode.A(nodeP1, r1, c1)) * x[1*M+c];
                        fr2 += (-12.0*h*ode.A(nodeP2, r1, c1)) * x[2*M+c];
                        fr3 += (-12.0*h*ode.A(nodeP3, r1, c1)) * x[3*M+c];
                        fr4 += (-12.0*h*ode.A(nodeP4, r1, c1)) * x[4*M+c];
                    }

                    for (unsigned int c=0; c<M; c++)
                    {
                        g[0*M+c] += 2.0 * (fr0 * (-25.0*E[r][c]-12.0*h*ode.A(nodeP0, r+1, c+1)));
                        g[1*M+c] += 2.0 * (fr0 * (+48.0*E[r][c]) + fr1 * (-25.0*E[r][c]-12.0*h*ode.A(nodeP1, r+1, c+1)));
                        g[2*M+c] += 2.0 * (fr0 * (-36.0*E[r][c]) + fr1 * (+48.0*E[r][c]) + fr2 * (-25.0*E[r][c]-12.0*h*ode.A(nodeP2, r+1, c+1)));
                        g[3*M+c] += 2.0 * (fr0 * (+16.0*E[r][c]) + fr1 * (-36.0*E[r][c]) + fr2 * (+48.0*E[r][c]) + fr3 * (-25.0*E[r][c]-12.0*h*ode.A(nodeP3, r+1, c+1)));
                        g[4*M+c] += 2.0 * (fr0 * ( -3.0*E[r][c]) + fr1 * (+16.0*E[r][c]) + fr2 * (-36.0*E[r][c]) + fr3 * (+48.0*E[r][c]) + fr4 * (-25.0*E[r][c]-12.0*h*ode.A(nodeP4, r+1, c+1)));
                    }
                }

                /********************************************************* loop from 5 to size-6 *********************************************************/
                for (unsigned int i=5; i<size-5; i++) // 5...
                {
                    PointNodeODE node0((i-4)*h, static_cast<int>(i-4));
                    PointNodeODE node1((i-3)*h, static_cast<int>(i-3));
                    PointNodeODE node2((i-2)*h, static_cast<int>(i-2));
                    PointNodeODE node3((i-1)*h, static_cast<int>(i-1));
                    PointNodeODE node4((i-0)*h, static_cast<int>(i-0));

                    for (unsigned int r=0; r<M; r++)
                    {
                        const unsigned int r1 = r+1;

                        double fr0 = -25.0 * x[(i-4)*M+r] + 48.0 * x[(i-3)*M+r] - 36.0 * x[(i-2)*M+r] + 16.0 * x[(i-1)*M+r] - 3.0 * x[(i+0)*M+r] - 12.0 * h*ode.B(node0, r1);
                        double fr1 = -25.0 * x[(i-3)*M+r] + 48.0 * x[(i-2)*M+r] - 36.0 * x[(i-1)*M+r] + 16.0 * x[(i+0)*M+r] - 3.0 * x[(i+1)*M+r] - 12.0 * h*ode.B(node1, r1);
                        double fr2 = -25.0 * x[(i-2)*M+r] + 48.0 * x[(i-1)*M+r] - 36.0 * x[(i+0)*M+r] + 16.0 * x[(i+1)*M+r] - 3.0 * x[(i+2)*M+r] - 12.0 * h*ode.B(node2, r1);
                        double fr3 = -25.0 * x[(i-1)*M+r] + 48.0 * x[(i+0)*M+r] - 36.0 * x[(i+1)*M+r] + 16.0 * x[(i+2)*M+r] - 3.0 * x[(i+3)*M+r] - 12.0 * h*ode.B(node3, r1);
                        double fr4 = -25.0 * x[(i-0)*M+r] + 48.0 * x[(i+1)*M+r] - 36.0 * x[(i+2)*M+r] + 16.0 * x[(i+3)*M+r] - 3.0 * x[(i+4)*M+r] - 12.0 * h*ode.B(node4, r1);

                        for (unsigned int c=0; c<M; c++)
                        {
                            const unsigned int c1 = c+1;

                            fr0 += (-12.0*h*ode.A(node0, r1, c1)) * x[(i-4)*M+c];
                            fr1 += (-12.0*h*ode.A(node1, r1, c1)) * x[(i-3)*M+c];
                            fr2 += (-12.0*h*ode.A(node2, r1, c1)) * x[(i-2)*M+c];
                            fr3 += (-12.0*h*ode.A(node3, r1, c1)) * x[(i-1)*M+c];
                            fr4 += (-12.0*h*ode.A(node4, r1, c1)) * x[(i-0)*M+c];
                        }

                        for (unsigned int c=0; c<M; c++)
                        {
                            g[i*M+c] += 2.0 * (fr0 * ( -3.0*E[r][c]) + fr1 * (+16.0*E[r][c]) + fr2 * (-36.0*E[r][c]) + fr3 * (+48.0*E[r][c]) + fr4 * (-25.0*E[r][c]-12.0*h*ode.A(node4, r+1, c+1)));
                        }
                    }
                }

                /********************************************************* loop from size-5 to size-1 *********************************************************/
                const unsigned int i=size-1; // 10
                const PointNodeODE nodeM8((i-8)*h, static_cast<int>(i-8));
                const PointNodeODE nodeM7((i-7)*h, static_cast<int>(i-7));
                const PointNodeODE nodeM6((i-6)*h, static_cast<int>(i-6));
                const PointNodeODE nodeM5((i-5)*h, static_cast<int>(i-5));
                const PointNodeODE nodeM4((i-4)*h, static_cast<int>(i-4));
                const PointNodeODE nodeM3((i-3)*h, static_cast<int>(i-3));
                const PointNodeODE nodeM2((i-2)*h, static_cast<int>(i-2));
                const PointNodeODE nodeM1((i-1)*h, static_cast<int>(i-1));

                for (unsigned int r=0; r<M; r++)
                {
                    const unsigned int r1 = r+1;

                    double fr7 = -25.0 * x[(i-8)*M+r] + 48.0 * x[(i-7)*M+r] - 36.0 * x[(i-6)*M+r] + 16.0 * x[(i-5)*M+r] - 3.0 * x[(i-4)*M+r] - 12.0 * h*ode.B(nodeM8, r1);
                    double fr6 = -25.0 * x[(i-7)*M+r] + 48.0 * x[(i-6)*M+r] - 36.0 * x[(i-5)*M+r] + 16.0 * x[(i-4)*M+r] - 3.0 * x[(i-3)*M+r] - 12.0 * h*ode.B(nodeM7, r1);
                    double fr5 = -25.0 * x[(i-6)*M+r] + 48.0 * x[(i-5)*M+r] - 36.0 * x[(i-4)*M+r] + 16.0 * x[(i-3)*M+r] - 3.0 * x[(i-2)*M+r] - 12.0 * h*ode.B(nodeM6, r1);
                    double fr4 = -25.0 * x[(i-5)*M+r] + 48.0 * x[(i-4)*M+r] - 36.0 * x[(i-3)*M+r] + 16.0 * x[(i-2)*M+r] - 3.0 * x[(i-1)*M+r] - 12.0 * h*ode.B(nodeM5, r1);

                    double fr3 = -25.0 * x[(i-4)*M+r] + 48.0 * x[(i-3)*M+r] - 36.0 * x[(i-2)*M+r] + 16.0 * x[(i-1)*M+r] - 3.0 * x[(i+0)*M+r] - 12.0 * h*ode.B(nodeM4, r1);
                    double fr2 =  -3.0 * x[(i-4)*M+r] - 10.0 * x[(i-3)*M+r] + 18.0 * x[(i-2)*M+r] -  6.0 * x[(i-1)*M+r] + 1.0 * x[(i+0)*M+r] - 12.0 * h*ode.B(nodeM3, r1);
                    double fr1 =  +1.0 * x[(i-4)*M+r] -  8.0 * x[(i-3)*M+r] +  0.0 * x[(i-2)*M+r] +  8.0 * x[(i-1)*M+r] - 1.0 * x[(i+0)*M+r] - 12.0 * h*ode.B(nodeM2, r1);
                    double fr0 =  -1.0 * x[(i-4)*M+r] +  6.0 * x[(i-3)*M+r] - 18.0 * x[(i-2)*M+r] + 10.0 * x[(i-1)*M+r] + 3.0 * x[(i+0)*M+r] - 12.0 * h*ode.B(nodeM1, r1);

                    for (unsigned int c=0; c<M; c++)
                    {
                        const unsigned int c1 = c+1;

                        fr7 += (-12.0*h*ode.A(nodeM8, r1, c1)) * x[(i-8)*M+c];
                        fr6 += (-12.0*h*ode.A(nodeM7, r1, c1)) * x[(i-7)*M+c];
                        fr5 += (-12.0*h*ode.A(nodeM6, r1, c1)) * x[(i-6)*M+c];
                        fr4 += (-12.0*h*ode.A(nodeM5, r1, c1)) * x[(i-5)*M+c];

                        fr3 += (-12.0*h*ode.A(nodeM4, r1, c1)) * x[(i-4)*M+c];
                        fr2 += (-12.0*h*ode.A(nodeM3, r1, c1)) * x[(i-3)*M+c];
                        fr1 += (-12.0*h*ode.A(nodeM2, r1, c1)) * x[(i-2)*M+c];
                        fr0 += (-12.0*h*ode.A(nodeM1, r1, c1)) * x[(i-1)*M+c];
                    }

                    for (unsigned int c=0; c<M; c++)
                    {
                        g[(i-4)*M+c] += 2.0 * (fr7 * ( -3.0*E[r][c]) + fr6 * (+16.0*E[r][c]) + fr5 * (-36.0*E[r][c]) + fr4 * (+48.0*E[r][c]) + fr3 * (-25.0*E[r][c] - 12.0*h*ode.A(nodeM4, r+1, c+1)) + fr2 * (-3.0*E[r][c]) + fr1 * (+1.0*E[r][c]) + fr0 * (-1.0*E[r][c]));
                        g[(i-3)*M+c] += 2.0 * (fr6 * ( -3.0*E[r][c]) + fr5 * (+16.0*E[r][c]) + fr4 * (-36.0*E[r][c]) + fr3 * (+48.0*E[r][c]) + fr2 * (-10.0*E[r][c] - 12.0*h*ode.A(nodeM3)) + fr1 * ( -8.0*E[r][c]) + fr0 * ( +6.0*E[r][c]));
                        g[(i-2)*M+c] += 2.0 * (fr5 * ( -3.0*E[r][c]) + fr4 * (+16.0*E[r][c]) + fr3 * (-36.0*E[r][c]) + fr2 * (+18.0*E[r][c]) + fr1 * (-12.0*h*ode.A(nodeM2)) + fr0 * (-18.0*E[r][c]));
                        g[(i-1)*M+c] += 2.0 * (fr4 * ( -3.0*E[r][c]) + fr3 * (+16.0*E[r][c]) + fr2 * ( -6.0*E[r][c]) + fr1 * ( +8.0*E[r][c]) + fr0 * (+10.0*E[r][c] - 12.0*h*ode.A(nodeM1)));
                        g[(i-0)*M+c] += 2.0 * (fr3 * ( -3.0*E[r][c]) + fr2 * ( +1.0*E[r][c]) + fr1 * ( -1.0*E[r][c]) + fr0 * ( +3.0*E[r][c]));
                    }
                }
            }

            /********************************************************* non-local conditions *********************************************************/
            for (unsigned int r=0; r<M; r++)
            {
                double fr = -d[r];
                for (unsigned int s=0; s<C.size(); s++)
                {
                    const NonLocalCondition &Cs = C[s];
                    const DoubleMatrix &m = Cs.m;
                    const unsigned int i = static_cast<unsigned int>(Cs.n.i);
                    for (unsigned int c=0; c<M; c++)
                    {
                        fr += m[r][c] * x[i*M+c];
                    }
                }
                for (unsigned int s=0; s<C.size(); s++)
                {
                    const NonLocalCondition &Cs = C[s];
                    const unsigned int i = static_cast<unsigned int>(Cs.n.i);

                    for (unsigned int c=0; c<M; c++)
                    {
                        g[i*M+c] += 2.0 * Cs.m[r][c] * fr;
                    }
                }
            }
        }

        virtual void print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double f, double, GradientMethod::MethodResult) const
        {
            //printf("%4d %14.8f | %14.8f %14.8f %14.8f %14.8f %14.8f \n", i, f, x[0], x[1], x[2], x[3], x[4]);
            //printf("%4d %14.8f | %14.8f %14.8f %14.8f %14.8f %14.8f \n", i, f, g[0], g[1], g[2], g[3], g[4]);
            //IPrinter::printVector(x);
        }
    };


    DoubleVector X;
    X.resize(size*M, 0.0);

    MyRnFunction func(*this, Cd, d, size, h, k, M);

    ConjugateGradient g;
    g.setFunction(&func);
    g.setGradient(&func);
    g.setPrinter(&func);
    //g.setProjection(&fw1);
    g.setOptimalityTolerance(DBL_EPSILON);
    g.setFunctionTolerance(DBL_EPSILON);
    g.setStepTolerance(DBL_EPSILON);
    g.setR1MinimizeEpsilon(1.0, 2.0*DBL_EPSILON);
    //g.setMaxIterationCount(200);
    g.setNormalize(false);
    g.setResetIteration(false);
    g.showExitMessage(false);
    g.calculate(X);

    for (unsigned int i=0; i<size; i++) x.push_back(X.mid(i*M, i*M+M));
}

void IFirstOrderLinearODEIVP::transferOfConditionS(const std::vector<NonLocalCondition> &co, const DoubleVector &d, std::vector<DoubleVector> &x, unsigned int k, unsigned int schema) const
{
    if (co.size() < 2) throw ExceptionODE(1);

    for (unsigned int i=0; i<co.size(); i++)
    {
        const DoubleMatrix &m = co[i].m;

        if (!m.squareMatrix()) throw ExceptionODE(2);

        if (m.rows() != count()) throw ExceptionODE(3);

        if (d.length() != count()) throw ExceptionODE(4);
    }

    const int min = dimension().min();
    const int max = dimension().max();
    const unsigned int size = static_cast<unsigned int>( dimension().size() );
    const unsigned int end = static_cast<unsigned int>( size-(k+1) );
    const double h = dimension().step();
    const unsigned int M = count();

    std::vector<NonLocalCondition> C = co;
    //discritize(co, C);

    double **betta = static_cast<double**>(malloc(sizeof(double*)*M));
    double *gamma = static_cast<double*>(malloc(sizeof(double)*M));
    for (unsigned int r=0; r<M; r++)
    {
        betta[r] = static_cast<double*>(malloc(sizeof(double)*size));
        for (unsigned int i=0; i<size; i++) betta[r][i] = 0.0;
        gamma[r] = d[r];
    }
    for (unsigned int s=0; s<C.size(); s++)
    {
        unsigned int i = static_cast<unsigned int>(C[s].n.i);
        for (unsigned int r=0; r<M; r++)
        {
            betta[r][i] += C[s].m[r][0];
        }
    }

    double **alpha = static_cast<double**>(malloc(sizeof(double*)*((k+1)*M)));
    for (unsigned int i=0; i<(k+1)*M; i++) alpha[i] = static_cast<double*>(malloc(sizeof(double)*(size*M)));
    double *zetta = static_cast<double*>(malloc(sizeof(double)*(size*M)));

    if (k==4)
    {
        if (schema == 1)
        {
            for (unsigned int i=0; i<size-5; i++)
            {
                const unsigned int i0 = i+0; PointNodeODE node0(i0*h, static_cast<int>(i0));
                const unsigned int i1 = i+1; PointNodeODE node1(i1*h, static_cast<int>(i1));
                const unsigned int i2 = i+2; PointNodeODE node2(i2*h, static_cast<int>(i2));
                const unsigned int i3 = i+3; PointNodeODE node3(i3*h, static_cast<int>(i3));
                const unsigned int i4 = i+4; PointNodeODE node4(i4*h, static_cast<int>(i4));

                alpha[0][i] = +1.0;
                alpha[1][i] = +0.8 - 2.4*h*A(node1, 1, 1);
                alpha[2][i] = +0.0 + 2.4*h*A(node2, 1, 1);
                alpha[3][i] = -0.8 - 2.4*h*A(node3, 1, 1);
                alpha[4][i] = +1.0;
                zetta[i] = -2.4*h*B(node1, 1) + 2.4*h*B(node2, 1) - 2.4*h*B(node3, 1);

                //alpha[0][i] = 1.0;
                //alpha[1][i] = (+3.0 - 13.2*h*A(node1, 1, 1))/(-5.3);
                //alpha[2][i] = (+1.8 + 12.0*h*A(node2, 1, 1))/(-5.3);
                //alpha[3][i] = (-4.6 - 12.0*h*A(node3, 1, 1))/(-5.3);
                //alpha[4][i] = (+5.1)/(-5.3);
                //zetta[i] = (13.2*h*B(node1, 1) - 12.0*h*B(node2, 1) + 12.0*h*B(node3, 1))/(-5.3);

                //alpha[0][i] = 1.0;
                //alpha[1][i] = -0.56603773584906 + 2.49056603773585*h*A(node1, 1, 1);
                //alpha[2][i] = -0.33962264150943 - 2.26415094339623*h*A(node2, 1, 1);
                //alpha[3][i] = +0.86792452830189 + 2.26415094339623*h*A(node3, 1, 1);
                //alpha[4][i] = -0.96226415094340;
                //zetta[i]    = -2.49056603773585*h*B(node1, 1) + 2.26415094339623*h*B(node2, 1) - 2.26415094339623*h*B(node3, 1);

                //alpha[0][i] = 1.0;
                //alpha[1][i] = (+4.0 - 12.0*h*A(node1, 1, 1))/(-5.0);
                //alpha[2][i] = (+0.0 + 12.0*h*A(node2, 1, 1))/(-5.0);
                //alpha[3][i] = (-4.0 - 12.0*h*A(node3, 1, 1))/(-5.0);
                //alpha[4][i] = (+5.0)/(-5.0);
                //zetta[i] = (12.0*h*B(node1, 1) - 12.0*h*B(node2, 1) + 12.0*h*B(node3, 1))/(-5.0);

                //alpha[0][i] = +1.0;
                //alpha[1][i] = -0.8 + 2.4*h*A(node1, 1, 1);
                //alpha[2][i] = +0.0 - 2.4*h*A(node2, 1, 1);
                //alpha[3][i] = +0.8 + 2.4*h*A(node3, 1, 1);
                //alpha[4][i] = -1.0;
                //zetta[i] = -2.4*h*B(node1, 1) + 2.4*h*B(node2, 1) - 2.4*h*B(node3, 1);
            }

            for (unsigned int i=0; i<size-5; i++)
            {
                //            if (i>99000)
                //            {
                //printf("%4d | %25.14f %25.14f %25.14f %25.14f %25.14f | %25.14f\n", i, alpha[0][i], alpha[1][i], alpha[2][i], alpha[3][i], alpha[4][i], zetta[i]);
                //printf("%4d | %25.14f %25.14f %25.14f %25.14f | %25.14f\n", i, betta[0][i], betta[0][i+1], betta[0][i+2], betta[0][i+3], gamma[0]);
                //            }

                betta[0][i+1] += betta[0][i]*alpha[1][i];
                betta[0][i+2] += betta[0][i]*alpha[2][i];
                betta[0][i+3] += betta[0][i]*alpha[3][i];
                betta[0][i+4] += betta[0][i]*alpha[4][i];
                gamma[0]      -= betta[0][i]*zetta[i];

                //if (fabs(betta[0][i+1])>10.0)
                //                {
                //                    betta[0][i+2] /= betta[0][i+1];
                //                    betta[0][i+3] /= betta[0][i+1];
                //                    betta[0][i+4] /= betta[0][i+1];
                //                    gamma[0]      /= betta[0][i+1];
                //                    betta[0][size-1] /= betta[0][i+1];
                //                    betta[0][i+1] /= betta[0][i+1];
                //                }
            }

            //printf("%4d | %25.14f %25.14f %25.14f %25.14f %25.14f | %25.14f\n", size-5, betta[0][size-5], betta[0][size-4], betta[0][size-3], betta[0][size-2], betta[0][size-1], gamma[0]);
            //printf("*********  %20.14f %20.14f\n", betta[0][size-5]*(0.96)*(0.96) + betta[0][size-4]*(0.97)*(0.97) + betta[0][size-3]*(0.98)*(0.98) + betta[0][size-2]*(0.99)*(0.99)
            //                + betta[0][size-1]*(1.00)*(1.00), gamma[0]);

            double t4 = (size-5)*h; PointNodeODE node4(t4, static_cast<int>((size-5)));
            double t3 = (size-4)*h; PointNodeODE node3(t3, static_cast<int>((size-4)));
            double t2 = (size-3)*h; PointNodeODE node2(t2, static_cast<int>((size-3)));
            double t1 = (size-2)*h; PointNodeODE node1(t1, static_cast<int>((size-2)));
            double t0 = (size-1)*h; PointNodeODE node0(t0, static_cast<int>((size-1)));

            DoubleMatrix F(5, 5);
            DoubleVector f(5);

            F[0][0] = betta[0][size-5];
            F[0][1] = betta[0][size-4];
            F[0][2] = betta[0][size-3];
            F[0][3] = betta[0][size-2];
            F[0][4] = betta[0][size-1];

            //        F[0][0] = -25.0 - 12.0*h*A(node4, 1, 1);
            //        F[0][1] = +48.0;
            //        F[0][2] = -36.0;
            //        F[0][3] = +16.0;
            //        F[0][4] = -3.0;

            F[1][0] = -3.0;
            F[1][1] = -10.0 - 12.0*h*A(node3, 1, 1);
            F[1][2] = +18.0;
            F[1][3] = -6.0;
            F[1][4] = +1.0;

            F[2][0] = +1.0;
            F[2][1] = -8.0;
            F[2][2] = +0.0 - 12.0*h*A(node2, 1, 1);
            F[2][3] = +8.0;
            F[2][4] = -1.0;

            F[3][0] = -1.0;
            F[3][1] = +6.0;
            F[3][2] = -18.0;
            F[3][3] = +10.0 - 12.0*h*A(node1, 1, 1);
            F[3][4] = +3.0;

            F[4][0] = +3.0;
            F[4][1] = -16.0;
            F[4][2] = +36.0;
            F[4][3] = -48.0;
            F[4][4] = +25.0 - 12.0*h*A(node0, 1, 1);


            f[0] = gamma[0];
            //        f[0] = +12.0*h*B(node4, 1);
            f[1] = +12.0*h*B(node3, 1);
            f[2] = +12.0*h*B(node2, 1);
            f[3] = +12.0*h*B(node1, 1);
            f[4] = +12.0*h*B(node0, 1);

            IPrinter::printSeperatorLine();
            IPrinter::print(F, F.rows(), F.cols(), 12, 6);
            IPrinter::printSeperatorLine();

            DoubleVector xf((k+1)*M);
            LinearEquation::GaussianElimination(F, f, xf);

            //    double a = 0.53582679*betta[size-5][0][0]+0.72896863*betta[size-4][0][0]+0.87630668*betta[size-3][0][0]+0.96858316*betta[size-2][0][0]+1.00000000*betta[size-1][0][0];
            //    double c =      xf[0]*betta[size-5][0][0]+     xf[1]*betta[size-4][0][0]+     xf[2]*betta[size-3][0][0]+     xf[3]*betta[size-2][0][0]+     xf[4]*betta[size-1][0][0];
            //    double b = gamma[0][0];
            //    printf("----- %14.8f %14.8f %14.8f\n", a, c, b);
            //    printf("***** %14.8f %14.8f %14.8f %14.8f %14.8f\n", xf[4], xf[3], xf[2], xf[1], xf[0]);


            F.clear();
            f.clear();

            for (unsigned int i=0; i<=k; i++)
            {
                printf("***** ");
                for (unsigned int r=1; r<=M; r++)
                {
                    printf("%4d %14.8f ", i*M+M-r, xf[i*M+M-r]);
                }
                printf("\n");
            }

            x.clear();
            x.resize(size); for (unsigned int n=0; n<size; n++) x[n].resize(M);

            unsigned int s = xf.length()-M;
            unsigned int e = xf.length()-1;

            for (unsigned int n=size-1; n>=size-k-1; n--)
            {
                x[n] = xf.mid(s, e);
                s -= M;
                e -= M;
            }
            xf.clear();

            if (k==4)
            {
                for (unsigned int i=size-(k+1); i!=0; i--)
                {
                    x[i-1][0] = alpha[1][i-1]*x[i+0][0] + alpha[2][i-1]*x[i+1][0] + alpha[3][i-1]*x[i+2][0] + alpha[4][i-1]*x[i+3][0] + zetta[i-1];

                    printf("************************  %4d %14.8f \n", i-1, x[i-1][0]);
                }
            }

        }

        if (schema == 11)
        {
            unsigned int sizeM = size*M;
            unsigned int sizeK = (k+1)*M;
            double **alpha1 = static_cast<double**>( malloc(sizeof(double*) * sizeM ) );
            for (unsigned int i=0; i<sizeM; i++) alpha1[i] = static_cast<double*>( malloc(sizeof(double)* sizeK) );
            double *zetta1 = static_cast<double*>( malloc(sizeof(double) * sizeM) );

            double **betta1 = static_cast<double**>(malloc(sizeof(double*)*M));
            double *gamma1 = static_cast<double*>(malloc(sizeof(double)*M));

            for (unsigned int r=0; r<M; r++)
            {
                betta1[r] = static_cast<double*>(malloc(sizeof(double) * sizeM ));
                for (unsigned int i=0; i<sizeM; i++) betta1[r][i] = 0.0;
                gamma1[r] = d[r];

                printf("%d %f %f\n", r, gamma1[r], d[r]);
            }

            for (unsigned int s=0; s<C.size(); s++)
            {
                unsigned int i = static_cast<unsigned int>(C[s].n.i);
                for (unsigned int r=0; r<M; r++)
                {
                    for (unsigned int c=0; c<M; c++)
                    {
                        betta1[r][i*M+c] += C[s].m[r][c];
                    }
                }
            }

            //            printf("%20.14f | %20.14f\n", betta1[0][0], betta1[0][100]);
            //            printf("%20.14f %20.14f %20.14f | %20.14f %20.14f %20.14f\n", betta1[0][0], betta1[0][1], betta1[0][2], betta1[0][300], betta1[0][301], betta1[0][302]);
            //            printf("%20.14f %20.14f %20.14f | %20.14f %20.14f %20.14f\n", betta1[1][0], betta1[1][1], betta1[1][2], betta1[1][300], betta1[1][301], betta1[1][302]);
            //            printf("%20.14f %20.14f %20.14f | %20.14f %20.14f %20.14f\n", betta1[2][0], betta1[2][1], betta1[2][2], betta1[2][300], betta1[2][301], betta1[2][302]);

            for (unsigned int i=0; i<=size-5; i++)
            {
                const unsigned int i0 = i+0; PointNodeODE node0(i0*h, static_cast<int>(i0));
                const unsigned int i1 = i+1; PointNodeODE node1(i1*h, static_cast<int>(i1));
                const unsigned int i2 = i+2; PointNodeODE node2(i2*h, static_cast<int>(i2));
                const unsigned int i3 = i+3; PointNodeODE node3(i3*h, static_cast<int>(i3));
                const unsigned int i4 = i+4; PointNodeODE node4(i4*h, static_cast<int>(i4));

                //                printf("------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");

                for (unsigned int r=0; r<M; r++)
                {
                    const unsigned int rw = i*M+r;

                    for (unsigned int c=0; c<M; c++)
                    {
                        alpha1[rw][0*M+c] = +0.0;
                        alpha1[rw][1*M+c] = -2.4*h*A(node1, r+1, c+1);
                        alpha1[rw][2*M+c] = +2.4*h*A(node2, r+1, c+1);
                        alpha1[rw][3*M+c] = -2.4*h*A(node3, r+1, c+1);
                        alpha1[rw][4*M+c] = +0.0;
                    }

                    alpha1[rw][0*M+r] += +1.0;
                    alpha1[rw][1*M+r] += +0.8;
                    alpha1[rw][2*M+r] += +0.0;
                    alpha1[rw][3*M+r] += -0.8;
                    alpha1[rw][4*M+r] += +1.0;
                    zetta1[rw] = -2.4*h*B(node1, r+1) + 2.4*h*B(node2, r+1) - 2.4*h*B(node3, r+1);


                    if (i==0)
                    {
                        //                        printf("alpha +++++\n");
                        //                        printf("%4d %4d | %10.6f %10.6f %10.6f | %10.6f %10.6f %10.6f | %10.6f %10.6f %10.6f | %10.6f %10.6f %10.6f | %10.6f %10.6f %10.6f | %10.6f\n", i, rw,
                        //                               alpha1[rw][0], alpha1[rw][1], alpha1[rw][2],
                        //                                alpha1[rw][3], alpha1[rw][4], alpha1[rw][5],
                        //                                alpha1[rw][6], alpha1[rw][7], alpha1[rw][8],
                        //                                alpha1[rw][9], alpha1[rw][10], alpha1[rw][11],
                        //                                alpha1[rw][12], alpha1[rw][13], alpha1[rw][14], zetta1[rw] );
                        //                        puts("-");
                    }

                    for (unsigned int r1=0; r1<M; r1++)
                    {
                        //printf("++++++++++++++++ %d ++++++++++++++\n", r1);
                        //                        printf("%4d %4d | %10.6f %10.6f %10.6f | %10.6f %10.6f %10.6f | %10.6f %10.6f %10.6f | %10.6f %10.6f %10.6f | %10.6f %10.6f %10.6f | %10.6f\n", i, r1,
                        //                               betta1[r1][i*M+0], betta1[r1][i*M+1], betta1[r1][i*M+2],
                        //                                betta1[r1][i*M+3], betta1[r1][i*M+4], betta1[r1][i*M+5],
                        //                                betta1[r1][i*M+6], betta1[r1][i*M+7], betta1[r1][i*M+8],
                        //                                betta1[r1][i*M+9], betta1[r1][i*M+10], betta1[r1][i*M+11],
                        //                                betta1[r1][i*M+12], betta1[r1][i*M+13], betta1[r1][i*M+14], gamma1[r1]);

                        double mm = 0.0;
                        for (unsigned int c=0; c<M; c++)
                        {
                            betta1[r1][(i+1)*M+c] += betta1[r1][rw]*alpha1[rw][1*M+c]; if (mm<fabs(betta1[r1][(i+1)*M+c])) { mm = fabs(betta1[r1][(i+1)*M+c]); }
                            betta1[r1][(i+2)*M+c] += betta1[r1][rw]*alpha1[rw][2*M+c]; if (mm<fabs(betta1[r1][(i+2)*M+c])) { mm = fabs(betta1[r1][(i+2)*M+c]); }
                            betta1[r1][(i+3)*M+c] += betta1[r1][rw]*alpha1[rw][3*M+c]; if (mm<fabs(betta1[r1][(i+3)*M+c])) { mm = fabs(betta1[r1][(i+3)*M+c]); }
                            betta1[r1][(i+4)*M+c] += betta1[r1][rw]*alpha1[rw][4*M+c]; if (mm<fabs(betta1[r1][(i+4)*M+c])) { mm = fabs(betta1[r1][(i+4)*M+c]); }
                        }
                        gamma1[r1] -= betta1[r1][rw]*zetta1[rw]; if (mm<fabs(gamma1[r1])) { mm = fabs(gamma1[r1]); }

                        //mm = 1.0;
                        if (i==size-5/* || mm > 10.0*/)
                        {
                            for (unsigned int c=0; c<M; c++)
                            {
                                betta1[r1][(i+0)*M+c] /= mm;
                                betta1[r1][(i+1)*M+c] /= mm;
                                betta1[r1][(i+2)*M+c] /= mm;
                                betta1[r1][(i+3)*M+c] /= mm;
                                betta1[r1][(i+4)*M+c] /= mm;

                                // betta1[r1][(size-1)*M+c] /= mm;
                            }
                            gamma1[r1] /= mm;
                        }
                        mm = 0.0;

                        //printf("********************* %f %f\n", zetta1[rw], gamma1[r1]);

                        //                        printf("%4d %4d | %10.6f %10.6f %10.6f | %10.6f %10.6f %10.6f | %10.6f %10.6f %10.6f | %10.6f %10.6f %10.6f | %10.6f %10.6f %10.6f | %10.6f\n", i, r1,
                        //                               betta1[r1][(i)*M+0], betta1[r1][(i)*M+1], betta1[r1][(i)*M+2],
                        //                                betta1[r1][(i)*M+3], betta1[r1][(i)*M+4], betta1[r1][(i)*M+5],
                        //                                betta1[r1][(i)*M+6], betta1[r1][(i)*M+7], betta1[r1][(i)*M+8],
                        //                                betta1[r1][(i)*M+9], betta1[r1][(i)*M+10], betta1[r1][(i)*M+11],
                        //                                betta1[r1][(i)*M+12], betta1[r1][(i)*M+13], betta1[r1][(i)*M+14], gamma1[r1]);
                    }
                    //                    puts("---");
                    // if (i==3) exit(0);
                }

                //                IPrinter::printSeperatorLine();
                //                printf("%20.14f %20.14f %20.14f | %20.14f %20.14f %20.14f\n", betta1[0][(i+1)*M+0], betta1[0][(i+1)*M+1], betta1[0][(i+1)*M+2], betta1[0][(i)*M+0], betta1[0][(i)*M+1], betta1[0][(i)*M+2]);
                //                printf("%20.14f %20.14f %20.14f | %20.14f %20.14f %20.14f\n", betta1[1][(i+1)*M+0], betta1[1][(i+1)*M+1], betta1[1][(i+1)*M+2], betta1[1][(i)*M+0], betta1[1][(i)*M+1], betta1[1][(i)*M+2]);
                //                printf("%20.14f %20.14f %20.14f | %20.14f %20.14f %20.14f\n", betta1[2][(i+1)*M+0], betta1[2][(i+1)*M+1], betta1[2][(i+1)*M+2], betta1[2][(i)*M+0], betta1[2][(i)*M+1], betta1[2][(i)*M+2]);

            }

            //            for (unsigned int i=0; i<size-5; i++)
            //            {
            //                //            if (i>99000)
            //                //            {
            //                //printf("%4d | %25.14f %25.14f %25.14f %25.14f %25.14f | %25.14f\n", i, alpha[0][i], alpha[1][i], alpha[2][i], alpha[3][i], alpha[4][i], zetta[i]);
            //                //printf("%4d | %25.14f %25.14f %25.14f %25.14f | %25.14f\n", i, betta[0][i], betta[0][i+1], betta[0][i+2], betta[0][i+3], gamma[0]);
            //                //            }

            //                betta[0][i+1] += betta[0][i]*alpha[1][i];
            //                betta[0][i+2] += betta[0][i]*alpha[2][i];
            //                betta[0][i+3] += betta[0][i]*alpha[3][i];
            //                betta[0][i+4] += betta[0][i]*alpha[4][i];
            //                gamma[0]      -= betta[0][i]*zetta[i];

            //                //if (fabs(betta[0][i+1])>10.0)
            //                //                {
            //                //                    betta[0][i+2] /= betta[0][i+1];
            //                //                    betta[0][i+3] /= betta[0][i+1];
            //                //                    betta[0][i+4] /= betta[0][i+1];
            //                //                    gamma[0]      /= betta[0][i+1];
            //                //                    betta[0][size-1] /= betta[0][i+1];
            //                //                    betta[0][i+1] /= betta[0][i+1];
            //                //                }
            //            }


            double t4 = (size-5)*h; PointNodeODE node4(t4, static_cast<int>((size-5)));
            double t3 = (size-4)*h; PointNodeODE node3(t3, static_cast<int>((size-4)));
            double t2 = (size-3)*h; PointNodeODE node2(t2, static_cast<int>((size-3)));
            double t1 = (size-2)*h; PointNodeODE node1(t1, static_cast<int>((size-2)));
            double t0 = (size-1)*h; PointNodeODE node0(t0, static_cast<int>((size-1)));

            printf("%20.14f %20.14f %20.14f %20.14f %20.14f\n", node4.x, node3.x, node2.x, node1.x, node0.x);

            DoubleMatrix F(5*M, 5*M);
            DoubleVector f(5*M);


            for (unsigned int r=0; r<M; r++)
            {
                for (unsigned int c=0; c<M; c++)
                {
                    //                    F[0*M+r][0*M+c] = -12.0*h*A(node4, r+1, c+1);
                    //                    F[0*M+r][1*M+c] = 0.0;
                    //                    F[0*M+r][2*M+c] = 0.0;
                    //                    F[0*M+r][3*M+c] = 0.0;
                    //                    F[0*M+r][4*M+c] = 0.0;

                    //                    F[0*M+r][0*M+c] = betta1[r][(size-5)*M+c];
                    //                    F[0*M+r][1*M+c] = betta1[r][(size-4)*M+c];
                    //                    F[0*M+r][2*M+c] = betta1[r][(size-3)*M+c];
                    //                    F[0*M+r][3*M+c] = betta1[r][(size-2)*M+c];
                    //                    F[0*M+r][4*M+c] = betta1[r][(size-1)*M+c];

                    F[0*M+r][0*M+c] = 0.0;
                    F[0*M+r][1*M+c] = betta1[r][(size-4)*M+c];
                    F[0*M+r][2*M+c] = betta1[r][(size-3)*M+c];
                    F[0*M+r][3*M+c] = betta1[r][(size-2)*M+c];
                    F[0*M+r][4*M+c] = betta1[r][(size-1)*M+c];

                    F[1*M+r][0*M+c] = 0.0;
                    F[1*M+r][1*M+c] = -12.0*h*A(node3, r+1, c+1);
                    F[1*M+r][2*M+c] = 0.0;
                    F[1*M+r][3*M+c] = 0.0;
                    F[1*M+r][4*M+c] = 0.0;

                    F[2*M+r][0*M+c] = 0.0;
                    F[2*M+r][1*M+c] = 0.0;
                    F[2*M+r][2*M+c] = -12.0*h*A(node2, r+1, c+1);
                    F[2*M+r][3*M+c] = 0.0;
                    F[2*M+r][4*M+c] = 0.0;

                    F[3*M+r][0*M+c] = 0.0;
                    F[3*M+r][1*M+c] = 0.0;
                    F[3*M+r][2*M+c] = 0.0;
                    F[3*M+r][3*M+c] = -12.0*h*A(node1, r+1, c+1);
                    F[3*M+r][4*M+c] = 0.0;

                    F[4*M+r][0*M+c] = 0.0;
                    F[4*M+r][1*M+c] = 0.0;
                    F[4*M+r][2*M+c] = 0.0;
                    F[4*M+r][3*M+c] = 0.0;
                    F[4*M+r][4*M+c] = -12.0*h*A(node0, r+1, c+1);
                }

                //                f[0*M+r] = +12.0*h*B(node4, r+1);
                f[0*M+r] = gamma1[r];
                f[1*M+r] = +12.0*h*B(node3, r+1);
                f[2*M+r] = +12.0*h*B(node2, r+1);
                f[3*M+r] = +12.0*h*B(node1, r+1);
                f[4*M+r] = +12.0*h*B(node0, r+1);

                //                F[0*M+r][0*M+r] += -25.0;
                //                F[0*M+r][1*M+r] += +48.0;
                //                F[0*M+r][2*M+r] += -36.0;
                //                F[0*M+r][3*M+r] += +16.0;
                //                F[0*M+r][4*M+r] += -3.0;

                F[1*M+r][0*M+r] += -3.0;
                F[1*M+r][1*M+r] += -10.0;
                F[1*M+r][2*M+r] += +18.0;
                F[1*M+r][3*M+r] += -6.0;
                F[1*M+r][4*M+r] += +1.0;

                F[2*M+r][0*M+r] += +1.0;
                F[2*M+r][1*M+r] += -8.0;
                F[2*M+r][2*M+r] += +0.0;
                F[2*M+r][3*M+r] += +8.0;
                F[2*M+r][4*M+r] += -1.0;

                F[3*M+r][0*M+r] += -1.0;
                F[3*M+r][1*M+r] += +6.0;
                F[3*M+r][2*M+r] += -18.0;
                F[3*M+r][3*M+r] += +10.0;
                F[3*M+r][4*M+r] += +3.0;

                F[4*M+r][0*M+r] += +3.0;
                F[4*M+r][1*M+r] += -16.0;
                F[4*M+r][2*M+r] += +36.0;
                F[4*M+r][3*M+r] += -48.0;
                F[4*M+r][4*M+r] += +25.0;
            }

            //            F[0][0] = betta[0][size-5];
            //            F[0][1] = betta[0][size-4];
            //            F[0][2] = betta[0][size-3];
            //            F[0][3] = betta[0][size-2];
            //            F[0][4] = betta[0][size-1];

            //        F[0][0] = -25.0 - 12.0*h*A(node4, 1, 1);
            //        F[0][1] = +48.0;
            //        F[0][2] = -36.0;
            //        F[0][3] = +16.0;
            //        F[0][4] = -3.0;

            //            F[1][0] = -3.0;
            //            F[1][1] = -10.0 - 12.0*h*A(node3, 1, 1);
            //            F[1][2] = +18.0;
            //            F[1][3] = -6.0;
            //            F[1][4] = +1.0;

            //            F[2][0] = +1.0;
            //            F[2][1] = -8.0;
            //            F[2][2] = +0.0 - 12.0*h*A(node2, 1, 1);
            //            F[2][3] = +8.0;
            //            F[2][4] = -1.0;

            //            F[3][0] = -1.0;
            //            F[3][1] = +6.0;
            //            F[3][2] = -18.0;
            //            F[3][3] = +10.0 - 12.0*h*A(node1, 1, 1);
            //            F[3][4] = +3.0;

            //            F[4][0] = +3.0;
            //            F[4][1] = -16.0;
            //            F[4][2] = +36.0;
            //            F[4][3] = -48.0;
            //            F[4][4] = +25.0 - 12.0*h*A(node0, 1, 1);


            //            f[0] = gamma[0];
            //            //        f[0] = +12.0*h*B(node4, 1);
            //            f[1] = +12.0*h*B(node3, 1);
            //            f[2] = +12.0*h*B(node2, 1);
            //            f[3] = +12.0*h*B(node1, 1);
            //            f[4] = +12.0*h*B(node0, 1);

            IPrinter::printSeperatorLine();
            IPrinter::print(F, F.rows(), F.cols(), 12, 6);
            IPrinter::printSeperatorLine();
            IPrinter::print(f, f.length(), 12, 6);
            IPrinter::printSeperatorLine();

            DoubleVector xf((k+1)*M);
            LinearEquation::GaussianElimination(F, f, xf);

            //    double a = 0.53582679*betta[size-5][0][0]+0.72896863*betta[size-4][0][0]+0.87630668*betta[size-3][0][0]+0.96858316*betta[size-2][0][0]+1.00000000*betta[size-1][0][0];
            //    double c =      xf[0]*betta[size-5][0][0]+     xf[1]*betta[size-4][0][0]+     xf[2]*betta[size-3][0][0]+     xf[3]*betta[size-2][0][0]+     xf[4]*betta[size-1][0][0];
            //    double b = gamma[0][0];
            //    printf("----- %14.8f %14.8f %14.8f\n", a, c, b);
            //    printf("***** %14.8f %14.8f %14.8f %14.8f %14.8f\n", xf[4], xf[3], xf[2], xf[1], xf[0]);


            F.clear();
            f.clear();

            IPrinter::printSeperatorLine();
            IPrinter::print(xf, xf.length());
            IPrinter::printSeperatorLine();

            for (unsigned int i=0; i<=k; i++)
            {
                printf("***** ");
                for (unsigned int r=1; r<=M; r++)
                {
                    printf("%4d %14.8f ", i*M+M-r, xf[i*M+M-r]);
                }
                printf("\n");
            }

            x.clear();
            x.resize(size); for (unsigned int n=0; n<size; n++) x[n].resize(M);

            unsigned int s = xf.length()-M;
            unsigned int e = xf.length()-1;

            IPrinter::printSeperatorLine();
            for (unsigned int n=size-1; n>=size-(k+1); n--)
            {
                x[n] = xf.mid(s, e);
                s -= M;
                e -= M;
                printf("%4d : ", n);
                IPrinter::print(x[n], x[n].length());
            }
            IPrinter::printSeperatorLine();
            xf.clear();

            if (k==4)
            {
                unsigned int rw;
                rw = 286;
                x[95][2] = zetta1[rw]
                        + alpha1[rw][3]*x[96][0]
                        + alpha1[rw][4]*x[96][1]
                        + alpha1[rw][5]*x[96][2]

                        + alpha1[rw][6]*x[97][0]
                        + alpha1[rw][7]*x[97][1]
                        + alpha1[rw][8]*x[97][2]

                        + alpha1[rw][9]*x[98][0]
                        + alpha1[rw][10]*x[98][1]
                        + alpha1[rw][11]*x[98][2]

                        + alpha1[rw][12]*x[99][0]
                        + alpha1[rw][13]*x[99][1]
                        + alpha1[rw][14]*x[99][2];
                rw = 284;
                x[95][1] = zetta1[rw]
                        + alpha1[rw][3]*x[95][2]

                        + alpha1[rw][4]*x[96][0]
                        + alpha1[rw][5]*x[96][1]
                        + alpha1[rw][6]*x[96][2]

                        + alpha1[rw][7]*x[97][0]
                        + alpha1[rw][8]*x[97][1]
                        + alpha1[rw][9]*x[97][2]

                        + alpha1[rw][10]*x[98][0]
                        + alpha1[rw][11]*x[98][1]
                        + alpha1[rw][12]*x[98][2]

                        + alpha1[rw][13]*x[99][0]
                        + alpha1[rw][14]*x[99][1];
                rw = 283;
                x[95][0] = zetta1[rw]
                        + alpha1[rw][3]*x[95][1]
                        + alpha1[rw][4]*x[95][2]

                        + alpha1[rw][5]*x[96][0]
                        + alpha1[rw][6]*x[96][1]
                        + alpha1[rw][7]*x[96][2]

                        + alpha1[rw][8]*x[97][0]
                        + alpha1[rw][9]*x[97][1]
                        + alpha1[rw][10]*x[97][2]

                        + alpha1[rw][11]*x[98][0]
                        + alpha1[rw][12]*x[98][1]
                        + alpha1[rw][13]*x[98][2]

                        + alpha1[rw][14]*x[99][0];

                IPrinter::print(x[95], x[95].length());


                for (unsigned int i=size-(k+1); i!=0; i--)
                {
                    unsigned int rw = (i-1)*M;

                    for (unsigned int r=0; r<M; r++)
                    {
                        x[i-1][r] = alpha1[(i-1)][1]*x[i+0][0] + alpha1[2][i-1]*x[i+1][0] + alpha1[3][i-1]*x[i+2][0] + alpha1[4][i-1]*x[i+3][0] + zetta1[i-1];
                    }


                    //printf("************************  %4d %14.8f \n", i-1, x[i-1][0]);
                }
            }

        } // schema == 11

        if (schema == 2)
        {
            for (unsigned int i=size-1; i>=5; i--)
            {
                const unsigned int i0 = i+0; PointNodeODE node0(i0*h, static_cast<int>(i0));
                const unsigned int i1 = i-1; PointNodeODE node1(i1*h, static_cast<int>(i1));
                const unsigned int i2 = i-2; PointNodeODE node2(i2*h, static_cast<int>(i2));
                const unsigned int i3 = i-3; PointNodeODE node3(i3*h, static_cast<int>(i3));
                const unsigned int i4 = i-4; PointNodeODE node4(i4*h, static_cast<int>(i4));

                alpha[0][i] = +1.0;
                alpha[1][i] = +0.8 + 2.4*h*A(node1, 1, 1);
                alpha[2][i] = +0.0 - 2.4*h*A(node2, 1, 1);
                alpha[3][i] = -0.8 + 2.4*h*A(node3, 1, 1);
                alpha[4][i] = +1.0;
                zetta[i] = +2.4*h*B(node1, 1) - 2.4*h*B(node2, 1) + 2.4*h*B(node3, 1);
            }

            for (unsigned int i=size-1; i>=5; i--)
            {
                //            if (i>99000)
                //            {
                //printf("%4d | %25.14f %25.14f %25.14f %25.14f %25.14f | %25.14f\n", i, alpha[4][i], alpha[3][i], alpha[2][i], alpha[1][i], alpha[0][i], zetta[i]);
                //printf("%4d | %25.14f %25.14f %25.14f %25.14f|%25.14f\n", i, betta[0][i-3], betta[0][i-2], betta[0][i-1], betta[0][i-0], gamma[0]);
                //            }

                betta[0][i-1] += betta[0][i]*alpha[1][i];
                betta[0][i-2] += betta[0][i]*alpha[2][i];
                betta[0][i-3] += betta[0][i]*alpha[3][i];
                betta[0][i-4] += betta[0][i]*alpha[4][i];
                gamma[0]      -= betta[0][i]*zetta[i];

                if (fabs(betta[0][i-1])>1000.0)
                {
                    betta[0][i-2] /= betta[0][i-1];
                    betta[0][i-3] /= betta[0][i-1];
                    betta[0][i-4] /= betta[0][i-1];
                    gamma[0]      /= betta[0][i-1];
                    betta[0][0] /= betta[0][i-1];
                    betta[0][i-1] /= betta[0][i-1];
                }
            }

            double t4 = (0)*h; PointNodeODE node4(t4, static_cast<int>((0)));
            double t3 = (1)*h; PointNodeODE node3(t3, static_cast<int>((1)));
            double t2 = (2)*h; PointNodeODE node2(t2, static_cast<int>((2)));
            double t1 = (3)*h; PointNodeODE node1(t1, static_cast<int>((3)));
            double t0 = (4)*h; PointNodeODE node0(t0, static_cast<int>((4)));

            DoubleMatrix F(5, 5);
            DoubleVector f(5);

            F[0][0] = betta[0][0];
            F[0][1] = betta[0][1];
            F[0][2] = betta[0][2];
            F[0][3] = betta[0][3];
            F[0][4] = betta[0][4];

            //F[0][0] = -25.0 - 12.0*h*A(node4, 1, 1);
            //F[0][1] = +48.0;
            //F[0][2] = -36.0;
            //F[0][3] = +16.0;
            //F[0][4] = -3.0;

            F[1][0] = -3.0;
            F[1][1] = -10.0 - 12.0*h*A(node3, 1, 1);
            F[1][2] = +18.0;
            F[1][3] = -6.0;
            F[1][4] = +1.0;

            F[2][0] = +1.0;
            F[2][1] = -8.0;
            F[2][2] = +0.0 - 12.0*h*A(node2, 1, 1);
            F[2][3] = +8.0;
            F[2][4] = -1.0;

            F[3][0] = -1.0;
            F[3][1] = +6.0;
            F[3][2] = -18.0;
            F[3][3] = +10.0 - 12.0*h*A(node1, 1, 1);
            F[3][4] = +3.0;

            F[4][0] = +3.0;
            F[4][1] = -16.0;
            F[4][2] = +36.0;
            F[4][3] = -48.0;
            F[4][4] = +25.0 - 12.0*h*A(node0, 1, 1);

            f[0] = gamma[0];
            //f[0] = +12.0*h*B(node4, 1);
            f[1] = +12.0*h*B(node3, 1);
            f[2] = +12.0*h*B(node2, 1);
            f[3] = +12.0*h*B(node1, 1);
            f[4] = +12.0*h*B(node0, 1);

            IPrinter::printSeperatorLine();
            IPrinter::print(F, F.rows(), F.cols(), 12, 6);
            IPrinter::printSeperatorLine();

            DoubleVector xf((k+1)*M);
            LinearEquation::GaussianElimination(F, f, xf);

            //    double a = 0.53582679*betta[size-5][0][0]+0.72896863*betta[size-4][0][0]+0.87630668*betta[size-3][0][0]+0.96858316*betta[size-2][0][0]+1.00000000*betta[size-1][0][0];
            //    double c =      xf[0]*betta[size-5][0][0]+     xf[1]*betta[size-4][0][0]+     xf[2]*betta[size-3][0][0]+     xf[3]*betta[size-2][0][0]+     xf[4]*betta[size-1][0][0];
            //    double b = gamma[0][0];
            //    printf("----- %14.8f %14.8f %14.8f\n", a, c, b);
            //    printf("***** %14.8f %14.8f %14.8f %14.8f %14.8f\n", xf[4], xf[3], xf[2], xf[1], xf[0]);


            F.clear();
            f.clear();

            for (unsigned int i=0; i<=k; i++)
            {
                printf("***** ");
                for (unsigned int r=1; r<=M; r++)
                {
                    printf("%4d %14.8f ", i*M+M-r, xf[i*M+M-r]);
                }
                printf("\n");
            }

            x.clear();
            x.resize(size); for (unsigned int n=0; n<size; n++) x[n].resize(M);

            unsigned int s = xf.length()-M;
            unsigned int e = xf.length()-1;

            for (unsigned int n=size-1; n>=size-k-1; n--)
            {
                x[n] = xf.mid(s, e);
                s -= M;
                e -= M;
            }
            xf.clear();

            if (k==4)
            {
                for (unsigned int i=5; i!=size-5; i++)
                {
                    x[i][0] = alpha[1][i]*x[i-1][0]
                            + alpha[2][i]*x[i-2][0]
                            + alpha[3][i]*x[i-4][0]
                            + alpha[4][i]*x[i-4][0] + zetta[i];
                    printf("************************  %4d %14.8f \n", i, x[i][0]);
                }
            }

        }
    }
}

void IFirstOrderLinearODEIVP::transferOfConditionP(const std::vector<NonLocalCondition> &co, const DoubleVector &d, std::vector<DoubleVector> &x, unsigned int k, unsigned int schema) const
{
    if (co.size() < 2) throw ExceptionODE(1);

    for (unsigned int i=0; i<co.size(); i++)
    {
        const DoubleMatrix &m = co[i].m;

        if (!m.squareMatrix()) throw ExceptionODE(2);

        if (m.rows() != count()) throw ExceptionODE(3);

        if (d.length() != count()) throw ExceptionODE(4);
    }

    const int min = dimension().min();
    const int max = dimension().max();
    const double h = dimension().step();
    const unsigned int N = static_cast<unsigned int>(max - min);
    const unsigned int size = static_cast<unsigned int>( dimension().size() );
    const unsigned int end = static_cast<unsigned int>( size-(k+1) );
    const unsigned int M = count();

    std::vector<NonLocalCondition> C = co;
    //discritize(co, C);

    unsigned int sizeM = size*M;
    unsigned int sizeK = (k+1)*M*sizeof(double);

    // begin allocating memory

    double **alpha = static_cast<double**>( malloc( sizeof(double*) * sizeM ) );
    for (unsigned int i=0; i<sizeM; i++)
    {
        alpha[i] = static_cast<double*>( malloc( sizeK ) );
    }
    double *betta  = static_cast<double*> ( malloc( sizeof(double)  * sizeM ) );

    double  *delta = static_cast< double*  >( malloc(sizeof(double)  * M) );
    double **gamma = static_cast< double** >( malloc(sizeof(double*) * M) );
    for (unsigned int r=0; r<M; r++)
    {
        gamma[r] = static_cast<double*> ( malloc( sizeof(double) * sizeM ) );
        memset (gamma[r],0,sizeof(double) * sizeM);
        delta[r] = d[r];
    }

    //--------------------------------------

    // setting non-separated conditions

    for (unsigned int s=0; s<C.size(); s++)
    {
        unsigned int i = static_cast<unsigned int>(C[s].n.i);
        unsigned int offset = i*M;
        for (unsigned int r=0; r<M; r++)
        {
            for (unsigned int c=0; c<M; c++)
            {
                gamma[r][offset+c] += C[s].m[r][c];
            }
        }
    }

    //--------------------------------------

    if (k==2)
    {
        const unsigned int Mx0=0, Mx1=M, Mx2=2*M;

        for (unsigned int i=0; i<end; i++)
        {
            const PointNodeODE node0((i+0)*h, static_cast<int>(i+0));
            const PointNodeODE node1((i+1)*h, static_cast<int>(i+1));
            const PointNodeODE node2((i+2)*h, static_cast<int>(i+2));

            const unsigned int Mi0=i*M, Mi1=Mi0+M, Mi2=Mi1+M;
            for (unsigned int r=0; r<M; r++)
            {
                const unsigned int rw = i*M+r;

                for (unsigned int c=0; c<M; c++)
                {
                    alpha[rw][Mx0+c] = +0.0;
                    alpha[rw][Mx1+c] = -1.6*h*A(node1, r+1, c+1);
                    alpha[rw][Mx2+c] = +0.4*h*A(node2, r+1, c+1);

                    //alpha[rw][Mx0+c] = +0.0;
                    //alpha[rw][Mx1+c] = -1.8*h*A(node1, r+1, c+1);
                    //alpha[rw][Mx2+c] = +0.2*h*A(node2, r+1, c+1);
                }

                alpha[rw][Mx0+r] += +1.0;
                alpha[rw][Mx1+r] += +0.8;
                alpha[rw][Mx2+r] += +0.2;
                betta[rw] = -0.4*h*(4.0*B(node1, r+1)-B(node2, r+1));

                //alpha[rw][Mx0+r] += +1.0;
                //alpha[rw][Mx1+r] += +0.4;
                //alpha[rw][Mx2+r] += +0.6;
                //betta[rw] = -0.2*h*(9.0*B(node1, r+1)-B(node2, r+1));

                for (unsigned int r1=0; r1<M; r1++)
                {
                    for (unsigned int c=0; c<M; c++)
                    {
                        gamma[r1][Mi1+c] += gamma[r1][rw]*alpha[rw][Mx1+c];
                        gamma[r1][Mi2+c] += gamma[r1][rw]*alpha[rw][Mx2+c];
                    }
                    delta[r1] -= gamma[r1][rw]*betta[rw];
                }
            }
        }

        const unsigned int N0=N, N1=N-1, N2=N-2;
        const PointNodeODE nodeN2(N2*h, static_cast<int>(N2));
        const PointNodeODE nodeN1(N1*h, static_cast<int>(N1));
        const PointNodeODE nodeN0(N0*h, static_cast<int>(N0));

        DoubleMatrix Mx((k+1)*M, (k+1)*M);
        DoubleVector f((k+1)*M);

        for (unsigned int r=0; r<M; r++)
        {
            for (unsigned int c=0; c<M; c++)
            {
                Mx[Mx0+r][Mx0+c] = -2.0*h*A(nodeN2, r+1, c+1);
                Mx[Mx0+r][Mx1+c] = +0.0;
                Mx[Mx0+r][Mx2+c] = +0.0;

                Mx[Mx1+r][Mx0+c] = 0.0;
                Mx[Mx1+r][Mx1+c] = -2.0*h*A(nodeN1, r+1, c+1);
                Mx[Mx1+r][Mx2+c] = 0.0;

                Mx[Mx2+r][Mx0+c] = gamma[r][N2*M+c];
                Mx[Mx2+r][Mx1+c] = gamma[r][N1*M+c];
                Mx[Mx2+r][Mx2+c] = gamma[r][N0*M+c];
            }

            f[Mx0+r] = +2.0*h*B(nodeN2, r+1);
            f[Mx1+r] = +2.0*h*B(nodeN1, r+1);
            f[Mx2+r] = delta[r];

            Mx[Mx0+r][Mx0+r] += -3.0; Mx[Mx0+r][Mx1+r] += +4.0; Mx[Mx0+r][Mx2+r] += -1.0;
            Mx[Mx1+r][Mx0+r] += -1.0; Mx[Mx1+r][Mx1+r] += +0.0; Mx[Mx1+r][Mx2+r] += +1.0;
        }




        ///////////////////////////////////////////////////////////////////////////
        IPrinter::printSeperatorLine("Mx");
        IPrinter::print(Mx, Mx.rows(), Mx.cols(), 12, 6);
        IPrinter::printSeperatorLine("f");
        IPrinter::print(f, f.length(), 12, 6);

        DoubleVector xf((k+1)*M);
        LinearEquation::GaussianElimination(Mx, f, xf);

        Mx.clear();
        f.clear();

        IPrinter::printSeperatorLine("x");
        IPrinter::print(xf, xf.length(), 12, 6);
        IPrinter::printSeperatorLine();

        for (unsigned int i=0; i<=k; i++)
        {
            printf("***** ");
            for (unsigned int r=1; r<=M; r++)
            {
                printf("%4d %20.16f ", i*M+M-r, xf[i*M+M-r]);
            }
            printf("\n");
        }

        x.clear();
        x.resize(size); for (unsigned int n=0; n<size; n++) x[n].resize(M);

        unsigned int s = xf.length()-M;
        unsigned int e = xf.length()-1;

        IPrinter::printSeperatorLine();
        for (unsigned int n=size-1; n>=end; n--)
        {
            x[n] = xf.mid(s, e);
            s -= M;
            e -= M;
            printf("%6d : ", n); IPrinter::print(x[n], x[n].length(), 28, 8);
        }
        IPrinter::printSeperatorLine();
        xf.clear();

        for (unsigned int j=0, i=N-(k+1); j<=N-(k+1); j++, i--)
        {
            for (unsigned int r=0; r<M; r++)
            {
                unsigned int rw = i*M+r;

                x[i][r] = betta[rw];
                for (unsigned int c=0; c<M; c++)
                {
                    x[i][r] += alpha[rw][Mx1+c]*x[i+1][c] + alpha[rw][Mx2+c]*x[i+2][c];
                }
            }
            if (i%(N/5)==0) { printf("%6d : ", i); IPrinter::print(x[i], x[i].length(), 28, 8); }
        }
        IPrinter::printSeperatorLine("=");
        ///////////////////////////////////////////////////////////////////////////
    }

    if (k==4)
    {
        const unsigned int Mx0=0, Mx1=M, Mx2=2*M, Mx3=3*M, Mx4=4*M;

        for (unsigned int i=0; i<end; i++)
        {
            const PointNodeODE node0((i+0)*h, static_cast<int>(i+0));
            const PointNodeODE node1((i+1)*h, static_cast<int>(i+1));
            const PointNodeODE node2((i+2)*h, static_cast<int>(i+2));
            const PointNodeODE node3((i+3)*h, static_cast<int>(i+3));
            const PointNodeODE node4((i+4)*h, static_cast<int>(i+4));

            const unsigned int Mi0=i*M, Mi1=Mi0+M, Mi2=Mi1+M, Mi3=Mi2+M, Mi4=Mi3+M;

            for (unsigned int r=0; r<M; r++)
            {
                const unsigned int rw = i*M+r;
                //double *pAlpha = alpha[rw];

                for (unsigned int c=0; c<M; c++)
                {
                    alpha[rw][Mx0+c] = +0.0;
                    alpha[rw][Mx1+c] = -2.4*h*A(node1, r+1, c+1);
                    alpha[rw][Mx2+c] = +2.4*h*A(node2, r+1, c+1);
                    alpha[rw][Mx3+c] = -2.4*h*A(node3, r+1, c+1);
                    alpha[rw][Mx4+c] = +0.0;

                    //alpha[rw][Mx0+c] = +0.0;
                    //alpha[rw][Mx1+c] = -2.6666666666666667*h*A(node1, r+1, c+1);
                    //alpha[rw][Mx2+c] = +1.3333333333333333*h*A(node2, r+1, c+1);
                    //alpha[rw][Mx3+c] = -2.6666666666666667*h*A(node3, r+1, c+1);
                    //alpha[rw][Mx4+c] = +0.0;

                    //alpha[rw][Mx0+c] = +0.0;
                    //alpha[rw][Mx1+c] = -2.6*h*A(node1, r+1, c+1);
                    //alpha[rw][Mx2+c] = +2.0*h*A(node2, r+1, c+1);
                    //alpha[rw][Mx3+c] = -2.2*h*A(node3, r+1, c+1);
                    //alpha[rw][Mx4+c] = +0.0;
                }

                alpha[rw][Mx0+r] += +1.0;
                alpha[rw][Mx1+r] += +0.8;
                alpha[rw][Mx2+r] += +0.0;
                alpha[rw][Mx3+r] += -0.8;
                alpha[rw][Mx4+r] += +1.0;
                betta[rw] = -2.4*h*(B(node1, r+1) - B(node2, r+1) + B(node3, r+1));

                //alpha[rw][Mx0+r] += +1.0;
                //alpha[rw][Mx1+r] += +0.0;
                //alpha[rw][Mx2+r] += +0.0;
                //alpha[rw][Mx3+r] += -0.0;
                //alpha[rw][Mx4+r] += +1.0;
                //betta[rw] = -1.3333333333333333*h*(2.0*B(node1, r+1)-B(node2, r+1)+2.0*B(node3, r+1));

                //alpha[rw][Mx0+r] += +1.0;
                //alpha[rw][Mx1+r] += +0.266666666666666667;
                //alpha[rw][Mx2+r] += +0.6;
                //alpha[rw][Mx3+r] += -0.8;
                //alpha[rw][Mx4+r] += +0.933333333333333333;
                //betta[rw] = -2.6*h*B(node1, r+1) + 2.0*h*B(node2, r+1) - 2.2*h*B(node3, r+1);

                for (unsigned int r1=0; r1<M; r1++)
                {
                    for (unsigned int c=0; c<M; c++)
                    {
                        gamma[r1][Mi1+c] += gamma[r1][rw]*alpha[rw][Mx1+c];
                        gamma[r1][Mi2+c] += gamma[r1][rw]*alpha[rw][Mx2+c];
                        gamma[r1][Mi3+c] += gamma[r1][rw]*alpha[rw][Mx3+c];
                        gamma[r1][Mi4+c] += gamma[r1][rw]*alpha[rw][Mx4+c];
                    }
                    delta[r1] -= gamma[r1][rw]*betta[rw];
                }

                // normalizing...

                //double mm = 0.0;
                //for (unsigned int r1=0; r1<M; r1++)
                //{
                //    for (unsigned int c=0; c<M; c++)
                //    {
                //        if (mm<fabs(gamma[r1][Mi1+c])) { mm = fabs(gamma[r1][Mi1+c]); }
                //        if (mm<fabs(gamma[r1][Mi2+c])) { mm = fabs(gamma[r1][Mi2+c]); }
                //        if (mm<fabs(gamma[r1][Mi3+c])) { mm = fabs(gamma[r1][Mi3+c]); }
                //        if (mm<fabs(gamma[r1][Mi4+c])) { mm = fabs(gamma[r1][Mi4+c]); }
                //    }
                //    if (mm<fabs(delta[r1])) { mm = fabs(delta[r1]); }
                //}

                //if ( mm > 1000.0 )
                //{
                //    for (unsigned int r1=0; r1<M; r1++)
                //    {
                //        for (unsigned int c=0; c<M; c++)
                //        {
                //            gamma[r1][Mi1+c] /= mm;
                //            gamma[r1][Mi2+c] /= mm;
                //            gamma[r1][Mi3+c] /= mm;
                //            gamma[r1][Mi4+c] /= mm;
                //        }
                //        delta[r1] /= mm;
                //        gamma[r1][N*M] /= mm;
                //    }
                //}
            }
        }

        const unsigned int N0=N, N1=N-1, N2=N-2, N3=N-3, N4=N-4;
        const PointNodeODE nodeN4(N4*h, static_cast<int>(N4));
        const PointNodeODE nodeN3(N3*h, static_cast<int>(N3));
        const PointNodeODE nodeN2(N2*h, static_cast<int>(N2));
        const PointNodeODE nodeN1(N1*h, static_cast<int>(N1));
        const PointNodeODE nodeN0(N0*h, static_cast<int>(N0));

        DoubleMatrix Mx((k+1)*M, (k+1)*M);
        DoubleVector f((k+1)*M);

        for (unsigned int r=0; r<M; r++)
        {
            for (unsigned int c=0; c<M; c++)
            {
                //Mx[Mx0+r][Mx0+c] = -25.0 - 12.0*h*A(nodeN4, r+1, c+1);;
                //Mx[Mx0+r][Mx1+c] = +48.0;
                //Mx[Mx0+r][Mx2+c] = -36.0;
                //Mx[Mx0+r][Mx3+c] = +16.0;
                //Mx[Mx0+r][Mx4+c] = -3.0;

                //Mx[Mx0+r][Mx0+c] = gamma[r][N0*M+c];
                //Mx[Mx0+r][Mx1+c] = gamma[r][N3*M+c];
                //Mx[Mx0+r][Mx2+c] = gamma[r][N2*M+c];
                //Mx[Mx0+r][Mx3+c] = gamma[r][N1*M+c];
                //Mx[Mx0+r][Mx4+c] = gamma[r][N0*M+c];

                Mx[Mx0+r][Mx0+c] = Mx[Mx0+r][Mx1+c] = Mx[Mx0+r][Mx2+c] = Mx[Mx0+r][Mx3+c] = 0.0; Mx[Mx0+r][Mx4+c] = -12.0*h*A(nodeN0, r+1, c+1);
                Mx[Mx1+r][Mx0+c] = Mx[Mx1+r][Mx1+c] = Mx[Mx1+r][Mx2+c] = Mx[Mx1+r][Mx4+c] = 0.0; Mx[Mx1+r][Mx3+c] = -12.0*h*A(nodeN1, r+1, c+1);
                Mx[Mx2+r][Mx0+c] = Mx[Mx2+r][Mx1+c] = Mx[Mx2+r][Mx3+c] = Mx[Mx2+r][Mx4+c] = 0.0; Mx[Mx2+r][Mx2+c] = -12.0*h*A(nodeN2, r+1, c+1);
                Mx[Mx3+r][Mx0+c] = Mx[Mx3+r][Mx2+c] = Mx[Mx3+r][Mx3+c] = Mx[Mx3+r][Mx4+c] = 0.0; Mx[Mx3+r][Mx1+c] = -12.0*h*A(nodeN3, r+1, c+1);

                Mx[Mx4+r][Mx0+c] = gamma[r][N4*M+c];
                Mx[Mx4+r][Mx1+c] = gamma[r][N3*M+c];
                Mx[Mx4+r][Mx2+c] = gamma[r][N2*M+c];
                Mx[Mx4+r][Mx3+c] = gamma[r][N1*M+c];
                Mx[Mx4+r][Mx4+c] = gamma[r][N0*M+c];
            }

            //f[0*M+r] = +12.0*h*B(nodeN4, r+1);
            f[Mx0+r] = +12.0*h*B(nodeN0, r+1);
            f[Mx1+r] = +12.0*h*B(nodeN1, r+1);
            f[Mx2+r] = +12.0*h*B(nodeN2, r+1);
            f[Mx3+r] = +12.0*h*B(nodeN3, r+1);
            f[Mx4+r] = delta[r];

            Mx[Mx0+r][Mx0+r] += +3.0; Mx[Mx0+r][Mx1+r] += -16.0; Mx[Mx0+r][Mx2+r] += +36.0; Mx[Mx0+r][Mx3+r] += -48.0; Mx[Mx0+r][Mx4+r] += +25.0;
            Mx[Mx1+r][Mx0+r] += -1.0; Mx[Mx1+r][Mx1+r] += +6.00; Mx[Mx1+r][Mx2+r] += -18.0; Mx[Mx1+r][Mx3+r] += +10.0; Mx[Mx1+r][Mx4+r] += +3.00;
            Mx[Mx2+r][Mx0+r] += +1.0; Mx[Mx2+r][Mx1+r] += -8.00; Mx[Mx2+r][Mx2+r] += +0.00; Mx[Mx2+r][Mx3+r] += +8.00; Mx[Mx2+r][Mx4+r] += -1.00;
            Mx[Mx3+r][Mx0+r] += -3.0; Mx[Mx3+r][Mx1+r] += -10.0; Mx[Mx3+r][Mx2+r] += +18.0; Mx[Mx3+r][Mx3+r] += -6.00; Mx[Mx3+r][Mx4+r] += +1.00;
        }


        ///////////////////////////////////////////////////////////////////////////
        IPrinter::printSeperatorLine("Mx");
        IPrinter::print(Mx, Mx.rows(), Mx.cols(), 12, 6);
        IPrinter::printSeperatorLine("f");
        IPrinter::print(f, f.length(), 12, 6);

        DoubleVector xf((k+1)*M);
        LinearEquation::GaussianElimination(Mx, f, xf);

        Mx.clear();
        f.clear();

        IPrinter::printSeperatorLine("x");
        IPrinter::print(xf, xf.length(), 12, 6);
        IPrinter::printSeperatorLine();

        for (unsigned int i=0; i<=k; i++)
        {
            printf("***** ");
            for (unsigned int r=1; r<=M; r++)
            {
                printf("%4d %20.16f ", i*M+M-r, xf[i*M+M-r]);
            }
            printf("\n");
        }

        x.clear();
        x.resize(size); for (unsigned int n=0; n<size; n++) x[n].resize(M);

        unsigned int s = xf.length()-M;
        unsigned int e = xf.length()-1;

        IPrinter::printSeperatorLine();
        for (unsigned int n=size-1; n>=end; n--)
        {
            x[n] = xf.mid(s, e);
            s -= M;
            e -= M;
            printf("%6d : ", n); IPrinter::print(x[n], x[n].length(), 28, 8);
        }
        IPrinter::printSeperatorLine();
        xf.clear();

        for (unsigned int j=0, i=N-(k+1); j<=N-(k+1); j++, i--)
        {
            for (unsigned int eq=0; eq<M; eq++)
            {
                unsigned int rw = i*M+eq;

                x[i][eq] = betta[rw];
                for (unsigned int c=0; c<M; c++)
                {
                    x[i][eq] += alpha[rw][Mx1+c]*x[i+1][c] + alpha[rw][Mx2+c]*x[i+2][c]
                            + alpha[rw][Mx3+c]*x[i+3][c] + alpha[rw][Mx4+c]*x[i+4][c];
                }
            }
            if (i%(N/5)==0) { printf("%6d : ", i); IPrinter::print(x[i], x[i].length(), 28, 8); }
        }
        IPrinter::printSeperatorLine("=");
        ///////////////////////////////////////////////////////////////////////////
    }

    if (k==6)
    {
        const unsigned int Mx0=0, Mx1=M, Mx2=2*M, Mx3=3*M, Mx4=4*M, Mx5=5*M, Mx6=6*M;

        for (unsigned int i=0; i<end; i++)
        {
            const PointNodeODE node0((i+0)*h, static_cast<int>(i+0));
            const PointNodeODE node1((i+1)*h, static_cast<int>(i+1));
            const PointNodeODE node2((i+2)*h, static_cast<int>(i+2));
            const PointNodeODE node3((i+3)*h, static_cast<int>(i+3));
            const PointNodeODE node4((i+4)*h, static_cast<int>(i+4));
            const PointNodeODE node5((i+5)*h, static_cast<int>(i+5));
            const PointNodeODE node6((i+6)*h, static_cast<int>(i+6));

            const unsigned int Mi0=i*M, Mi1=Mi0+M, Mi2=Mi1+M, Mi3=Mi2+M, Mi4=Mi3+M, Mi5=Mi4+M, Mi6=Mi5+M;

            for (unsigned int r=0; r<M; r++)
            {
                const unsigned int rw = i*M+r;
                //double *pAlpha = alpha[rw];

                for (unsigned int c=0; c<M; c++)
                {
                    alpha[rw][Mx0+c] = +0.0;
                    alpha[rw][Mx1+c] = -2.4*h*A(node1, r+1, c+1);
                    alpha[rw][Mx2+c] = +2.4*h*A(node2, r+1, c+1);
                    alpha[rw][Mx3+c] = -2.4*h*A(node3, r+1, c+1);
                    alpha[rw][Mx4+c] = +0.0;
                    alpha[rw][Mx5+c] = -2.4*h*A(node3, r+1, c+1);
                    alpha[rw][Mx6+c] = +0.0;
                }

                alpha[rw][Mx0+r] += +1.0;
                alpha[rw][Mx1+r] += +0.8;
                alpha[rw][Mx2+r] += +0.0;
                alpha[rw][Mx3+r] += -0.8;
                alpha[rw][Mx4+r] += +1.0;
                alpha[rw][Mx5+r] += -0.8;
                alpha[rw][Mx6+r] += +1.0;
                betta[rw] = -2.4*h*(B(node1, r+1) - B(node2, r+1) + B(node3, r+1));

                for (unsigned int r1=0; r1<M; r1++)
                {
                    for (unsigned int c=0; c<M; c++)
                    {
                        gamma[r1][Mi1+c] += gamma[r1][rw]*alpha[rw][Mx1+c];
                        gamma[r1][Mi2+c] += gamma[r1][rw]*alpha[rw][Mx2+c];
                        gamma[r1][Mi3+c] += gamma[r1][rw]*alpha[rw][Mx3+c];
                        gamma[r1][Mi4+c] += gamma[r1][rw]*alpha[rw][Mx4+c];
                        gamma[r1][Mi5+c] += gamma[r1][rw]*alpha[rw][Mx5+c];
                        gamma[r1][Mi6+c] += gamma[r1][rw]*alpha[rw][Mx6+c];
                    }
                    delta[r1] -= gamma[r1][rw]*betta[rw];
                }

                // normalizing...

                //double mm = 0.0;
                //for (unsigned int r1=0; r1<M; r1++)
                //{
                //    for (unsigned int c=0; c<M; c++)
                //    {
                //        if (mm<fabs(gamma[r1][Mi1+c])) { mm = fabs(gamma[r1][Mi1+c]); }
                //        if (mm<fabs(gamma[r1][Mi2+c])) { mm = fabs(gamma[r1][Mi2+c]); }
                //        if (mm<fabs(gamma[r1][Mi3+c])) { mm = fabs(gamma[r1][Mi3+c]); }
                //        if (mm<fabs(gamma[r1][Mi4+c])) { mm = fabs(gamma[r1][Mi4+c]); }
                //    }
                //    if (mm<fabs(delta[r1])) { mm = fabs(delta[r1]); }
                //}

                //if ( mm > 1000.0 )
                //{
                //    for (unsigned int r1=0; r1<M; r1++)
                //    {
                //        for (unsigned int c=0; c<M; c++)
                //        {
                //            gamma[r1][Mi1+c] /= mm;
                //            gamma[r1][Mi2+c] /= mm;
                //            gamma[r1][Mi3+c] /= mm;
                //            gamma[r1][Mi4+c] /= mm;
                //        }
                //        delta[r1] /= mm;
                //        gamma[r1][N*M] /= mm;
                //    }
                //}
            }
        }

        const unsigned int N0=N, N1=N-1, N2=N-2, N3=N-3, N4=N-4;
        const PointNodeODE nodeN4(N4*h, static_cast<int>(N4));
        const PointNodeODE nodeN3(N3*h, static_cast<int>(N3));
        const PointNodeODE nodeN2(N2*h, static_cast<int>(N2));
        const PointNodeODE nodeN1(N1*h, static_cast<int>(N1));
        const PointNodeODE nodeN0(N0*h, static_cast<int>(N0));

        DoubleMatrix Mx((k+1)*M, (k+1));
        DoubleVector f((k+1));

        for (unsigned int r=0; r<M; r++)
        {
            for (unsigned int c=0; c<M; c++)
            {
                //Mx[Mx0+r][Mx0+c] = -25.0 - 12.0*h*A(nodeN4, r+1, c+1);;
                //Mx[Mx0+r][Mx1+c] = +48.0;
                //Mx[Mx0+r][Mx2+c] = -36.0;
                //Mx[Mx0+r][Mx3+c] = +16.0;
                //Mx[Mx0+r][Mx4+c] = -3.0;

                //Mx[Mx0+r][Mx0+c] = gamma[r][N0*M+c];
                //Mx[Mx0+r][Mx1+c] = gamma[r][N3*M+c];
                //Mx[Mx0+r][Mx2+c] = gamma[r][N2*M+c];
                //Mx[Mx0+r][Mx3+c] = gamma[r][N1*M+c];
                //Mx[Mx0+r][Mx4+c] = gamma[r][N0*M+c];

                Mx[Mx0+r][Mx0+c] = 0.0;
                Mx[Mx0+r][Mx1+c] = 0.0;
                Mx[Mx0+r][Mx2+c] = 0.0;
                Mx[Mx0+r][Mx3+c] = 0.0;
                Mx[Mx0+r][Mx4+c] = -12.0*h*A(nodeN0, r+1, c+1);

                Mx[Mx1+r][Mx0+c] = 0.0;
                Mx[Mx1+r][Mx1+c] = 0.0;
                Mx[Mx1+r][Mx2+c] = 0.0;
                Mx[Mx1+r][Mx3+c] = -12.0*h*A(nodeN1, r+1, c+1);
                Mx[Mx1+r][Mx4+c] = 0.0;

                Mx[Mx2+r][Mx0+c] = 0.0;
                Mx[Mx2+r][Mx1+c] = 0.0;
                Mx[Mx2+r][Mx2+c] = -12.0*h*A(nodeN2, r+1, c+1);
                Mx[Mx2+r][Mx3+c] = 0.0;
                Mx[Mx2+r][Mx4+c] = 0.0;

                Mx[Mx3+r][Mx0+c] = 0.0;
                Mx[Mx3+r][Mx1+c] = -12.0*h*A(nodeN3, r+1, c+1);
                Mx[Mx3+r][Mx2+c] = 0.0;
                Mx[Mx3+r][Mx3+c] = 0.0;
                Mx[Mx3+r][Mx4+c] = 0.0;

                //Mx[Mx4+r][Mx0+c] = 0.0;
                Mx[Mx4+r][Mx0+c] = gamma[r][N4*M+c];
                Mx[Mx4+r][Mx1+c] = gamma[r][N3*M+c];
                Mx[Mx4+r][Mx2+c] = gamma[r][N2*M+c];
                Mx[Mx4+r][Mx3+c] = gamma[r][N1*M+c];
                Mx[Mx4+r][Mx4+c] = gamma[r][N0*M+c];
            }

            //f[0*M+r] = +12.0*h*B(nodeN4, r+1);
            f[Mx0+r] = +12.0*h*B(nodeN0, r+1);
            f[Mx1+r] = +12.0*h*B(nodeN1, r+1);
            f[Mx2+r] = +12.0*h*B(nodeN2, r+1);
            f[Mx3+r] = +12.0*h*B(nodeN3, r+1);
            f[Mx4+r] = delta[r];

            Mx[Mx0+r][Mx0+r] += +3.0; Mx[Mx0+r][Mx1+r] += -16.0; Mx[Mx0+r][Mx2+r] += +36.0; Mx[Mx0+r][Mx3+r] += -48.0; Mx[Mx0+r][Mx4+r] += +25.0;
            Mx[Mx1+r][Mx0+r] += -1.0; Mx[Mx1+r][Mx1+r] += +6.00; Mx[Mx1+r][Mx2+r] += -18.0; Mx[Mx1+r][Mx3+r] += +10.0; Mx[Mx1+r][Mx4+r] += +3.00;
            Mx[Mx2+r][Mx0+r] += +1.0; Mx[Mx2+r][Mx1+r] += -8.00; Mx[Mx2+r][Mx2+r] += +0.00; Mx[Mx2+r][Mx3+r] += +8.00; Mx[Mx2+r][Mx4+r] += -1.00;
            Mx[Mx3+r][Mx0+r] += -3.0; Mx[Mx3+r][Mx1+r] += -10.0; Mx[Mx3+r][Mx2+r] += +18.0; Mx[Mx3+r][Mx3+r] += -6.00; Mx[Mx3+r][Mx4+r] += +1.00;
        }

        ///////////////////////////////////////////////////////////////////////////
        IPrinter::printSeperatorLine("Mx");
        IPrinter::print(Mx, Mx.rows(), Mx.cols(), 12, 6);
        IPrinter::printSeperatorLine("f");
        IPrinter::print(f, f.length(), 12, 6);

        DoubleVector xf((k+1)*M);
        LinearEquation::GaussianElimination(Mx, f, xf);

        Mx.clear();
        f.clear();

        IPrinter::printSeperatorLine("x");
        IPrinter::print(xf, xf.length(), 12, 6);
        IPrinter::printSeperatorLine();

        for (unsigned int i=0; i<=k; i++)
        {
            printf("***** ");
            for (unsigned int r=1; r<=M; r++)
            {
                printf("%4d %20.16f ", i*M+M-r, xf[i*M+M-r]);
            }
            printf("\n");
        }

        x.clear();
        x.resize(size); for (unsigned int n=0; n<size; n++) x[n].resize(M);

        unsigned int s = xf.length()-M;
        unsigned int e = xf.length()-1;

        IPrinter::printSeperatorLine();
        for (unsigned int n=size-1; n>=end; n--)
        {
            x[n] = xf.mid(s, e);
            s -= M;
            e -= M;
            printf("%6d : ", n); IPrinter::print(x[n], x[n].length(), 20, 16);
        }
        IPrinter::printSeperatorLine();
        xf.clear();

        for (unsigned int j=0, i=N-(k+1); j<=N-(k+1); j++, i--)
        {
            for (unsigned int r=0; r<M; r++)
            {
                unsigned int rw = i*M+r;

                x[i][r] = betta[rw];
                for (unsigned int c=0; c<M; c++)
                {
                    x[i][r] += alpha[rw][Mx1+c]*x[i+1][c] + alpha[rw][Mx2+c]*x[i+2][c]
                            + alpha[rw][Mx3+c]*x[i+3][c] + alpha[rw][Mx4+c]*x[i+4][c]
                            + alpha[rw][Mx5+c]*x[i+5][c] + alpha[rw][Mx6+c]*x[i+6][c];
                }
                printf("%6d : %20.16f\n", i, x[i][r]);
            }
            //printf("%4d : ", i); IPrinter::print(x[i], x[i].length());
        }
        IPrinter::printSeperatorLine("=");
    }

    // begin cleanig memory

    free(delta);
    for (unsigned int r=0; r<M; r++) { free(gamma[r]); }
    free(gamma);
    free(betta);
    for (unsigned int i=0; i<sizeM; i++) { free(alpha[i]); }
    free(alpha);

    // end cleanig memory
}

void IFirstOrderLinearODEIVP::solveInitialValueProblem(DoubleVector &rv) const
{
    if (count() != 1) throw ExceptionODE(5);

    const Dimension &dim = dimension();
    const int min = dim.min();
    const int max = dim.max();
    const double h = dim.step();
    const unsigned int size = static_cast<unsigned int>(max-min);

    rv.resize(size+1);

    double value = initial(InitialCondition::InitialValue);

    rv[0] = value;
    unsigned int j = 1;
    for (int i=min; i<max; i++, j++)
    {
        PointNodeODE node(static_cast<double>(i*h), i);
        double an = A(node)*h + 1.0;
        double bn = B(node)*h;
        rv[j] = an*rv[j-1] + bn;
    }
}

void IFirstOrderLinearODEIVP::solveInitialValueProblem(std::vector<DoubleVector> &rv, ODESolverMethod method) const
{
    switch (method)
    {
    case ODESolverMethod::RUNGE_KUTTA_2: solveInitialValueProblemRK2(rv); break;
    case ODESolverMethod::RUNGE_KUTTA_4: solveInitialValueProblemRK4(rv); break;
    case ODESolverMethod::RUNGE_KUTTA_6: solveInitialValueProblemRK4(rv); break;
    case ODESolverMethod::EULER: solveInitialValueProblemEuler(rv); break;
    case ODESolverMethod::EULER_MOD: solveInitialValueProblemEulerMod(rv); break;
    }
}

void IFirstOrderLinearODEIVP::solveInitialValueProblem(ODESolverMethod method) const
{
    switch (method)
    {
    case ODESolverMethod::RUNGE_KUTTA_2: solveInitialValueProblemRK2(); break;
    case ODESolverMethod::RUNGE_KUTTA_4: solveInitialValueProblemRK4(); break;
    case ODESolverMethod::RUNGE_KUTTA_6: solveInitialValueProblemRK4(); break;
    case ODESolverMethod::EULER: solveInitialValueProblemEuler(); break;
    case ODESolverMethod::EULER_MOD: solveInitialValueProblemEulerMod(); break;
    }
}

void IFirstOrderLinearODEIVP::solveInitialValueProblemRK2(std::vector<DoubleVector> &rv) const
{
    const Dimension &dim = dimension();
    const int min = dim.min();
    const int max = dim.max();
    const double h1 = dim.step();
    const double h2 = 0.5*h1;
    const unsigned int m = count();
    const unsigned int size = static_cast<unsigned int>(max-min);

    rv.resize(size+1); for (unsigned int i=0; i<=size; i++) rv[i].resize(m);

    double *k1 = new double[m];
    double *k2 = new double[m];

    for (unsigned int r=1; r<=m; r++)
    {
        rv[0][r-1] = initial(InitialCondition::InitialValue, r);
    }
    iterationInfo(rv[0], PointNodeODE(min*h1, min));

    unsigned int i=0;
    PointNodeODE node1, node2;
    double *v1 = new double[m];
    for (int n=min; n<max; n++, i++)
    {
        node1.i = n+0; node1.x = node1.i*h1;
        node2.i = n+1; node2.x = node2.i*h1;

        DoubleVector v0 = rv[i];

        // k1
        for (unsigned int row=1, j=0; row<=m; row++, j++)
        {
            double sum = B(node1, row);
            for (unsigned int col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node1, row, col)*v0[j];
            }
            k1[row-1] = sum;
            v1[j]=v0[j]+h1*k1[j];
        }

        // k2
        for (unsigned int row=1; row<=m; row++)
        {
            double sum = B(node2, row);
            for (unsigned int col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node2, row, col)*v1[j];
            }
            k2[row-1] = sum;
        }

        for (unsigned int j=0; j<m; j++)
        {
            rv[i+1][j] = rv[i][j] +  h2 * (k1[j] + k2[j]);
        }
        iterationInfo(rv[i+1], node2);
    }

    delete [] v1;

    delete [] k2;
    delete [] k1;
}

void IFirstOrderLinearODEIVP::solveInitialValueProblemRK4(std::vector<DoubleVector> &rv) const
{
    const Dimension &dim = dimension();
    const int min = dim.min();
    const int max = dim.max();
    const double h1 = dim.step();
    const double h2 = h1/2.0;
    const double h6 = h1/6.0;
    const unsigned int size = static_cast<unsigned int>(max-min);
    const unsigned int m = count();

    rv.resize(size+1); for (unsigned int i=0; i<=size; i++) rv[i].resize(size+1);

    double *k1 = new double[m];
    double *k2 = new double[m];
    double *k3 = new double[m];
    double *k4 = new double[m];

    for (unsigned int r=1; r<=m; r++)
    {
        rv[0][r-1] = initial(InitialCondition::InitialValue, r);
    }
    iterationInfo(rv[0], PointNodeODE(min*h1, min));

    unsigned int i=0;
    PointNodeODE node1, node2, node3, node4;
    double *v1 = new double[m];
    for (int n=min; n<max; n++, i++)
    {
        node1.i = n;   node1.x = node1.i*h1;
        node2.i = n;   node2.x = n*h1+h2;
        node3.i = n;   node3.x = n*h1+h2;
        node4.i = n+1; node4.x = node4.i*h1;

        DoubleVector v0 = rv[i];

        // k1
        for (unsigned int j=0; j<m; v1[j]=v0[j], j++);
        for (unsigned int row=1; row<=m; row++)
        {
            double sum = B(node1, row);
            for (unsigned int col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node1, row, col)*v1[j];
            }
            k1[row-1] = sum;
        }

        // k2
        for (unsigned int j=0; j<m; v1[j]=v0[j]+h2*k1[j], j++);
        for (unsigned int row=1; row<=m; row++)
        {
            double sum = B(node2, row);
            for (unsigned int col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node2, row, col)*v1[j];
            }
            k2[row-1] = sum;
        }

        // k3
        for (unsigned int j=0; j<m; v1[j]=v0[j]+h2*k2[j], j++);
        for (unsigned int row=1; row<=m; row++)
        {
            double sum = B(node3, row);
            for (unsigned int col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node3, row, col)*v1[j];
            }
            k3[row-1] = sum;
        }

        // k4
        for (unsigned int j=0; j<m; v1[j]=v0[j]+h1*k3[j], j++);
        for (unsigned int row=1; row<=m; row++)
        {
            double sum = B(node4, row);
            for (unsigned int col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node4, row, col)*v1[j];
            }
            k4[row-1] = sum;
        }

        for (unsigned int j=0; j<m; j++)
        {
            rv[i+1][j] = rv[i][j] + h6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
        }
        iterationInfo(rv[i+1], node4);
    }

    delete [] v1;

    delete [] k4;
    delete [] k3;
    delete [] k2;
    delete [] k1;
}

void IFirstOrderLinearODEIVP::solveInitialValueProblemRK6(std::vector<DoubleVector> &) const {}

void IFirstOrderLinearODEIVP::solveInitialValueProblemEuler(std::vector<DoubleVector> &rv) const
{
    const Dimension &dim = dimension();
    const int min = dim.min();
    const int max = dim.max();
    const double h = dim.step();
    const size_t size = static_cast<size_t>(max-min);
    const size_t m = count();

    rv.resize(size+1); for (unsigned int i=0; i<=size; i++) rv[i].resize(m);

    DoubleVector &rv0 = rv[0];
    for (unsigned int row=1; row<=m; row++)
    {
        rv0[row-1] = initial(InitialCondition::InitialValue, row);
    }

    unsigned int i=0;
    PointNodeODE node;
    for (int n=min; n<max; n++, i++)
    {
        node.i = n; node.x = n*h;

        DoubleVector v0 = rv[i];

        for (size_t row=1; row<=m; row++)
        {
            double sum = h*B(node, row);
            for (size_t col=1; col<=m; col++)
            {
                sum += h*A(node, row, col)*rv[i][col-1];
            }
            rv[i+1][row-1] = rv[i][row-1] + sum;
        }
    }
}

void IFirstOrderLinearODEIVP::solveInitialValueProblemEulerMod(std::vector<DoubleVector> &) const {}

/**
 * @brief Second-Order Runge-Kutta Methods
 */
void IFirstOrderLinearODEIVP::solveInitialValueProblemRK2() const
{
    const Dimension &dim = dimension();
    const int min = dim.min(), max = dim.max();
    const double h1 = dim.step(), h2 = dim.step()*0.5;
    const size_t m = count(), m2 = count()/2;

    double *k1 = new double[m];
    double *k2 = new double[m];

    double *v0 = new double[m];
    double *v1 = new double[m];

    for (unsigned int r=1; r<=m; r++)
    {
        v0[r-1] = initial(InitialCondition::InitialValue, r);
    }
    iterationInfo(DoubleVector(v0, m2), PointNodeODE(min*h1, min));

    PointNodeODE node1, node2;
    for (int n=min; n<max; n++)
    {
        node1.i = n+0; node1.x = node1.i*h1;
        node2.i = n+1; node2.x = node2.i*h1;

        // k1
        for (size_t row=1, j=0; row<=m; row++, j++)
        {
            double sum = B(node1, row);
            for (size_t col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node1, row, col)*v0[j];
            }
            k1[row-1] = sum;
            v1[j]=v0[j]+h1*k1[j];
        }

        // k2
        for (size_t row=1; row<=m; row++)
        {
            double sum = B(node2, row);
            for (size_t col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node2, row, col)*v1[j];
            }
            k2[row-1] = sum;
        }

        for (size_t j=0; j<m; j++)
        {
            v0[j] = v0[j] + h2 * (k1[j] + k2[j]);
        }
        iterationInfo(DoubleVector(v0, m2), node2);
    }

    delete [] v0;
    delete [] v1;

    delete [] k2;
    delete [] k1;
}

/**
 * @brief Fourth-Order Runge-Kutta Methods
 */
void IFirstOrderLinearODEIVP::solveInitialValueProblemRK4() const
{
    const Dimension &dim = dimension();
    const int min = dim.min(), max = dim.max();
    const double h1 = dim.step(), h2 = dim.step()*0.5, h6 = dim.step()/6.0;
    const size_t m = count(), m2 = count()/2;

    double *k1 = new double[m];
    double *k2 = new double[m];
    double *k3 = new double[m];
    double *k4 = new double[m];

    double *v0 = new double[m];
    double *v1 = new double[m];

    for (unsigned int r=1; r<=m; r++)
    {
        v0[r-1] = initial(InitialCondition::InitialValue, r);
    }
    iterationInfo(DoubleVector(v0, m2), PointNodeODE(min*h1, min));

    unsigned int i = 0;
    PointNodeODE node1, node2, node3, node4;
    for (int n=min; n<max; n++, i++)
    {
        node1.i = n;   node1.x = node1.i*h1;
        node2.i = n;   node2.x = node2.i*h1+h2;
        node3.i = n;   node3.x = node3.i*h1+h2;
        node4.i = n+1; node4.x = node4.i*h1;

        memcpy(v1, v0, sizeof (double) * m);

        // k1
        for (size_t row=1; row<=m; row++)
        {
            double sum = B(node1, row);
            for (size_t col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node1, row, col)*v1[j];
            }
            k1[row-1] = sum;
        }

        // k2
        for (size_t j=0; j<m; v1[j]=v0[j]+h2*k1[j], j++);
        for (size_t row=1; row<=m; row++)
        {
            double sum = B(node2, row);
            for (size_t col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node2, row, col)*v1[j];
            }
            k2[row-1] = sum;
        }

        // k3
        for (size_t j=0; j<m; v1[j]=v0[j]+h2*k2[j], j++);
        for (size_t row=1; row<=m; row++)
        {
            double sum = B(node3, row);
            for (size_t col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node3, row, col)*v1[j];
            }
            k3[row-1] = sum;
        }

        // k4
        for (size_t j=0; j<m; v1[j]=v0[j]+h1*k3[j], j++);
        for (size_t row=1; row<=m; row++)
        {
            double sum = B(node4, row);
            for (size_t col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node4, row, col)*v1[j];
            }
            k4[row-1] = sum;
        }

        for (size_t j=0; j<m; j++)
        {
            v0[j] = v0[j] + h6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
        }
        iterationInfo(DoubleVector(v0, m2), node4);
    }

    delete [] v0;
    delete [] v1;

    delete [] k4;
    delete [] k3;
    delete [] k2;
    delete [] k1;
}

/**
 * @brief Sixth-Order Runge-Kutta Methods
 */
void IFirstOrderLinearODEIVP::solveInitialValueProblemRK6() const {}

/**
 * @brief IFirstOrderLinearODEIVP::solveInitialValueProblemEuler
 */
void IFirstOrderLinearODEIVP::solveInitialValueProblemEuler() const
{
    const Dimension &dim = dimension();
    const double h = dim.step();
    const int min = dim.min(), max = dim.max();
    const size_t m = count(), m2 = count()/2;

    double *v0 = new double[m];
    double *v1 = new double[m];

    for (size_t row=1; row<=m; row++)
    {
        v0[row-1] = initial(InitialCondition::InitialValue, row);
    }
    iterationInfo(DoubleVector(v0, m2), PointNodeODE(min*h, min));

    size_t i = 0;
    PointNodeODE node;
    for (int n=min; n<max; n++, i++)
    {
        node.i = n; node.x = n*h;

        for (size_t row=1; row<=m; row++)
        {
            double sum = h*B(node, row);
            for (size_t col=1; col<=m; col++)
            {
                sum += h*A(node, row, col)*v0[col-1];
            }
            v1[row-1] = v0[row-1] + sum;
        }

        memcpy(v0, v1, sizeof (double)*m);

        iterationInfo(DoubleVector(v0, m2), PointNodeODE((n+1)*h, (n+1)));
    }

    delete [] v1;
    delete [] v0;
}

void IFirstOrderLinearODEIVP::solveInitialValueProblemEulerMod() const {}

void IFirstOrderLinearODEIVP::start(DoubleVector &x, PointNodeODE &n) const
{
    const Dimension &dim = dimension();
    const size_t m = count();
    const int min = dim.min();
    const double h = dim.step();
    n.i = min; n.x = n.i*h;

    if (x.length() != m) x.resize(m);

    for (size_t row=1; row<=m; row++)
    {
        x[row-1] = initial(InitialCondition::InitialValue, row);
    }
    iterationInfo(x, n);
}

void IFirstOrderLinearODEIVP::next(const DoubleVector &x0, const PointNodeODE &n0, DoubleVector &x1, PointNodeODE &n1, ODESolverMethod method) const
{
    const Dimension &dim = dimension();
    const double h = dim.step();
    const size_t m = count();
    n1.i = n0.i + 1; n1.x = n1.i*h;

    if (x1.length() != m) x1.resize(m);

    if (method == ODESolverMethod::EULER)
    {
        for (size_t row=1; row<=m; row++)
        {
            double sum = 0.0;
            for (size_t col=1; col<=m; col++)
            {
                sum += h*A(n0, row, col)*x0[col-1];
            }
            x1[row-1] = x0[row-1] + sum + h*B(n0, row);
        }
    }
/*
    if (method == ODESolverMethod::EULER_MOD)
    {
    }

    if (method == ODESolverMethod::RUNGE_KUTTA_2)
    {
        const double h2 = h/2.0;
        const double h6 = h/6.0;

        double *k1 = new double[m];
        double *k2 = new double[m];

        double *v0 = new double[m];
        double *v1 = new double[m];

        PointNodeODE node1, node2;
        node1.i = n0.i;   node1.x = node1.i*h;
        node2.i = n0.i+1; node2.x = node2.i*h;

        memcpy(v0, x0.data(), sizeof (double) * m);

        // k1
        for (size_t row=1, j=0; row<=m; row++, j++)
        {
            double sum = B(node1, row);
            for (size_t col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node1, row, col)*v0[j];
            }
            k1[row-1] = sum;
            v1[j]=v0[j]+h*k1[j];
        }

        // k2
        for (size_t row=1; row<=m; row++)
        {
            double sum = B(node2, row);
            for (size_t col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node2, row, col)*v1[j];
            }
            k2[row-1] = sum;
        }

        for (size_t j=0; j<m; j++)
        {
            x1[j] = x0[j] + h2 * (k1[j] + k2[j]);
        }

        delete [] v0;
        delete [] v1;

        delete [] k2;
        delete [] k1;
    }

    if (method == ODESolverMethod::RUNGE_KUTTA_4)
    {
        const double h2 = h/2.0;
        const double h6 = h/6.0;

        double *k1 = new double[m];
        double *k2 = new double[m];
        double *k3 = new double[m];
        double *k4 = new double[m];

        double *v0 = new double[m];
        double *v1 = new double[m];

        PointNodeODE node1, node2, node3, node4;
        node1.i = n0.i;   node1.x = node1.i*h;
        node2.i = n0.i;   node2.x = node2.i*h+h2;
        node3.i = n0.i;   node3.x = node3.i*h+h2;
        node4.i = n0.i+1; node4.x = node4.i*h;

        memcpy(v0, x0.data(), sizeof (double) * m);

        // k1
        for (size_t row=1; row<=m; row++)
        {
            double sum = B(node1, row);
            for (size_t col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node1, row, col)*v0[j];
            }
            k1[row-1] = sum;
        }

        // k2
        for (size_t j=0; j<m; v1[j]=v0[j]+h2*k1[j], j++);
        for (size_t row=1; row<=m; row++)
        {
            double sum = B(node2, row);
            for (size_t col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node2, row, col)*v1[j];
            }
            k2[row-1] = sum;
        }

        // k3
        for (size_t j=0; j<m; v1[j]=v0[j]+h2*k2[j], j++);
        for (size_t row=1; row<=m; row++)
        {
            double sum = B(node3, row);
            for (size_t col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node3, row, col)*v1[j];
            }
            k3[row-1] = sum;
        }

        // k4
        for (size_t j=0; j<m; v1[j]=v0[j]+h*k3[j], j++);
        for (size_t row=1; row<=m; row++)
        {
            double sum = B(node4, row);
            for (size_t col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node4, row, col)*v1[j];
            }
            k4[row-1] = sum;
        }

        for (size_t j=0; j<m; j++)
        {
            x1[j] = x0[j] + h6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
        }

        delete [] v0;
        delete [] v1;

        delete [] k4;
        delete [] k3;
        delete [] k2;
        delete [] k1;
    }

    if (method == ODESolverMethod::RUNGE_KUTTA_6)
    {
    }
*/
    iterationInfo(x1, n1);
}

//void LinearODE1stOrder::calculate(const std::vector<Condition> &nscs, const DoubleVector &bt, std::vector<DoubleVector> &x)
//{
//    double h = dimension().step();

//    std::vector<Condition> cs = nscs;
//    DoubleVector beta = bt;

//    unsigned int L = cs.size();
//    unsigned int n = 0;
//    if (L!=0) n = cs.at(0).mtrx.rows();

//    DoubleVector x0(n+2);
//    DoubleVector rx(n+2);
//    for (unsigned int row=0; row<n; row++)
//    {
//        for (unsigned int s=0; s<L-1; s++)
//        {
//            Condition &sc = cs.at(s);
//            Condition &ec = cs.at(s+1);

//            for (unsigned int i=0; i<n; i++) x0[i] = sc.mtrx[row][i];
//            x0[n] = beta[row];
//            x0[n+1] = 1.0;

//            struct HelperB : public NonLinearODE1stOrder
//            {
//                LinearODE1stOrder *p;
//                unsigned int n;

//                virtual double f(double, double, unsigned int) const { return NAN; }

//                virtual double f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const
//                {
//                    //unsigned int n = p->systemOrder();
//                    double _SO = S0(t,x,k);

//                    PointNodeODE node(t, k);

//                    if (i<n)
//                    {
//                        double res = _SO*x[i];
//                        for (unsigned int j=0; j<n; j++) res -= p->A(node,j,i)*x[j];
//                        return res;
//                    }
//                    else if (i==n)
//                    {
//                        double res = _SO*x[n];
//                        for (unsigned int j=0; j<n; j++) res += p->B(node,j)*x[j];
//                        return res;
//                    }
//                    else
//                    {
//                        return _SO*x[n+1];
//                    }

//                    return NAN;
//                }

//                double S0(double t, const DoubleVector &x, unsigned int k) const
//                {
//                    //unsigned int n = p.systemOrder();
//                    double btr = x[n];

//                    PointNodeODE node(t, k);

//                    double s1 = 0.0;
//                    for (unsigned int i=0; i<n; i++)
//                    {
//                        double aa = 0.0;
//                        for (unsigned int j=0; j<n; j++) aa += x[j]*p->A(node,j,i);
//                        s1 += aa*x[i];
//                    }

//                    double s2 = 0.0;
//                    for (unsigned int i=0; i<n; i++)
//                    {
//                        s2 += x[i]*p->B(node,i);
//                    }
//                    s2 *= btr;


//                    double m1 = 0.0;
//                    for (unsigned int i=0; i<n; i++)
//                    {
//                        m1 += x[i]*x[i];
//                    }
//                    m1 += btr*btr;

//                    return (s1-s2)/m1;
//                }
//            };

//            HelperB helper;
//            helper.p = this;
//            helper.n = n;
//            helper.setDimension(Dimension(h, sc.nmbr, ec.nmbr));
//            helper.cauchyProblem(sc.time, x0, rx, NonLinearODE1stOrder::RK4);

//            for (unsigned int i=0; i<n; i++) sc.mtrx[row][i] = rx[i];
//            beta[row] = rx[n];
//            double M = rx[n+1];

//            for (unsigned int j=s+1; j<L; j++)
//            {
//                Condition &cc = cs.at(j);
//                for (unsigned int i=0; i<n; i++) cc.mtrx[row][i] *= M;
//            }

//            for (unsigned int i=0; i<n; i++) ec.mtrx[row][i] += rx[i];
//        }
//    }
//    x.clear();
//    rx.clear();

//    // Separated conditions in right side

//    DoubleMatrix A(n, n);
//    DoubleVector b(n);
//    DoubleVector x1(n);

//    Condition c0 = cs.at(L-1);
//    for (unsigned int row=0; row<n; row++)
//    {
//        for (unsigned int col=0; col<n; col++)
//        {
//            A[row][col] = c0.mtrx[row][col];
//        }
//        b[row] = beta[row];
//    }

//    LinearEquation::GaussianElimination(A, b, x1);

//    struct HelperB : public NonLinearODE1stOrder
//    {
//        LinearODE1stOrder *p;
//        unsigned int n;
//    protected:
//        virtual double f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const
//        {
//            PointNodeODE node(t,k);

//            double res = p->B(node,i);
//            for (unsigned int j=0; j<n; j++) res += p->A(node,i,j)*x[j];
//            return res;
//        }
//    };

//    HelperB helper;
//    helper.p = this;
//    helper.n = n;
//    helper.setDimension(dimension());
//    helper.cauchyProblem(nscs.back().time, x1, x, HelperB::RK4, HelperB::R2L);
//}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//void discretizationL2(const std::vector<LinearODE1stOrder::Condition> &cs, double *b, double h, unsigned int N)
//{
//    const unsigned int cnd_size = static_cast<unsigned int>( cs.size() );

//    for (unsigned int s=0; s<cnd_size; s++)
//    {
//        const LinearODE1stOrder::Condition &cnd = cs[s];

//        double alpha = cnd.mtrx.at(0,0);
//        double time  = cnd.time;

//        for (unsigned int n=0; n<=N; n++)
//        {
//            double dh = fabs(time - n*h);
//            if (dh <= h) b[n] += (1.0 - dh/h)*alpha;
//        }
//    }
//}

//void discretizationL3(const std::vector<LinearODE1stOrder::Condition> &cs, double *b, double h, unsigned int N)
//{
//    const unsigned int cnd_size = static_cast<unsigned int>( cs.size() );

//    for (unsigned int s=0; s<cnd_size; s++)
//    {
//        const LinearODE1stOrder::Condition &cnd = cs[s];

//        double alpha = cnd.mtrx.at(0,0);
//        double time  = cnd.time;

//        double h2 = h*h;
//        double h21 = (1.0/h2) * alpha;
//        double h22 = (1.0/(2.0*h2)) * alpha;
//        for (unsigned int n=0; n<=N; n++)
//        {
//            double curt = n*h;
//            double dh = fabs(time - curt);

//            if (0.0 < time && time < h/* && n < 4*/)
//            {
//                if (n==0) b[n] += ((h-dh)*(2.0*h-dh)) * h22;
//                if (n==1) b[n] += ((h-dh)*(h+dh)) * h21;
//                if (n==2) b[n] += ((2.0*h-dh)*(h-dh)) * h22;
//            }
//            else if ((N-1)*h < time && time < N*h/* && n > N-4*/)
//            {
//                if (n==N-2) b[n] += ((h-dh)*(2.0*h-dh)) * h22;
//                if (n==N-1) b[n] += ((h-dh)*(h+dh)) * h21;
//                if (n==N-0) b[n] += ((2.0*h-dh)*(h-dh)) * h22;
//            }
//            else
//            {
//                if (dh <= h && n*h <= time)               b[n] += ((2.0*h-dh)*(h-dh)) * h22;
//                if (dh <= h && n*h >= time)               b[n] += ((h-dh)*(h+dh)) * h21;
//                if (dh > h && dh <= 2.0*h && n*h >= time) b[n] += ((2.0*h-dh)*(h-dh)) * h22;
//            }
//        }
//    }
//}

//void discretizationL4(const std::vector<LinearODE1stOrder::Condition> &cs, double *b, double h, unsigned int N)
//{
//    const unsigned int cnd_size = static_cast<unsigned int>( cs.size() );

//    for (unsigned int s=0; s<cnd_size; s++)
//    {
//        const LinearODE1stOrder::Condition &cnd = cs[s];

//        double alpha = cnd.mtrx.at(0,0);
//        double time  = cnd.time;

//        double h3 = h*h*h;
//        double h32 = (1.0/(2.0*h3)) * alpha;
//        double h36 = (1.0/(6.0*h3)) * alpha;
//        for (unsigned int n=0; n<=N; n++)
//        {
//            double curt = n*h;
//            double dh = fabs(time - curt);

//            if (0.0 < time && time < h/* && n < 4*/)
//            {
//                if (n==0) b[n] += ((2.0*h-dh)*(h-dh)*(3.0*h-dh)) * h36;
//                if (n==1) b[n] += ((2.0*h+dh)*(h-dh)*(h+dh)) * h32;
//                if (n==2) b[n] += ((2.0*h-dh)*(h-dh)*(h+dh)) * h32;
//                if (n==3) b[n] += ((2.0*h-dh)*(h-dh)*(3.0*h-dh)) * h36;
//            }
//            else if ((N-1)*h < time && time < N*h/* && n > N-4*/)
//            {
//                if (n==N-3) b[n] += ((2.0*h-dh)*(h-dh)*(3.0*h-dh)) * h36;
//                if (n==N-2) b[n] += ((2.0*h-dh)*(h-dh)*(h+dh)) * h32;
//                if (n==N-1) b[n] += ((2.0*h+dh)*(h-dh)*(h+dh)) * h32;
//                if (n==N-0) b[n] += ((2.0*h-dh)*(h-dh)*(3.0*h-dh)) * h36;
//            }
//            else
//            {
//                if (dh <= h)               b[n] += ((2.0*h-dh)*(h-dh)*(h+dh)) * h32;
//                if (dh > h && dh <= 2.0*h) b[n] += ((2.0*h-dh)*(h-dh)*(3.0*h-dh)) * h36;
//            }
//        }
//    }
//}

//void discretization(const std::vector<LinearODE1stOrder::Condition> &cs, double *b, double h, unsigned int N)
//{
//    discretizationL4(cs, b, h, N);
//}

//void discretizationL2(const std::vector<LinearODE1stOrder::Condition> &cs, DoubleMatrix* b, double h, unsigned int N)
//{
//    const unsigned int cnd_size = static_cast<unsigned int>( cs.size() );

//    for (unsigned int s=0; s<cnd_size; s++)
//    {
//        const LinearODE1stOrder::Condition &cnd = cs[s];

//        const DoubleMatrix &alpha = cnd.mtrx;
//        double time  = cnd.time;

//        for (unsigned int n=0; n<=N; n++)
//        {
//            double dh = fabs(time - n*h);
//            if (dh <= h) b[n] += (1.0 - dh/h)*alpha;
//        }
//    }
//}

//void discretizationL4(const std::vector<LinearODE1stOrder::Condition> &cs, DoubleMatrix* b, double h, unsigned int N)
//{
//    const unsigned int cnd_size = static_cast<unsigned int>( cs.size() );

//    for (unsigned int s=0; s<cnd_size; s++)
//    {
//        const LinearODE1stOrder::Condition &cnd = cs[s];

//        const DoubleMatrix &alpha = cnd.mtrx;
//        double time  = cnd.time;

//        double h3 = h*h*h;
//        DoubleMatrix h32 = (1.0/(2.0*h3)) * alpha;
//        DoubleMatrix h36 = (1.0/(6.0*h3)) * alpha;
//        for (unsigned int n=0; n<=N; n++)
//        {
//            double curt = n*h;
//            double dh = fabs(time - curt);

//            if (0.0 < time && time < h && n < 4)
//            {
//                if (n==0) b[n] += ((2.0*h-dh)*(h-dh)*(3.0*h-dh)) * h36;
//                if (n==1) b[n] += ((2.0*h+dh)*(h-dh)*(h+dh)) * h32;
//                if (n==2) b[n] += ((2.0*h-dh)*(h-dh)*(h+dh)) * h32;
//                if (n==3) b[n] += ((2.0*h-dh)*(h-dh)*(3.0*h-dh)) * h36;
//            }
//            else if ((N-1)*h < time && time < N*h && n > N-4)
//            {
//                if (n==N-3) b[n] += ((2.0*h-dh)*(h-dh)*(3.0*h-dh)) * h36;
//                if (n==N-2) b[n] += ((2.0*h-dh)*(h-dh)*(h+dh)) * h32;
//                if (n==N-1) b[n] += ((2.0*h+dh)*(h-dh)*(h+dh)) * h32;
//                if (n==N-0) b[n] += ((2.0*h-dh)*(h-dh)*(3.0*h-dh)) * h36;
//            }
//            else
//            {
//                if (dh <= h)               b[n] += ((2.0*h-dh)*(h-dh)*(h+dh)) * h32;
//                if (dh > h && dh <= 2.0*h) b[n] += ((2.0*h-dh)*(h-dh)*(3.0*h-dh)) * h36;
//            }
//        }
//    }
//}

//void discretization(const std::vector<LinearODE1stOrder::Condition> &cs, DoubleMatrix* b, double h, unsigned int N)
//{
//    discretizationL4(cs, b, h, N);
//}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//void LinearODE1stOrder::solveHighOderAccuracy(const std::vector<Condition> &cs, const DoubleVector &rs, std::vector<DoubleVector> &x, unsigned int k, Direction direction UNUSED_PARAM)
//{
//    switch (k)
//    {
//    case 2: highOder2Accuracy(cs, rs, x, direction); break;
//    case 4: highOder4Accuracy(cs, rs, x, direction); break;
//    case 6: highOder6Accuracy(cs, rs, x, direction); break;
//    default: break;
//    }
//}

//void LinearODE1stOrder::highOder2Accuracy(const std::vector<Condition> &cs, const DoubleVector & rs, std::vector<DoubleVector> &x, Direction direction UNUSED_PARAM)
//{
//    unsigned int en = equationsNumber();

//    double h = dimension().step();
//    int N = dimension().size();

//    if (en == 1)
//    {
//        double *p = (double*) malloc(sizeof(double)*N);
//        double *q = (double*) malloc(sizeof(double)*N);
//        double *r = (double*) malloc(sizeof(double)*N);

//        std::vector<unsigned int> ind;
//        std::vector<double>       ems;

//        /**********************************************************************
//         *                  Discretization using trapesium rule
//         *********************************************************************/
//        double *b = (double*) malloc(sizeof(double)*(N+1));
//        for (unsigned int n=0; n<=N; n++) b[n] = 0.0;

//        /********************* discretisation ********************************/
//        discretization(cs, b, h, N);
//        /*********************************************************************/

//        p[0] = b[0];
//        q[0] = b[1];
//        for (unsigned int n=2; n<=N; n++)
//        {
//            if (fabs(b[n]) > DBL_EPSILON)
//            {
//                ind.push_back(n);
//                ems.push_back(b[n]);
//            }
//        }
//        free(b);
//        r[0] = rs[0];

//        unsigned int ind_size = ind.size();
//        unsigned int ems_size = ems.size();

//        /**********************************************************************
//         *                          End of discretization
//         *********************************************************************/

//        /**********************************************************************
//         *                      Finding function at end of grid
//         *********************************************************************/

//        for (unsigned int n=0; n<=N-2; n++)
//        {
//            PointNodeODE node(n*h, n);

//            double m = 1.0/(2.0*h*A(node)+3.0);
//            double alpha1 = +4.0*m;
//            double alpha2 = -1.0*m;
//            double alpha0 = -2.0*h*B(node)*m;

//            r[n+1] = r[n] - p[n]*alpha0;
//            q[n+1] = p[n]*alpha2;
//            p[n+1] = p[n]*alpha1 + q[n];


//            for (unsigned int i=0; i<ind_size; i++)
//            {
//                if (n+2 == ind[i]) q[n+1] += ems[i];
//            }
//        }

//        DoubleMatrix m(3,3);
//        DoubleVector c(3);
//        DoubleVector xT(3);

//        m[0][0] = 0.0;
//        m[0][1] = p[N-1];
//        m[0][2] = q[N-1];
//        c[0] = r[N-1];

//        PointNodeODE nodeN1((N-1)*h, N-1);
//        m[1][0] = -1.0;
//        m[1][1] = -2.0*h*A(nodeN1);
//        m[1][2] = +1.0;
//        c[1] = 2.0*h*B(nodeN1);

//        PointNodeODE nodeN(N*h, N);
//        m[2][0] = +1.0;
//        m[2][1] = -4.0;
//        m[2][2] = 3.0-2.0*h*A(nodeN);
//        c[2] = 2.0*h*B(nodeN);

//        LinearEquation::GaussianElimination(m, c, xT);
//        m.clear();
//        c.clear();

//        /**********************************************************************
//         *                      Finding function at end of grid
//         *********************************************************************/

//        DoubleVector x0(N+1);
//        x0[N-0] = xT[2];
//        x0[N-1] = xT[1];
//        x0[N-2] = xT[0];
//        xT.clear();


//        for (unsigned int i=N-2; i!=0; i--)
//        {
//            x0[i-1] = -q[i-1]*x0[i]+r[i-1];
//            for (unsigned int s=0; s<ems_size; s++)
//                if (ind[s] > i) x0[i-1] -= ems[s]*x0[ind[s]];
//            x0[i-1] /= p[i-1];
//        }

//        x.clear();
//        x.push_back(x0);

//        free(p);
//        free(q);
//        free(r);
//        ems.clear();
//        ind.clear();
//    }
//    else
//    {
//        DoubleMatrix* p = new DoubleMatrix[N]; for (unsigned int i=0; i<N; i++) p[i].resize(en, en, 0.0);
//        DoubleMatrix* q = new DoubleMatrix[N]; for (unsigned int i=0; i<N; i++) q[i].resize(en, en, 0.0);
//        DoubleMatrix* r = new DoubleMatrix[N]; for (unsigned int i=0; i<N; i++) r[i].resize(en, 1, 0.0);

//        std::vector<unsigned int> ind;
//        std::vector<DoubleMatrix> ems;

//        /**********************************************************************
//         *                  Discretization using trapesium rule
//         *********************************************************************/
//        DoubleMatrix* b = new DoubleMatrix[N+1];
//        for (unsigned int n=0; n<=N; n++) b[n].resize(en, en, 0.0);

//        /********************* discretisation ********************************/
//        discretization(cs, b, h, N);
//        /*********************************************************************/

//        p[0] = b[0];
//        q[0] = b[1];
//        for (unsigned int n=2; n<=N; n++)
//        {
//            if (!b[n].zeroMatrix())
//            {
//                ind.push_back(n);
//                ems.push_back(b[n]);
//            }
//        }
//        free(b);
//        r[0] = rs;

//        unsigned int ind_size = ind.size();
//        unsigned int ems_size = ems.size();

//        /**********************************************************************
//         *                          End of discretization
//         *********************************************************************/

//        /**********************************************************************
//         *                      Finding function at end of grid
//         *********************************************************************/

//        DoubleMatrix m(en, en);
//        DoubleMatrix alpha1(en, en, 0.0);
//        DoubleMatrix alpha2(en, en, 0.0);
//        DoubleMatrix alpha0(en, 1, 0.0);

//        for (unsigned int n=0; n<=N-2; n++)
//        {
//            PointNodeODE node(n*h, n);

//            for (unsigned int row=0; row<en; row++)
//            {
//                for (unsigned int col=0; col<en; col++)
//                {
//                    m[row][col] = 2.0*h*A(node,row,col);
//                    if (row==col) m[row][col] += 3.0;

//                    alpha1[row][col] = 0.0;
//                    if (row==col) alpha1[row][col] = 4.0;

//                    alpha2[row][col] = 0.0;
//                    if (row==col) alpha2[row][col] = -1.0;
//                }
//                alpha0[row][0] = -2.0*h*B(node,row);
//            }

//            m.inverse();

//            alpha1 = m*alpha1;
//            alpha2 = m*alpha2;
//            alpha0 = m*alpha0;

//            r[n+1] = r[n] - p[n]*alpha0;
//            q[n+1] = p[n]*alpha2;
//            p[n+1] = p[n]*alpha1 + q[n];

//            for (unsigned int i=0; i<ind_size; i++)
//            {
//                if (n+2 == ind[i]) q[n+1] += ems[i];
//            }
//        }

//        m.clear();
//        alpha0.clear();
//        alpha1.clear();
//        alpha2.clear();

//        DoubleMatrix M(3*en,3*en);
//        DoubleVector C(3*en);
//        DoubleVector xT(3*en);

//        for (unsigned int row=0*en; row<1*en; row++)
//        {
//            for (unsigned int col=0*en; col<1*en; col++) M[row][col] = 0.0;
//            for (unsigned int col=1*en; col<2*en; col++) M[row][col] = p[N-1][row%en][col%en];
//            for (unsigned int col=2*en; col<3*en; col++) M[row][col] = q[N-1][row%en][col%en];
//            C[row] = r[N-1][row%en][0];
//        }

//        PointNodeODE nodeN1((N-1)*h, N-1);
//        for (unsigned int row=1*en; row<2*en; row++)
//        {
//            for (unsigned int col=0*en; col<1*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -1.0; }
//            for (unsigned int col=1*en; col<2*en; col++) { M[row][col] = -2.0*h*A(nodeN1, row%en, col%en); }
//            for (unsigned int col=2*en; col<3*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +1.0; }
//            C[row] = 2.0*h*B(nodeN1, row%en);
//        }

//        PointNodeODE nodeN(N*h, N);
//        for (unsigned int row=2*en; row<3*en; row++)
//        {
//            for (unsigned int col=0*en; col<1*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +1.0; }
//            for (unsigned int col=1*en; col<2*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -4.0; }
//            for (unsigned int col=2*en; col<3*en; col++) { M[row][col] = -2.0*h*A(nodeN, row%en, col%en); if (row%en == col%en) M[row][col] += +3.0; }
//            C[row] = 2.0*h*B(nodeN, row%en);
//        }

//        LinearEquation::GaussianElimination(M, C, xT);
//        M.clear();
//        C.clear();

//        std::vector<DoubleMatrix> x0(N+1);
//        for (unsigned int i=0; i<=N; i++) x0[i].resize(en,1);

//        for (unsigned int row=0; row<en; row++)
//        {
//            x0[N-2][row][0] = xT[row+0*en];
//            x0[N-1][row][0] = xT[row+1*en];
//            x0[N-0][row][0] = xT[row+2*en];
//        }
//        xT.clear();

//        for (unsigned int i=N-2; i!=0; i--)
//        {
//            x0[i-1] = r[i-1]-q[i-1]*x0[i];
//            for (unsigned int s=0; s<ems_size; s++)
//                if (ind[s] > i) x0[i-1] -= (ems[s]*x0[ind[s]]);
//            p[i-1].inverse();
//            x0[i-1] = p[i-1] * x0[i-1];
//        }

//        x.resize(en);
//        for (unsigned int i=0; i<en; i++) x[i].resize(N+1);

//        for (unsigned int i=0; i<=N; i++)
//        {
//            for (unsigned int row=0; row<en; row++)
//            {
//                x[row][i] = x0[i][row][0];
//            }
//        }

//        for (unsigned int i=0; i<=N; i++)
//        {
//            x0[i].clear();
//        }
//        x0.clear();

//        ind.clear();
//        ems.clear();

//        for (unsigned int i=0; i<N; i++)
//        {
//            r[i].clear();
//            q[i].clear();
//            p[i].clear();
//        }
//        delete [] r;
//        delete [] q;
//        delete [] p;
//    }
//}

//void LinearODE1stOrder::highOder4Accuracy(const std::vector<Condition> &cs, const DoubleVector& rs, std::vector<DoubleVector> &x, Direction direction UNUSED_PARAM)
//{
//    unsigned int en = equationsNumber();

//    double h = dimension().step();
//    unsigned int N = dimension().size();

//    if (en == 1)
//    {
//        double *p = (double*) malloc(sizeof(double)*(N-2));
//        double *q = (double*) malloc(sizeof(double)*(N-2));
//        double *v = (double*) malloc(sizeof(double)*(N-2));
//        double *u = (double*) malloc(sizeof(double)*(N-2));
//        double *r = (double*) malloc(sizeof(double)*(N-2));

//        std::vector<unsigned int> ind;
//        std::vector<double>       ems;

//        /**********************************************************************
//         *                  Discretization using trapesium rule
//         *********************************************************************/
//        double *b = (double*) malloc(sizeof(double)*(N+1));
//        for (unsigned int m=0; m<=N; m++) b[m] = 0.0;

//        /********************* discretisation ********************************/
//        discretization(cs, b, h, N);
//        /*********************************************************************/

//        p[0] = b[0];
//        q[0] = b[1];
//        v[0] = b[2];
//        u[0] = b[3];
//        for (unsigned int n=4; n<=N; n++)
//        {
//            if (fabs(b[n]) > DBL_EPSILON)
//            {
//                ind.push_back(n);
//                ems.push_back(b[n]);
//            }
//        }
//        free(b);
//        r[0] = rs[0];

//        unsigned int ind_size = ind.size();
//        unsigned int ems_size = ems.size();

//        /**********************************************************************
//         *                          End of discretization
//         *********************************************************************/

//        /**********************************************************************
//         *                      Finding function at end of grid
//         *********************************************************************/

//        for (unsigned int n=0; n<=N-4; n++)
//        {
//            PointNodeODE node(n*h, n);

//            double m = +1.0/(-12.0*h*A(node)-25.0);
//            double alpha1 = -48.0*m;
//            double alpha2 = +36.0*m;
//            double alpha3 = -16.0*m;
//            double alpha4 = +3.0*m;
//            double alpha0 = +12.0*h*B(node)*m;

//            r[n+1] = r[n] - p[n]*alpha0;
//            u[n+1] = p[n]*alpha4;
//            v[n+1] = p[n]*alpha3 + u[n];
//            q[n+1] = p[n]*alpha2 + v[n];
//            p[n+1] = p[n]*alpha1 + q[n];

//            for (unsigned int i=0; i<ind_size; i++)
//            {
//                if (n+4 == ind[i]) u[n+1] += ems[i];
//            }
//        }

//        DoubleMatrix m(5,5);
//        DoubleVector c(5);
//        DoubleVector xT(5);

//        m[0][0] = 0.0;
//        m[0][1] = p[N-3];
//        m[0][2] = q[N-3];
//        m[0][3] = v[N-3];
//        m[0][4] = u[N-3];
//        c[0] = r[N-3];

//        PointNodeODE nodeN3((N-3)*h, N-3);
//        m[1][0] = -3.0;
//        m[1][1] = -10.0 - 12.0*h*A(nodeN3);
//        m[1][2] = +18.0;
//        m[1][3] = -6.0;
//        m[1][4] = +1.0;
//        c[1] = 12.0*h*B(nodeN3);

//        PointNodeODE nodeN2((N-2)*h, N-2);
//        m[2][0] = +1.0;
//        m[2][1] = -8.0;
//        m[2][2] = +0.0 -12.0*h*A(nodeN2);
//        m[2][3] = +8.0;
//        m[2][4] = -1.0;
//        c[2] = 12.0*h*B(nodeN2);

//        PointNodeODE nodeN1((N-1)*h, N-1);
//        m[3][0] = -1.0;
//        m[3][1] = +6.0;
//        m[3][2] = -18.0;
//        m[3][3] = +10.0 - 12.0*h*A(nodeN1);
//        m[3][4] = +3.0;
//        c[3] = 12.0*h*B(nodeN1);

//        PointNodeODE nodeN0(N*h, N);
//        m[4][0] = +3.0;
//        m[4][1] = -16.0;
//        m[4][2] = +36.0;
//        m[4][3] = -48.0;
//        m[4][4] = +25.0 - 12.0*h*A(nodeN0);
//        c[4] = 12.0*h*B(nodeN0);

//        LinearEquation::GaussianElimination(m, c, xT);
//        m.clear();
//        c.clear();

//        /**********************************************************************
//         *                      Finding function at end of grid
//         *********************************************************************/

//        DoubleVector x0(N+1);
//        x0[N-0] = xT[4];
//        x0[N-1] = xT[3];
//        x0[N-2] = xT[2];
//        x0[N-3] = xT[1];
//        x0[N-4] = xT[0];
//        xT.clear();

//        for (unsigned int i=N-4; i!=0; i--)
//        {
//            x0[i-1] = -q[i-1]*x0[i]-v[i-1]*x0[i+1]-u[i-1]*x0[i+2]+r[i-1];
//            for (unsigned int s=0; s<ems_size; s++)
//                if (ind[s]-2 > i) x0[i-1] -= ems[s]*x0[ind[s]];
//            x0[i-1] /= p[i-1];
//        }

//        x.clear();
//        x.push_back(x0);

//        free(p);
//        free(q);
//        free(v);
//        free(u);
//        free(r);
//        ems.clear();
//        ind.clear();

//    }
//    else
//    {
//        DoubleMatrix* p = new DoubleMatrix[N-2]; for (unsigned int i=0; i<N-2; i++) p[i].resize(en, en, 0.0);
//        DoubleMatrix* q = new DoubleMatrix[N-2]; for (unsigned int i=0; i<N-2; i++) q[i].resize(en, en, 0.0);
//        DoubleMatrix* v = new DoubleMatrix[N-2]; for (unsigned int i=0; i<N-2; i++) v[i].resize(en, en, 0.0);
//        DoubleMatrix* u = new DoubleMatrix[N-2]; for (unsigned int i=0; i<N-2; i++) u[i].resize(en, en, 0.0);
//        DoubleMatrix* r = new DoubleMatrix[N-2]; for (unsigned int i=0; i<N-2; i++) r[i].resize(en, 1, 0.0);

//        std::vector<unsigned int> ind;
//        std::vector<DoubleMatrix> ems;

//        /**********************************************************************
//         *                  Discretization using trapesium rule
//         *********************************************************************/
//        DoubleMatrix* b = new DoubleMatrix[N+1];
//        for (unsigned int n=0; n<=N; n++) b[n].resize(en, en, 0.0);

//        /********************* discretisation ********************************/
//        discretization(cs, b, h, N);
//        /*********************************************************************/

//        p[0] = b[0];
//        q[0] = b[1];
//        v[0] = b[2];
//        u[0] = b[3];
//        for (unsigned int n=4; n<=N; n++)
//        {
//            if (!b[n].zeroMatrix())
//            {
//                ind.push_back(n);
//                ems.push_back(b[n]);
//            }
//        }
//        free(b);
//        r[0] = rs;

//        unsigned int ind_size = ind.size();
//        unsigned int ems_size = ems.size();

//        /**********************************************************************
//         *                          End of discretization
//         *********************************************************************/

//        /**********************************************************************
//         *                      Finding function at end of grid
//         *********************************************************************/

//        DoubleMatrix m(en, en);
//        DoubleMatrix alpha1(en, en, 0.0);
//        DoubleMatrix alpha2(en, en, 0.0);
//        DoubleMatrix alpha3(en, en, 0.0);
//        DoubleMatrix alpha4(en, en, 0.0);
//        DoubleMatrix alpha0(en, 1, 0.0);

//        for (unsigned int n=0; n<=N-4; n++)
//        {
//            PointNodeODE node(n*h, n);

//            for (unsigned int row=0; row<en; row++)
//            {
//                for (unsigned int col=0; col<en; col++)
//                {
//                    m[row][col] = -12.0*h*A(node,row,col);
//                    if (row==col) m[row][col] -= 25.0;

//                    alpha1[row][col] = 0.0;
//                    if (row==col) alpha1[row][col] = -48.0;

//                    alpha2[row][col] = 0.0;
//                    if (row==col) alpha2[row][col] = +36.0;

//                    alpha3[row][col] = 0.0;
//                    if (row==col) alpha3[row][col] = -16.0;

//                    alpha4[row][col] = 0.0;
//                    if (row==col) alpha4[row][col] = +3.0;

//                }
//                alpha0[row][0] = +12.0*h*B(node,row);
//            }

//            m.inverse();

//            alpha1 = m*alpha1;
//            alpha2 = m*alpha2;
//            alpha3 = m*alpha3;
//            alpha4 = m*alpha4;
//            alpha0 = m*alpha0;

//            r[n+1] = r[n] - p[n]*alpha0;
//            u[n+1] = p[n]*alpha4;
//            v[n+1] = p[n]*alpha3 + u[n];
//            q[n+1] = p[n]*alpha2 + v[n];
//            p[n+1] = p[n]*alpha1 + q[n];

//            for (unsigned int i=0; i<ind_size; i++)
//            {
//                if (n+4 == ind[i]) u[n+1] += ems[i];
//            }
//        }

//        m.clear();
//        alpha0.clear();
//        alpha1.clear();
//        alpha2.clear();
//        alpha3.clear();
//        alpha4.clear();

//        DoubleMatrix M(5*en,5*en);
//        DoubleVector C(5*en);
//        DoubleVector xT(5*en);

//        for (unsigned int row=0*en; row<1*en; row++)
//        {
//            for (unsigned int col=0*en; col<1*en; col++) M[row][col] = 0.0;
//            for (unsigned int col=1*en; col<2*en; col++) M[row][col] = p[N-3][row%en][col%en];
//            for (unsigned int col=2*en; col<3*en; col++) M[row][col] = q[N-3][row%en][col%en];
//            for (unsigned int col=3*en; col<4*en; col++) M[row][col] = v[N-3][row%en][col%en];
//            for (unsigned int col=4*en; col<5*en; col++) M[row][col] = u[N-3][row%en][col%en];
//            C[row] = r[N-3][row%en][0];
//        }

//        PointNodeODE nodeN3((N-3)*h, N-3);
//        for (unsigned int row=1*en; row<2*en; row++)
//        {
//            for (unsigned int col=0*en; col<1*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -3.0; }
//            for (unsigned int col=1*en; col<2*en; col++) { M[row][col] = -12.0*h*A(nodeN3, row%en, col%en); if (row%en == col%en) M[row][col] += -10.0; }
//            for (unsigned int col=2*en; col<3*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +18.0; }
//            for (unsigned int col=3*en; col<4*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -6.0; }
//            for (unsigned int col=4*en; col<5*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +1.0; }
//            C[row] = 12.0*h*B(nodeN3, row%en);
//        }

//        PointNodeODE nodeN2((N-2)*h, N-2);
//        for (unsigned int row=2*en; row<3*en; row++)
//        {
//            for (unsigned int col=0*en; col<1*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +1.0; }
//            for (unsigned int col=1*en; col<2*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -8.0; }
//            for (unsigned int col=2*en; col<3*en; col++) { M[row][col] = -12.0*h*A(nodeN2, row%en, col%en); }
//            for (unsigned int col=3*en; col<4*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +8.0; }
//            for (unsigned int col=4*en; col<5*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -1.0; }
//            C[row] = 12.0*h*B(nodeN2, row%en);
//        }

//        PointNodeODE nodeN1((N-1)*h, N-1);
//        for (unsigned int row=3*en; row<4*en; row++)
//        {
//            for (unsigned int col=0*en; col<1*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -1.0; }
//            for (unsigned int col=1*en; col<2*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +6.0; }
//            for (unsigned int col=2*en; col<3*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -18.0; }
//            for (unsigned int col=3*en; col<4*en; col++) { M[row][col] = -12.0*h*A(nodeN1, row%en, col%en); if (row%en == col%en) M[row][col] += +10.0; }
//            for (unsigned int col=4*en; col<5*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +3.0; }
//            C[row] = 12.0*h*B(nodeN1, row%en);
//        }

//        PointNodeODE nodeN0(N*h, N);
//        for (unsigned int row=4*en; row<5*en; row++)
//        {
//            for (unsigned int col=0*en; col<1*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +3.0; }
//            for (unsigned int col=1*en; col<2*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -16.0; }
//            for (unsigned int col=2*en; col<3*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +36.0; }
//            for (unsigned int col=3*en; col<4*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -48.0; }
//            for (unsigned int col=4*en; col<5*en; col++) { M[row][col] = -12.0*h*A(nodeN0, row%en, col%en); if (row%en == col%en) M[row][col] += +25.0;}
//            C[row] = 12.0*h*B(nodeN0, row%en);
//        }

//        LinearEquation::GaussianElimination(M, C, xT);
//        M.clear();
//        C.clear();

//        std::vector<DoubleMatrix> x0(N+1);
//        for (unsigned int i=0; i<=N; i++) x0[i].resize(en,1);

//        for (unsigned int row=0; row<en; row++)
//        {
//            x0[N-4][row][0] = xT[row+0*en];
//            x0[N-3][row][0] = xT[row+1*en];
//            x0[N-2][row][0] = xT[row+2*en];
//            x0[N-1][row][0] = xT[row+3*en];
//            x0[N-0][row][0] = xT[row+4*en];
//        }
//        xT.clear();

//        for (unsigned int i=N-4; i!=0; i--)
//        {
//            x0[i-1] = r[i-1] - q[i-1]*x0[i] - v[i-1]*x0[i+1] - u[i-1]*x0[i+2];
//            for (unsigned int s=0; s<ems_size; s++)
//                if (ind[s]-2 > i) x0[i-1] += -1.0*(ems[s]*x0[ind[s]]);
//            p[i-1].inverse();
//            x0[i-1] = p[i-1] * x0[i-1];
//        }

//        x.resize(en);
//        for (unsigned int i=0; i<en; i++) x[i].resize(N+1);

//        for (unsigned int i=0; i<=N; i++)
//        {
//            for (unsigned int row=0; row<en; row++)
//            {
//                x[row][i] = x0[i][row][0];
//            }
//        }

//        for (unsigned int i=0; i<=N; i++)
//        {
//            x0[i].clear();
//        }
//        x0.clear();

//        ind.clear();
//        ems.clear();

//        for (unsigned int i=0; i<N-2; i++)
//        {
//            r[i].clear();
//            u[i].clear();
//            v[i].clear();
//            q[i].clear();
//            p[i].clear();
//        }
//        delete [] u;
//        delete [] v;
//        delete [] r;
//        delete [] q;
//        delete [] p;
//    }
//}

//void LinearODE1stOrder::highOder6Accuracy(const std::vector<Condition> &cs, const DoubleVector& rs, std::vector<DoubleVector> &x, Direction direction UNUSED_PARAM)
//{
//    unsigned int en = equationsNumber();

//    double h = dimension().step();
//    unsigned int N = dimension().size();

//    if (en == 1)
//    {
//        double *p = (double*) malloc(sizeof(double)*(N-4));
//        double *q = (double*) malloc(sizeof(double)*(N-4));
//        double *v = (double*) malloc(sizeof(double)*(N-4));
//        double *u = (double*) malloc(sizeof(double)*(N-4));
//        double *w = (double*) malloc(sizeof(double)*(N-4));
//        double *z = (double*) malloc(sizeof(double)*(N-4));
//        double *r = (double*) malloc(sizeof(double)*(N-4));

//        std::vector<unsigned int> ind;
//        std::vector<double>       ems;

//        /**********************************************************************
//         *                  Discretization using trapesium rule
//         *********************************************************************/
//        double *b = (double*) malloc(sizeof(double)*(N+1));
//        for (unsigned int m=0; m<=N; m++) b[m] = 0.0;

//        /********************* discretisation ********************************/
//        discretization(cs, b, h, N);
//        /*********************************************************************/

//        p[0] = b[0];
//        q[0] = b[1];
//        v[0] = b[2];
//        u[0] = b[3];
//        w[0] = b[4];
//        z[0] = b[5];
//        for (unsigned int n=6; n<=N; n++)
//        {
//            if (fabs(b[n]) > DBL_EPSILON)
//            {
//                ind.push_back(n);
//                ems.push_back(b[n]);
//            }
//        }
//        free(b);
//        r[0] = rs[0];

//        unsigned int ind_size = ind.size();
//        unsigned int ems_size = ems.size();

//        /**********************************************************************
//         *                          End of discretization
//         *********************************************************************/

//        /**********************************************************************
//         *                      Finding function at end of grid
//         *********************************************************************/


//        for (unsigned int n=0; n<=N-6; n++)
//        {
//            PointNodeODE node(n*h, n);

//            double m = +1.0/(-60.0*h*A(node)-147.0);
//            double alpha1 = -360.0*m;
//            double alpha2 = +450.0*m;
//            double alpha3 = -400.0*m;
//            double alpha4 = +225.0*m;
//            double alpha5 = -72.0*m;
//            double alpha6 = +10.0*m;
//            double alpha0 = +60.0*h*B(node)*m;

//            r[n+1] = r[n] - p[n]*alpha0;
//            z[n+1] = p[n]*alpha6;
//            w[n+1] = p[n]*alpha5 + z[n];
//            u[n+1] = p[n]*alpha4 + w[n];
//            v[n+1] = p[n]*alpha3 + u[n];
//            q[n+1] = p[n]*alpha2 + v[n];
//            p[n+1] = p[n]*alpha1 + q[n];

//            for (unsigned int i=0; i<ind_size; i++)
//            {
//                if (n+6 == ind[i]) z[n+1] += ems[i];
//            }
//        }

//        DoubleMatrix m(7,7);
//        DoubleVector c(7);
//        DoubleVector xT(7);

//        m[0][0] = 0.0;
//        m[0][1] = p[N-5];
//        m[0][2] = q[N-5];
//        m[0][3] = v[N-5];
//        m[0][4] = u[N-5];
//        m[0][5] = w[N-5];
//        m[0][6] = z[N-5];
//        c[0] = r[N-5];

//        PointNodeODE nodeN5((N-5)*h, N-5);
//        m[1][0] = -10.0;
//        m[1][1] = -77.0 - 60.0*h*A(nodeN5);
//        m[1][2] = +150.0;
//        m[1][3] = -100.0;
//        m[1][4] = +50.0;
//        m[1][5] = -15.0;
//        m[1][6] = +2.0;
//        c[1] = 60.0*h*B(nodeN5);

//        PointNodeODE nodeN4((N-4)*h, N-4);
//        m[2][0] = +2.0;
//        m[2][1] = -24.0;
//        m[2][2] = -35.0 - 60.0*h*A(nodeN4);
//        m[2][3] = +80.0;
//        m[2][4] = -30.0;
//        m[2][5] = +8.0;
//        m[2][6] = -1.0;
//        c[2] = 60.0*h*B(nodeN4);

//        PointNodeODE nodeN3((N-3)*h, N-3);
//        m[3][0] = -1.0;
//        m[3][1] = +9.0;
//        m[3][2] = -45.0;
//        m[3][3] = -60.0*h*A(nodeN3);
//        m[3][4] = +45.0;
//        m[3][5] = -9.0;
//        m[3][6] = +1.0;
//        c[3] = 60.0*h*B(nodeN3);

//        PointNodeODE nodeN2((N-2)*h, N-2);
//        m[4][0] = +1.0;
//        m[4][1] = -8.0;
//        m[4][2] = +30.0;
//        m[4][3] = -80.0;
//        m[4][4] = +35.0 - 60.0*h*A(nodeN2);
//        m[4][5] = +24.0;
//        m[4][6] = -2.0;
//        c[4] = 60.0*h*B(nodeN2);

//        PointNodeODE nodeN1((N-1)*h, N-1);
//        m[5][0] = -2.0;
//        m[5][1] = +15.0;
//        m[5][2] = -50.0;
//        m[5][3] = +100.0;
//        m[5][4] = -150.0;
//        m[5][5] = +77.0 - 60.0*h*A(nodeN1);
//        m[5][6] = +10.0;
//        c[5] = 60.0*h*B(nodeN1);

//        PointNodeODE nodeN0(N*h, N);
//        m[6][0] = +10.0;
//        m[6][1] = -72.0;
//        m[6][2] = +225.0;
//        m[6][3] = -400.0;
//        m[6][4] = +450.0;
//        m[6][5] = -360.0;
//        m[6][6] = +147.0 - 60.0*h*A(nodeN0);
//        c[6] = 60.0*h*B(nodeN0);

//        LinearEquation::GaussianElimination(m, c, xT);
//        m.clear();
//        c.clear();

//        /**********************************************************************
//         *                      Finding function at end of grid
//         *********************************************************************/
//        DoubleVector x0(N+1);
//        x0[N-0] = xT[6];
//        x0[N-1] = xT[5];
//        x0[N-2] = xT[4];
//        x0[N-3] = xT[3];
//        x0[N-4] = xT[2];
//        x0[N-5] = xT[1];
//        x0[N-6] = xT[0];
//        xT.clear();

//        for (unsigned int i=N-6; i!=0; i--)
//        {
//            x0[i-1] = -q[i-1]*x0[i]-v[i-1]*x0[i+1]-u[i-1]*x0[i+2]-w[i-1]*x0[i+3]-z[i-1]*x0[i+4]+r[i-1];
//            for (unsigned int s=0; s<ems_size; s++)
//                if (ind[s]-4 > i) x0[i-1] -= ems[s]*x0[ind[s]];
//            x0[i-1] /= p[i-1];
//        }

//        x.clear();
//        x.push_back(x0);

//        free(p);
//        free(q);
//        free(v);
//        free(u);
//        free(r);
//        free(w);
//        free(z);
//        ems.clear();
//        ind.clear();

//    }
//    else
//    {
//        DoubleMatrix* p = new DoubleMatrix[N-4]; for (unsigned int i=0; i<N-4; i++) p[i].resize(en, en, 0.0);
//        DoubleMatrix* q = new DoubleMatrix[N-4]; for (unsigned int i=0; i<N-4; i++) q[i].resize(en, en, 0.0);
//        DoubleMatrix* v = new DoubleMatrix[N-4]; for (unsigned int i=0; i<N-4; i++) v[i].resize(en, en, 0.0);
//        DoubleMatrix* u = new DoubleMatrix[N-4]; for (unsigned int i=0; i<N-4; i++) u[i].resize(en, en, 0.0);
//        DoubleMatrix* w = new DoubleMatrix[N-4]; for (unsigned int i=0; i<N-4; i++) w[i].resize(en, en, 0.0);
//        DoubleMatrix* z = new DoubleMatrix[N-4]; for (unsigned int i=0; i<N-4; i++) z[i].resize(en, en, 0.0);
//        DoubleMatrix* r = new DoubleMatrix[N-4]; for (unsigned int i=0; i<N-4; i++) r[i].resize(en, 1, 0.0);

//        std::vector<unsigned int> ind;
//        std::vector<DoubleMatrix> ems;

//        /**********************************************************************
//         *                  Discretization using trapesium rule
//         *********************************************************************/
//        DoubleMatrix* b = new DoubleMatrix[N+1];
//        for (unsigned int n=0; n<=N; n++) b[n].resize(en, en, 0.0);

//        /********************* discretisation ********************************/
//        discretization(cs, b, h, N);
//        /*********************************************************************/

//        p[0] = b[0];
//        q[0] = b[1];
//        v[0] = b[2];
//        u[0] = b[3];
//        w[0] = b[4];
//        z[0] = b[5];
//        for (unsigned int n=6; n<=N; n++)
//        {
//            if (!b[n].zeroMatrix())
//            {
//                ind.push_back(n);
//                ems.push_back(b[n]);
//            }
//        }
//        free(b);
//        r[0] = rs;

//        unsigned int ind_size = ind.size();
//        unsigned int ems_size = ems.size();

//        /**********************************************************************
//         *                          End of discretization
//         *********************************************************************/

//        /**********************************************************************
//         *                      Finding function at end of grid
//         *********************************************************************/

//        DoubleMatrix m(en, en);
//        DoubleMatrix alpha1(en, en, 0.0);
//        DoubleMatrix alpha2(en, en, 0.0);
//        DoubleMatrix alpha3(en, en, 0.0);
//        DoubleMatrix alpha4(en, en, 0.0);
//        DoubleMatrix alpha5(en, en, 0.0);
//        DoubleMatrix alpha6(en, en, 0.0);
//        DoubleMatrix alpha0(en, 1, 0.0);

//        for (unsigned int n=0; n<=N-6; n++)
//        {
//            PointNodeODE node(n*h, n);

//            for (unsigned int row=0; row<en; row++)
//            {
//                for (unsigned int col=0; col<en; col++)
//                {
//                    m[row][col] = -60.0*h*A(node,row,col);
//                    if (row==col) m[row][col] -= 147.0;

//                    alpha1[row][col] = 0.0;
//                    if (row==col) alpha1[row][col] = -360.0;

//                    alpha2[row][col] = 0.0;
//                    if (row==col) alpha2[row][col] = +450.0;

//                    alpha3[row][col] = 0.0;
//                    if (row==col) alpha3[row][col] = -400.0;

//                    alpha4[row][col] = 0.0;
//                    if (row==col) alpha4[row][col] = +225.0;

//                    alpha5[row][col] = 0.0;
//                    if (row==col) alpha5[row][col] = -72.0;

//                    alpha6[row][col] = 0.0;
//                    if (row==col) alpha6[row][col] = +10.0;
//                }
//                alpha0[row][0] = +60.0*h*B(node,row);
//            }

//            m.inverse();

//            alpha1 = m*alpha1;
//            alpha2 = m*alpha2;
//            alpha3 = m*alpha3;
//            alpha4 = m*alpha4;
//            alpha5 = m*alpha5;
//            alpha6 = m*alpha6;
//            alpha0 = m*alpha0;

//            r[n+1] = r[n] - p[n]*alpha0;
//            z[n+1] = p[n]*alpha6;
//            w[n+1] = p[n]*alpha5 + z[n];
//            u[n+1] = p[n]*alpha4 + w[n];
//            v[n+1] = p[n]*alpha3 + u[n];
//            q[n+1] = p[n]*alpha2 + v[n];
//            p[n+1] = p[n]*alpha1 + q[n];

//            for (unsigned int i=0; i<ind_size; i++)
//            {
//                if (n+6 == ind[i]) z[n+1] += ems[i];
//            }
//        }

//        m.clear();
//        alpha0.clear();
//        alpha1.clear();
//        alpha2.clear();
//        alpha3.clear();
//        alpha4.clear();
//        alpha5.clear();
//        alpha6.clear();

//        DoubleMatrix M(7*en,7*en);
//        DoubleVector C(7*en);
//        DoubleVector xT(7*en);

//        for (unsigned int row=0*en; row<1*en; row++)
//        {
//            for (unsigned int col=0*en; col<1*en; col++) M[row][col] = 0.0;
//            for (unsigned int col=1*en; col<2*en; col++) M[row][col] = p[N-5][row%en][col%en];
//            for (unsigned int col=2*en; col<3*en; col++) M[row][col] = q[N-5][row%en][col%en];
//            for (unsigned int col=3*en; col<4*en; col++) M[row][col] = v[N-5][row%en][col%en];
//            for (unsigned int col=4*en; col<5*en; col++) M[row][col] = u[N-5][row%en][col%en];
//            for (unsigned int col=5*en; col<6*en; col++) M[row][col] = w[N-5][row%en][col%en];
//            for (unsigned int col=6*en; col<7*en; col++) M[row][col] = z[N-5][row%en][col%en];
//            C[row] = r[N-5][row%en][0];
//        }

//        PointNodeODE nodeN5((N-5)*h, N-5);
//        for (unsigned int row=1*en; row<2*en; row++)
//        {
//            for (unsigned int col=0*en; col<1*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -10.0; }
//            for (unsigned int col=1*en; col<2*en; col++) { M[row][col] = -60.0*h*A(nodeN5, row%en, col%en); if (row%en == col%en) M[row][col] += -77.0; }
//            for (unsigned int col=2*en; col<3*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +150.0; }
//            for (unsigned int col=3*en; col<4*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -100.0; }
//            for (unsigned int col=4*en; col<5*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +50.0; }
//            for (unsigned int col=5*en; col<6*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -15.0; }
//            for (unsigned int col=6*en; col<7*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +2.0; }
//            C[row] = 60.0*h*B(nodeN5, row%en);
//        }

//        PointNodeODE nodeN4((N-4)*h, N-4);
//        for (unsigned int row=2*en; row<3*en; row++)
//        {
//            for (unsigned int col=0*en; col<1*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +2.0; }
//            for (unsigned int col=1*en; col<2*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -24.0; }
//            for (unsigned int col=2*en; col<3*en; col++) { M[row][col] = -60.0*h*A(nodeN4, row%en, col%en); if (row%en == col%en) M[row][col] += -35.0; }
//            for (unsigned int col=3*en; col<4*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +80.0; }
//            for (unsigned int col=4*en; col<5*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -30.0; }
//            for (unsigned int col=5*en; col<6*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +8.0; }
//            for (unsigned int col=6*en; col<7*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -1.0; }
//            C[row] = 60.0*h*B(nodeN4, row%en);
//        }

//        PointNodeODE nodeN3((N-3)*h, N-3);
//        for (unsigned int row=3*en; row<4*en; row++)
//        {
//            for (unsigned int col=0*en; col<1*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -1.0; }
//            for (unsigned int col=1*en; col<2*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +9.0; }
//            for (unsigned int col=2*en; col<3*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -45.0; }
//            for (unsigned int col=3*en; col<4*en; col++) { M[row][col] = -60.0*h*A(nodeN3, row%en, col%en);}
//            for (unsigned int col=4*en; col<5*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +45.0; }
//            for (unsigned int col=5*en; col<6*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -9.0; }
//            for (unsigned int col=6*en; col<7*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +1.0; }
//            C[row] = 60.0*h*B(nodeN3, row%en);
//        }

//        PointNodeODE nodeN2((N-2)*h, N-2);
//        for (unsigned int row=4*en; row<5*en; row++)
//        {
//            for (unsigned int col=0*en; col<1*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +1.0; }
//            for (unsigned int col=1*en; col<2*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -8.0; }
//            for (unsigned int col=2*en; col<3*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +30.0; }
//            for (unsigned int col=3*en; col<4*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -80.0; }
//            for (unsigned int col=4*en; col<5*en; col++) { M[row][col] = -60.0*h*A(nodeN2, row%en, col%en); if (row%en == col%en) M[row][col] += +35.0; }
//            for (unsigned int col=5*en; col<6*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +24.0; }
//            for (unsigned int col=6*en; col<7*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -2.0; }
//            C[row] = 60.0*h*B(nodeN2, row%en);
//        }

//        PointNodeODE nodeN1((N-1)*h, N-1);
//        for (unsigned int row=5*en; row<6*en; row++)
//        {
//            for (unsigned int col=0*en; col<1*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -2.0; }
//            for (unsigned int col=1*en; col<2*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +15.0; }
//            for (unsigned int col=2*en; col<3*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -50.0; }
//            for (unsigned int col=3*en; col<4*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +100.0; }
//            for (unsigned int col=4*en; col<5*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -150.0; }
//            for (unsigned int col=5*en; col<6*en; col++) { M[row][col] = -60.0*h*A(nodeN1, row%en, col%en); if (row%en == col%en) M[row][col] += +77.0; }
//            for (unsigned int col=6*en; col<7*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +10.0; }
//            C[row] = 60.0*h*B(nodeN1, row%en);
//        }

//        PointNodeODE nodeN0(N*h, N);
//        for (unsigned int row=6*en; row<7*en; row++)
//        {
//            for (unsigned int col=0*en; col<1*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +10.0; }
//            for (unsigned int col=1*en; col<2*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -72.0; }
//            for (unsigned int col=2*en; col<3*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +225.0; }
//            for (unsigned int col=3*en; col<4*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -400.0; }
//            for (unsigned int col=4*en; col<5*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += +450.0; }
//            for (unsigned int col=5*en; col<6*en; col++) { M[row][col] = +0.0; if (row%en == col%en) M[row][col] += -360.0; }
//            for (unsigned int col=6*en; col<7*en; col++) { M[row][col] = -60.0*h*A(nodeN0, row%en, col%en); if (row%en == col%en) M[row][col] += 147.0; }
//            C[row] = 60.0*h*B(nodeN0, row%en);
//        }

//        LinearEquation::GaussianElimination(M, C, xT);
//        M.clear();
//        C.clear();

//        std::vector<DoubleMatrix> x0(N+1);
//        for (unsigned int i=0; i<=N; i++) x0[i].resize(en,1);

//        for (unsigned int row=0; row<en; row++)
//        {
//            x0[N-6][row][0] = xT[row+0*en];
//            x0[N-5][row][0] = xT[row+1*en];
//            x0[N-4][row][0] = xT[row+2*en];
//            x0[N-3][row][0] = xT[row+3*en];
//            x0[N-2][row][0] = xT[row+4*en];
//            x0[N-1][row][0] = xT[row+5*en];
//            x0[N-0][row][0] = xT[row+6*en];
//        }
//        xT.clear();

//        for (unsigned int i=N-6; i!=0; i--)
//        {
//            x0[i-1] = r[i-1] - q[i-1]*x0[i] - v[i-1]*x0[i+1] - u[i-1]*x0[i+2] - w[i-1]*x0[i+3] - z[i-1]*x0[i+4];
//            for (unsigned int s=0; s<ems_size; s++)
//                if (ind[s]-4 > i) x0[i-1] += -1.0*(ems[s]*x0[ind[s]]);
//            p[i-1].inverse();
//            x0[i-1] = p[i-1] * x0[i-1];
//        }

//        x.resize(en);
//        for (unsigned int i=0; i<en; i++) x[i].resize(N+1);

//        for (unsigned int i=0; i<=N; i++)
//        {
//            for (unsigned int row=0; row<en; row++)
//            {
//                x[row][i] = x0[i][row][0];
//            }
//        }

//        for (unsigned int i=0; i<=N; i++)
//        {
//            x0[i].clear();
//        }
//        x0.clear();

//        ind.clear();
//        ems.clear();

//        for (unsigned int i=0; i<N-4; i++)
//        {
//            r[i].clear();
//            z[i].clear();
//            w[i].clear();
//            u[i].clear();
//            v[i].clear();
//            q[i].clear();
//            p[i].clear();
//        }
//        delete [] z;
//        delete [] w;
//        delete [] u;
//        delete [] v;
//        delete [] r;
//        delete [] q;
//        delete [] p;
//    }
//}

IFirstOrderLinearODEFVP::IFirstOrderLinearODEFVP() {}

IFirstOrderLinearODEFVP::IFirstOrderLinearODEFVP(const IFirstOrderLinearODEFVP &) {}

IFirstOrderLinearODEFVP& IFirstOrderLinearODEFVP::operator=(const IFirstOrderLinearODEFVP &other)
{
    if (this == &other) { return *this; }
    return *this;
}

IFirstOrderLinearODEFVP::~IFirstOrderLinearODEFVP() {}

void IFirstOrderLinearODEFVP::solveFinalValueProblem(std::vector<DoubleVector> &rv, ODESolverMethod method) const
{
    switch (method)
    {
    case ODESolverMethod::RUNGE_KUTTA_2: solveFinalValueProblemRK2(rv); break;
    case ODESolverMethod::RUNGE_KUTTA_4: solveFinalValueProblemRK4(rv); break;
    case ODESolverMethod::RUNGE_KUTTA_6: solveFinalValueProblemRK6(rv); break;
    case ODESolverMethod::EULER: solveFinalValueProblemEuler(rv); break;
    case ODESolverMethod::EULER_MOD: solveFinalValueProblemEulerMod(rv); break;
    }
}

void IFirstOrderLinearODEFVP::solveFinalValueProblemRK2(std::vector<DoubleVector> &rv) const
{
    const Dimension &dim = dimension();
    const int min = dim.min();
    const int max = dim.max();
    const double h1 = dim.step();
    const double h2 = 0.5*h1;
    const unsigned int m = count();
    const unsigned int size = static_cast<unsigned int>(max-min);

    rv.resize(size+1); for (unsigned int i=0; i<=size; i++) rv[i].resize(m);

    double *k1 = new double[m];
    double *k2 = new double[m];

    for (unsigned int r=1; r<=m; r++)
    {
        rv[size][r-1] = final(FinalCondition::FinalValue, r);
    }
    iterationInfo(rv[size], PointNodeODE(max*h1, max));

    unsigned int i=size;
    PointNodeODE node1, node2;
    double *v1 = new double[m];
    for (int n=max; n>min; n--, i--)
    {
        node1.i = n+0; node1.x = node1.i*h1;
        node2.i = n-1; node2.x = node2.i*h1;

        const DoubleVector &v0 = rv[i];

        // k1
        for (unsigned int row=1, j=0; row<=m; row++, j++)
        {
            double sum = B(node1, row);
            for (unsigned int col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node1, row, col)*v0[j];
            }
            k1[row-1] = sum;
            v1[j]=v0[j]-h1*k1[j];
        }

        // k2
        for (unsigned int row=1; row<=m; row++)
        {
            double sum = B(node2, row);
            for (unsigned int col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node2, row, col)*v1[j];
            }
            k2[row-1] = sum;
        }

        for (unsigned int j=0; j<m; j++)
        {
            rv[i-1][j] = rv[i][j] -  h2 * (k1[j] + k2[j]);
        }
        iterationInfo(rv[i-1], node2);
    }

    delete [] v1;

    delete [] k2;
    delete [] k1;
}

void IFirstOrderLinearODEFVP::solveFinalValueProblemRK4(std::vector<DoubleVector> &rv) const
{
    const Dimension &dim = dimension();
    const int min = dim.min();
    const int max = dim.max();
    const double h1 = dim.step();
    const double h2 = h1/2.0;
    const double h6 = h1/6.0;
    const size_t m = count();
    const size_t size = static_cast<size_t>(max-min);

    rv.resize(size+1); for (unsigned int i=0; i<=size; i++) rv[i].resize(size+1);

    double *k1 = new double[m];
    double *k2 = new double[m];
    double *k3 = new double[m];
    double *k4 = new double[m];

    for (unsigned int r=1; r<=m; r++)
    {
        rv[size][r-1] = final(FinalCondition::FinalValue, r);
    }
    iterationInfo(rv[size], PointNodeODE(max*h1, max));

    size_t i=size;
    PointNodeODE node1, node2, node3, node4;
    double *v1 = new double[m];
    for (int n=max; n>min; n--, i--)
    {
        node1.i = n;   node1.x = n*h1;
        node2.i = n;   node2.x = n*h1-h2;
        node3.i = n;   node3.x = n*h1-h2;
        node4.i = n-1; node4.x = node4.i*h1;

        const DoubleVector &v0 = rv[i];

        // k1
        for (unsigned int j=0; j<m; v1[j]=v0[j], j++);
        for (unsigned int row=1; row<=m; row++)
        {
            double sum = B(node1, row);
            for (unsigned int col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node1, row, col)*v1[j];
            }
            k1[row-1] = sum;
        }

        // k2
        for (unsigned int j=0; j<m; v1[j]=v0[j]-h2*k1[j], j++);
        for (unsigned int row=1; row<=m; row++)
        {
            double sum = B(node2, row);
            for (unsigned int col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node2, row, col)*v1[j];
            }
            k2[row-1] = sum;
        }

        // k3
        for (unsigned int j=0; j<m; v1[j]=v0[j]-h2*k2[j], j++);
        for (unsigned int row=1; row<=m; row++)
        {
            double sum = B(node3, row);
            for (unsigned int col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node3, row, col)*v1[j];
            }
            k3[row-1] = sum;
        }

        // k4
        for (unsigned int j=0; j<m; v1[j]=v0[j]-h1*k3[j], j++);
        for (unsigned int row=1; row<=m; row++)
        {
            double sum = B(node4, row);
            for (unsigned int col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node4, row, col)*v1[j];
            }
            k4[row-1] = sum;
        }

        for (unsigned int j=0; j<m; j++)
        {
            rv[i-1][j] = rv[i][j] - h6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
        }
        iterationInfo(rv[i-1], node4);
    }

    delete [] v1;

    delete [] k4;
    delete [] k3;
    delete [] k2;
    delete [] k1;
}

void IFirstOrderLinearODEFVP::solveFinalValueProblemRK6(std::vector<DoubleVector> &rv) const { C_UNUSED(rv); }

void IFirstOrderLinearODEFVP::solveFinalValueProblemEuler(std::vector<DoubleVector> &rv) const
{
    const Dimension &dim = dimension();
    const int min = dim.min();
    const int max = dim.max();
    const double h = dim.step();
    const size_t size = static_cast<size_t>(max-min);
    const size_t m = count();

    rv.resize(size+1); for (unsigned int i=0; i<=size; i++) rv[i].resize(m);

    DoubleVector &rv0 = rv[size];
    for (size_t row=1; row<=m; row++)
    {
        rv0[row-1] = final(FinalCondition::FinalValue, row);
    }

    size_t i=size;//, j=0;
    PointNodeODE node;
    for (int n=max; n>min; n--, i--)
    {
        node.i = n; node.x = n*h;

        DoubleVector v0 = rv[i];

        for (size_t row=1; row<=m; row++)
        {
            double sum = h*B(node, row);
            for (size_t col=1; col<=m; col++)
            {
                sum += h*A(node, row, col)*rv[i][col-1];
            }
            rv[i-1][row-1] = rv[i][row-1] - sum;
        }
    }
}

void IFirstOrderLinearODEFVP::start(DoubleVector &x, PointNodeODE &n) const
{
    const Dimension &dim = dimension();
    const int max = dim.max();
    const double h = dim.step();
    const size_t m = count();
    n.i = max;
    n.x = n.i*h;
    x.clear();
    x.resize(m);
    for (size_t row=1; row<=m; row++)
    {
        x[row-1] = final(FinalCondition::FinalValue, row);
    }
    iterationInfo(x, n);
}

void IFirstOrderLinearODEFVP::next(const DoubleVector &x0, const PointNodeODE &n0, DoubleVector &x1, PointNodeODE &n1, ODESolverMethod method) const
{
    const Dimension &dim = dimension();
    const double h = dim.step();
    const size_t m = count();

    n1.i = n0.i - 1; n1.x = n1.i*h;

    if (method == ODESolverMethod::EULER)
    {
        for (size_t row=1; row<=m; row++)
        {
            double sum = h*B(n0, row);
            for (size_t col=1; col<=m; col++)
            {
                sum += h*A(n0, row, col)*x0[col-1];
            }
            x1[row-1] = x0[row-1] - sum;
        }
    }

    if (method == ODESolverMethod::EULER_MOD)
    {
    }

    if (method == ODESolverMethod::RUNGE_KUTTA_2)
    {
    }

    if (method == ODESolverMethod::RUNGE_KUTTA_4)
    {
        const double h2 = h/2.0;
        const double h6 = h/6.0;

        double *k1 = new double[m];
        double *k2 = new double[m];
        double *k3 = new double[m];
        double *k4 = new double[m];

        double *v0 = new double[m];
        double *v1 = new double[m];

        PointNodeODE node1, node2, node3, node4;
        node1.i = n0.i;   node1.x = node1.i*h;
        node2.i = n0.i;   node2.x = node2.i*h+h2;
        node3.i = n0.i;   node3.x = node3.i*h+h2;
        node4.i = n0.i-1; node4.x = node4.i*h;

        memcpy(v0, x0.data(), sizeof (double) * m);

        // k1
        for (unsigned int row=1; row<=m; row++)
        {
            double sum = B(node1, row);
            for (unsigned int col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node1, row, col)*v0[j];
            }
            k1[row-1] = sum;
        }

        // k2
        for (unsigned int j=0; j<m; v1[j]=v0[j]-h2*k1[j], j++);
        for (unsigned int row=1; row<=m; row++)
        {
            double sum = B(node2, row);
            for (unsigned int col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node2, row, col)*v1[j];
            }
            k2[row-1] = sum;
        }

        // k3
        for (unsigned int j=0; j<m; v1[j]=v0[j]-h2*k2[j], j++);
        for (unsigned int row=1; row<=m; row++)
        {
            double sum = B(node3, row);
            for (unsigned int col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node3, row, col)*v1[j];
            }
            k3[row-1] = sum;
        }

        // k4
        for (unsigned int j=0; j<m; v1[j]=v0[j]-h*k3[j], j++);
        for (unsigned int row=1; row<=m; row++)
        {
            double sum = B(node4, row);
            for (unsigned int col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node4, row, col)*v1[j];
            }
            k4[row-1] = sum;
        }

        for (unsigned int j=0; j<m; j++)
        {
            x1[j] = x0[j] - h6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
        }

        delete [] v0;
        delete [] v1;

        delete [] k4;
        delete [] k3;
        delete [] k2;
        delete [] k1;
    }

    if (method == ODESolverMethod::RUNGE_KUTTA_6)
    {
    }

    iterationInfo(x1, n1);
}

void IFirstOrderLinearODEFVP::solveFinalValueProblemEulerMod(std::vector<DoubleVector> &rv) const { C_UNUSED(rv); }

void IFirstOrderLinearODEFVP::solveFinalValueProblem(ODESolverMethod method) const
{
    switch (method)
    {
    case ODESolverMethod::RUNGE_KUTTA_2: solveFinalValueProblemRK2(); break;
    case ODESolverMethod::RUNGE_KUTTA_4: solveFinalValueProblemRK4(); break;
    case ODESolverMethod::RUNGE_KUTTA_6: solveFinalValueProblemRK6(); break;
    case ODESolverMethod::EULER: solveFinalValueProblemEuler(); break;
    case ODESolverMethod::EULER_MOD: solveFinalValueProblemEulerMod(); break;
    }
}

void IFirstOrderLinearODEFVP::solveFinalValueProblemRK2() const
{
    const Dimension &dim = dimension();
    const int min = dim.min();
    const int max = dim.max();
    const double h1 = dim.step();
    const double h2 = 0.5*h1;
    const unsigned int m = count();

    double *k1 = new double[m];
    double *k2 = new double[m];

    double *v0 = new double[m];
    double *v1 = new double[m];

    for (unsigned int r=1; r<=m; r++)
    {
        v0[r-1] = final(FinalCondition::FinalValue, r);
    }
    iterationInfo(DoubleVector(v0, m), PointNodeODE(max*h1, max));

    PointNodeODE node1, node2;
    for (int n=max; n>min; n--)
    {
        node1.i = n+0; node1.x = node1.i*h1;
        node2.i = n-1; node2.x = node2.i*h1;

        // k1
        for (unsigned int row=1, j=0; row<=m; row++, j++)
        {
            double sum = B(node1, row);
            for (unsigned int col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node1, row, col)*v0[j];
            }
            k1[row-1] = sum;
            v1[j]=v0[j]-h1*k1[j];
        }

        // k2
        for (unsigned int row=1; row<=m; row++)
        {
            double sum = B(node2, row);
            for (unsigned int col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node2, row, col)*v1[j];
            }
            k2[row-1] = sum;
        }

        for (unsigned int j=0; j<m; j++)
        {
            v0[j] = v0[j] - h2 * (k1[j] + k2[j]);
        }
        iterationInfo(DoubleVector(v0, m), node2);
    }

    delete [] v0;
    delete [] v1;

    delete [] k2;
    delete [] k1;
}

void IFirstOrderLinearODEFVP::solveFinalValueProblemRK4() const
{
    const Dimension &dim = dimension();
    const int min = dim.min();
    const int max = dim.max();
    const double h1 = dim.step();
    const double h2 = h1/2.0;
    const double h6 = h1/6.0;
    const unsigned int m = count();

    double *k1 = new double[m];
    double *k2 = new double[m];
    double *k3 = new double[m];
    double *k4 = new double[m];

    double *v0 = new double[m];
    double *v1 = new double[m];

    for (unsigned int r=1; r<=m; r++)
    {
        v0[r-1] = final(FinalCondition::FinalValue, r);
    }
    iterationInfo(DoubleVector(v0, m), PointNodeODE(max*h1, max));

    unsigned int i = static_cast<unsigned int>(max-min+1);
    PointNodeODE node1, node2, node3, node4;
    for (int n=max; n>min; n--, i--)
    {
        node1.i = n;   node1.x = node1.i*h1;
        node2.i = n;   node2.x = node2.i*h1-h2;
        node3.i = n;   node3.x = node3.i*h1-h2;
        node4.i = n-1; node4.x = node4.i*h1;

        memcpy(v1, v0, sizeof (double) * m);

        // k1
        for (unsigned int row=1; row<=m; row++)
        {
            double sum = B(node1, row);
            for (unsigned int col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node1, row, col)*v1[j];
            }
            k1[row-1] = sum;
        }

        // k2
        for (unsigned int j=0; j<m; v1[j]=v0[j]-h2*k1[j], j++);
        for (unsigned int row=1; row<=m; row++)
        {
            double sum = B(node2, row);
            for (unsigned int col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node2, row, col)*v1[j];
            }
            k2[row-1] = sum;
        }

        // k3
        for (unsigned int j=0; j<m; v1[j]=v0[j]-h2*k2[j], j++);
        for (unsigned int row=1; row<=m; row++)
        {
            double sum = B(node3, row);
            for (unsigned int col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node3, row, col)*v1[j];
            }
            k3[row-1] = sum;
        }

        // k4
        for (unsigned int j=0; j<m; v1[j]=v0[j]-h1*k3[j], j++);
        for (unsigned int row=1; row<=m; row++)
        {
            double sum = B(node4, row);
            for (unsigned int col=1, j=0; col<=m; col++, j++)
            {
                sum += A(node4, row, col)*v1[j];
            }
            k4[row-1] = sum;
        }

        for (unsigned int j=0; j<m; j++)
        {
            v0[j] = v0[j] - h6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
        }
        iterationInfo(DoubleVector(v0, m), node4);
    }

    delete [] v0;
    delete [] v1;

    delete [] k4;
    delete [] k3;
    delete [] k2;
    delete [] k1;
}

void IFirstOrderLinearODEFVP::solveFinalValueProblemRK6() const { }

void IFirstOrderLinearODEFVP::solveFinalValueProblemEuler() const
{
    const Dimension &dim = dimension();
    const int min = dim.min();
    const int max = dim.max();
    const double h = dim.step();
    const unsigned int m = count();

    double *v0 = new double[m];
    double *v1 = new double[m];

    for (unsigned int row=1; row<=m; row++)
    {
        v0[row-1] = final(FinalCondition::FinalValue, row);
    }
    iterationInfo(DoubleVector(v0, m), PointNodeODE(min*h, max));

    unsigned int i = static_cast<unsigned int>(max-min);
    PointNodeODE node;
    for (int n=max; n>min; n--, i--)
    {
        node.i = n; node.x = n*h;

        for (unsigned int row=1; row<=m; row++)
        {
            double sum = h*B(node, row);
            for (unsigned int col=1; col<=m; col++)
            {
                sum += h*A(node, row, col)*v0[col-1];
            }
            v1[row-1] = v0[row-1] - sum;
        }

        memcpy(v0, v1, sizeof (double)*m);

        iterationInfo(DoubleVector(v0, m), PointNodeODE((n-1)*h, (n-1)));
    }

    delete [] v1;
    delete [] v0;
}

void IFirstOrderLinearODEFVP::solveFinalValueProblemEulerMod() const { }


