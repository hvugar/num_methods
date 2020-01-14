#include "lode1o.h"
#include "nlode1o.h"
#include "../matrix2d.h"
#include "../printer.h"
#include <math.h>
#include <cmethods.h>
#include <float.h>

NonLocalCondition::NonLocalCondition() {}

NonLocalCondition::NonLocalCondition(unsigned int i, const PointNodeODE &node, const DoubleMatrix &m) : i(i), n(node), m(m)  {}

NonLocalCondition::~NonLocalCondition()
{
    m.clear();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

FirstOrderLinearODE::~FirstOrderLinearODE() {}

void FirstOrderLinearODE::discritize(const std::vector<NonLocalCondition> &co, std::vector<NonLocalCondition> &cn, unsigned int k) const
{
    const auto cnd_size = static_cast<unsigned int>( co.size() );
    const auto h = dimension().step();
    const auto min = dimension().min();
    const auto max = dimension().max();
    //const auto sze = dimension().size();

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

void FirstOrderLinearODE::transferOfCondition(const std::vector<NonLocalCondition> &co, const DoubleVector &d, std::vector<DoubleVector> &x, unsigned int k) const
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
    const unsigned int M = count();

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

            if (i%10==0)
            {
                //if (k==2) printf("%d %f %f\n", i, betta[i][1][0][0], betta[i][2][0][0]);
                if (k==4) printf("%4d %10.6f %10.6f %10.6f %10.6f\n", i, betta[i][1][0][0], betta[i][2][0][0], betta[i][3][0][0], betta[i][4][0][0]);
                //if (k==6) printf("%d %f %f %f %f %f %f\n", i, betta[i][1][0][0], betta[i][2][0][0], betta[i][3][0][0], betta[i][4][0][0], betta[i][5][0][0], betta[i][6][0][0]);
            }

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

void FirstOrderLinearODE::transferOfCondition1(const std::vector<NonLocalCondition> &C, const DoubleVector &d, std::vector<DoubleVector> &x, unsigned int k) const
{
    const int min = dimension().min();
    const int max = dimension().max();
    const double h = dimension().step();

    const unsigned int size = static_cast<unsigned int>(max - min + 1);
    const unsigned int end = static_cast<unsigned int>( size-k + 1 );
    const unsigned int sze = static_cast<unsigned int>( size - k );
    const unsigned int M = count();

    DoubleMatrix **coof = new DoubleMatrix*[end+1];
    DoubleMatrix *alpha = new DoubleMatrix[k+1];
    DoubleMatrix mx(M, M);

    DoubleMatrix * const D = new DoubleMatrix[size+1];
    for (unsigned int i=0; i<=size; i++) D[i].resize(M,M);

    for (unsigned int i=0; i<C.size(); i++)
    {
        unsigned int index = static_cast<unsigned int>(C[i].n.i);
        D[index] += C[i].m;
    }

    if (k==2)
    {
        //betta[i][0] = d;
        //for (unsigned int j=1; j<=k; j++) betta[i][j] = D[j-1];
        //

        coof[0][0] = d;
        coof[0][1] = D[0];
        coof[0][2] = D[1];

        for (int i=min-1; i<=end; i++)
        {
            PointNodeODE node((i-1)*h,static_cast<int>(i-1));

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

            coof[i][0] = coof[i-1][0] - coof[i-1][1]*alpha[0];
            coof[i][1] = coof[i-1][1]*alpha[1] + coof[i-1][2];
            coof[i][2] = coof[i-1][1]*alpha[2] + D[k+(i-1)];
        }
    }

    if (k==4)
    {

    }

    if (k==6)
    {

    }

    mx.clear();


    unsigned int n = 0;
    for (int i=min; i<=max; i++, n++)
    {
    }


    for (unsigned int i=0; i<size; i++)
    {
        for (unsigned int j=0; j<k; j++) coof[i][j].clear();
        delete [] coof[i];
    }
    delete [] coof;
}

void FirstOrderLinearODE::solveInitialValueProblem(DoubleVector &rv) const
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
    unsigned int ai = 1; // array index
    for (int i=min; i<max; i++, ai++)
    {
        PointNodeODE node(static_cast<double>(i*h), i);
        double an = A(node)*h + 1.0;
        double bn = B(node)*h;
        rv[ai] = an*rv[ai-1] + bn;
    }
}

void FirstOrderLinearODE::solveInitialValueProblem(std::vector<DoubleVector> &rv, ODESolverMethod method) const
{
    switch (method)
    {
    case ODESolverMethod::RUNGE_KUTTA_2:
        solveInitialValueProblemRK2(rv);
        break;
    case ODESolverMethod::RUNGE_KUTTA_4:
        solveInitialValueProblemRK4(rv);
        break;
    case ODESolverMethod::RUNGE_KUTTA_6:
        solveInitialValueProblemRK4(rv);
        break;
    case ODESolverMethod::EULER:
        solveInitialValueProblemEuler(rv);
        break;
    case ODESolverMethod::EULER_MOD:
        solveInitialValueProblemEulerMod(rv);
        break;
    }
}

void FirstOrderLinearODE::solveInitialValueProblemEuler(std::vector<DoubleVector> &rv) const
{
    const Dimension &dim = dimension();
    const int min = dim.min();
    const int max = dim.max();
    const double h = dim.step();
    const unsigned int size = static_cast<unsigned int>(max-min);
    const unsigned int k = count();

    rv.resize(size+1); for (unsigned int i=0; i<=size; i++) rv[i].resize(k);

    DoubleVector &rv0 = rv[0];
    for (unsigned int row=1; row<=k; row++)
    {
        rv0[row-1] = initial(InitialCondition::InitialValue, row);
    }

    unsigned int ai = 1; // array index
    DoubleMatrix _A(k, k, 0.0);
    DoubleMatrix _B(k, 1, 0.0);
    DoubleMatrix _I = DoubleMatrix::IdentityMatrix(k);
    for (int i=min; i<max; i++, ai++)
    {
        PointNodeODE node(static_cast<double>(i*h), i);

        for (unsigned int row=1; row<=k; row++)
        {
            for (unsigned int col=1; col<=k; col++)
            {
                _A[row-1][col-1] = A(node, row, col);
            }
            _B[row-1][0] = B(node, row);
        }

        rv[ai] = (h*_A + _I)*rv[ai-1] + h*_B;
    }
}

void FirstOrderLinearODE::solveInitialValueProblemRK2(std::vector<DoubleVector> &rv) const
{
    const Dimension &dim = dimension();
    const int min = dim.min();
    const int max = dim.max();
    //const double h = dim.step();
    const unsigned int size = static_cast<unsigned int>(max-min);
    const unsigned int k = count();

    rv.resize(size+1); for (unsigned int i=0; i<=size; i++) rv[i].resize(k);

    DoubleVector &rv0 = rv[0];
    for (unsigned int row=1; row<=k; row++)
    {
        rv0[row-1] = initial(InitialCondition::InitialValue, row);
    }


    //    double *k1 = static_cast<double*>( malloc(sizeof(double)*m) );
    //    double *k2 = static_cast<double*>( malloc(sizeof(double)*m) );

    //    if (direction == L2R)
    //    {
    //        double h2 = h/2.0;
    //        ry[0] = y0;

    ////        double xn = x0;

    //        DoubleVector yn(m);

    //        /* initializing */
    ////        for (unsigned int j=0; j<n; j++) ry[j][0] = y0[j];

    //        for (int i=min+1; i<=max; i++)
    //        {
    //            const unsigned int mi = static_cast<unsigned int>(i-min);
    //            const PointNodeODE node1((i-1)*h, i-1);
    //            const PointNodeODE node2((i-1)*h+h2, i-1);

    //            // k1
    //            //for (unsigned int j=0; j<n; j++) yn[j] = ry[j][(i-1)-min];
    //            for (unsigned int j=0; j<m; j++) k1[j] = f(node1, ry[mi-1], j+1);

    //            // k2
    //            for (unsigned int j=0; j<m; j++) yn[j] = ry[mi-1][j]+h2*k1[j];
    //            for (unsigned int j=0; j<m; j++) k2[j] = f(node2, yn, j+1);

    //            for (unsigned int j=0; j<m; j++) ry[mi][j] = ry[mi-1][j] + h * k2[j];
    //        }
    //    }

    //    const unsigned int n = y0.length();

    //    if (direction == R2L)
    //    {
    //        double h2 = h/2.0;

    //        double xn = x0;

    //        DoubleVector yn(n);

    //        /* initializing */
    //        for (unsigned int j=0; j<n; j++) ry[j][N] = y0[j];

    //        for (unsigned int i=N-1; i!=UINT32_MAX; i--)
    //        {
    //            // k1
    //            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][i+1];
    //            for (unsigned int j=0; j<n; j++) k1[j] = f(xn, yn, i+1, j);

    //            // k2
    //            for (unsigned int j=0; j<n; j++) yn[j] = ry[j][i+1]-h2*k1[j];
    //            for (unsigned int j=0; j<n; j++) k2[j] = f(xn-h2, yn, i+1, j);

    //            for (unsigned int j=0; j<n; j++) ry[j][i] = ry[j][i+1] - h * k2[j];

    //            xn -= h;
    //        }
    //    }

    //    free(k2);
    //    free(k1);
}

void FirstOrderLinearODE::solveInitialValueProblemRK4(std::vector<DoubleVector> &) const {}

void FirstOrderLinearODE::solveInitialValueProblemRK6(std::vector<DoubleVector> &) const {}

void FirstOrderLinearODE::solveInitialValueProblemEulerMod(std::vector<DoubleVector> &) const {}

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
