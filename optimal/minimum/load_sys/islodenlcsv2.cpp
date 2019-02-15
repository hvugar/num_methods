#include "islodenlcsv2.h"

void ISystemLinearODENonLocalContionsV2::calculateBackward()
{}

void ISystemLinearODENonLocalContionsV2::addCondition(const Condition &nlsc)
{
    nlscs.push_back(nlsc);
}

void ISystemLinearODENonLocalContionsV2::addLoadPoint(const LoadPoint &lpnt)
{
    lpnts.push_back(lpnt);
}

void ISystemLinearODENonLocalContionsV2::setRightSize(const DoubleVector &gamma)
{
    this->gamma = gamma;
}

void ISystemLinearODENonLocalContionsV2::calculateCauchyProblem(const Condition &sc, const Condition &ec,
                                                                const DoubleVector &x0, std::vector<DoubleVector> &rx,
                                                                double h)
{
    class CauchyProblemM1stOrderA1 : public NonLinearODE1stOrder
    {
    public:
        CauchyProblemM1stOrderA1(ISystemLinearODENonLocalContionsV2 &parent) : p(parent) {}
    protected:
        virtual double f(double t, const DoubleVector &x, unsigned int k, unsigned int i) const
        {
            unsigned int n = 3;//p.systemOrder();
            unsigned int k1 = p.lpnts.size();
            double _SO = S0(t,x,k);

            // Alpha
            if (i<n)
            {
                double res = _SO*x[i];
                TimeNodePDE node;
                node.t = t;
                node.i = k;
                for (unsigned int row=0; row<n; row++) res -= p.A(node,row,i)*x[row];
                return res;
            }
            // Bettas
            if (i>=n && i<(k1+1)*n)
            {
                unsigned int s = (i-n)/n;
                unsigned int col = i%n;
                double res = _SO*x[i];
                TimeNodePDE node;
                node.t = t;
                node.i = k;
                for (unsigned int row=0; row<n; row++) res -= p.B(node,s,row,col)*x[row];
                return res;
            }
            else if (i==(k1+1)*n)
            {
                double res = _SO*x[i];
                TimeNodePDE node;
                node.t = t;
                node.i = k;
                for (unsigned int row=0; row<n; row++) res += p.C(node,row)*x[row];
                return res;
            }
            else
            {
                return _SO*x[(k1+1)*n+1];
            }
            return 0.0;//NAN;
        }

        double S0(double t, const DoubleVector &x, unsigned int k) const
        {
            unsigned int n = 3;//p.systemOrder();
            unsigned int k1 = p.lpnts.size();

            // Calculating alpha
            double s1 = 0.0;
            for (unsigned int col=0; col<n; col++)
            {
                double aa = 0.0;
                TimeNodePDE node;
                node.t = t;
                node.i = k;
                for (unsigned int row=0; row<n; row++) aa += x[row]*p.A(node,row,col);
                s1 += aa*x[col];
            }

            // Calculating bettas
            double s2 = 0.0;
            for (unsigned int s=1; s<=k1; s++)
            {
                double ss = 0.0;
                for (unsigned int col=0; col<n; col++)
                {
                    double aa = 0.0;
                    TimeNodePDE node;
                    node.t = t;
                    node.i = k;
                    for (unsigned int row=0; row<n; row++) aa += x[row]*p.B(node,s-1,row,col);
                    ss += aa*x[s*n+col];
                }
                s2 += ss;
            }

            // Calculating gamma
            double s3 = 0.0;
            for (unsigned int col=0; col<n; col++)
            {
                TimeNodePDE node;
                node.t = t;
                node.i = k;
                s3 += x[col]*p.C(node,col);
            }
            s3 *= x[(k1+1)*n];

            double R = 0.0;
            // Calculating alpha
            for (unsigned int i=0; i<n; i++) R += x[i]*x[i];
            // Calculateing bettas
            for (unsigned int s=1; s<=k1; s++)
                for (unsigned int i=0; i<n; i++) R += x[s*n+i];
            // Calculating gamma
            R += x[(k1+1)*n]*x[(k1+1)*n];

            return (s1+s2-s3)/R;
        }
    private:
        ISystemLinearODENonLocalContionsV2 &p;
    };

    CauchyProblemM1stOrderA1 cpa(*this);
    cpa.setDimension(Dimension(h, ec.nmbr, sc.nmbr));
    cpa.cauchyProblem(sc.time, x0, rx, CauchyProblemM1stOrderA1::RK4);
}
