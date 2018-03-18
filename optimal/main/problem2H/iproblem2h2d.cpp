#include "iproblem2h2d.h"

#include "iproblem2hforward2d.h"
#include "iproblem2hbackward2d.h"

void IProblem2H2D::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    IProblem2HForward2D p;
    p.addSpaceDimension(Dimension(0.01, 0, 100));
    p.addSpaceDimension(Dimension(0.01, 0, 100));
    p.setTimeDimension(Dimension(0.001, 0, 1000));

    IProblem2H2D::Parameter prm;
    prm.Ns = 2;
    prm.q.resize(prm.Ns);
    prm.q[0] = 1000.1; prm.q[1] = 2000.2;
    prm.theta.resize(prm.Ns);
    prm.Nc = 2;
    prm.No = 2;

    prm.xi.resize(prm.No);
    prm.xi[0].x = 0.442; prm.xi[0].y = 0.443;
    prm.xi[1].x = 0.603; prm.xi[1].y = 0.608;

    prm.eta.resize(prm.Nc);
    prm.eta[0].x = 0.312; prm.eta[0].y = 0.723;
    prm.eta[1].x = 0.756; prm.eta[1].y = 0.344;

    prm.k.resize(prm.Nc, prm.No, 0.1);
    prm.z.resize(prm.Nc, prm.No, 0.0);

    prm.a = 1.0;
    prm.lambda = 0.01;

    p.mParameter = prm;

    DoubleMatrix u;
    p.calculateMVD(u);
}

void IProblem2H2D::distributeDelta(const SpacePoint &pt, std::vector<IProblem2H2D::ExtendedSpacePointNode> &nodes, unsigned int id,
                                   const Dimension &xd, const Dimension &yd)
{
    double hx = xd.step();
    double hy = yd.step();
    unsigned int Nx = xd.sizeN();
    unsigned int Ny = yd.sizeN();

    double sigmaX = hx;
    double sigmaY = hy;

    unsigned int rx = (unsigned int)(round(pt.x*Nx));
    unsigned int ry = (unsigned int)(round(pt.y*Ny));

    unsigned int k=3;

    double sumX = 0.0;
    for (unsigned int n=rx-k; n<=rx+k; n++)
    {
        sumX += exp(-((n*hx-pt.x)*(n*hx-pt.x))/(2.0*sigmaX*sigmaX));
    }
    sumX *= hx;

    double sumY = 0.0;
    for (unsigned int m=ry-k; m<=ry+k; m++)
    {
        sumY += exp(-((m*hy-pt.y)*(m*hy-pt.y))/(2.0*sigmaY*sigmaY));
    }
    sumY *= hy;

    double sigma = (sumX*sumY) / (2.0*M_PI);
    double factor = 1.0/((2.0*M_PI)*sigma);

    for (unsigned int n=rx-k; n<=rx+k; n++)
    {
        for (unsigned int m=ry-k; m<=ry+k; m++)
        {
            IProblem2H2D::ExtendedSpacePointNode node;
            node.i = n; node.x = n*hx; node.j = m; node.y = m*hy; node.pt = pt; node.id = id;
            node.w = factor*exp(-0.5*(((node.x-pt.x)*(node.x-pt.x))/(sigmaX*sigmaX)+((node.y-pt.y)*(node.y-pt.y))/(sigmaY*sigmaY)));
            nodes.push_back(node);
        }
    }
}
