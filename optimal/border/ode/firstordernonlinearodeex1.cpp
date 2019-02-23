#include "firstordernonlinearodeex1.h"

void FirstOrderNonLinearODErEx1::Main(int argc, char **argv)
{
    FirstOrderNonLinearODErEx1 nnl;

    unsigned int N = 100;
    unsigned int M = nnl.count();
    double h = 0.01;
    nnl.setDimension(Dimension(h, 0, static_cast<int>(N)));

    DoubleVector y0(M);
    for (unsigned int r=0; r<M; r++)
    {
        y0[r] = nnl.x(PointNodeODE(0.0, 0), r+1);
    }

    std::vector<DoubleVector> x;
    nnl.cauchyProblem(0.0, y0, x, OdeSolverMethod::RK4);

    for (unsigned int m=0; m<M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            if (n%(N/10)==0) printf("%14.10f ", x[n][m]);
        }
        puts("");
    }
    IPrinter::printSeperatorLine();

    for (unsigned int m=0; m<M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            if (n%(N/10)==0) printf("%14.10f ", nnl.x(n*h, m+1));
        }
        puts("");
    }
    IPrinter::printSeperatorLine();

    double norm[] = {0.0, 0.0, 0.0};
    for (unsigned int m=0; m<M; m++)
    {
        norm[m] = 0.0;
        for (unsigned int n=0; n<=N; n++)
        {
            norm[m] += (x[n][m]-nnl.x(n*h, m+1))*(x[n][m]-nnl.x(n*h, m+1));
        }
        norm[m] = sqrt(norm[m]);
    }
    printf("Norms: %.10f %.10f %.10f\n", norm[0], norm[1], norm[2]);
}

unsigned int FirstOrderNonLinearODErEx1::count() const
{
    return 3;
}

double FirstOrderNonLinearODErEx1::x(const PointNodeODE &node, unsigned int r UNUSED_PARAM) const
{
    double t =  node.x;

    if (r == 1) return t*t;
    if (r == 2) return t*t+t;
    if (r == 3) return t;

    throw std::exception();
}

double FirstOrderNonLinearODErEx1::f(const PointNodeODE &node, const DoubleVector &y, unsigned int r) const
{
    if (r == 1) return 3.0*y[0] + 2.0*y[1] - 1.0*y[2] - ( 3.0*x(node,1) + 2.0*x(node,2) - 1.0*x(node,3) ) + 2.0*node.x;
    if (r == 2) return 5.0*y[0] - 3.0*y[1] + 2.0*y[2] - ( 5.0*x(node,1) - 3.0*x(node,2) + 2.0*x(node,3) ) + 2.0*node.x+1.0;
    if (r == 3) return 2.0*y[0] + 8.0*y[1] - 5.0*y[2] - ( 2.0*x(node,1) + 8.0*x(node,2) - 5.0*x(node,3) ) + 1.0;

    throw std::exception();
}

double FirstOrderNonLinearODEEx2::f(double x, double y, unsigned int k) const
{
    double sigma = 0.01;
//    return  3.0*x*x;// + 4.0*(1.0/(sqrt(2.0*M_PI*sigma*sigma))) * exp(((x-0.5)*(x-0.5))/(-2.0*sigma*sigma));
//    return  4.0*x*x*x*exp(pow(x,4.0));// + 4.0*(1.0/(sqrt(2.0*M_PI*sigma*sigma))) * exp(((x-0.8)*(x-0.8))/(-2.0*sigma*sigma));

    double w = 0.0;
    if (k==80) w = 100.0;
    return  4.0*x*x*x*exp(pow(x,4.0)) + 4.0*w;
}
