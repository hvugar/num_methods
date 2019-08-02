#include "problem2h_solver.h"
#include <time.h>

SpacePointInfo::SpacePointInfo(unsigned int length)
    : point(SpacePoint(0.0, 0.0)), length(length), vl(nullptr), dx(nullptr), dy(nullptr)
{
    if (length == 0) return;
    vl = new double[length];
    dx = new double[length];
    dy = new double[length];
}

SpacePointInfo::SpacePointInfo(const SpacePoint &point, unsigned int length)
    : point(point), length(length), vl(nullptr), dx(nullptr), dy(nullptr)
{
    if (length == 0) return;
    vl = new double[length];
    dx = new double[length];
    dy = new double[length];
}

SpacePointInfo::SpacePointInfo(const SpacePointInfo& spi)
{
    if (length != 0) { this->~SpacePointInfo(); }

    if (spi.length == 0) return;

    point = spi.point;
    length = spi.length;
    vl = new double[length]; memcpy(vl, spi.vl, sizeof (double)*length);
    dx = new double[length]; memcpy(dx, spi.dx, sizeof (double)*length);
    dy = new double[length]; memcpy(dy, spi.dy, sizeof (double)*length);
}

void SpacePointInfo::clear()
{
    point.x = 0.0;
    point.y = 0.0;
    delete [] dy; dy = nullptr;
    delete [] dx; dx = nullptr;
    delete [] vl; vl = nullptr;
    length = 0;
}

SpacePointInfo& SpacePointInfo::operator=(const SpacePointInfo& other)
{
    if (this == &other) return *this;

    if (length != 0) { this->clear(); }

    if (other.length == 0) return *this;

    length = other.length;
    point = other.point;
    vl = new double[length]; memcpy(vl, other.vl, sizeof (double)*length);
    dx = new double[length]; memcpy(dx, other.dx, sizeof (double)*length);
    dy = new double[length]; memcpy(dy, other.dy, sizeof (double)*length);

    return *this;
}

SpacePointInfo::~SpacePointInfo()
{
    clear();
}

/**************************************************************************************************/

void Problem2HSolver::Main(int argc UNUSED_PARAM, char* argv[] UNUSED_PARAM)
{
    example2();
}

void Problem2HSolver::example1()
{
    Problem2HSolver ps;
    ps.L = 300;
    ps.D = 30;
    ps.setDimensions(Dimension(0.01, 0, 100), Dimension(0.01, 0, 100), Dimension(0.01, 0, static_cast<int>(ps.L+ps.D)));
    ps.setEquationParameters(1.0, 0.01);
    ps.u_list.resize(2*ps.D+1);

    /****************************************************************************************************************/

    ps.Nt = 10;
    ps.times = new TimeNodePDE[ps.Nt];
    ps.times[0] = TimeNodePDE(30,  0.3);
    ps.times[1] = TimeNodePDE(60,  0.6);
    ps.times[2] = TimeNodePDE(90,  0.9);
    ps.times[3] = TimeNodePDE(120, 1.2);
    ps.times[4] = TimeNodePDE(150, 1.5);
    ps.times[5] = TimeNodePDE(180, 1.8);
    ps.times[6] = TimeNodePDE(210, 2.1);
    ps.times[7] = TimeNodePDE(240, 2.4);
    ps.times[8] = TimeNodePDE(270, 2.7);
    ps.times[9] = TimeNodePDE(300, 3.0);

    /****************************************************************************************************************/

    const unsigned int initialPulsesCount = 2;
    InitialPulse *initialPulses = new InitialPulse[initialPulsesCount];
    initialPulses[0] = { SpacePoint(0.25, 0.25), 0.052 };
    initialPulses[1] = { SpacePoint(0.75, 0.75), 0.048 };
    ps.Problem2HWaveEquationIBVP::setInitialConditionMatrix(initialPulses, initialPulsesCount);
    delete [] initialPulses;
    /****************************************************************************************************************/

    const unsigned int Nc = 2;
    const unsigned int No = 3;
    const unsigned int length = 2*Nc*No + 2*(Nc+No);
    ps.setParameters(Nc, No, ps.Problem2HWaveEquationIBVP::spaceDimensionX(), ps.Problem2HWaveEquationIBVP::spaceDimensionY(),
                     ps.Problem2HWaveEquationIBVP::timeDimension());

    // Optimization Vector
    double ox[] = { -0.5412, -0.8412, +0.5745, -0.8259, +0.8482, +0.3751,
                    -0.0084, -0.0075, +0.0086, -0.0035, +0.0022, +0.0031,
                    +0.1418, +0.2914, +0.8724, +0.1008, +0.6332, +0.6028,
                    +0.7608, +0.2631, +0.3725, +0.6045 };

    // Regularization Vector
    double rx[] = { +0.3068, +0.3486, -0.2102, +0.3355, +0.3874, -0.2490,
                    +0.0028, +0.0075, -0.0026, +0.0030, +0.0065, -0.0027,
                    +0.2059, +0.3967, +0.7639, +0.1825, +0.5215, +0.4707,
                    +0.6734, +0.2997, +0.2968, +0.6893};

    ps.setOptimizationVector(ox);
    const double epsilon = 0.2;
    ps.setRegularizationVector(rx, epsilon, epsilon, epsilon, epsilon);
    ps.noise = 0.0;

    for (unsigned int i=0; i<Nc; i++) { ps.vmin[i] = -0.05; ps.vmax[i] = +0.05; }
    ps.r = 0.1;

    //checkGradient3(ps, ox, length);

    //SteepestDescentGradient g;
    ConjugateGradient g;
    g.setFunction(&ps);
    g.setGradient(&ps);
    g.setPrinter(&ps);
    g.setProjection(&ps);
    g.setOptimalityTolerance(0.0);
    g.setFunctionTolerance(0.0);
    g.setStepTolerance(0.0);
    g.setR1MinimizeEpsilon(1.0, 0.01);
    g.setMaxIterationCount(20);
    g.setNormalize(false);
    g.showExitMessage(true);

    IPrinter::printSeperatorLine();

    DoubleVector x(ox, length);
    g.calculate(x);

    puts("Optimization is finished...");

    DoubleVector rv(rx, length);
    //example1_1(ps, rv);
    example1_2(ps, rv);

    puts("Finished");

    delete [] ps.times;
}

void Problem2HSolver::example1_1(Problem2HSolver& ps, DoubleVector &x)
{
    ps.L = 1000;
    ps.setDimensions(Dimension(0.01, 0, 100), Dimension(0.01, 0, 100), Dimension(0.01, 0, static_cast<int>(ps.L+ps.D)));
    ps.setParameters(ps.Nc, ps.No, ps.Problem2HWaveEquationIBVP::spaceDimensionX(),
                     ps.Problem2HWaveEquationIBVP::spaceDimensionY(),
                     ps.Problem2HWaveEquationIBVP::timeDimension());
    ps.setOptimizationVector(x.data());
    ps.setRegularizationVector(x.data(), 0.0, 0.0, 0.0, 0.0);
    ps.r = 0.0;
    ps.noise = 0.0;
    ps.Problem2HWaveEquationIBVP::setWaveDissipation(0.0);
    ps.Problem2HWaveEquationIBVP::save = true;
    ps.Problem2HWaveEquationIBVP::save_u.resize(2061);
    ps.Problem2HWaveEquationIBVP::implicit_calculate_D2V1();
    std::vector<DoubleMatrix> mm(61);
    for (unsigned int i=0; i<=2000; i+=2)
    {
        for (unsigned int j=0; j<=60; j++) mm[j] = ps.save_u[j+i];
        double fx = ps.integral(mm);
        //printf("%4d %10.8f %10.8f\n", i/2, i*0.005, fx);
        printf("%14.10f\n", fx);
    }

    //    double f1 = ps.fx(x);
    //    x[0] = x[1] = x[2] = x[3] = x[4] = x[5] = 0.0;
    //    double f2 = ps.fx(x);

    //    printf("%f %f\n", f1, f2);
}

void Problem2HSolver::example1_2(Problem2HSolver& ps, DoubleVector &x)
{
    const double epsilon = 0.0;
    ps.setRegularizationVector(x.data(), epsilon, epsilon, epsilon, epsilon);
    ps.noise = 0.0;
    ps.r = 0.0;

    double **result = static_cast<double**>(malloc(sizeof(double*)*22));

//    for (unsigned int j=0; j<12; j++)
//    {
//        DoubleVector tx = x;
//        result[j] = static_cast<double*>(malloc(sizeof(double)*201));
//        for (unsigned int i=0; i<=200; i++)
//        {
//            tx[j] = i*0.01 - 1.0;
//            double fx = ps.fx(tx);
//            result[j][i] = fx;
//        }
//        printf("%4d finished.\n", j);
//    }

//    for (unsigned int i=0; i<=200; i++)
//    {
//        double a = i*0.01 - 1.0;
//        printf("%4d %8.4f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f\n", (i-100), a,
//               result[0][i], result[1][i], result[2][i], result[3][i], result[4][i], result[5][i],
//               result[6][i], result[7][i], result[8][i], result[9][i], result[10][i], result[11][i]);
//    }

    for (unsigned int j=12; j<22; j++)
    {
        DoubleVector tx = x;
        result[j] = static_cast<double*>(malloc(sizeof(double)*101));
        result[j][0] = result[j][1] = result[j][2] = result[j][3] = result[j][4] = 0.0;
        for (unsigned int i=5; i<=95; i++)
        {
            tx[j] = i*0.01;
            double fx = ps.fx(tx);
            result[j][i] = fx;
        }
        result[j][96] = result[j][97] = result[j][98] = result[j][99] = result[j][100] = 0.0;
        printf("%4d finished.\n", j);
    }

    for (unsigned int i=0; i<=100; i++)
    {
        double a = i*0.01;
        printf("%4d %8.4f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f\n", i, a,
               result[12][i], result[13][i], result[14][i], result[15][i], result[16][i], result[17][i],
               result[18][i], result[19][i], result[20][i], result[21][i] );
    }
}

void Problem2HSolver::example2()
{
    Problem2HSolver ps;
    ps.L = 300;
    ps.D = 30;
    ps.setDimensions(Dimension(0.01, 0, 100), Dimension(0.01, 0, 100), Dimension(0.01, 0, static_cast<int>(ps.L+ps.D)));
    ps.setEquationParameters(1.0, 0.01);
    ps.u_list.resize(2*ps.D+1);

    /****************************************************************************************************************/

    ps.Nt = 10;
    ps.times = new TimeNodePDE[ps.Nt];
    ps.times[0] = TimeNodePDE(30,  0.3);
    ps.times[1] = TimeNodePDE(60,  0.6);
    ps.times[2] = TimeNodePDE(90,  0.9);
    ps.times[3] = TimeNodePDE(120, 1.2);
    ps.times[4] = TimeNodePDE(150, 1.5);
    ps.times[5] = TimeNodePDE(180, 1.8);
    ps.times[6] = TimeNodePDE(210, 2.1);
    ps.times[7] = TimeNodePDE(240, 2.4);
    ps.times[8] = TimeNodePDE(270, 2.7);
    ps.times[9] = TimeNodePDE(300, 3.0);

    /****************************************************************************************************************/

    const unsigned int initialPulsesCount = 2;
    InitialPulse *initialPulses = new InitialPulse[initialPulsesCount];
    initialPulses[0] = { SpacePoint(0.25, 0.25), 0.052 };
    initialPulses[1] = { SpacePoint(0.75, 0.75), 0.048 };
    ps.Problem2HWaveEquationIBVP::setInitialConditionMatrix(initialPulses, initialPulsesCount);
    delete [] initialPulses;

    /****************************************************************************************************************/

    const unsigned int Nc = 2;
    const unsigned int No = 3;
    const unsigned int length = 2*Nc*No + 2*(Nc+No);
    ps.setParameters(Nc, No, ps.Problem2HWaveEquationIBVP::spaceDimensionX(), ps.Problem2HWaveEquationIBVP::spaceDimensionY(),
                     ps.Problem2HWaveEquationIBVP::timeDimension());

    // Optimization Vector
//    double ox[] = { -0.0412, -0.4182, +0.7455, -0.2589, -0.7521, +0.3848,
//                    +0.4008, +0.2805, +0.3006, +0.8128, -0.2070, +0.3001,
//                    +0.3105, +0.2914, +0.7248, +0.2153, +0.4687, +0.6215,
//                    +0.7524, +0.4135, +0.2015, +0.7549 };
//    ox[0] = ox[1] = ox[2] = ox[3] = ox[4] = ox[5] = 0.0;

    // 0.00009911
    // Regularization Vector
    double rx[] = { -0.5259, -0.1830, +0.2802, -0.4922, -0.1949, +0.3130,
                    +0.0083, +0.0167, +0.0217, +0.0415, -0.0264, +0.0457,
                    +0.3388, +0.2948, +0.8400, +0.0840, +0.6400, +0.4634,
                    +0.6012, +0.3929, +0.3740, +0.6324 };

    // Regularization Vector Example1 and optimal 0.000090
    //double rx[] = { +0.3068, +0.3486, -0.2102, +0.3355, +0.3874, -0.2490,
    //                +0.0028, +0.0075, -0.0026, +0.0030, +0.0065, -0.0027,
    //                +0.2059, +0.3967, +0.7639, +0.1825, +0.5215, +0.4707,
    //                +0.6734, +0.2997, +0.2968, +0.6893};
    //

//    double rx[] = {  0.2013,  0.4131, -0.2228,  0.2543,  0.3740, -0.1758,
//                     -0.0140,  0.0193,  0.0164,  0.0218,  0.0016,  0.0290,
//                     0.3181,  0.3256,  0.8085,  0.1588,  0.6400,  0.4933,
//                     0.6325,  0.2997,  0.3480,  0.7467 };

    ps.setOptimizationVector(rx);
    const double epsilon = 0.0;
    ps.setRegularizationVector(rx, epsilon, epsilon, epsilon, epsilon);
    ps.noise = 0.0;

    for (unsigned int i=0; i<Nc; i++) { ps.vmin[i] = -0.05; ps.vmax[i] = +0.05; }
    ps.r = 0.0;

    //checkGradient3(ps, ox, length);

    //SteepestDescentGradient g;
    ConjugateGradient g;
    //ConstStepGradient g;
    g.setFunction(&ps);
    g.setGradient(&ps);
    g.setPrinter(&ps);
    g.setProjection(&ps);
    g.setOptimalityTolerance(-0.01);
    g.setFunctionTolerance(-0.01);
    g.setStepTolerance(-0.01);
    g.setR1MinimizeEpsilon(1.0, 0.01);
    g.setMaxIterationCount(2000);
    g.setNormalize(false);
    g.showExitMessage(true);

    IPrinter::printSeperatorLine();
    DoubleVector x(rx, length);
    g.calculate(x);

    puts("Optimization is finished...");

    delete [] ps.times;
}

void Problem2HSolver::gradient(const DoubleVector &pv, DoubleVector &g) const
{
    g.clear();
    g.resize(pv.length(), 0.0);
    Problem2HSolver* solver = const_cast<Problem2HSolver*>(this);
    //solver->OptimalParameterFromVector(pv);
    solver->setOptimizationVector(pv.data());
    //solver->Problem2HWaveEquationIBVP::setInitialConditionMatrix(initialPulses, initialPulsesCount);
    gradient_one(pv, g, solver);

    //    for (unsigned int i=0; i<g.length(); i++)
    //    {
    //        if (i<12) continue;
    //        g[i] = 0.0;
    //    }
}

void Problem2HSolver::checkGradient3(const Problem2HSolver &prob, const double *data, unsigned int length)
{
    IPrinter::printSeperatorLine();
    DoubleVector x(data, length);
    printf("ok: "); IPrinter::print(x.mid(0,  5), x.mid(0,  5).length(), 9, 4);
    printf("oz: "); IPrinter::print(x.mid(6, 11), x.mid(6, 11).length(), 9, 4);
    printf("xy: "); IPrinter::print(x.mid(12,21), x.mid(12,21).length(), 9, 4);
    IPrinter::printSeperatorLine();

    //    DoubleVector rv;
    //    prob.equaPrm.RegularParameterToVector(rv);
    //    printf("rk: "); IPrinter::print(rv.mid(0,  3), rv.mid(0,  3).length(), 9, 6);
    //    printf("rz: "); IPrinter::print(rv.mid(4,  7), rv.mid(4,  7).length(), 9, 6);
    //    printf("xy: "); IPrinter::print(rv.mid(8, 15), rv.mid(8, 15).length(), 9, 6);
    //    IPrinter::printSeperatorLine();

    DoubleVector ag(x.length());
    double functional = prob.fx(x);
    printf("Functional: %f\n", functional);
    puts("Calculating gradients....");
    prob.gradient(x, ag);
    puts("Gradients are calculated.");
    //return;

    DoubleVector ng1(x.length(), 0.0);
    DoubleVector ng2(x.length(), 0.0);

    const unsigned int Nc = prob.Nc;
    const unsigned int No = prob.No;
    const unsigned int offset = Nc*No;

    puts("Calculating numerical gradients.... dh=0.01");
    printf("*** Calculating numerical gradients for k...... dh=0.01 ");
    IGradient::Gradient(&prob, 0.01, x, ng1, 0*offset, 1*offset-1);
    printf("Calculated.\n");
    printf("*** Calculating numerical gradients for z...... dh=0.01 ");
    IGradient::Gradient(&prob, 0.01, x, ng1, 1*offset, 2*offset-1);
    printf("Calculated.\n");
    printf("*** Calculating numerical gradients for xi..... dh=0.01 ");
    IGradient::Gradient(&prob, 0.01, x, ng1, 2*offset+0*No, 2*offset+2*No-1);
    printf("Calculated.\n");
    printf("*** Calculating numerical gradients for eta.... dh=0.01 ");
    IGradient::Gradient(&prob, 0.01, x, ng1, 2*offset+2*No, 2*offset+2*No+2*Nc-1);
    printf("Calculated.\n");

    puts("Calculating numerical gradients.... dh=0.001");
    printf("*** Calculating numerical gradients for k...... dh=0.001 ");
    IGradient::Gradient(&prob, 0.001, x, ng2, 0*offset, 1*offset-1);
    printf("Calculated.\n");
    printf("*** Calculating numerical gradients for z...... dh=0.001 ");
    IGradient::Gradient(&prob, 0.001, x, ng2, 1*offset, 2*offset-1);
    printf("Calculated.\n");
    printf("*** Calculating numerical gradients for xi..... dh=0.001 ");
    IGradient::Gradient(&prob, 0.001, x, ng2, 2*offset+0*No, 2*offset+2*No-1);
    printf("Calculated.\n");
    printf("*** Calculating numerical gradients for eta.... dh=0.001 ");
    IGradient::Gradient(&prob, 0.001, x, ng2, 2*offset+2*No, 2*offset+2*No+2*Nc-1);
    printf("Calculated.\n");

    const unsigned int N = 20;
    const unsigned int W = 9;
    const unsigned int P = 6;
    //k------------------------------------------------------//
    {
        IPrinter::printSeperatorLine("k");
        DoubleVector pk0 = x.mid(0, offset-1);
        DoubleVector ak0 = ag.mid(0, offset-1);
        DoubleVector nk1 = ng1.mid(0, offset-1);
        DoubleVector nk2 = ng2.mid(0, offset-1);

        IPrinter::print(pk0,pk0.length(),W,P);
        IPrinter::print(ak0,ak0.length(),W,P); ak0.L2Normalize();
        IPrinter::print(nk1,nk1.length(),W,P); nk1.L2Normalize();
        IPrinter::print(nk2,nk2.length(),W,P); nk2.L2Normalize();
        IPrinter::print(ak0,ak0.length(),W,P);
        IPrinter::print(nk1,nk1.length(),W,P);
        IPrinter::print(nk2,nk2.length(),W,P);
    }
    //z------------------------------------------------------//
    {
        IPrinter::printSeperatorLine("z");
        DoubleVector pz0 = x.mid(offset, 2*offset-1);
        DoubleVector az0 = ag.mid(offset, 2*offset-1);
        DoubleVector nz1 = ng1.mid(offset, 2*offset-1);
        DoubleVector nz2 = ng2.mid(offset, 2*offset-1);

        IPrinter::print(pz0,pz0.length(),W,P);
        IPrinter::print(az0,az0.length(),W,P); az0.L2Normalize();
        IPrinter::print(nz1,nz1.length(),W,P); nz1.L2Normalize();
        IPrinter::print(nz2,nz2.length(),W,P); nz2.L2Normalize();
        IPrinter::print(az0,az0.length(),W,P);
        IPrinter::print(nz1,nz1.length(),W,P);
        IPrinter::print(nz2,nz2.length(),W,P);
    }
    //xi------------------------------------------------------//
    {
        IPrinter::printSeperatorLine("xi");
        DoubleVector pe0 = x.mid(2*offset, 2*offset+2*No-1);
        DoubleVector ae0 = ag.mid(2*offset, 2*offset+2*No-1);
        DoubleVector ne1 = ng1.mid(2*offset, 2*offset+2*No-1);
        DoubleVector ne2 = ng2.mid(2*offset, 2*offset+2*No-1);

        IPrinter::print(pe0,pe0.length(),W,P);
        IPrinter::print(ae0,ae0.length(),W,P); ae0.L2Normalize();
        IPrinter::print(ne1,ne1.length(),W,P); ne1.L2Normalize();
        IPrinter::print(ne2,ne2.length(),W,P); ne2.L2Normalize();
        IPrinter::print(ae0,ae0.length(),W,P);
        IPrinter::print(ne1,ne1.length(),W,P);
        IPrinter::print(ne2,ne2.length(),W,P);
    }

    //eta------------------------------------------------------//
    {
        IPrinter::printSeperatorLine("eta");
        DoubleVector px0 = x.mid(2*offset+2*No, 2*offset+2*No+2*Nc-1);
        DoubleVector ax0 = ag.mid(2*offset+2*No, 2*offset+2*No+2*Nc-1);
        DoubleVector nx1 = ng1.mid(2*offset+2*No, 2*offset+2*No+2*Nc-1);
        DoubleVector nx2 = ng2.mid(2*offset+2*No, 2*offset+2*No+2*Nc-1);

        IPrinter::print(px0,px0.length(),W,P);
        IPrinter::print(ax0,ax0.length(),W,P); ax0.L2Normalize();
        IPrinter::print(nx1,nx1.length(),W,P); nx1.L2Normalize();
        IPrinter::print(nx2,nx2.length(),W,P); nx2.L2Normalize();
        IPrinter::print(ax0,ax0.length(),W,P);
        IPrinter::print(nx1,nx1.length(),W,P);
        IPrinter::print(nx2,nx2.length(),W,P);
        IPrinter::printSeperatorLine();
    }
}

double Problem2HSolver::fx(const DoubleVector &x) const
{
    double SUM = 0.0;
    Problem2HSolver* solver = const_cast<Problem2HSolver*>(this);
    //solver->OptimalParameterFromVector(x);
    solver->setOptimizationVector(x.data());
    //solver->Problem2HWaveEquationIBVP::setInitialConditionMatrix(this->initialPulses, this->initialPulsesCount);
    SUM += fx_one(x, solver);
    return SUM;
}

double Problem2HSolver::fx_one(const DoubleVector &, Problem2HSolver *solver) const
{
    solver->u_list.clear();
    solver->Problem2HWaveEquationIBVP::implicit_calculate_D2V1();

    double sum = integral(u_list);
#ifdef USE_NORM
    sum += fx_norm();
#endif
#ifdef USE_PENALTY
    sum += r*penalty();
#endif
    for (unsigned int i=0; i<u_list.size(); i++)
        solver->u_list[i].clear();
    solver->u_list.clear();
    return sum;
}

double Problem2HSolver::fx_norm() const
{
    double norm = 0.0;

    double norm_k = 0.0;
    double norm_z = 0.0;
    double norm_o = 0.0;
    double norm_c = 0.0;
    for (unsigned int i=0; i<Nc; i++)
    {
        norm_c += (eta[i].x - r_eta[i].x)*(eta[i].x - r_eta[i].x);
        norm_c += (eta[i].y - r_eta[i].y)*(eta[i].y - r_eta[i].y);
        for (unsigned int j=0; j<No; j++)
        {
            if (i==0)
            {
                norm_o += (ksi[j].x - r_ksi[j].x)*(ksi[j].x - r_ksi[j].x);
                norm_o += (ksi[j].y - r_ksi[j].y)*(ksi[j].y - r_ksi[j].y);
            }
            norm_k += (k[i][j] - r_k[i][j])*(k[i][j] - r_k[i][j]);
            norm_z += (z[i][j] - r_z[i][j])*(z[i][j] - r_z[i][j]);
        }
    }

    norm += regEpsilon1*norm_k;
    norm += regEpsilon2*norm_z;
    norm += regEpsilon3*norm_o;
    norm += regEpsilon4*norm_c;

    return norm;
}

double Problem2HSolver::integral(const std::vector<DoubleMatrix> &vu) const
{
    const double ht = Problem2HWaveEquationIBVP::timeDimension().step();

    double sum = 0.0;
    sum += 0.5*integralU(vu[0]);

    for (unsigned int ln=2; ln<=2*(D-1); ln+=2)
    {
        sum += integralU(vu[ln]);
    }

    sum += 0.5*integralU(vu[2*D]);

    return sum*ht;
}

double Problem2HSolver::integralU(const DoubleMatrix &u) const
{
    const unsigned int N = static_cast<unsigned int> ( Problem2HWaveEquationIBVP::spaceDimensionX().size() );
    const unsigned int M = static_cast<unsigned int> ( Problem2HWaveEquationIBVP::spaceDimensionY().size() );
    const double hx = Problem2HWaveEquationIBVP::spaceDimensionX().step();
    const double hy = Problem2HWaveEquationIBVP::spaceDimensionY().step();

    double udiff = 0.0;
    double usum = 0.0;

    udiff = u[0][0]; usum += 0.25 * udiff * udiff;// * mu(0, 0);
    udiff = u[0][N]; usum += 0.25 * udiff * udiff;// * mu(N, 0);
    udiff = u[M][0]; usum += 0.25 * udiff * udiff;// * mu(0, M);
    udiff = u[M][N]; usum += 0.25 * udiff * udiff;// * mu(N, M);

    for (unsigned int n=1; n<=N-1; n++)
    {
        udiff = u[0][n]; usum += 0.5 * udiff * udiff;// * mu(n, 0);
        udiff = u[M][n]; usum += 0.5 * udiff * udiff;// * mu(n, M);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        udiff = u[m][0]; usum += 0.5 * udiff * udiff;// * mu(0, m);
        udiff = u[m][N]; usum += 0.5 * udiff * udiff;// * mu(N, m);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        for (unsigned int n=1; n<=N-1; n++)
        {
            udiff = u[m][n]; usum += udiff * udiff;// * mu(n, m);
        }
    }

    return usum*(hx*hy);
}

void Problem2HSolver::gradient_one(const DoubleVector &, DoubleVector &g, Problem2HSolver* solver) const
{
    solver->u_list.clear();
    solver->Problem2HWaveEquationIBVP::implicit_calculate_D2V1();
    solver->Problem2HConjugateWaveEquationIBVP::implicit_calculate_D2V1();
    solver->u_list.clear();

    unsigned int gi = 0;

    // k
    if (optimizeK)
    {
        //puts("Calculating k gradients...");
        for (unsigned int i=0; i<Nc; i++)
        {
            const SpacePointInfo &eta_spi = eta_info[i];
            for (unsigned int j=0; j<No; j++)
            {
                const SpacePointInfo &ksi_spi = ksi_info[j];

                double zij = z[i][j];

                double grad_Kij = 0.0;
                for (unsigned int s=0; s<Nt; s++)
                {
                    const unsigned int ln = 2*times[s].i;
                    grad_Kij += -(ksi_spi.vl[ln] - zij) * eta_spi.vl[ln];
#ifdef USE_PENALTY
                    grad_Kij += -(ksi_spi.vl[ln] - zij) * 2.0*r*gpi(i,ln)*sgn(g0i(i,ln));
#endif
                }

#ifdef USE_NORM
                grad_Kij += +2.0*regEpsilon1*(k[i][j] - r_k[i][j]);
#endif
                g[gi++] += grad_Kij;
            }
        }
    }
    else
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            for (unsigned int j=0; j<No; j++)
            {
                g[gi++] = 0.0;
            }
        }
    }

    // z
    if (optimizeZ)
    {
        //puts("Calculating z gradients...");
        for (unsigned int i=0; i<Nc; i++)
        {
            const SpacePointInfo &eta_spi = eta_info[i];
            for (unsigned int j=0; j<No; j++)
            {
                double grad_Zij = 0.0;

                double kij = k[i][j];
                for (unsigned int s=0; s<Nt; s++)
                {
                    const unsigned int ln = 2*times[s].i;
                    grad_Zij += kij * eta_spi.vl[ln];
#ifdef USE_PENALTY
                    grad_Zij += kij * 2.0*r*gpi(i,ln)*sgn(g0i(i,ln));
#endif
                }
#ifdef USE_NORM
                grad_Zij += +2.0*regEpsilon2*(z[i][j] - r_z[i][j]);
#endif
                g[gi++] = grad_Zij;
            }
        }
    }
    else
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            for (unsigned int j=0; j<No; j++)
            {
                g[gi++] = 0.0;
            }
        }
    }

    // xi
    if (optimizeO)
    {
        //puts("Calculating o gradients...");
        for (unsigned int j=0; j<No; j++)
        {
            const SpacePointInfo &ksi_spi = ksi_info[j];

            double gradXijX = 0.0;
            double gradXijY = 0.0;

            for (unsigned int s=0; s<Nt; s++)
            {
                const unsigned int ln = 2*times[s].i;
                for (unsigned int i=0; i<Nc; i++)
                {
                    const SpacePointInfo &eta_spi = eta_info[i];

                    double kij = k[i][j];

                    gradXijX += -kij * ksi_spi.dx[ln] * eta_spi.vl[ln];
                    gradXijY += -kij * ksi_spi.dy[ln] * eta_spi.vl[ln];
#ifdef USE_PENALTY
                    gradXijX += -kij * ksi_spi.dx[ln] * 2.0*r*gpi(i,ln)*sgn(g0i(i,ln));
                    gradXijY += -kij * ksi_spi.dy[ln] * 2.0*r*gpi(i,ln)*sgn(g0i(i,ln));
#endif
                }
            }
#ifdef USE_NORM
            gradXijX += 2.0*regEpsilon3*(ksi[j].x - r_ksi[j].x);
            gradXijY += 2.0*regEpsilon3*(ksi[j].y - r_ksi[j].y);
#endif

            g[gi++] += gradXijX;
            g[gi++] += gradXijY;
        }
    }
    else
    {
        for (unsigned int j=0; j<No; j++)
        {
            g[gi++] = 0.0;
            g[gi++] = 0.0;
        }
    }

    // eta
    if (optimizeC)
    {
        //puts("Calculating c gradients...");
        for (unsigned int i=0; i<Nc; i++)
        {
            const SpacePointInfo &eta_spi = eta_info[i];

            double gradEtaiX = 0.0;
            double gradEtaiY = 0.0;

            for (unsigned int s=0; s<Nt; s++)
            {
                const unsigned int ln = 2*times[s].i;
                for (unsigned int j=0; j<No; j++)
                {
                    const SpacePointInfo &ksi_spi = ksi_info[j];

                    double kij = k[i][j];
                    double zij = z[i][j];

                    gradEtaiX += -kij * eta_spi.dx[ln] * (ksi_spi.vl[ln] - zij);
                    gradEtaiY += -kij * eta_spi.dy[ln] * (ksi_spi.vl[ln] - zij);
                }
            }
#ifdef USE_NORM
            gradEtaiX += 2.0*regEpsilon4*(eta[i].x - r_eta[i].x);
            gradEtaiY += 2.0*regEpsilon4*(eta[i].y - r_eta[i].y);
#endif
            g[gi++] += gradEtaiX;
            g[gi++] += gradEtaiY;
        }
    }
    else
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            g[gi++] = 0.0;
            g[gi++] = 0.0;
        }
    }
}

void Problem2HSolver::setDimensions(const Dimension &dimensionX, const Dimension &dimensionY, const Dimension &timeDimension)
{
    Problem2HWaveEquationIBVP::setTimeDimension(timeDimension);
    Problem2HWaveEquationIBVP::setSpaceDimensions(dimensionX, dimensionY);

    Problem2HConjugateWaveEquationIBVP::setTimeDimension(timeDimension);
    Problem2HConjugateWaveEquationIBVP::setSpaceDimensions(dimensionX, dimensionY);
}

void Problem2HSolver::setEquationParameters(double waveSpeed, double waveDissipation)
{
    Problem2HWaveEquationIBVP::setWaveSpeed(waveSpeed);
    Problem2HWaveEquationIBVP::setWaveDissipation(waveDissipation);

    Problem2HConjugateWaveEquationIBVP::setWaveSpeed(waveSpeed);
    Problem2HConjugateWaveEquationIBVP::setWaveDissipation(waveDissipation);
}

void Problem2HSolver::project(DoubleVector &x) const
{
    unsigned int start = 2*Nc*No;
    unsigned int end = 2*Nc*No + 2*No + 2*Nc - 1;

    //    for (unsigned int index = start; index <= end; index++)
    //    {
    //        if (x[index] <= 0.05) x[index] = 0.05;
    //        if (x[index] >= 0.95) x[index] = 0.95;
    //    }


    if (x[start+0] <= 0.12) x[start+0] = 0.12; if (x[start+0] >= 0.32) x[start+0] = 0.32;
    if (x[start+1] <= 0.28) x[start+1] = 0.28; if (x[start+1] >= 0.48) x[start+1] = 0.48;
    if (x[start+2] <= 0.71) x[start+2] = 0.71; if (x[start+2] >= 0.89) x[start+2] = 0.89;
    if (x[start+3] <= 0.05) x[start+3] = 0.05; if (x[start+3] >= 0.22) x[start+3] = 0.22;
    if (x[start+4] <= 0.44) x[start+4] = 0.44; if (x[start+4] >= 0.64) x[start+4] = 0.64;
    if (x[start+5] <= 0.45) x[start+5] = 0.45; if (x[start+5] >= 0.63) x[start+5] = 0.63;
    if (x[start+6] <= 0.58) x[start+6] = 0.58; if (x[start+6] >= 0.78) x[start+6] = 0.78;
    if (x[start+7] <= 0.25) x[start+7] = 0.25; if (x[start+7] >= 0.42) x[start+7] = 0.42;
    if (x[start+8] <= 0.19) x[start+8] = 0.19; if (x[start+8] >= 0.39) x[start+8] = 0.39;
    if (x[start+9] <= 0.59) x[start+9] = 0.59; if (x[start+9] >= 0.78) x[start+9] = 0.78;


}

void Problem2HSolver::project(DoubleVector &, unsigned int) {}

void Problem2HSolver::print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientMethod::MethodResult result) const
{
    Problem2HSolver* solver = const_cast<Problem2HSolver*>(this);
    solver->setOptimizationVector(x.data());
    solver->u_list.clear();
    solver->Problem2HWaveEquationIBVP::implicit_calculate_D2V1();
    double integr = solver->integral(u_list);
    double nrm = solver->fx_norm();
    double plty = solver->penalty();
    for (unsigned int i=0; i<u_list.size(); i++) solver->u_list[i].clear(); solver->u_list.clear();

    printf("I[%3d] %.8f %7.4f %.8f %.8f %.8f  ", iteration, f, alpha, integr, nrm, plty);
    IPrinter::print(x, x.length(), 7, 4);

    //printf("k: "); IPrinter::print(x.mid(0,  5), x.mid(0,  5).length(), 9, 4);
    //printf("z: "); IPrinter::print(x.mid(6, 11), x.mid(6, 11).length(), 9, 4);
    //printf("o: "); IPrinter::print(x.mid(12,21), x.mid(12,21).length(), 9, 4);
    //IPrinter::printSeperatorLine();

    //    if (iteration == 10)
    //    if (fabs(alpha) <= DBL_EPSILON)
    //    if (iteration > 0 && iteration % 9 == 0)
    //    {
    //        Problem2HSolver *ps = const_cast<Problem2HSolver*>(this);
    //        if (ps->optimizeK) { ps->optimizeZ = true; ps->optimizeK = false; } else
    //        if (ps->optimizeZ) { ps->optimizeK = true; ps->optimizeZ = false; }
    //        //ps->regEpsilon1 = ps->regEpsilon2 = ps->regEpsilon3 = ps->regEpsilon4 = 1.0;
    //    }
}

auto Problem2HCommon::penalty() const -> double
{
    double pnlt = 0.0;

    for (unsigned int s=0; s<Nt; s++)
    {
        for (unsigned int i=0; i<Nc; i++)
        {
            double _gpi_s = gpi(i, 2*times[s].i);
            pnlt += _gpi_s*_gpi_s;
        }
    }

    return pnlt;
}

auto Problem2HCommon::gpi(unsigned int i, unsigned int ln) const -> double
{
    double gpi_ln = fabs( g0i(i, ln) ) - (vmax[i]-vmin[i])/2.0;
    return gpi_ln > 0.0 ? gpi_ln : 0.0;
}

auto Problem2HCommon::g0i(unsigned int i, unsigned int ln) const -> double
{
    double vi = 0.0;

    for (unsigned int j=0; j<No; j++) { vi += k[i][j] * ( ksi_info[j].vl[ln] - z[i][j] ); }

    return (vmax[i]+vmin[i])/2.0 - vi;
}

/*********************************** Problem2HCommon********** ***********************************************************/

Problem2HCommon::~Problem2HCommon() {}

void Problem2HCommon::setParameters(unsigned int Nc, unsigned int No, const Dimension &dimensionX,
                                    const Dimension &dimensionY, const Dimension &timeDimension)
{
    const unsigned int length = 2*static_cast<unsigned int>(timeDimension.size())+1;

    k.clear(); k.resize(Nc, No, 0.0);
    z.clear(); z.resize(Nc, No, 0.0);

    if (ksi != nullptr) { delete [] ksi; } ksi = new SpacePoint[No];
    if (eta != nullptr) { delete [] eta; } eta = new SpacePoint[Nc];

    if (ksi_info != nullptr)
    {
        for (unsigned int j=0; j<this->No; j++) ksi_info[j].deltaGrid.cleanGrid();
        delete [] ksi_info;
    }
    ksi_info = new SpacePointInfo[No];
    for (unsigned int j=0; j<No; j++)
    {
        ksi_info[j] = SpacePointInfo(length);
        ksi_info[j].deltaGrid.initGrid(dimensionX, dimensionY);
    }

    if (eta_info != nullptr)
    {
        for (unsigned int i=0; i<this->Nc; i++) eta_info[i].deltaGrid.cleanGrid();
        delete [] eta_info;
    }
    eta_info = new SpacePointInfo[Nc];
    for (unsigned int i=0; i<Nc; i++)
    {
        eta_info[i] = SpacePointInfo(length);
        eta_info[i].deltaGrid.initGrid(dimensionX, dimensionY);
    }

    this->Nc = Nc;
    this->No = No;

    r_k.clear(); r_k.resize(Nc, No, 0.0);
    r_z.clear(); r_z.resize(Nc, No, 0.0);
    if (r_ksi != nullptr) { delete [] r_ksi; } r_ksi = new SpacePoint[No];
    if (r_eta != nullptr) { delete [] r_eta; } r_eta = new SpacePoint[Nc];

    if (vmin != nullptr) { delete [] vmin; } vmin = new double[Nc];
    if (vmax != nullptr) { delete [] vmax; } vmax = new double[Nc];
    for (unsigned int i=0; i<Nc; i++) { vmin[i] = 0.0; vmax[i] = 0.0; }
}

void Problem2HCommon::setOptimizationVector(const double *x)
{
    unsigned int index = 0;
    unsigned int offset = Nc*No;

    for (unsigned int i=0; i<Nc; i++)
    {
        for (unsigned int j=0; j<No; j++)
        {
            k[i][j] = x[index];
            z[i][j] = x[index+offset];
            index++;
        }
    }

    index = 2*Nc*No;

    for (unsigned int j=0; j<No; j++)
    {
        ksi[j].x = x[index];
        ksi[j].y = x[index+1];
        index+=2;
        ksi_info[j].point = ksi[j];
        ksi_info[j].deltaGrid.distributeGauss(ksi[j], 1, 1);

    }
    for (unsigned int i=0; i<Nc; i++)
    {
        eta[i].x = x[index];
        eta[i].y = x[index+1];
        index+=2;
        eta_info[i].point = eta[i];
        eta_info[i].deltaGrid.distributeGauss(eta[i], 1, 1);
    }
}

void Problem2HCommon::setRegularizationVector(const double *x, double epsilon1, double epsilon2, double epsilon3, double epsilon4)
{
    this->regEpsilon1 = epsilon1;
    this->regEpsilon2 = epsilon2;
    this->regEpsilon3 = epsilon3;
    this->regEpsilon4 = epsilon4;

    unsigned int index = 0;
    unsigned int offset = Nc*No;

    for (unsigned int i=0; i<Nc; i++)
    {
        for (unsigned int j=0; j<No; j++)
        {
            r_k[i][j] = x[index];
            r_z[i][j] = x[index+offset];
            index++;
        }
    }
    index = 2*Nc*No;

    for (unsigned int j=0; j<No; j++) { r_ksi[j].x = x[index]; r_ksi[j].y = x[index+1]; index+=2; }
    for (unsigned int i=0; i<Nc; i++) { r_eta[i].x = x[index]; r_eta[i].y = x[index+1]; index+=2; }
}

/*********************************** Problem2HWaveEquationIBVP ***********************************************************/

double Problem2HWaveEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition condition) const
{
    if (condition == InitialCondition::FirstDerivative)
        return f_initialMatrix[static_cast<uint32_t>(sn.j)][static_cast<uint32_t>(sn.i)];
    else
        return 0.0;
}

double Problem2HWaveEquationIBVP::boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryConditionPDE &condition) const
{
    condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet);
    return 0.0;
}

double Problem2HWaveEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &) const
{
    if (f_return_zero) return 0.0;
    return f_crLayerMatrix[static_cast<unsigned int>(sn.j)][static_cast<unsigned int>(sn.i)];
}

void Problem2HWaveEquationIBVP::setSpaceDimensions(const Dimension &dimensionX, const Dimension &dimensionY)
{
    IWaveEquationIBVP::setSpaceDimensions(dimensionX, dimensionY);
    f_initialMatrix.resize(static_cast<unsigned int>(dimensionY.size())+1, static_cast<unsigned int>(dimensionX.size())+1, 0.0);
    f_crLayerMatrix.resize(static_cast<unsigned int>(dimensionY.size())+1, static_cast<unsigned int>(dimensionX.size())+1, 0.0);
}

void Problem2HWaveEquationIBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    const_cast<Problem2HWaveEquationIBVP*>(this)->layerInfoPrepareLayerMatrix(u, tn);
    //if (tn.i == 600) { layerInfoSave2TextFile(u, tn); }

    //    if (tn.i == 0 || tn.i == 1 || tn.i == 2 || tn.i == 3)
    //    {
    //        IPrinter::printMatrix(u);
    //        IPrinter::printSeperatorLine();
    //    }
}

void Problem2HWaveEquationIBVP::layerInfoSave2TextFile(const DoubleMatrix &u UNUSED_PARAM, const TimeNodePDE &tn) const
{
    static double MIN = +100000.0;
    static double MAX = -100000.0;

    std::string txt_number = std::to_string(tn.i);
    std::string filename = std::string("data/problem2H/f/txt/f_") +
            std::string(4 - txt_number.length(), '0') + txt_number + std::string(".txt");
    IPrinter::print(u, filename.c_str());
    if (MIN > u.min()) MIN = u.min();
    if (MAX < u.max()) MAX = u.max();
    printf("Forward: %4d %0.3f %10.8f %10.8f %10.8f %10.8f %4d %4d\n", tn.i, tn.t, u.min(), u.max(), MIN, MAX, 0, 0);
}

void Problem2HWaveEquationIBVP::setInitialConditionMatrix(InitialPulse *initialPulses, unsigned int initialPulsesCount)
{
    if (this->initialPulses == initialPulses) return;

    if (this->initialPulses != nullptr)
    {
        delete [] this->initialPulses;
        this->initialPulses = nullptr;
    }

    this->initialPulsesCount = initialPulsesCount;
    this->initialPulses = new InitialPulse[initialPulsesCount];

    for (unsigned int s=0; s<initialPulsesCount; s++)
    {
        this->initialPulses[s].point = initialPulses[s].point;
        this->initialPulses[s].blow = initialPulses[s].blow;
    }

    f_initialMatrix.reset();

    for (unsigned int s=0; s<initialPulsesCount; s++)
    {
        const InitialPulse &initialPulse = this->initialPulses[s];

        DeltaGrid2D deltaGrid;
        deltaGrid.initGrid(_spaceDimensionX, _spaceDimensionY);
        deltaGrid.distributeGauss(initialPulse.point, 5, 5);
        const unsigned minX = deltaGrid.minX();
        const unsigned maxX = deltaGrid.maxX();
        const unsigned minY = deltaGrid.minY();
        const unsigned maxY = deltaGrid.maxY();

        for (unsigned int m=minY; m<=maxY; m++)
        {
            for (unsigned int n=minX; n<=maxX; n++)
            {
                f_initialMatrix[m][n] += initialPulse.blow * deltaGrid.weight(n, m);
            }
        }

        deltaGrid.cleanGrid();
    }
}

void Problem2HWaveEquationIBVP::clrInitialConditionMatrix()
{
    f_initialMatrix.clear();
}

void Problem2HWaveEquationIBVP::layerInfoPrepareLayerMatrix(const DoubleMatrix &u, const TimeNodePDE& tn)
{
    //if (save) { save_u[tn.i] = u; }

    if (tn.i >= 2*L) u_list[tn.i-2*L] = u;

    const double ht = timeDimension().step();

    for (unsigned int j=0; j<No; j++)
    {
        double u_vl, u_dx, u_dy;
        u_vl = ksi_info[j].deltaGrid.lumpPointGauss(u, u_dx, u_dy);
        ksi_info[j].vl[tn.i] = u_vl * (1.0 + noise * ( rand()%2 == 0 ? -1.0 : +1.0));
        ksi_info[j].dx[tn.i] = u_dx;
        ksi_info[j].dy[tn.i] = u_dy;
    }

    f_crLayerMatrix.reset();

    f_return_zero = false;

    for (unsigned int s=0; s<Nt; s++)
    {
        double wt = 0.0;
        if (tn.i == 2*times[s].i) { wt = 2.0/ht; } else { continue; }

        //if (f_return_zero) { f_crLayerMatrix.reset(); f_return_zero = false; }

        double* v = new double[Nc];

        for (unsigned int i=0; i<Nc; i++)
        {
            v[i] = 0.0;
            for (unsigned int j=0; j<No; j++)
            {
                v[i] += k.at(i,j) * (ksi_info[j].vl[tn.i]-z[i][j]);
            }
        }

        for (unsigned int i=0; i<Nc; i++)
        {
            const unsigned minX = eta_info[i].deltaGrid.minX();
            const unsigned maxX = eta_info[i].deltaGrid.maxX();
            const unsigned minY = eta_info[i].deltaGrid.minY();
            const unsigned maxY = eta_info[i].deltaGrid.maxY();

            for (unsigned int m=minY; m<=maxY; m++)
            {
                for (unsigned int n=minX; n<=maxX; n++)
                {
                    if (i==0) f_crLayerMatrix[m][n] = 0.0;
                    f_crLayerMatrix[m][n] += v[i] * eta_info[i].deltaGrid.weight(n, m) * wt;
                }
            }
        }

        delete [] v;
    }
}

/*********************************** Problem2HWaveEquationIBVP ***********************************************************/

/*********************************** Problem2HConjugateWaveEquationIBVP***************************************************/

double Problem2HConjugateWaveEquationIBVP::final(const SpaceNodePDE &, FinalCondition) const
{
    return 0.0;
}

double Problem2HConjugateWaveEquationIBVP::boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryConditionPDE &condition) const
{
    condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet);
    return 0.0;
}

double Problem2HConjugateWaveEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &) const
{
    return b_crLayerMatrix[static_cast<unsigned int>(sn.j)][static_cast<unsigned int>(sn.i)];
}

void Problem2HConjugateWaveEquationIBVP::setSpaceDimensions(const Dimension &dimensionX, const Dimension &dimensionY)
{
    IFinalWaveEquationIBVP::setSpaceDimensions(dimensionX, dimensionY);
    //b_initialMatrix.resize(static_cast<unsigned int>(dimensionY.size())+1, static_cast<unsigned int>(dimensionX.size())+1, 0.0);
    b_crLayerMatrix.resize(static_cast<unsigned int>(dimensionY.size())+1, static_cast<unsigned int>(dimensionX.size())+1, 0.0);
}

void Problem2HConjugateWaveEquationIBVP::layerInfo(const DoubleMatrix &p, const TimeNodePDE &tn) const
{
    const_cast<Problem2HConjugateWaveEquationIBVP*>(this)->layerInfoPrepareLayerMatrix(p, tn);
}

void Problem2HConjugateWaveEquationIBVP::layerInfoPrepareLayerMatrix(const DoubleMatrix &p UNUSED_PARAM, const TimeNodePDE &tn)
{
    const double ht = timeDimension().step();

    const unsigned int N = static_cast<const unsigned int> ( spaceDimensionX().size() );
    const unsigned int M = static_cast<const unsigned int> ( spaceDimensionY().size() );

    for (unsigned int i=0; i<Nc; i++)
    {
        double p_vl, p_dx, p_dy;
        p_vl = eta_info[i].deltaGrid.lumpPointGauss(p, p_dx, p_dy);
        eta_info[i].vl[tn.i] = p_vl;
        eta_info[i].dx[tn.i] = p_dx;
        eta_info[i].dy[tn.i] = p_dy;
    }

    if (tn.i >= 600)
    {
        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                const DoubleMatrix &u = u_list[tn.i-600];
                b_crLayerMatrix[m][n] = -2.0*u[m][n];
            }
        }
    }
    else
    {
        b_crLayerMatrix.reset();
    }

    for (unsigned int s=0; s<Nt; s++)
    {
        double wt = 0.0;
        if (tn.i == 2*times[s].i) { wt = 2.0/ht; } else { continue; }

        double* w = new double[No];

        for (unsigned int j=0; j<No; j++)
        {
            w[j] = 0.0;
            for (unsigned int i=0; i<Nc; i++)
            {
                w[j] += k[i][j] * ( eta_info[i].vl[tn.i] + 2.0*r*gpi(i,2*times[s].i)*sgn(g0i(i,2*times[s].i)) );
            }
        }

        for (unsigned int j=0; j<No; j++)
        {
            const unsigned minX = ksi_info[j].deltaGrid.minX();
            const unsigned maxX = ksi_info[j].deltaGrid.maxX();
            const unsigned minY = ksi_info[j].deltaGrid.minY();
            const unsigned maxY = ksi_info[j].deltaGrid.maxY();

            for (unsigned int m=minY; m<=maxY; m++)
            {
                for (unsigned int n=minX; n<=maxX; n++)
                {
                    if (j==0) b_crLayerMatrix[m][n] = 0.0;
                    b_crLayerMatrix[m][n] += w[j] * ksi_info[j].deltaGrid.weight(n, m) * wt;
                }
            }
        }

        delete [] w;
    }
}

void Problem2HConjugateWaveEquationIBVP::layerInfoSave2TextFile(const DoubleMatrix &u, const TimeNodePDE & tn) const {}

/*********************************** Problem2HConjugateWaveEquationIBVP***************************************************/

//void Problem2HSolver::OptimalParameterFromVector(const DoubleVector &x)
//{
//    unsigned int index = 0;

//    k.clear();   k.resize(Nc, No);
//    z.clear();   z.resize(Nc, No);

//    if (ksi != nullptr) delete [] ksi;
//    if (eta != nullptr) delete [] eta;

//    ksi = new SpacePoint[No];
//    eta = new SpacePoint[Nc];

//    for (unsigned int i=0; i<Nc; i++)
//    {
//        for (unsigned int j=0; j<No; j++)
//        {
//            k[i][j] = x[index++];
//        }
//    }

//    for (unsigned int i=0; i<Nc; i++)
//    {
//        for (unsigned int j=0; j<No; j++)
//        {
//            z[i][j] = x[index++];
//        }
//    }

//    for (unsigned int j=0; j<No; j++)
//    {
//        ksi[j].x = x[index++];
//        ksi[j].y = x[index++];
//    }

//    for (unsigned int i=0; i<Nc; i++)
//    {
//        eta[i].x = x[index++];
//        eta[i].y = x[index++];
//    }
//}

//void Problem2HSolver::OptimalParameterToVector(DoubleVector &x) const
//{
//    x.clear();
//    x.resize(2*Nc*No+2*No+2*Nc);

//    unsigned int index = 0;
//    for (unsigned int i=0; i<Nc; i++)
//    {
//        for (unsigned int j=0; j<No; j++)
//        {
//            x[index++] = k[i][j];
//        }
//    }

//    for (unsigned int i=0; i<Nc; i++)
//    {
//        for (unsigned int j=0; j<No; j++)
//        {
//            x[index++] = z[i][j];
//        }
//    }

//    for (unsigned int j=0; j<No; j++)
//    {
//        x[index++] = ksi[j].x;
//        x[index++] = ksi[j].y;
//    }

//    for (unsigned int i=0; i<Nc; i++)
//    {
//        x[index++] = eta[i].x;
//        x[index++] = eta[i].y;
//    }
//}
