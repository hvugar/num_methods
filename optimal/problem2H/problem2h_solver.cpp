#include "problem2h_solver.h"

void Problem2HSolver::Main(int argc UNUSED_PARAM, char* argv[] UNUSED_PARAM)
{

    Problem2HSolver ps;
    ps.L = 300;
    ps.D = 30;
    ps.setDimensions(Dimension(0.01, 0, 100), Dimension(0.01, 0, 100), Dimension(0.01, 0, static_cast<int>(ps.L+ps.D)));
    ps.setEquationParameters(1.0, 0.0);
    ps.u_list.resize(61);

    ps.Nt = 10;
    ps.times.resize(ps.Nt);
    ps.times[0] = TimeNodePDE(30, 0.3);
    ps.times[1] = TimeNodePDE(60, 0.6);
    ps.times[2] = TimeNodePDE(90, 0.9);
    ps.times[3] = TimeNodePDE(120, 1.2);
    ps.times[4] = TimeNodePDE(150, 1.5);
    ps.times[5] = TimeNodePDE(180, 1.8);
    ps.times[6] = TimeNodePDE(210, 2.1);
    ps.times[7] = TimeNodePDE(240, 2.4);
    ps.times[8] = TimeNodePDE(270, 2.7);
    ps.times[9] = TimeNodePDE(300, 3.0);

    const unsigned int initialPulsesCount = 2;
    InitialPulse *initialPulses = new InitialPulse[initialPulsesCount];
    initialPulses[0] = { SpacePoint(0.25, 0.25), 0.05 };
    initialPulses[1] = { SpacePoint(0.75, 0.75), 0.05 };
    ps.Problem2HWaveEquationIBVP::setInitialConditionMatrix(initialPulses, initialPulsesCount);
    //delete [] initialPulses;

    std::vector<SpacePoint> eta(2);
    eta[0] = SpacePoint(0.65, 0.34);
    eta[1] = SpacePoint(0.25, 0.75);

    std::vector<SpacePoint> ksi(2);
    ksi[0] = SpacePoint(0.22, 0.54);
    ksi[1] = SpacePoint(0.82, 0.27);

    DoubleMatrix k(2, 2, -0.01);
    DoubleMatrix z(2, 2, +0.00);

    ps.setParameterCounts(2, 2, ps.Problem2HWaveEquationIBVP::spaceDimensionX(),
                          ps.Problem2HWaveEquationIBVP::spaceDimensionY(),
                          ps.Problem2HWaveEquationIBVP::timeDimension());
    ps.setOptimizedParameters(k, z, eta, ksi);

    checkGradient3(ps);

    //ps.distributeControlDeltaGrid();
    //ps.distributeMeasurementDeltaGrid();

    //ps.Problem2HWaveEquationIBVP::implicit_calculate_D2V1();

    puts("Finished");
}

void Problem2HSolver::checkGradient3(const Problem2HSolver &prob)
{
    IPrinter::printSeperatorLine();
    DoubleVector pv;
    prob.OptimalParameterToVector(pv);
    printf("ok: "); IPrinter::print(pv.mid(0,  3), pv.mid(0,  3).length(), 9, 6);
    printf("oz: "); IPrinter::print(pv.mid(4,  7), pv.mid(4,  7).length(), 9, 6);
    printf("xy: "); IPrinter::print(pv.mid(8, 15), pv.mid(8, 15).length(), 9, 6);
    IPrinter::printSeperatorLine();

//    DoubleVector rv;
//    prob.equaPrm.RegularParameterToVector(rv);
//    printf("rk: "); IPrinter::print(rv.mid(0,  3), rv.mid(0,  3).length(), 9, 6);
//    printf("rz: "); IPrinter::print(rv.mid(4,  7), rv.mid(4,  7).length(), 9, 6);
//    printf("xy: "); IPrinter::print(rv.mid(8, 15), rv.mid(8, 15).length(), 9, 6);
//    IPrinter::printSeperatorLine();

    DoubleVector ag(pv.length());
    double functional = prob.fx(pv);
    printf("Functional: %f\n", functional);
    puts("Calculating gradients....");
    //const_cast<Problem2HDirichletDelta&>(prob).printLayers = false;
    prob.gradient(pv, ag);
    //const_cast<Problem2HDirichletDelta&>(prob).printLayers = false;
    puts("Gradients are calculated.");
    //    return;

    DoubleVector ng1(pv.length(), 0.0);
    DoubleVector ng2(pv.length(), 0.0);

    const unsigned int Nc = prob.Nc;
    const unsigned int No = prob.No;
    const unsigned int offset = Nc*No;

    puts("Calculating numerical gradients.... dh=0.01");
    printf("*** Calculating numerical gradients for k...... dh=0.01 ");
    IGradient::Gradient(&prob, 0.01, pv, ng1, 0*offset, 1*offset-1);
    printf("Calculated.\n");
    printf("*** Calculating numerical gradients for z...... dh=0.01 ");
    IGradient::Gradient(&prob, 0.01, pv, ng1, 1*offset, 2*offset-1);
    printf("Calculated.\n");
    printf("*** Calculating numerical gradients for xi..... dh=0.01 ");
    IGradient::Gradient(&prob, 0.01, pv, ng1, 2*offset+0*No, 2*offset+2*No-1);
    printf("Calculated.\n");
    printf("*** Calculating numerical gradients for eta.... dh=0.01 ");
    IGradient::Gradient(&prob, 0.01, pv, ng1, 2*offset+2*No, 2*offset+2*No+2*Nc-1);
    printf("Calculated.\n");

    puts("Calculating numerical gradients.... dh=0.001");
    printf("*** Calculating numerical gradients for k...... dh=0.001 ");
    IGradient::Gradient(&prob, 0.001, pv, ng2, 0*offset, 1*offset-1);
    printf("Calculated.\n");
    printf("*** Calculating numerical gradients for z...... dh=0.001 ");
    IGradient::Gradient(&prob, 0.001, pv, ng2, 1*offset, 2*offset-1);
    printf("Calculated.\n");
    printf("*** Calculating numerical gradients for xi..... dh=0.001 ");
    IGradient::Gradient(&prob, 0.001, pv, ng2, 2*offset+0*No, 2*offset+2*No-1);
    printf("Calculated.\n");
    printf("*** Calculating numerical gradients for eta.... dh=0.001 ");
    IGradient::Gradient(&prob, 0.001, pv, ng2, 2*offset+2*No, 2*offset+2*No+2*Nc-1);
    printf("Calculated.\n");

    const unsigned int N = 20;
    const unsigned int W = 9;
    const unsigned int P = 6;
    //k------------------------------------------------------//
    {
        IPrinter::printSeperatorLine("k");
        DoubleVector pk0 = pv.mid(0, offset-1);
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
        DoubleVector pz0 = pv.mid(offset, 2*offset-1);
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
        DoubleVector pe0 = pv.mid(2*offset, 2*offset+2*No-1);
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
        DoubleVector px0 = pv.mid(2*offset+2*No, 2*offset+2*No+2*Nc-1);
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
    Problem2HSolver* solver = const_cast<Problem2HSolver*>(this);
    solver->OptimalParameterFromVector(x);

    //const DoubleVector &Q1 = funcPrm.Q1;
    //const DoubleVector &Q2 = funcPrm.Q2;
    //std::vector<InitialPulse2D> &pulses = prob->equaPrm.pulses;

    //const double regEpsilon = funcPrm.regEpsilon;
    //const double r = funcPrm.r;

    double SUM = 0.0;
    //for (unsigned int q1=0; q1<Q1.length(); q1++)
    {
        //pulses[0].q = Q1[q1];
        //for (unsigned int q2=0; q2<Q2.length(); q2++)
        {
            //pulses[1].q = Q2[q2];

            solver->u_list.clear();
            solver->Problem2HWaveEquationIBVP::setInitialConditionMatrix(initialPulses, initialPulsesCount);
            solver->Problem2HWaveEquationIBVP::distributeControlDeltaGrid();
            solver->Problem2HWaveEquationIBVP::distributeMeasurementDeltaGrid();
            solver->Problem2HWaveEquationIBVP::implicit_calculate_D2V1();

            double intgrl = integral(u_list);

            //        double pnt = 0.0;
            //#ifdef USE_PENALTY
            //            pnt = penalty(u_info, equaPrm);
            //#endif
            //            double nrm = 0.0;
            //#ifdef USE_NORM
            //            nrm = norm(equaPrm);
            //#endif
            double sum = intgrl;// + regEpsilon*nrm + r*pnt;

            for (unsigned int i=0; i<u_list.size(); i++) solver->u_list[i].clear(); solver->u_list.clear();

            SUM += sum;// * (1.0/(double(Q1.length())*double(Q2.length())));
        }
    }
    return SUM;
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

void Problem2HSolver::gradient(const DoubleVector &pv, DoubleVector &g) const
{
    g.clear();
    g.resize(pv.length(), 0.0);

    //const double regEpsilon = funcPrm.regEpsilon;
    //const double r = funcPrm.r;

    Problem2HSolver* solver = const_cast<Problem2HSolver*>(this);
    solver->OptimalParameterFromVector(pv);

    //    const DoubleVector &Q1 = funcPrm.Q1;
    //    const DoubleVector &Q2 = funcPrm.Q2;

    //for (unsigned int q1=0; q1<Q1.length(); q1++)
    {
        //pulses[0].q = Q1[q1];
        //for (unsigned int q2=0; q2<Q2.length(); q2++)
        {
            //pulses[1].q = Q2[q2];

            solver->u_list.clear();
            solver->Problem2HWaveEquationIBVP::setInitialConditionMatrix(initialPulses, initialPulsesCount);
            solver->Problem2HWaveEquationIBVP::distributeControlDeltaGrid();
            solver->Problem2HWaveEquationIBVP::distributeMeasurementDeltaGrid();
            solver->Problem2HWaveEquationIBVP::implicit_calculate_D2V1();

            solver->Problem2HConjugateWaveEquationIBVP::implicit_calculate_D2V1();

            unsigned int gi = 0;

            // k
            if (true)
            {
                //puts("Calculating k gradients...");
                for (unsigned int i=0; i<Nc; i++)
                {
                    for (unsigned int j=0; j<No; j++)
                    {
                        double zij = z[i][j];

                        double grad_Kij = 0.0;
                        for (unsigned int s=0; s<Nt; s++)
                        {
                            const unsigned int ln = 2*solver->times[s].i;
                            grad_Kij += -(u_info[j].vl[ln] - zij) * p_info[i].vl[ln];
#ifdef USE_PENALTY
                            grad_Kij += -(uj.vl[ln] - zij) * 2.0*r*gpi(i,s,u_info,equaPrm)*sgn(g0i(i,s,u_info,equaPrm));
#endif
                        }


#ifdef USE_NORM
                        grad_Kij += +2.0*regEpsilon*(equaPrm.opt.k[s][i][j] - equaPrm.reg.k[s][i][j]);
#endif
                        g[gi++] += grad_Kij;// * (1.0/(double(Q1.length())*double(Q2.length())));
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
            if (true)
            {
                //puts("Calculating z gradients...");
                for (unsigned int i=0; i<Nc; i++)
                {
                    for (unsigned int j=0; j<No; j++)
                    {
                        double grad_Zij = 0.0;

                        double kij = k[i][j];
                        for (unsigned int s=0; s<Nt; s++)
                        {
                            const unsigned int ln = 2*times[s].i;
                            grad_Zij += kij * p_info[i].vl[ln];
#ifdef USE_PENALTY
                            grad_Zij += kij * 2.0*r*gpi(i,s,u_info,equaPrm)*sgn(g0i(i,s,u_info,equaPrm));
#endif
                        }
#ifdef USE_NORM
                        grad_Zij += +2.0*regEpsilon*(equaPrm.opt.z[s][i][j] - equaPrm.opt.z[s][i][j]);
#endif
                        g[gi++] = grad_Zij;// * (1.0/(double(Q1.length())*double(Q2.length())));
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
            if (true)
            {
                //puts("Calculating o gradients...");
                for (unsigned int j=0; j<No; j++)
                {
                    double gradXijX = 0.0;
                    double gradXijY = 0.0;

                    for (unsigned int s=0; s<Nt; s++)
                    {
                        const unsigned int ln = 2*solver->times[s].i;
                        for (unsigned int i=0; i<Nc; i++)
                        {
                            double kij = k[i][j];

                            gradXijX += kij * u_info[j].dx[ln] * p_info[i].vl[ln];
                            gradXijY += kij * u_info[j].dy[ln] * p_info[i].vl[ln];
#ifdef USE_PENALTY
                            gradXijX += kij * uj.dx[ln] * 2.0*r*gpi(i,s,u_info,equaPrm)*sgn(g0i(i,s,u_info,equaPrm));
                            gradXijY += kij * uj.dy[ln] * 2.0*r*gpi(i,s,u_info,equaPrm)*sgn(g0i(i,s,u_info,equaPrm));
#endif
                        }
                    }
#ifdef USE_NORM
                    gradXijX += 2.0*regEpsilon*(equaPrm.opt.ksi[j].x - equaPrm.reg.ksi[j].x);
                    gradXijY += 2.0*regEpsilon*(equaPrm.opt.ksi[j].y - equaPrm.reg.ksi[j].y);
#endif

                    g[gi++] += gradXijX;// * (1.0/(double(Q1.length())*double(Q2.length())));
                    g[gi++] += gradXijY;// * (1.0/(double(Q1.length())*double(Q2.length())));
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
            if (true)
            {
                //puts("Calculating c gradients...");
                for (unsigned int i=0; i<Nc; i++)
                {
                    double gradEtaiX = 0.0;
                    double gradEtaiY = 0.0;

                    for (unsigned int s=0; s<Nt; s++)
                    {
                        const unsigned int ln = 2*times[s].i;
                        for (unsigned int j=0; j<No; j++)
                        {
                            double kij = k[i][j];
                            double zij = z[i][j];

                            gradEtaiX += p_info[i].dx[ln] * kij * (u_info[j].vl[ln] - zij);
                            gradEtaiY += p_info[i].dy[ln] * kij * (u_info[j].vl[ln] - zij);
                        }
                    }
#ifdef USE_NORM
                    gradEtaiX += 2.0*regEpsilon*(equaPrm.opt.eta[i].x - equaPrm.reg.eta[i].x);
                    gradEtaiY += 2.0*regEpsilon*(equaPrm.opt.eta[i].y - equaPrm.reg.eta[i].y);
#endif
                    g[gi++] += gradEtaiX;// * (1.0/(double(Q1.length())*double(Q2.length())));
                    g[gi++] += gradEtaiY;// * (1.0/(double(Q1.length())*double(Q2.length())));
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

            solver->u_list.clear();
        }
    }
}

void Problem2HSolver::OptimalParameterFromVector(const DoubleVector &x)
{
    unsigned int index = 0;

    k.clear();   k.resize(Nc, No);
    z.clear();   z.resize(Nc, No);
    ksi.clear(); ksi.resize(No);
    eta.clear(); eta.resize(Nc);

    for (unsigned int i=0; i<Nc; i++)
    {
        for (unsigned int j=0; j<No; j++)
        {
            k[i][j] = x[index++];
        }
    }

    for (unsigned int i=0; i<Nc; i++)
    {
        for (unsigned int j=0; j<No; j++)
        {
            z[i][j] = x[index++];
        }
    }

    for (unsigned int j=0; j<No; j++)
    {
        ksi[j].x = x[index++];
        ksi[j].y = x[index++];
    }

    for (unsigned int i=0; i<Nc; i++)
    {
        eta[i].x = x[index++];
        eta[i].y = x[index++];
    }
}

void Problem2HSolver::OptimalParameterToVector(DoubleVector &x) const
{
    x.clear();
    x.resize(2*Nc*No+2*No+2*Nc);

    unsigned int index = 0;
    for (unsigned int i=0; i<Nc; i++)
    {
        for (unsigned int j=0; j<No; j++)
        {
            x[index++] = k[i][j];
        }
    }

    for (unsigned int i=0; i<Nc; i++)
    {
        for (unsigned int j=0; j<No; j++)
        {
            x[index++] = z[i][j];
        }
    }

    for (unsigned int j=0; j<No; j++)
    {
        x[index++] = ksi[j].x;
        x[index++] = ksi[j].y;
    }

    for (unsigned int i=0; i<Nc; i++)
    {
        x[index++] = eta[i].x;
        x[index++] = eta[i].y;
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

Problem2HCommon::~Problem2HCommon() {}

void Problem2HCommon::setParameterCounts(unsigned int Nc, unsigned int No, const Dimension &dimensionX, const Dimension &dimensionY, const Dimension &timeDimension)
{
    this->Nc = Nc;
    this->No = No;

    const unsigned int N = static_cast<unsigned int>(dimensionX.size());
    const unsigned int M = static_cast<unsigned int>(dimensionY.size());
    const unsigned int L = static_cast<unsigned int>(timeDimension.size());
    const double hx = dimensionX.step();
    const double hy = dimensionY.step();
    const unsigned int length = 2*L+1;

    _deltaGridControl = new DeltaGrid2D[Nc];
    p_info = new SpacePointInfo[Nc];
    for (unsigned int i=0; i<Nc; i++)
    {
        _deltaGridControl[i].initGrid(N, hx, M, hy);
        p_info[i].vl.resize(length);
        p_info[i].dx.resize(length);
        p_info[i].dy.resize(length);
    }

    _deltaGridMeasurement = new DeltaGrid2D[No];
    u_info = new SpacePointInfo[No];
    for (unsigned int j=0; j<No; j++)
    {
        _deltaGridMeasurement[j].initGrid(N, hx, M, hy);
        u_info[j].vl.resize(length);
        u_info[j].dx.resize(length);
        u_info[j].dy.resize(length);
    }
}

void Problem2HCommon::setOptimizedParameters(const DoubleMatrix &k, const DoubleMatrix &z, const std::vector<SpacePoint> &eta, const std::vector<SpacePoint> &ksi)
{
    this->k = k;
    this->z = z;
    this->eta = eta;
    this->ksi = ksi;
}

void Problem2HCommon::distributeControlDeltaGrid()
{
    for (unsigned int i=0; i<Nc; i++)
    {
        _deltaGridControl[i].resetGrid();
        _deltaGridControl[i].distributeGauss(eta[i], 1, 1);
    }
}

void Problem2HCommon::distributeMeasurementDeltaGrid()
{
    for (unsigned int j=0; j<No; j++)
    {
        _deltaGridMeasurement[j].resetGrid();
        _deltaGridMeasurement[j].distributeGauss(ksi[j], 1, 1);
    }
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

//    if (tn.i == 0 || tn.i == 1 || tn.i == 2 || tn.i == 3)
//    {
//        IPrinter::printMatrix(u);
//        IPrinter::printSeperatorLine();
//    }

    //layerInfoSave2TextFile(u, tn);
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
    const unsigned int N = static_cast<unsigned int>(_spaceDimensionX.size());
    const unsigned int M = static_cast<unsigned int>(_spaceDimensionY.size());
    const double hx = _spaceDimensionX.step();
    const double hy = _spaceDimensionY.step();

    this->initialPulses = initialPulses;
    this->initialPulsesCount = initialPulsesCount;

    for (unsigned int s=0; s<initialPulsesCount; s++)
    {
        InitialPulse &initialPulse = initialPulses[s];

        DeltaGrid2D deltaGrid;
        deltaGrid.initGrid(N, hx, M, hy);
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
    if (tn.i >= 600) u_list[tn.i-600] = u;

    const double ht = timeDimension().step();

    for (unsigned int j=0; j<No; j++)
    {
        double u_vl, u_dx, u_dy;
        u_vl = _deltaGridMeasurement[j].lumpPointGauss(u, u_dx, u_dy);
        u_info[j].vl[tn.i] = u_vl;
        u_info[j].dx[tn.i] = u_dx;
        u_info[j].dy[tn.i] = u_dy;
    }

    for (unsigned int s=0; s<Nt; s++)
    {
        double wt = 0.0;
        if (tn.i == 2*times[s].i) { wt = 2.0/ht; } else { continue; }

        double* v = new double[Nc];

        for (unsigned int i=0; i<Nc; i++)
        {
            v[i] = 0.0;
            for (unsigned int j=0; j<No; j++)
            {
                v[i] += k.at(i,j) * (u_info[j].vl[tn.i]-z.at(i,j));
            }
        }

        for (unsigned int i=0; i<Nc; i++)
        {
            const unsigned minX = _deltaGridControl[i].minX();
            const unsigned maxX = _deltaGridControl[i].maxX();
            const unsigned minY = _deltaGridControl[i].minY();
            const unsigned maxY = _deltaGridControl[i].maxY();

            for (unsigned int m=minY; m<=maxY; m++)
            {
                for (unsigned int n=minX; n<=maxX; n++)
                {
                    if (i==0) f_crLayerMatrix[m][n] = 0.0;
                    f_crLayerMatrix[m][n] += v[i] * _deltaGridControl[i].weight(n, m) * wt;
                }
            }
        }

        delete [] v;
    }
}

/*********************************** Problem2HWaveEquationIBVP ***********************************************************/

/*********************************** Problem2HConjugateWaveEquationIBVP***************************************************/

double Problem2HConjugateWaveEquationIBVP::initial(const SpaceNodePDE &, InitialCondition) const
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
    IConjugateWaveEquationIBVP::setSpaceDimensions(dimensionX, dimensionY);
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
        p_vl = _deltaGridControl[i].lumpPointGauss(p, p_dx, p_dy);
        p_info[i].vl[tn.i] = p_vl;
        p_info[i].dx[tn.i] = p_dx;
        p_info[i].dy[tn.i] = p_dy;
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
        for (unsigned int m=0; m<=M; m++)
        {
            for (unsigned int n=0; n<=N; n++)
            {
                b_crLayerMatrix[m][n] = 0.0;
            }
        }
    }

    for (unsigned int s=0; s<Nt; s++)
    {
        double wt = 0.0;
        if (tn.i == 2*times[s].i) { wt = 2.0/ht; } else { continue; }

        double* w = new double[No];

        for (unsigned int i=0; i<No; i++)
        {
            w[i] = 0.0;
            for (unsigned int j=0; j<No; j++)
            {
                w[i] += k.at(i,j) * p_info[j].vl[tn.i];
            }
        }

        for (unsigned int j=0; j<Nc; j++)
        {
            const unsigned minX = _deltaGridMeasurement[j].minX();
            const unsigned maxX = _deltaGridMeasurement[j].maxX();
            const unsigned minY = _deltaGridMeasurement[j].minY();
            const unsigned maxY = _deltaGridMeasurement[j].maxY();

            for (unsigned int m=minY; m<=maxY; m++)
            {
                for (unsigned int n=minX; n<=maxX; n++)
                {
                    b_crLayerMatrix[m][n] += w[j] * _deltaGridMeasurement[j].weight(n, m) * wt;
                }
            }
        }

        delete [] w;
    }
}

void Problem2HConjugateWaveEquationIBVP::layerInfoSave2TextFile(const DoubleMatrix &u, const TimeNodePDE & tn) const {}

/*********************************** Problem2HConjugateWaveEquationIBVP***************************************************/
