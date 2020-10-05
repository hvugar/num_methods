#include "solver2.h"

using namespace p3p0;

void Functional::Main(int /*argc*/, char **/*argv*/)
{
    const double _thermalDiffusivity  = 0.01;
    const double _thermalConductivity = 0.00;
    const double _thermalConvection   = 0.0001;
    const size_t _heatSourceNumber = 2;
    const size_t _measrPointNumber = 4;

    const auto time_min = 0;
    const auto time_max = 200;
    const auto time_step = 0.005;
    const auto time_size = 201;
    const Dimension timeDimension(time_step, time_min, time_max);

    const auto dimX_min = 0;
    const auto dimX_max = 100;
    const auto dimX_step = 0.01;
    const auto dimX_size = 101;
    const Dimension spaceDimensionX(dimX_step, dimX_min, dimX_max);

    Functional functional(_thermalDiffusivity, _thermalConductivity, _thermalConvection);
    functional.setTimeDimension(timeDimension);
    functional.setSpaceDimensionX(spaceDimensionX);
    functional.setControlSize(_heatSourceNumber, _measrPointNumber);

#ifdef OPTIMIZE_Y
    const auto vector_size = _measrPointNumber * (3*_heatSourceNumber + 1);
#else
    const auto vector_size = time_size;
    const unsigned int start[] = { 0*vector_size+0, 1*vector_size+0 };
    const unsigned int finsh[] = { 1*vector_size-1, 2*vector_size-1 };
#endif


    /******************************************************************
     *                     Calculating V(x) function
     ******************************************************************/
    {
#ifdef OPTIMIZE_Y
        DoubleVector x(vector_size, 0.0);
        for (size_t j=0; j<_measrPointNumber; j++)
        {
            for (size_t i=0; i<_heatSourceNumber; i++)
            {
                x.at(0*_heatSourceNumber*_measrPointNumber + i*_measrPointNumber + j) = -0.5*sin(2.0*(30.0*i+1)+3.0*(12.0*j+1));
                x.at(1*_measrPointNumber*_heatSourceNumber + i*_measrPointNumber + j) = +0.4*sin(5.0*(10.0*i+1)+5.0*(15.0*j+1));
                x.at(2*_measrPointNumber*_heatSourceNumber + i*_measrPointNumber + j) = +4.0*cos(2.0*(i+1.0)+3.0*(j+2.0));
            }
            x.at(3*_measrPointNumber*_heatSourceNumber + j) = 0.2*(j+1);
        }

        puts("---");
        printf("a:");IPrinter::print(x.mid(0x00, 0x07), 8, functional._w, functional._p);
        printf("b:");IPrinter::print(x.mid(0x08, 0x0F), 8, functional._w, functional._p);
        printf("o:");IPrinter::print(x.mid(0x10, 0x17), 8, functional._w, functional._p);
        printf("e:");IPrinter::print(x.mid(0x18, 0x1B), 4, functional._w, functional._p);

        functional.convertFromVector(x);
        functional.forward.implicit_calculate_D1V1();

        //printf("q1:"); for (size_t ln=0; ln<time_size; ln++) { printf("%12.6f ", functional.mq[ln][0]); } printf("\n");
        //printf("q2:"); for (size_t ln=0; ln<time_size; ln++) { printf("%12.6f ", functional.mq[ln][1]); } printf("\n");

#else
        DoubleVector x(vector_size*functional.heatSourceNumber);
        for (unsigned int ln=0; ln<vector_size; ln++)
        {
            DoubleVector _qNorm = functional.qNorm1(ln * time_step);
            x[0*vector_size + ln] = _qNorm[0];
            x[1*vector_size + ln] = _qNorm[1];
        }
        IPrinter::printVector(functional._w, functional._p, x.mid(start[0], finsh[0]), "q1:");
        IPrinter::printVector(functional._w, functional._p, x.mid(start[1], finsh[1]), "q2:");
        functional.convertFromVector(x);
        functional.forward.implicit_calculate_D1V1();
#endif
        functional.V = functional.U;
        functional.V.clear();
        functional.V.resize(dimX_size, 5.0);
        IPrinter::printVector(functional._w, functional._p, functional.U, "V: ");
    }
    /******************************************************************
     *                     Calculating V(x) function
     ******************************************************************/

#ifdef OPTIMIZE_Y
    DoubleVector x(vector_size, 0.0);
    for (size_t j=0; j<_measrPointNumber; j++)
    {
        for (size_t i=0; i<_heatSourceNumber; i++)
        {
            x.at(0*_heatSourceNumber*_measrPointNumber + i*_measrPointNumber + j) = +0.5*sin(2.0*(30.0*i+1)+3.0*(12.0*j+1));
            x.at(1*_measrPointNumber*_heatSourceNumber + i*_measrPointNumber + j) = +0.6*cos(5.0*(12.0*i+1)+5.0*(10.0*j+1));
            x.at(2*_measrPointNumber*_heatSourceNumber + i*_measrPointNumber + j) = 4.0*cos(2.0*(i+1.0)+3.0*(j+2.0));
        }
        x.at(3*_measrPointNumber*_heatSourceNumber + j) = 0.2*(j+1);
    }
    //x[0x18] = 0.24; x[0x19] = 0.36; x[0x1A] = 0.64; x[0x1B] = 0.76;

    x[0] = -1.0000; x[2] = -1.0000; x[2] = -1.0000; x[3] = -1.0000; x[4] = -1.0000; x[5] = -1.0000; x[6] = -1.0000; x[7] = -1.0000;
    x[16] = 5.0000; x[17] = 5.0000; x[18] = 5.0000; x[19] = 5.0000; x[20] = 5.0000; x[21] = 5.0000; x[22] = 5.0000; x[23] = 5.0000;
    x[24] = 0.1634; x[25] = 0.3895; x[26] = 0.6836; x[27] = 0.8595;
    //x[24] = 0.2094; x[25] = 0.4734; x[26] = 0.5592; x[27] = 0.7868;

    x.clear(); /*0.00007882*/
    x << -1.0545 <<  0.6535 << -1.1357 << -1.1477 << -0.7377 << 0.9597 << -1.1885 << -0.7221
      << -0.0478 << -0.4230 << -0.8102 << -0.9891 << -0.0093 << 0.2945 <<  0.5956 <<  0.7514
      <<  5.2754 <<  5.0960 <<  5.2618 <<  5.1903 <<  4.9103 << 5.0278 <<  4.7885 <<  4.8992
      <<  0.1634 <<  0.3895 <<  0.6836 <<  0.8595;


    //    x.clear(); /*0.00003079*/
    //    x << 0.9616 << -1.9615 <<  0.1521 << -1.2175 << -0.0817 << -1.5069 << -0.5075 << -1.2325
    //      << 0.3984 <<  0.1484 << -0.3091 << -0.4197 <<  0.2276 <<  0.4323 <<  0.7706 <<  0.9026
    //      << 5.2333 <<  5.0948 <<  5.4212 <<  5.3341 <<  4.9449 <<  4.8119 <<  4.9440 <<  4.8672
    //      << 0.1078 <<  0.2612 <<  0.7232 << 0.9125;

    DoubleVector g(x.length());
    DoubleVector g1(x.length());
    puts("---");
    printf("a:");IPrinter::print(x.mid(0x00, 0x07), 8, functional._w, functional._p);
    printf("b:");IPrinter::print(x.mid(0x08, 0x0F), 8, functional._w, functional._p);
    printf("o:");IPrinter::print(x.mid(0x10, 0x17), 8, functional._w, functional._p);
    printf("e:");IPrinter::print(x.mid(0x18, 0x1B), 4, functional._w, functional._p);
#else
    DoubleVector x(vector_size*functional.heatSourceNumber, 2.0);
    //x[0*vctr_size] = 0.0;
    //x[1*vctr_size] = 0.0;

    //    for (unsigned int ln=0; ln<vctr_size; ln++)
    //    {
    //        double t = ln * time_step;
    //        DoubleVector _qNorm = s.qNorm1(t);
    //        x[0*vctr_size + ln] = _qNorm[0] * 1.1;//( rand()%2==0 ? 0.8 : 1.2 );
    //        x[1*vctr_size + ln] = _qNorm[1] * 1.1;//( rand()%2==0 ? 0.8 : 1.2 );
    //    }

    DoubleVector g(x.length());
    DoubleVector g1(x.length());

    puts("---");
    IPrinter::printVector(functional._w, functional._p, x.mid(start[0], finsh[0]), "q1:");
    IPrinter::printVector(functional._w, functional._p, x.mid(start[1], finsh[1]), "q2:");
#endif

    puts("---");
#ifdef OPTIMIZE_Y
    functional.gradient(x, g);
    printf("a:");IPrinter::print(g.mid(0x00, 0x07).L2Normalize(), 8, functional._w, functional._p);
    printf("b:");IPrinter::print(g.mid(0x08, 0x0F).L2Normalize(), 8, functional._w, functional._p);
    printf("o:");IPrinter::print(g.mid(0x10, 0x17).L2Normalize(), 8, functional._w, functional._p);
    printf("e:");IPrinter::print(g.mid(0x18, 0x1B).L2Normalize(), 4, functional._w, functional._p);

    puts("---");
    IGradient::Gradient(&functional, 0.001, x, g1, static_cast<size_t>(0x00), static_cast<size_t>(0x07));
    IGradient::Gradient(&functional, 0.001, x, g1, static_cast<size_t>(0x08), static_cast<size_t>(0x0F));
    IGradient::Gradient(&functional, 0.001, x, g1, static_cast<size_t>(0x10), static_cast<size_t>(0x17));
    IGradient::Gradient(&functional, 0.001, x, g1, static_cast<size_t>(0x18), static_cast<size_t>(0x1B));
    printf("a:");IPrinter::print(g1.mid(0x00, 0x07).L2Normalize(), 0x08, functional._w, functional._p);
    printf("b:");IPrinter::print(g1.mid(0x08, 0x0F).L2Normalize(), 0x08, functional._w, functional._p);
    printf("o:");IPrinter::print(g1.mid(0x10, 0x17).L2Normalize(), 0x08, functional._w, functional._p);
    printf("e:");IPrinter::print(g1.mid(0x18, 0x1B).L2Normalize(), 0x04, functional._w, functional._p);
#else
    functional.gradient(x, g);
    IPrinter::printVector(functional._w, functional._p, g.mid(start[0], finsh[0]).L2Normalize(), "g1:");
    IPrinter::printVector(functional._w, functional._p, g.mid(start[1], finsh[1]).L2Normalize(), "g2:");

    puts("---");
    IGradient::Gradient(&functional, 0.01, x, g1);
    IPrinter::printVector(functional._w, functional._p, g1.mid(start[0], finsh[0]).L2Normalize(), "g1:");
    IPrinter::printVector(functional._w, functional._p, g1.mid(start[1], finsh[1]).L2Normalize(), "g2:");
#endif

    IPrinter::printSeperatorLine();
    puts("Starting optimization...");

    // while (functional.epsilon > 0.00001)
    {

        // while (functional.R < 2.0)
        double step = 0.5;
        double epsl = 0.0001;
        int i=0;
        while (true)
        {
            i++;
            GradientBasedMethod *gm;
            gm = new ConjugateGradient;
            //gm = new SteepestDescentGradient;
            //gm = new ConstStepGradient;
            gm->setNormalize(true);

            gm->setFunction(&functional);
            gm->setGradient(&functional);
            gm->setPrinter(&functional);
            gm->setProjection(&functional);
            //gm.setGradientNormalizer(&prob);
            gm->setOptimalityTolerance(0.000001);
            gm->setStepTolerance(0.000001);
            gm->setFunctionTolerance(0.000001);
            gm->setR1MinimizeEpsilon(step, epsl);
            //gm.setMaxIterationCount(10);
            //gm.setNormalize(false);
            //gm->setNormalize(i%2==0);
            gm->showExitMessage(true);
            gm->calculate(x);

            delete gm;

            //            functional.R *= 10.0;

            //step *= 0.9;
            //epsl *= 0.9;

            IPrinter::print(x, x.length(), functional._w, functional._p);
        }

        //        functional.R = 1.0;
        //        functional.epsilon *= 0.1;
    }

    puts("Optimization is finished.");
    IPrinter::print(x, x.length(), functional._w, functional._p);


    //    for (size_t ln=0; ln<=time_max; ln++)
    //    {
    //        printf("%4d %10.8f %10.8f\n", ln, functional.mq[ln][0], functional.mq[ln][1]);
    //    }

    //    for (int ln=time_max; ln<=10*time_size; ln++)
    //    {
    //        functional.setTimeDimension(Dimension(time_step, 0, ln));

    //        double t = ln*time_step;
    //        double fx = functional.fx(x);
    //        double a1 = functional.mq[ln][0];
    //        double a2 = functional.mq[ln][1];

    //        printf("time: %6.4f fx: %16.10f q1: %16.10f q2: %16.10f\n", t, fx, a1, a2);
    //    }
}

/**************************************************************************************************************************/

auto Functional::print(unsigned int it, const DoubleVector &x, const DoubleVector &g, double /*f*/, double _alpha, GradientBasedMethod::MethodResult result) const -> void
{
    //const size_t time_size = timeDimension().size();
    //const auto vector_size = time_size;
    //IPrinter::printSeperatorLine();

    const auto _fx = fx(x);
    const auto _integral = integral(x);
    const auto _norm = norm(x);
    const auto _penalty = penalty(x);
    printf_s("I[%4d] fx: %10.8f int: %10.8f nrm: %10.8f penalty: %10.8f eps: %10.8f R: %10.8f alpha: %10.8f res: %d\n", it, _fx, _integral, _norm, _penalty, epsilon, R, _alpha, result);
#ifdef OPTIMIZE_Y

    if (
            result == GradientBasedMethod::MethodResult::BREAK_STEP_TOLERANCE ||
            result == GradientBasedMethod::MethodResult::BREAK_FUNCTION_TOLERANCE ||
            result == GradientBasedMethod::MethodResult::BREAK_OPTIMALITY_TOLERANCE
            )
    {
        //        IPrinter::print(x, x.length(), _w, _p);
        //        IPrinter::print(g, g.length(), _w, _p);

        //printf("a:");IPrinter::print(x.mid(0x00, 0x07), 8, _w, _p);
        //printf("b:");IPrinter::print(x.mid(0x08, 0x0F), 8, _w, _p);
        //printf("o:");IPrinter::print(x.mid(0x10, 0x17), 8, _w, _p);
        //printf("e:");IPrinter::print(x.mid(0x18, 0x1B), 4, _w, _p);
        //    puts("---");Q
        //    printf("a:");IPrinter::print(g.mid(0x00, 0x07).L2Normalize(), 8, _w, _p);
        //    printf("b:");IPrinter::print(g.mid(0x08, 0x0F).L2Normalize(), 8, _w, _p);
        //    printf("o:");IPrinter::print(g.mid(0x10, 0x17).L2Normalize(), 8, _w, _p);
        //    printf("e:");IPrinter::print(g.mid(0x18, 0x1B).L2Normalize(), 4, _w, _p);
        //puts("---");
    }
#else
    IPrinter::printVector(_w, _p, x.mid(0*vector_size, 1*vector_size-1), "q1");
    IPrinter::printVector(_w, _p, x.mid(1*vector_size, 2*vector_size-1), "q2");
    IPrinter::printVector(_w, _p, g.mid(0*vector_size, 1*vector_size-1), "g1");
    IPrinter::printVector(_w, _p, g.mid(1*vector_size, 2*vector_size-1), "g2");
#endif
    //    IPrinter::printVector(_w, _p, U, "U:");
    //    IPrinter::printVector(_w, _p, V, "V:");
}

auto CommonParameter::setTimeDimension(const Dimension &timeDimension) -> void
{
    const auto last_size = _timeDimension.size();

    _timeDimension = timeDimension;
    const auto time_size = _timeDimension.size();
    const auto time_step = _timeDimension.step();

    if (mq != nullptr)
    {
        for (size_t ln=0; ln<last_size; ln++) { mq[ln].clear(); }
        delete [] mq;
        mq = nullptr;
    }
    mq = new DoubleVector[time_size];
    for (size_t ln=0; ln<time_size; ln++) { mq[ln].resize(heatSourceNumber); }

    if (mz != nullptr)
    {
        for (size_t ln=0; ln<last_size; ln++) { mz[ln].clear(); }
        delete [] mz;
        mz = nullptr;
    }
    mz = new DoubleVector[time_size];
    for (size_t ln=0; ln<time_size; ln++)
    {
        mz[ln].resize(heatSourceNumber);

        double t = ln * time_step;
        mz[ln][0] = 0.05 + 0.90*fabs(sin(M_PI*t));
        mz[ln][1] = 0.95 - 0.90*fabs(sin(M_PI*t));
        //mz[ln][0] = 0.30;
        //mz[ln][1] = 0.70;
    }

    if (mp != nullptr)
    {
        for (size_t ln=0; ln<last_size; ln++) { mp[ln].clear(); }
        delete [] mp;
        mp = nullptr;
    }
    mp = new DoubleVector[time_size];
    for (size_t ln=0; ln<time_size; ln++) { mp[ln].resize(heatSourceNumber); }

#ifdef OPTIMIZE_Y
    if (uv != nullptr)
    {
        for (size_t ln=0; ln<last_size; ln++) { uv[ln].clear(); ud[ln].clear(); }
        delete [] uv;
        delete [] ud;
        uv = nullptr;
        ud = nullptr;
    }
    uv = new DoubleVector[time_size];
    ud = new DoubleVector[time_size];
    for (size_t ln=0; ln<time_size; ln++) { uv[ln].resize(measrPointNumber); ud[ln].resize(measrPointNumber); }
#endif

    if (qMin == nullptr || qMax == nullptr)
    {
        qMin = new DoubleVector[time_size];
        qMax = new DoubleVector[time_size];
        for (size_t ln=0; ln<time_size; ln++)
        {
            qMin[ln] << -0.1 << -0.1;
            qMax[ln] << +0.1 << +0.1;
        }
    }
    else
    {
        for (size_t ln=0; ln<last_size; ln++) { qMin[ln].clear(); qMax[ln].clear();  }
        delete [] qMin;
        delete [] qMax;
        qMin = nullptr;
        qMax = nullptr;
    }
}

auto CommonParameter::setSpaceDimensionX(const Dimension &spaceDimensionX) -> void
{
    _spaceDimensionX = spaceDimensionX;
    const auto spaceX_size = _spaceDimensionX.size();
    U.resize(spaceX_size, 0.0);
    V.resize(spaceX_size, 0.0);
}

auto CommonParameter::setControlSize(size_t heatSourceNumber, size_t measrPointNumber) -> void
{
    this->heatSourceNumber = heatSourceNumber;

#ifdef OPTIMIZE_Y
    this->measrPointNumber = measrPointNumber;

    alpha.resize(heatSourceNumber, measrPointNumber, 0.0);
    betta.resize(heatSourceNumber, measrPointNumber, 0.0);
    omega.resize(heatSourceNumber, measrPointNumber, 0.0);
    mPnts.resize(measrPointNumber, 0.0);

    alphaN.resize(heatSourceNumber, measrPointNumber, 0.0);
    bettaN.resize(heatSourceNumber, measrPointNumber, 0.0);
    omegaN.resize(heatSourceNumber, measrPointNumber, 0.0);
    mPntsN.resize(measrPointNumber, 0.0);

    alphaN[0][0] = -0.7873; alphaN[0][1] = -1.0482; alphaN[0][2] = +0.1091; alphaN[0][3] = +1.7487;
    alphaN[1][0] = +0.2138; alphaN[1][1] = -1.6883; alphaN[1][2] = -0.0406; alphaN[1][3] = +0.2155;
    bettaN[0][0] = -1.1766; bettaN[0][1] = -1.1083; bettaN[0][2] = -0.9647; bettaN[0][3] = -0.7437;
    bettaN[1][0] = +0.1263; alphaN[1][1] =  0.0839; bettaN[1][2] =  0.0227; bettaN[1][3] = -0.0627;
    omegaN[0][0] = -0.5157; omegaN[0][1] = -0.0206; omegaN[0][2] =  0.3726; omegaN[0][3] = -1.0293;
    omegaN[1][0] = -3.3899; omegaN[1][1] =  3.7042; omegaN[1][2] = -3.7916; omegaN[1][3] =  3.9368;
    mPntsN[0]    =  0.2000; mPntsN[1]    =  0.4000; mPntsN[2]    =  0.6000; mPntsN[3]    =  0.8000;
#endif
}

auto CommonParameter::q(const TimeNodePDE &tn) const -> DoubleVector { return mq[tn.i]; }

auto CommonParameter::z(const TimeNodePDE &tn) const -> DoubleVector { return mz[tn.i]; }

auto CommonParameter::g0(size_t i, size_t ln) const -> double
{
    return 0.5*(qMax[ln][i] + qMin[ln][i]) - mq[ln][i];
}

auto CommonParameter::gi(size_t i, size_t ln) const -> double
{
    double _g0 = g0(i, ln);
    return fabs(_g0) - 0.5*(qMax[ln][i] - qMin[ln][i]);
}

auto CommonParameter::gp(size_t i, size_t ln) const -> double
{
    double _gi = gi(i, ln);
    return _gi > 0.0 ? _gi : 0.0;
}

/**************************************************************************************************************************/

HeatEquationIBVP::HeatEquationIBVP() {}

HeatEquationIBVP::HeatEquationIBVP(const HeatEquationIBVP &h) : IHeatEquationIBVP (h) { *this = *this; }

HeatEquationIBVP & HeatEquationIBVP::operator =(const HeatEquationIBVP &) { return *this; }

HeatEquationIBVP::~HeatEquationIBVP() {}

auto HeatEquationIBVP::initial(const SpaceNodePDE &/*sn*/, InitialCondition /*cn*/) const -> double
{
    return common->initialTemperature;
}

auto HeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &/*tn*/, BoundaryConditionPDE &bc) const -> double
{
    const double _lmbd1 = common->lambda1();
    const double _theta = common->theta();
    if (sn.i == spaceDimensionX().min()) { bc = BoundaryConditionPDE::Robin(_lmbd1, -1.0, _lmbd1); return _theta; }
    if (sn.i == spaceDimensionX().max()) { bc = BoundaryConditionPDE::Robin(_lmbd1, +1.0, _lmbd1); return _theta; }
    //if (sn.i == spaceDimensionX().min()) { bc = BoundaryConditionPDE::Neumann(+1.0, 0.0); return 0.0; }
    //if (sn.i == spaceDimensionX().max()) { bc = BoundaryConditionPDE::Neumann(+1.0, 0.0); return 0.0; }
    throw std::exception();
}

auto HeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    const auto sigma = 0.01;
    DoubleVector _z = common->z(tn);
    DoubleVector _q = common->q(tn);
    auto fx = 0.0;
    fx += _q[0] * DeltaFunction::gaussian(sn.x, _z[0], sigma);
    fx += _q[1] * DeltaFunction::gaussian(sn.x, _z[1], sigma);
    return fx - thermalConvection()*common->theta();
}

auto HeatEquationIBVP::layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const -> void
{
    if (static_cast<int>(tn.i) == timeDimension().max()) { common->U = u; }

#ifdef OPTIMIZE_Y
    const size_t heatSourceNumber = common->heatSourceNumber;
    const size_t measrPointNumber = common->measrPointNumber;
    const auto hx = spaceDimensionX().step();
    const auto Nx = spaceDimensionX().size();

    for (size_t j=0; j<measrPointNumber; j++)
    {
        const double mp = common->mPnts[j];
        common->uv[tn.i][j] = DeltaFunction::lumpedPoint4(u, mp, hx, Nx, common->ud[tn.i][j]);

        //printf("%4d %12.6f %12.6f ", j, common->uv[tn.i][j], common->ud[tn.i][j]);
    }
    //puts("");

    if (tn.i < static_cast<unsigned int>(timeDimension().max()))
    {
        for (size_t i=0; i<heatSourceNumber; i++)
        {
            common->mq[tn.i+1].at(i) = 0.0;
            const auto zi = common->mz[tn.i].at(i);

            for (size_t j=0; j<measrPointNumber; j++)
            {
                const double mp = common->mPnts[j];
                common->mq[tn.i+1].at(i) += common->alpha.at(i,j) * (common->uv[tn.i].at(j) - common->omega.at(i,j));
                common->mq[tn.i+1].at(i) += common->betta.at(i,j) * (zi - mp) * (zi - mp);
            }
        }
    }

    if (tn.i == 0)
    {
        for (size_t i=0; i<heatSourceNumber; i++) { common->mq[0].at(i) = common->mq[1].at(i); }
    }
#endif
}

auto HeatEquationIBVP::timeDimension() const -> Dimension { return common->timeDimension(); }

auto HeatEquationIBVP::spaceDimensionX() const -> Dimension { return common->spaceDimensionX(); }

auto HeatEquationIBVP::spaceDimensionY() const -> Dimension { return spaceDimensionX(); }

auto HeatEquationIBVP::spaceDimensionZ() const -> Dimension { return spaceDimensionX(); }

/**************************************************************************************************************************/

HeatEquationFBVP::HeatEquationFBVP() {}

HeatEquationFBVP::HeatEquationFBVP(const HeatEquationFBVP &h) : IHeatEquationFBVP (h) {}

HeatEquationFBVP & HeatEquationFBVP::operator =(const HeatEquationFBVP &) { return *this; }

HeatEquationFBVP::~HeatEquationFBVP() {}

auto HeatEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition /*ic*/) const -> double
{
    size_t n = static_cast<size_t>(sn.i);
    return -2.0 * common->mu(sn) * ( common->U[n] - common->V[n] );
}

auto HeatEquationFBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &/*tn*/, BoundaryConditionPDE &bc) const -> double
{
    const double _lmbd1 = common->lambda1();
    if (sn.i == spaceDimensionX().min()) { bc = BoundaryConditionPDE::Robin(_lmbd1, -1.0, 0.0); return 0.0; }
    if (sn.i == spaceDimensionX().max()) { bc = BoundaryConditionPDE::Robin(_lmbd1, +1.0, 0.0); return 0.0; }
    //if (sn.i == spaceDimensionX().min()) { bc = BoundaryConditionPDE::Neumann(+1.0, 0.0); return 0.0; }
    //if (sn.i == spaceDimensionX().max()) { bc = BoundaryConditionPDE::Neumann(+1.0, 0.0); return 0.0; }
    throw std::exception();
}

auto HeatEquationFBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    double fx = 0.0;

#ifdef OPTIMIZE_Y
    const auto sigma = 0.01;
    for (size_t j=0; j<common->measrPointNumber; j++)
    {
        double fxj = 0.0;
        for (size_t i=0; i<common->heatSourceNumber; i++)
        {
            fxj += common->alpha[i][j] * ( common->mp[tn.i+1][i]
                    + 2.0 * common->R * common->gp(i,tn.i) * sgn( common->g0(i,tn.i) ));
        }
        fx -= fxj * DeltaFunction::gaussian(sn.x, common->mPnts[j], sigma);
    }
#endif
    return fx;
}

auto HeatEquationFBVP::HeatEquationFBVP::layerInfo(const DoubleVector &psi, const TimeNodePDE &tn) const -> void
{
    const auto hx = spaceDimensionX().step();
    const auto Nx = spaceDimensionX().size();
    DoubleVector _z = common->z(tn);
    common->mp[tn.i][0] = DeltaFunction::lumpedPoint4(psi, _z[0], hx, Nx);
    common->mp[tn.i][1] = DeltaFunction::lumpedPoint4(psi, _z[1], hx, Nx);
}

auto HeatEquationFBVP::timeDimension() const -> Dimension { return common->timeDimension(); }

auto HeatEquationFBVP::spaceDimensionX() const -> Dimension { return common->spaceDimensionX(); }

auto HeatEquationFBVP::spaceDimensionY() const -> Dimension { return spaceDimensionX(); }

auto HeatEquationFBVP::spaceDimensionZ() const -> Dimension { return spaceDimensionX(); }

/**************************************************************************************************************************/

CommonParameter::~CommonParameter() {}

auto CommonParameter::convertFromVector(const DoubleVector &x) -> void
{
#ifdef OPTIMIZE_Y
    for (unsigned int j=0; j<measrPointNumber; j++)
    {
        for (unsigned int i=0; i<heatSourceNumber; i++)
        {
            alpha.at(i,j) = x.at((0*heatSourceNumber*measrPointNumber) + i*measrPointNumber + j);
            betta.at(i,j) = x.at((1*heatSourceNumber*measrPointNumber) + i*measrPointNumber + j);
            omega.at(i,j) = x.at((2*heatSourceNumber*measrPointNumber) + i*measrPointNumber + j);
        }
        mPnts.at(j) = x.at((3*heatSourceNumber*measrPointNumber) + j);
    }

    for (unsigned int i=0; i<heatSourceNumber; i++) { mq[0].at(i) = 0.0; }

#else
    const auto time_size = _timeDimension.size();

    for (size_t ln=0; ln<time_size; ln++)
    {
        for (size_t i=0; i<heatSourceNumber; i++)
        {
            mq[ln][i] = x[i*time_size+ln];
        }
    }
#endif
}

auto CommonParameter::convertToVector(DoubleVector &x) const -> void
{
    const auto time_size = _timeDimension.size();

    x.clear();
    x.resize(heatSourceNumber*time_size);

    if (mq != nullptr)
    {
        for (size_t ln=0; ln<time_size; ln++)
        {
            for (size_t i=0; i<heatSourceNumber; i++)
            {
                x[i*time_size+ln] = mq[ln][i];
            }
        }
    }
}


auto CommonParameter::qNorm1(double t) const -> DoubleVector
{
    DoubleVector q(2, 0.0);
    q[0] = (t+2.0)*(t+2.0);
    q[1] = (t+2.0)*(t+2.0)*(t+2.0);
    return q;
}

/**************************************************************************************************************************/

Functional::Functional(double thermalDiffusivity, double thermalConductivity, double thermalConvection) : CommonParameter ()
{
    forward.setThermalDiffusivity(thermalDiffusivity);
    forward.setThermalConductivity(thermalConductivity);
    forward.setThermalConvection(-thermalConvection);

    backward.setThermalDiffusivity(-thermalDiffusivity);
    backward.setThermalConductivity(thermalConductivity);
    backward.setThermalConvection(+thermalConvection);

    forward.common = this;
    backward.common = this;
    forward._userHalfValues = backward._userHalfValues = false;
}

auto Functional::gradient(const DoubleVector &x, DoubleVector &g) const -> void
{
    const auto time_size = timeDimension().size();
    const auto time_step = timeDimension().step();
    const_cast<Functional*>(this)->convertFromVector(x);

    forward.implicit_calculate_D1V1();
    backward.implicit_calculate_D1V1();

    g.resize(x.length(), 0.0);

#ifdef OPTIMIZE_Y

    for (size_t j=0; j<measrPointNumber; j++)
    {
        for (size_t i=0; i<heatSourceNumber; i++)
        {
            double g1 = 0.0;
            double g2 = 0.0;
            double g3 = 0.0;
            double g4 = 0.0;

            size_t ln = 0;
            double _gpln = 2.0 * R * gp(i,ln) * sgn( g0(i,ln) );
            g1 += 0.5*(mp[ln][i] + _gpln) * (uv[ln][j]-omega[i][j]);
            g2 += 0.5*(mp[ln][i] + _gpln) * (mz[ln][i]-mPnts[j])*(mz[ln][i]-mPnts[j]);
            g3 += 0.5*(mp[ln][i] + _gpln) * (alpha[i][j]);
            g4 += 0.5*(mp[ln][i] + _gpln) * (alpha[i][j]*ud[ln][j] - 2.0*betta[i][j]*(mz[ln][j]-mPnts[j]));
            for (ln=1; ln<time_size-1; ln++)
            {
                _gpln = 2.0 * R * gp(i,ln) * sgn( g0(i,ln) );
                g1 += (mp[ln][i] + _gpln) * (uv[ln][j]-omega[i][j]);
                g2 += (mp[ln][i] + _gpln) * (mz[ln][i]-mPnts[j])*(mz[ln][i]-mPnts[j]);
                g3 += (mp[ln][i] + _gpln) * (alpha[i][j]);
                g4 += (mp[ln][i] + _gpln) * (alpha[i][j]*ud[ln][j] - 2.0*betta[i][j]*(mz[ln][i]-mPnts[j]));
            }
            ln = time_size-1;
            _gpln = 2.0 * R * gp(i,ln) * sgn( g0(i,ln) );
            g1 += 0.5*(mp[ln][i] + _gpln) * (uv[ln][j]-omega[i][j]);
            g2 += 0.5*(mp[ln][i] + _gpln) * (mz[ln][i]-mPnts[j])*(mz[ln][i]-mPnts[j]);
            g3 += 0.5*(mp[ln][i] + _gpln) * (alpha[i][j]);
            g4 += 0.5*(mp[ln][i] + _gpln) * (alpha[i][j]*ud[ln][j] - 2.0*betta[i][j]*(mz[ln][i]-mPnts[j]));

            g1 *= -time_step;
            g2 *= -time_step;
            g3 *= +time_step;
            g4 *= -time_step;

            g.at(0*heatSourceNumber*measrPointNumber + i*measrPointNumber + j) = g1 + 2.0*epsilon*(alpha[i][j] - alphaN[i][j]*no_norm);
            g.at(1*heatSourceNumber*measrPointNumber + i*measrPointNumber + j) = g2 + 2.0*epsilon*(betta[i][j] - bettaN[i][j]*no_norm);
            g.at(2*heatSourceNumber*measrPointNumber + i*measrPointNumber + j) = g3 + 2.0*epsilon*(omega[i][j] - omegaN[i][j]*no_norm);
#ifdef OPTIMIZE_ETA
            g.at(3*heatSourceNumber*measrPointNumber + j) += g4;
#endif
        }
#ifdef OPTIMIZE_ETA
        g.at(3*heatSourceNumber*measrPointNumber + j) += 2.0*epsilon*(mPnts[j]-mPntsN[j]*no_norm);
#endif
    }


#else
    for (size_t ln=0; ln<time_size; ln++)
    {
        double t = ln * time_step;
        DoubleVector _qNorm = qNorm1(t);
        for (size_t i=0; i<heatSourceNumber; i++)
        {
            g[i*time_size+ln] = -mp[ln][i] + 2.0*epsilon*(mq[ln][i] - _qNorm[i]*no_norm);
        }
    }
    //for (size_t i=0; i<heatSourceNumber; i++) { g[i*time_size] = 0.0; }
#endif
}

auto Functional::fx(const DoubleVector &x) const -> double
{
    const double _integral = integral(x);
    return _integral;// + epsilon * norm(x) + R * penalty(x);
}

auto Functional::integral(const DoubleVector &x) const -> double
{
    const_cast<Functional*>(this)->convertFromVector(x);
    forward.implicit_calculate_D1V1();
    const auto N = spaceDimensionX().size()-1;

    double sum = 0.0;

    sum += 0.5*(U[0]-V[0])*(U[0]-V[0]);
    for (unsigned int i=1; i<N; i++)
    {
        sum += (U[i]-V[i])*(U[i]-V[i]);
    }
    sum += 0.5*(U[N]-V[N])*(U[N]-V[N]);

    return sum*spaceDimensionX().step();
}

auto Functional::norm(const DoubleVector &x) const -> double
{
    const_cast<Functional*>(this)->convertFromVector(x);

    double _norm = 0.0;
#ifdef OPTIMIZE_Y
    for (size_t j=0; j<measrPointNumber; j++)
    {
        for (size_t i=0; i<heatSourceNumber; i++)
        {
            _norm += sqr(alpha[i][j] - alphaN[i][j]*no_norm);
            _norm += sqr(betta[i][j] - bettaN[i][j]*no_norm);
            _norm += sqr(omega[i][j] - omegaN[i][j]*no_norm);
        }
#ifdef OPTIMIZE_ETA
        _norm += sqr(mPnts[j] - mPntsN[j]*no_norm);
#endif
    }
#else
    const auto time_size = timeDimension().size();
    const auto time_step = timeDimension().step();

    size_t ln = 0;
    DoubleVector _qNorm = qNorm1(ln * time_step);
    for (size_t i=0; i<heatSourceNumber; i++) { _norm += 0.5*sqr(mq[ln][i]-_qNorm[i]*no_norm); }
    for (ln=1; ln<time_size-1; ln++)
    {
        _qNorm = qNorm1(ln * time_step);
        for (unsigned int i=0; i<heatSourceNumber; i++) { _norm += sqr(mq[ln][i]-_qNorm[i]*no_norm); }
    }
    ln = time_size-1;
    _qNorm = qNorm1(ln * time_step);
    for (size_t i=0; i<heatSourceNumber; i++) { _norm += 0.5*sqr(mq[ln][i]-_qNorm[i]*no_norm); }

    _norm *= time_step;
#endif
    return _norm;
}

auto Functional::penalty(const DoubleVector &x) const -> double
{
    double _penalty = 0.0;
    const auto time_size = timeDimension().size();
    const auto time_step = timeDimension().step();

    for (size_t i=0; i<heatSourceNumber; i++)
    {
        double pnlt = 0.0;
        pnlt += 0.5*sqr(gp(i, 0));
        for (size_t ln=1; ln<time_size-1; ln++) { pnlt += sqr(gp(i, ln)); }
        pnlt += 0.5*sqr(gp(i, time_size-1));
        pnlt *= time_step;

        _penalty += pnlt;
    }

    return _penalty;
}

/*virtual*/ void Functional::project(DoubleVector &x, size_t i)
{
    if (i < 3*heatSourceNumber*measrPointNumber) return;

    const auto n = 3*heatSourceNumber*measrPointNumber;
    if (i == (n+0)) { if (x[i] < 0.05) x[i] = 0.05; if (x[i] > 0.25) x[i] = 0.25; }
    if (i == (n+1)) { if (x[i] < 0.25) x[i] = 0.25; if (x[i] > 0.50) x[i] = 0.50; }
    if (i == (n+2)) { if (x[i] < 0.50) x[i] = 0.50; if (x[i] > 0.75) x[i] = 0.75; }
    if (i == (n+3)) { if (x[i] < 0.75) x[i] = 0.75; if (x[i] > 0.95) x[i] = 0.95; }

    //    if (i >= 3*heatSourceNumber*measrPointNumber)
    //    {
    //        if (x[i] < 0.05) x[i] = 0.05;
    //        if (x[i] > 0.95) x[i] = 0.95;
    //    }
}

/**************************************************************************************************************************/
