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
    functional.R = 0.0;
    functional.epsilon = 0.0;
    functional.error = 0.0;
    functional.withError = false;

#ifdef OPTIMIZE_Y
    size_t vector_size = 2*_measrPointNumber*_heatSourceNumber + _heatSourceNumber + _measrPointNumber;
#else
    const auto vector_size = time_size;
    const unsigned int start[] = { 0*vector_size+0, 1*vector_size+0 };
    const unsigned int finsh[] = { 1*vector_size-1, 2*vector_size-1 };
#endif


    /******************************************************************
     *                     Calculating V(x) function
     ******************************************************************/
    functional.V.clear();
    functional.V.resize(dimX_size, 5.0);
    IPrinter::printVector(functional._w, functional._p, functional.V, "V: ");
    /******************************************************************
     *                     Calculating V(x) function
     ******************************************************************/

#ifdef OPTIMIZE_Y
    DoubleVector x(functional.VCTR_1, vector_size);
    IPrinter::printSeperatorLine("x vector");
    functional.printVectorY(x);
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

#ifdef OPTIMIZE_Y

    DoubleVector g1(x.length());
    functional.gradient(x, g1);
    std::string msg = std::string("gradient vector: R:") + std::to_string(functional.R) + std::string(", epsilon: ") + std::to_string(functional.epsilon);
    IPrinter::printSeperatorLine(msg.data());
    functional.printVectorY(g1, true);

    DoubleVector g2(x.length());
    IPrinter::printSeperatorLine("numerical gradient vector");
#if defined (OPTIMIZE_ALPHA)
    IGradient::Gradient(&functional, 0.001, x, g2, static_cast<size_t>(0x00), static_cast<size_t>(0x07));
#endif
#if defined (OPTIMIZE_BETTA)
    IGradient::Gradient(&functional, 0.001, x, g2, static_cast<size_t>(0x08), static_cast<size_t>(0x0F));
#endif
#if defined (OPTIMIZE_OMEGA)
    IGradient::Gradient(&functional, 0.001, x, g2, static_cast<size_t>(0x10), static_cast<size_t>(0x11));
#endif
#if defined (OPTIMIZE_ETA)
    IGradient::Gradient(&functional, 0.001, x, g2, static_cast<size_t>(0x12), static_cast<size_t>(0x15));
#endif
    functional.printVectorY(g2, true);
#else
    functional.gradient(x, g);
    IPrinter::printVector(functional._w, functional._p, g.mid(start[0], finsh[0]).L2Normalize(), "g1:");
    IPrinter::printVector(functional._w, functional._p, g.mid(start[1], finsh[1]).L2Normalize(), "g2:");

    puts("---");
    IGradient::Gradient(&functional, 0.01, x, g1);
    IPrinter::printVector(functional._w, functional._p, g1.mid(start[0], finsh[0]).L2Normalize(), "g1:");
    IPrinter::printVector(functional._w, functional._p, g1.mid(start[1], finsh[1]).L2Normalize(), "g2:");
#endif

    IPrinter::printSeperatorLine("Starting optimization...");
    double step = 0.5;
    double epsl = 0.0001;
    functional.R = 0.0;
    functional.epsilon = 0.0;
    functional.no_norm = 0.0;
    functional.error = 0.0;
    functional.withError = false;
//    while (functional.epsilon > 0.000001)
    {
//        while (functional.R <= 100000.0)
        {
            GradientBasedMethod *gm;
            //gm = new ConjugateGradient; gm->setNormalize(false);
            gm = new SteepestDescentGradient; gm->setNormalize(true);
            gm->setFunction(&functional);
            gm->setGradient(&functional);
            gm->setPrinter(&functional);
            gm->setProjection(&functional);
            //gm.setGradientNormalizer(&prob);
            gm->setOptimalityTolerance(0.000000001);
            gm->setStepTolerance(0.000000001);
            gm->setFunctionTolerance(0.000000001);
            gm->setR1MinimizeEpsilon(step, epsl);
            //gm->setMaxIterationCount(100);
            gm->showExitMessage(false);
            gm->calculate(x);
            delete gm;

            //if (functional.R <= 100000.0) functional.R *= 10.0;
            //IPrinter::print(x, x.length(), functional._w, functional._p);
        }
//        functional.R = 1.0;
//        functional.epsilon *= 0.1;
    }
    IPrinter::printSeperatorLine("Optimization is finished.");
    IPrinter::printVector(functional._w, functional._p, functional.V, "V: ");

    IPrinter::printSeperatorLine();
    DoubleVector x1(functional.NORM_1, vector_size);
    x1= x;
    IPrinter::print(x1, x1.length(), functional._w, functional._p);
    functional.R = 0.0;
    functional.epsilon = 0.0;
    functional.no_norm = 1.0;
    functional.error = 0.03;
    functional.withError = true;
    //functional.setTimeDimension(Dimension(time_step, 0, static_cast<int>(5.0*(time_size-1))));
    functional.convertFromVector(x1);
    functional.forward.implicit_calculate_D1V1();

    IPrinter::printVector(functional._w, functional._p, functional.U, "U: ", functional.U.length());

//    for (size_t ln=1; ln<1.1*(time_size-1); ln++)
//    {
//        TimeNodePDE tn; tn.i = static_cast<unsigned int>(ln); tn.t = tn.i*time_step;
//        functional.setTimeDimension(Dimension(time_step, 0, static_cast<int>(ln)));
//        functional.convertFromVector(x1);

//        double t = ln*time_step;
//        double fx = functional.fx(x1);
//        double in = functional.integral(functional.U);
//        double nr = functional.norm(x1);
//        double pn = functional.penalty(x1);
//        double q1 = functional.mq.at(ln,0);
//        double q2 = functional.mq.at(ln,1);
//        double qMin1 = functional.qMin.at(ln,0);
//        double qMin2 = functional.qMin.at(ln,1);
//        double qMax1 = functional.qMax.at(ln,0);
//        double qMax2 = functional.qMax.at(ln,1);

//        DoubleVector _z = functional.z(tn);
//        printf("time: %6.4f fx: %16.10f in: %16.10f nr: %16.10f pn: %16.10f | q1: %16.10f %16.10f %16.10f | q2: %16.10f %16.10f %16.10f | %16.10f %16.10f | %0.4f %0.4f\n", t, fx, in, nr, pn,
//               qMin1, q1, qMax1,
//               qMin2, q2, qMax2,
//               functional.R, functional.epsilon,
//               _z[0], _z[1]);
//    }
}

/**************************************************************************************************************************/

auto Functional::print(unsigned int it, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientBasedMethod::MethodResult result) const -> void
{
#ifdef OPTIMIZE_Y
    printY(it, x, g, f, alpha, result);
#else
    printQ(it, x, g, f, alpha, result);
#endif
}

void CommonParameter::printVectorY(const DoubleVector &x, bool normolize) const
{
    const size_t offset = measrPointNumber*heatSourceNumber;
    size_t start = 0;

    if (normolize) {
        printf("a:");IPrinter::print(x.mid(start,start+offset-1).EuclideanNormalize(), offset, _w, _p);
    } else {
        printf("a:");IPrinter::print(x.mid(start,start+offset-1), offset, _w, _p);
    }
    start += offset;

    if (normolize) {
        printf("b:");IPrinter::print(x.mid(start,start+offset-1).EuclideanNormalize(), offset, _w, _p);
    } else {
        printf("b:");IPrinter::print(x.mid(start,start+offset-1), offset, _w, _p);
    }
    start += offset;

    if (normolize) {
        printf("o:");IPrinter::print(x.mid(start,start+heatSourceNumber-1).EuclideanNormalize(), heatSourceNumber, _w, _p);
    } else {
        printf("o:");IPrinter::print(x.mid(start,start+heatSourceNumber-1), heatSourceNumber, _w, _p);
    }
    start += heatSourceNumber;

    if (normolize) {
        printf("e:");IPrinter::print(x.mid(start,start+measrPointNumber-1).EuclideanNormalize(), measrPointNumber, _w, _p);
    } else {
        printf("e:");IPrinter::print(x.mid(start,start+measrPointNumber-1), measrPointNumber, _w, _p);
    }
}

auto CommonParameter::setTimeDimension(const Dimension &timeDimension) -> void
{
    const auto last_size = _timeDimension.size();

    _timeDimension = timeDimension;
    const auto time_size = _timeDimension.size();
    const auto time_step = _timeDimension.step();

    //    if (mq != nullptr)
    //    {
    //        for (size_t ln=0; ln<last_size; ln++) { mq[ln].clear(); }
    //        delete [] mq;
    //        mq = nullptr;
    //    }
    //    mq = new DoubleVector[time_size];
    //    for (size_t ln=0; ln<time_size; ln++) { mq[ln].resize(heatSourceNumber); }
    if (!mq.empty()) { mq.clear(); } mq.resize(time_size, heatSourceNumber);


    //    if (mz != nullptr)
    //    {
    //        for (size_t ln=0; ln<last_size; ln++) { mz[ln].clear(); }
    //        delete [] mz;
    //        mz = nullptr;
    //    }
    //    mz = new DoubleVector[time_size];
    //    for (size_t ln=0; ln<time_size; ln++)
    //    {
    //        mz[ln].resize(heatSourceNumber);

    //        double t = ln * time_step;
    //        mz[ln][0] = 0.05 + 0.90*fabs(sin(M_PI*t));
    //        mz[ln][1] = 0.95 - 0.90*fabs(sin(M_PI*t));
    //        //mz[ln][0] = 0.30;
    //        //mz[ln][1] = 0.70;
    //    }
    if (!mz.empty()) { mz.clear(); } mz.resize(time_size, heatSourceNumber);
    for (size_t ln=0; ln<time_size; ln++)
    {
        double t = ln * time_step;
        mz.at(ln,0) = 0.05 + 0.90*fabs(sin(M_PI*t));
        mz.at(ln,1) = 0.95 - 0.90*fabs(sin(M_PI*t));
    }


    //    if (mp != nullptr)
    //    {
    //        for (size_t ln=0; ln<last_size; ln++) { mp[ln].clear(); }
    //        delete [] mp;
    //        mp = nullptr;
    //    }
    //    mp = new DoubleVector[time_size];
    //    for (size_t ln=0; ln<time_size; ln++) { mp[ln].resize(heatSourceNumber); }
    if (!mp.empty()) { mp.clear(); } mp.resize(time_size, heatSourceNumber);


#ifdef OPTIMIZE_Y
    //    if (uv != nullptr)
    //    {
    //        for (size_t ln=0; ln<last_size; ln++) { uv[ln].clear(); ud[ln].clear(); }
    //        delete [] uv;
    //        delete [] ud;
    //        uv = nullptr;
    //        ud = nullptr;
    //    }
    //    uv = new DoubleVector[time_size];
    //    ud = new DoubleVector[time_size];
    //    for (size_t ln=0; ln<time_size; ln++) { uv[ln].resize(measrPointNumber); ud[ln].resize(measrPointNumber); }
    if (!uv.empty()) { uv.clear(); } uv.resize(time_size, measrPointNumber);
    if (!ud.empty()) { ud.clear(); } ud.resize(time_size, measrPointNumber);
#endif

    //    if (qMin != nullptr && qMax != nullptr)
    //    {
    //        for (size_t ln=0; ln<last_size; ln++) { qMin[ln].clear(); qMax[ln].clear();  }
    //        delete [] qMin;
    //        delete [] qMax;
    //        qMin = nullptr;
    //        qMax = nullptr;
    //    }
    //    qMin = new DoubleVector[time_size];
    //    qMax = new DoubleVector[time_size];
    //    for (size_t ln=0; ln<time_size; ln++)
    //    {
    //        qMin[ln] <<  +0.00 <<  +0.00;
    //        qMax[ln] << +10.00 << +10.00;
    //    }

    if (!qMin.empty()) { qMin.clear(); } qMin.resize(time_size, heatSourceNumber);
    if (!qMax.empty()) { qMax.clear(); } qMax.resize(time_size, heatSourceNumber);
    for (size_t ln=0; ln<time_size; ln++)
    {
        qMin.at(ln,0) =  0.00; qMin.at(ln,1) =  0.0;
        qMax.at(ln,0) = 10.00; qMax.at(ln,1) = 10.0;
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
    omega.resize(heatSourceNumber, 0.0);
    mPnts.resize(measrPointNumber, 0.0);

    alphaN.resize(heatSourceNumber, measrPointNumber, 0.0);
    bettaN.resize(heatSourceNumber, measrPointNumber, 0.0);
    omegaN.resize(heatSourceNumber, 0.0);
    mPntsN.resize(measrPointNumber, 0.0);
#endif
}

auto CommonParameter::q(const TimeNodePDE &tn) const -> DoubleVector
{
    DoubleVector _q; _q << mq.at(tn.i, 0) << mq.at(tn.i, 1); return _q;
    //return mq[tn.i];
}

auto CommonParameter::z(const TimeNodePDE &tn) const -> DoubleVector
{
    DoubleVector _z; _z << mz.at(tn.i, 0) << mz.at(tn.i, 1); return _z;
    //return mz[tn.i];
}

auto CommonParameter::g0(size_t i, size_t ln) const -> double
{
    double _qMax = qMax.at(ln,i);
    double _qMin = qMin.at(ln,i);

    return 0.5*(_qMax + _qMin) - mq[ln][i];
}

auto CommonParameter::gi(size_t i, size_t ln) const -> double
{
    double _qMax = qMax.at(ln,i);
    double _qMin = qMin.at(ln,i);

    double _g0 = g0(i, ln);
    return fabs(_g0) - 0.5*(_qMax - _qMin);
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
    return common->_initialTemperature;
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
    const size_t heatSourceNumber = common->heatSourceNumber;
    for (size_t i=0; i<heatSourceNumber; i++)
    {
        fx += _q[i] * DeltaFunction::gaussian(sn.x, _z[i], sigma);
    }
    return fx - thermalConvection()*common->theta();
}

auto HeatEquationIBVP::layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const -> void
{
    if (static_cast<int>(tn.i) == timeDimension().max()) { common->U = u; }

//    if (common->withError) {
//        printf("%6.3f %14.10f u: ", tn.t, common->integral(u));
//        IPrinter::printVector(u);
//    }

#ifdef OPTIMIZE_Y
    const size_t heatSourceNumber = common->heatSourceNumber;
    const size_t measrPointNumber = common->measrPointNumber;
    const auto hx = spaceDimensionX().step();
    const auto Nx = spaceDimensionX().size();

    for (size_t j=0; j<measrPointNumber; j++)
    {
        //const double mp = common->mPnts[j];
        //common->uv[tn.i][j] = DeltaFunction::lumpedPoint4(u, mp, hx, Nx, common->ud[tn.i][j]);
        common->uv.at(tn.i,j) = DeltaFunction::lumpedPoint4(u, common->mPnts.at(j), hx, Nx, common->ud.at(tn.i,j));
    }

    if (tn.i < static_cast<unsigned int>(timeDimension().max()))
    {
        for (size_t i=0; i<heatSourceNumber; i++)
        {
            //common->mq[tn.i+1].at(i) = 0.0;
            //const auto zi = common->mz[tn.i].at(i);

            const double mzlni = common->mz.at(tn.i, i);
            double &mqlni = common->mq.at(tn.i+1,i);

            mqlni = 0.0;

            for (size_t j=0; j<measrPointNumber; j++)
            {
                const double mp = common->mPnts[j];
                double err = rand() % 2==0 ? 1.0+common->error : 1.0-common->error;
                const double uv = common->uv.at(tn.i,j) * err;
                const double alpha = common->alpha.at(i,j);
                const double betta = common->betta.at(i,j);
                const double omega = common->omega.at(i);
                mqlni += alpha*uv - omega + betta * ( mzlni - mp ) * (mzlni - mp);
            }
        }
    }

    if (tn.i == 0)
    {
        //for (size_t i=0; i<heatSourceNumber; i++) { common->mq[0].at(i) = common->mq[1].at(i); }
        for (size_t i=0; i<heatSourceNumber; i++) { common->mq.at(0,i) = common->mq.at(1,i); }
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
    //printf(">>> %d %d rows: %d cols: %d\n", sn.i, tn.i, common->mp.rows(), common->mp.cols());
#ifdef OPTIMIZE_Y
    const auto sigma = 0.01;
    for (size_t j=0; j<common->measrPointNumber; j++)
    {
        double fxj = 0.0;
        for (size_t i=0; i<common->heatSourceNumber; i++)
        {
            fxj += common->alpha.at(i,j) * ( common->mp.at(tn.i, i)
                                             + 2.0 * common->R * common->gp(i,tn.i) * sgn( common->g0(i,tn.i) ));
        }
        fx -= fxj * DeltaFunction::gaussian(sn.x, common->mPnts.at(j), sigma);
    }
#endif
    //printf("<<< %d %d rows: %d cols: %d\n", sn.i, tn.i, common->mp.rows(), common->mp.cols());
    return fx;
}

auto HeatEquationFBVP::HeatEquationFBVP::layerInfo(const DoubleVector &psi, const TimeNodePDE &tn) const -> void
{
    const auto hx = spaceDimensionX().step();
    const auto Nx = spaceDimensionX().size();
    const auto heatSourceNumber = common->heatSourceNumber;
    DoubleVector _z = common->z(tn);
    //    common->mp[tn.i][0] = DeltaFunction::lumpedPoint4(psi, _z[0], hx, Nx);
    //    common->mp[tn.i][1] = DeltaFunction::lumpedPoint4(psi, _z[1], hx, Nx);

    if (tn.i > timeDimension().min())
    {
        for (size_t i=0; i<heatSourceNumber; i++)
        {
            common->mp.at(tn.i-1,i) = DeltaFunction::lumpedPoint4(psi, _z[i], hx, Nx);
        }
    }

    if (tn.i == timeDimension().max())
    {
        for (size_t i=0; i<heatSourceNumber; i++) { common->mp.at(tn.i,i) = common->mp.at(tn.i-1,i); }
    }
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

            alphaN.at(i,j) = NORM_1[(0*heatSourceNumber*measrPointNumber) + i*measrPointNumber + j];
            bettaN.at(i,j) = NORM_1[(1*heatSourceNumber*measrPointNumber) + i*measrPointNumber + j];

        }
    }

    for (unsigned int i=0; i<heatSourceNumber; i++)
    {
        omega.at(i)  = x.at((2*heatSourceNumber*measrPointNumber) + i);
        omegaN.at(i) = NORM_1[(2*heatSourceNumber*measrPointNumber) + i];
    }

    for (unsigned int j=0; j<measrPointNumber; j++)
    {
        mPnts.at(j) = x.at((2*heatSourceNumber*measrPointNumber) + heatSourceNumber + j);
        mPntsN.at(j) = NORM_1[(2*heatSourceNumber*measrPointNumber) + heatSourceNumber + j];
    }

    //for (unsigned int i=0; i<heatSourceNumber; i++) { mq[0].at(i) = 0.0; }
    for (unsigned int i=0; i<heatSourceNumber; i++) { mq.at(0,i) = 0.0; }

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

//auto CommonParameter::convertToVector(DoubleVector &x) const -> void
//{
//    const auto time_size = _timeDimension.size();

//    x.clear();
//    x.resize(heatSourceNumber*time_size);

//    if (mq != nullptr)
//    {
//        for (size_t ln=0; ln<time_size; ln++)
//        {
//            for (size_t i=0; i<heatSourceNumber; i++)
//            {
//                x[i*time_size+ln] = mq[ln][i];
//            }
//        }
//    }
//}


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

/*virtual*/ auto Functional::gradient(const DoubleVector &x, DoubleVector &g) const -> void
{
#ifdef OPTIMIZE_Y
    gradientY(x, g);
#else
    gradientQ(x, g);
#endif
}

/*virtual*/ auto Functional::fx(const DoubleVector &x) const -> double
{
    double _fx = 0.0;
    const_cast<Functional*>(this)->convertFromVector(x);

    for (unsigned int k1=0; k1<initialTemperature_nmb; k1++)
    {
        const_cast<Functional*>(this)->_initialTemperature = initialTemperature_vls[k1];

        for (unsigned int k2=0; k2<thetaTemperature_nmb; k2++)
        {
            const_cast<Functional*>(this)->_theta = thetaTemperature_vls[k2];

            double rho = initialTemperature_rho[k1] * thetaTemperature_rho[k2];

            forward.implicit_calculate_D1V1();

            const double _integral = integral(U);
            double __fx = (_integral + epsilon * norm(x) + R * penalty(x)) * rho;

            _fx += __fx;
        }
    }
    return _fx;
}

auto CommonParameter::integral(const DoubleVector &u) const -> double
{
    const auto N = spaceDimensionX().size()-1;

    double sum = 0.0;
    sum += 0.5*(u[0]-V[0])*(u[0]-V[0]);
    for (unsigned int i=1; i<N; i++)
    {
        sum += (u[i]-V[i])*(u[i]-V[i]);
    }
    sum += 0.5*(u[N]-V[N])*(u[N]-V[N]);

    return sum*spaceDimensionX().step();
}

auto Functional::norm(const DoubleVector &x) const -> double
{
#ifdef OPTIMIZE_Y
    return normY(x, DoubleVector(x.length(), 0.0));
#else
    return normY(x, DoubleVector(x.length(), 0.0));
#endif
}

auto Functional::penalty(const DoubleVector &x) const -> double
{
#ifdef OPTIMIZE_Y
    return penaltyY(x);
#else
    return penaltyQ(x);
#endif

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

auto Functional::normY(const DoubleVector &x, const DoubleVector &/*n*/) const -> double
{
    const_cast<Functional*>(this)->convertFromVector(x);

    double norm = 0.0;
    for (size_t j=0; j<measrPointNumber; j++)
    {
        for (size_t i=0; i<heatSourceNumber; i++)
        {
            norm += sqr(alpha.at(i,j) - alphaN.at(i,j)*no_norm);
            norm += sqr(betta.at(i,j) - bettaN.at(i,j)*no_norm);
        }
    }

    for (size_t i=0; i<heatSourceNumber; i++) { norm += sqr(omega.at(i) - omegaN.at(i)*no_norm); }

    for (size_t j=0; j<measrPointNumber; j++) { norm += sqr(mPnts.at(j) - mPntsN.at(j)*no_norm); }

    return norm;
}

auto Functional::normQ(const DoubleVector &x, const DoubleVector &/*n*/) const -> double
{
    const_cast<Functional*>(this)->convertFromVector(x);

    const auto time_size = timeDimension().size();
    const auto time_step = timeDimension().step();

    double norm = 0.0;
    size_t ln = 0;
    DoubleVector _qNorm = qNorm1(ln * time_step);
    for (size_t i=0; i<heatSourceNumber; i++) { norm += 0.5*sqr(mq[ln][i]-_qNorm[i]*no_norm); }
    for (ln=1; ln<time_size-1; ln++)
    {
        _qNorm = qNorm1(ln * time_step);
        for (unsigned int i=0; i<heatSourceNumber; i++) { norm += sqr(mq[ln][i]-_qNorm[i]*no_norm); }
    }
    ln = time_size-1;
    _qNorm = qNorm1(ln * time_step);
    for (size_t i=0; i<heatSourceNumber; i++) { norm += 0.5*sqr(mq[ln][i]-_qNorm[i]*no_norm); }

    norm *= time_step;
    return norm;
}

auto Functional::penaltyY(const DoubleVector &/*x*/) const -> double
{
    const auto time_size = timeDimension().size()-1;
    const auto time_step = timeDimension().step();

    double penalty = 0.0;
    for (size_t i=0; i<heatSourceNumber; i++)
    {
        double penaltyi = 0.0;
        penaltyi += 0.5 * sqr(gp(i, 0));
        for (size_t ln=1; ln<time_size; ln++) { penaltyi += sqr( gp(i, ln) ); }
        penaltyi += 0.5 * sqr( gp(i, time_size) );
        penaltyi *= time_step;
        penalty += penaltyi;
    }

    return penalty;
}

auto Functional::penaltyQ(const DoubleVector &/*x*/) const -> double
{
    return 0.0;
}

auto Functional::gradientY(const DoubleVector &x, DoubleVector &g) const -> void
{
    const auto time_size = timeDimension().size();
    const auto time_step = timeDimension().step();
    const_cast<Functional*>(this)->convertFromVector(x);

    g.resize(x.length(), 0.0);

    for (unsigned int k1=0; k1<initialTemperature_nmb; k1++)
    {
        const_cast<Functional*>(this)->_initialTemperature = initialTemperature_vls[k1];

        for (unsigned int k2=0; k2<thetaTemperature_nmb; k2++)
        {
            const_cast<Functional*>(this)->_theta = thetaTemperature_vls[k2];

            double rho = initialTemperature_rho[k1] * thetaTemperature_rho[k2];

            forward.implicit_calculate_D1V1();
            backward.implicit_calculate_D1V1();

            for (size_t j=0; j<measrPointNumber; j++)
            {
                const auto mPntsj = mPnts.at(j);
                const auto mPntsNj = mPntsN.at(j)*no_norm;

                for (size_t i=0; i<heatSourceNumber; i++)
                {
                    const auto alpha_ij = alpha.at(i,j);
                    const auto betta_ij = betta.at(i,j);

                    const auto alphaN_ij = alphaN.at(i,j)*no_norm;
                    const auto bettaN_ij = bettaN.at(i,j)*no_norm;

                    double g1 = 0.0;
                    double g2 = 0.0;
                    double g4 = 0.0;

                    size_t ln = 0;
                    double _mpln = mp.at(ln,i) + 2.0*R*gp(i,ln)*sgn( g0(i,ln) );
                    double _uvlnj = uv.at(ln,j);
                    double _mzlni = mz.at(ln,i);
                    double _udlnj = ud.at(ln,j);

                    g1 += 0.5*_mpln * (_uvlnj /*- omega_ij*/);
                    g2 += 0.5*_mpln * (_mzlni - mPntsj)*(_mzlni - mPntsj);
                    g4 += 0.5*_mpln * (alpha_ij*_udlnj - 2.0*betta_ij*(_mzlni - mPntsj));
                    for (ln=1; ln<time_size-1; ln++)
                    {
                        _mpln = mp.at(ln,i) + 2.0*R*gp(i,ln)*sgn( g0(i,ln) );
                        _uvlnj = uv.at(ln,j);
                        _mzlni = mz.at(ln,i);
                        _udlnj = ud.at(ln,j);

                        g1 += _mpln * (_uvlnj /*- omega_ij*/);
                        g2 += _mpln * (_mzlni - mPntsj)*(_mzlni - mPntsj);
                        g4 += _mpln * (alpha_ij*_udlnj - 2.0*betta_ij*(_mzlni - mPntsj));
                    }
                    ln = time_size-1;
                    _mpln = mp.at(ln,i) + 2.0*R*gp(i,ln)*sgn( g0(i,ln) );
                    _uvlnj = uv.at(ln,j);
                    _mzlni = mz.at(ln,i);
                    _udlnj = ud.at(ln,j);

                    g1 += 0.5*_mpln * (_uvlnj /*- omega_ij*/);
                    g2 += 0.5*_mpln * (_mzlni - mPntsj)*(_mzlni - mPntsj);
                    g4 += 0.5*_mpln * (alpha_ij*_udlnj - 2.0*betta_ij*(_mzlni - mPntsj));

                    g1 *= -time_step;
                    g2 *= -time_step;
                    g4 *= -time_step;

#ifdef OPTIMIZE_ALPHA
                    g.at(0*heatSourceNumber*measrPointNumber + i*measrPointNumber + j) += (g1 + 2.0*epsilon*(alpha_ij - alphaN_ij)) * rho;
#endif
#ifdef OPTIMIZE_BETTA
                    g.at(1*heatSourceNumber*measrPointNumber + i*measrPointNumber + j) += (g2 + 2.0*epsilon*(betta_ij - bettaN_ij)) * rho;
#endif
#ifdef OPTIMIZE_ETA
                    g.at(2*heatSourceNumber*measrPointNumber + heatSourceNumber + j) += g4 * rho;
#endif
                }
#ifdef OPTIMIZE_ETA
                g.at(2*heatSourceNumber*measrPointNumber + heatSourceNumber + j) += (2.0*epsilon*(mPntsj - mPntsNj)) * rho;
#endif
            }

#ifdef OPTIMIZE_OMEGA
            for (size_t i=0; i<heatSourceNumber; i++)
            {
                const auto omega_i = omega.at(i);
                const auto omegaN_i = omegaN.at(i)*no_norm;

                double g3 = 0.0;

                size_t ln = 0;
                double _mpln = mp.at(ln,i) + 2.0*R*gp(i,ln)*sgn( g0(i,ln) );
                g3 += 0.5*_mpln;
                for (ln=1; ln<time_size-1; ln++)
                {
                    _mpln = mp.at(ln,i) + 2.0*R*gp(i,ln)*sgn( g0(i,ln) );
                    g3 += _mpln;
                }
                ln = time_size-1;
                _mpln = mp.at(ln,i) + 2.0*R*gp(i,ln)*sgn( g0(i,ln) );
                g3 += 0.5*_mpln;
                g3 *= +time_step;

                g.at(2*heatSourceNumber*measrPointNumber + i) += (g3 + 2.0*epsilon*(omega_i - omegaN_i)) * rho;
            }
#endif

        }
    }
}

auto Functional::gradientQ(const DoubleVector &x, DoubleVector &g) const -> void
{
    const auto time_size = timeDimension().size();
    const auto time_step = timeDimension().step();
    const_cast<Functional*>(this)->convertFromVector(x);

    forward.implicit_calculate_D1V1();
    backward.implicit_calculate_D1V1();

    g.resize(x.length(), 0.0);

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
}

auto Functional::printY(unsigned int it, const DoubleVector &x, const DoubleVector &/*g*/, double /*f*/, double /*alpha*/, GradientBasedMethod::MethodResult result) const -> void
{
    //const size_t time_size = timeDimension().size();
    //const auto vector_size = time_size;
    //IPrinter::printSeperatorLine();

    const auto _fx = fx(x);
    const auto _integral = integral(U);
    const auto _norm = norm(x);
    const auto _penalty = penalty(x);
    printf_s("I[%4d] fx:%10.8f int: %8.5f nrm: %8.5f penalty: %8.5f eps: %8.6f R:%8.1f ", it, _fx, _integral, _norm, _penalty, epsilon, R);
    IPrinter::print(x, x.length(), _w, _p);
    //    IPrinter::print(g, g.length(), _w, _p);
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
    //    IPrinter::printVector(_w, _p, U, "U:");
    //    IPrinter::printVector(_w, _p, V, "V:");
}

auto Functional::printQ(unsigned int it, const DoubleVector &x, const DoubleVector &g, double /*f*/, double /*alpha*/, GradientBasedMethod::MethodResult /*result*/) const -> void
{
    const size_t time_size = timeDimension().size();
    const auto vector_size = time_size;
    //IPrinter::printSeperatorLine();

    const auto _fx = fx(x);
    const auto _integral = integral(x);
    const auto _norm = norm(x);
    const auto _penalty = penalty(x);
    printf_s("I[%4d] fx:%10.8f int: %8.5f nrm: %8.5f penalty: %8.5f eps: %8.6f R:%8.1f ", it, _fx, _integral, _norm, _penalty, epsilon, R);
    IPrinter::print(x, x.length(), _w, _p);
    //    IPrinter::print(g, g.length(), _w, _p);

    IPrinter::printVector(_w, _p, x.mid(0*vector_size, 1*vector_size-1), "q1");
    IPrinter::printVector(_w, _p, x.mid(1*vector_size, 2*vector_size-1), "q2");
    IPrinter::printVector(_w, _p, g.mid(0*vector_size, 1*vector_size-1), "g1");
    IPrinter::printVector(_w, _p, g.mid(1*vector_size, 2*vector_size-1), "g2");
    //    IPrinter::printVector(_w, _p, U, "U:");
    //    IPrinter::printVector(_w, _p, V, "V:");
}


/**************************************************************************************************************************/
