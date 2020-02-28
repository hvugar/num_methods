#include "wave_equation_ibvp.h"

//#define __DIRICHLET_LEFT_
//#define __DIRICHLET_RIGHT_
//#define __NEUMANN__
#define __ROBIN_LEFT__
#define __ROBIN_RIGHT__
#define x2_t2

const double a = +0.2;//+1.2;
const double b = +0.0;//-0.5;
const double c = +0.0;//+0.3;
const double d = +0.0;//+0.4;

void WaveEquationIBVP::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
#ifdef USE_LIB_IMAGING
    QGuiApplication app(argc, argv);
#endif

    WaveEquationIBVP w;
    w.setTimeDimension(Dimension(0.01, 0, 10000));
    w.setSpaceDimensionX(Dimension(0.01, 100, 200));
    //w.setSpaceDimensionY(Dimension(0.01, 0, 100));

    w.setWaveSpeed(a);//a
    w.setWaveDissipation(d);//d
    w.setUnknownB(b);//b
    w.setRestoration(c);//c

    Benchmark bm;
    bm.tick();
    w.implicit_calculate_D1V1();
    bm.tock();
    bm.printDuration();
}

double WaveEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition condition) const
{
    TimeNodePDE tn; tn.t = 0.0;
    if (condition == InitialCondition::InitialValue) return U(sn, tn);
    if (condition == InitialCondition::InitialFirstDerivative) return Udt1(sn, tn);
    throw std::exception();
}

double WaveEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{
    if (sn.i == spaceDimensionX().min())
    {
#if defined(__DIRICHLET_LEFT_)
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
        return U(sn, tn)*(condition.alpha()/condition.gamma());
#endif
#if defined(__ROBIN_LEFT__)
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, +4.0, +2.0, +1.0);
        return (condition.alpha()*U(sn, tn)+condition.beta()*Udx1(sn, tn))/condition.gamma();
#endif
    }

    if (sn.i == spaceDimensionX().max())
    {
#if defined(__DIRICHLET_RIGHT_)
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
        return U(sn, tn)*(condition.alpha()/condition.gamma());
#endif
#if defined(__ROBIN_RIGHT__)
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, +4.0, +2.0, +1.0);
        return (condition.alpha()*U(sn, tn)+condition.beta()*Udx1(sn, tn))/condition.gamma();
#endif
    }

    throw std::exception();
}

double WaveEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    const double a = waveSpeed();//ws
    const double b = unknownB();//wd
    const double c = restoration();//rs
    const double d = waveDissipation();//wd

#if defined (x1_t1) || defined (x2_t1) || defined (x3_t1) || defined (x1_t2) || defined (x2_t2) || defined (x3_t2) || defined (x1_t3) || defined (x2_t3) || defined (x3_t3)
    return Udt2(sn,tn)+d*Udt1(sn,tn)-a*a*Udx2(sn,tn)-b*Udx1(sn,tn)-c*U(sn,tn);
#endif
}

void WaveEquationIBVP::layerInfo(const DoubleVector& u, const TimeNodePDE& tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);

    if (tn.i % (timeDimension().size() / 10) == 0)
    IPrinter::printVector(16, 8, u);

//    if (tn.i==200)
//    {
//        double norm = 0.0;
//        double max = 0.0;
//        TimeNodePDE tn; tn.t = 1.0;
//        SpaceNodePDE sn;
//        for (unsigned int i=0; i<=100; i++)
//        {
//            sn.x = i*0.01;
//            double k = 1.0; if (i==0 || i== 100) k = 0.5;
//            norm += 0.01*k*(u[i]-U(sn, tn))*(u[i]-U(sn, tn));

//            if (max < fabs(u[i]-U(sn, tn))) max = fabs(u[i]-U(sn, tn));
//        }
//        printf("norm: %.10f max: %.10f\n", sqrt(norm), max);
//    }

    //    if (tn.i == 0 || tn.i == 1 || tn.i == 2 || tn.i == 2000) IPrinter::printVector(u);

#ifdef USE_LIB_IMAGING1
    QPixmap pxm;
    visualString(u, -0.30, +0.30, 500, 200, pxm, Qt::yellow, Qt::red, QString("d:/img/1/%1.png").arg(tn.i, 4, 10, QChar('0')));

    //    if (method == 1)
    //    {
    //        //visualString(u, -0.32, +0.32, 1000, 200, pxm, Qt::transparent, Qt::red, QString("d:/img/1/%1.png").arg(ln));
    //        visualString(u, -0.60, +0.60, 1000, 200, pxm, Qt::transparent, Qt::red, QString("d:/img/1/%1.png").arg(tn.i));
    //    }
    //    if (method == 2)
    //    {
    //        //visualString(u, -0.32, +0.32, 1000, 200, pxm, Qt::transparent, Qt::green, QString("d:/img/2/%1.png").arg(ln));
    //        visualString(u, -0.60, +0.60, 1000, 200, pxm, Qt::transparent, Qt::green, QString("d:/img/2/%1.png").arg(tn.i));
    //    }
    //    if (method == 3)
    //    {
    //        //visualString(u, -0.32, +0.32, 1000, 200, pxm, Qt::transparent, Qt::blue, QString("d:/img/3/%1.png").arg(ln));
    //        visualString(u, -0.60, +0.60, 1000, 200, pxm, Qt::transparent, Qt::blue, QString("d:/img/3/%1.png").arg(tn.i));

    //        QPixmap pxm(1000, 600);
    //        pxm.fill(Qt::yellow);
    //        QPainter painter(&pxm);
    //        QPixmap img1(QString("d:/img/1/%1.png").arg(tn.i), "PNG");
    //        QPixmap img2(QString("d:/img/2/%1.png").arg(tn.i), "PNG");
    //        QPixmap img3(QString("d:/img/3/%1.png").arg(tn.i), "PNG");
    //        painter.drawPixmap(5, 200, img1);
    //        painter.drawPixmap(5, 200, img2);
    //        painter.drawPixmap(2, 200, img3);
    //        pxm.save(QString("d:/img/0/%1.png").arg(tn.i));
    //    }

    //    //    QPixmap pxm(1000, 600);
    //    //    pxm.fill(Qt::yellow);
    //    //    QPainter painter(&pxm);
    //    //    QPixmap img1(QString("d:/img/1/%1.png").arg(ln), "PNG");
    //    //    QPixmap latex(QString("d:/img/latex_1_200.png"), "PNG");
    //    //    painter.drawPixmap(5, 200, img1);
    //    //    painter.drawPixmap(1000-540, 10, latex);
    //    //    pxm.save(QString("d:/img/4/%1.png").arg(ln));

    //    //printf("\rSaved: %d\n", ln);

    //    //IPrinter::printVector(u);
    //    //double w = u.min(); if (w < 0.0) {} else { w = u.max(); } printf("%14.10f\n", w);
    //    return;

    //    const static double ht = timeDimension().step();
    //    const static unsigned int N = static_cast<unsigned int> ( spaceDimension(Dimension::DimensionX).size() );

    //    DoubleVector &_u0 = const_cast<WaveEquationIBVP*>(this)->vu0;
    //    DoubleVector &_u1 = const_cast<WaveEquationIBVP*>(this)->vu1;
    //    DoubleVector &_u2 = const_cast<WaveEquationIBVP*>(this)->vu2;
    //    DoubleVector &_ut = const_cast<WaveEquationIBVP*>(this)->vut;

    //    _ut.resize(N+1, 0.0);

    //    if (tn.i == 0) { const_cast<WaveEquationIBVP*>(this)->vu0 = u; }
    //    if (tn.i == 1) { const_cast<WaveEquationIBVP*>(this)->vu1 = u; }
    //    if (tn.i == 2) { const_cast<WaveEquationIBVP*>(this)->vu2 = u; }

    //    if (tn.i >= 3)
    //    {
    //        _u0 = _u1; _u1 = _u2; _u2 = u;
    //    }

    //    if (tn.i >= 2)
    //    {
    //        for (unsigned int n=0; n<=N; n++)
    //        {
    //            _ut[n] = (_u2[n] - _u0[n])/(2.0*ht);
    //        }
    //    }

    //    //if (ln%20000==0) IPrinter::printVector(u, nullptr, u.length());
    //    //double k = 9.8681840001001186485325310599663;
    //    //double k = 1.0/523.3314251378401924542377492593;
    //    double k = 23.788030395009524606979129597455;
    //    if (tn.i > 1) printf("%5d %.10f %.10f %.10f\n", tn.i, k*integralUP(u), integralUK(u), k*integralUP(u)+integralUK(u));
#endif
}

void WaveEquationIBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    if (tn.i==200 || tn.i==199 || tn.i==198 || /*tn.i==397 || tn.i==396 ||*/ tn.i==2 || tn.i==1 || tn.i==0)
    {
        IPrinter::printMatrix(u);
        IPrinter::printSeperatorLine();
    }
    return;


    //if (ln%1==0) { double w = u.min(); if (w < 0.0) {} else { w = u.max(); } printf("%14.10f\n", w); }

    //    if (ln%2==0)
    {
        const double ht = _timeDimension.step();
        const unsigned int N = static_cast<unsigned int> ( _spaceDimensionX.size() );
        const unsigned int M = static_cast<unsigned int> ( _spaceDimensionY.size() );
        const double hx = _spaceDimensionX.step();
        const double hy = _spaceDimensionY.step();

        DoubleMatrix &_u0 = const_cast<WaveEquationIBVP*>(this)->mu0;
        DoubleMatrix &_u1 = const_cast<WaveEquationIBVP*>(this)->mu1;
        DoubleMatrix &_u2 = const_cast<WaveEquationIBVP*>(this)->mu2;

        DoubleMatrix &_ut = const_cast<WaveEquationIBVP*>(this)->mut;

        if (tn.i == 0) { _u0 = u; }
        if (tn.i == 1) { _u1 = u; }
        if (tn.i == 2) { _u2 = u; }
        if (tn.i >= 3) { _u0 = _u1; _u1 = _u2; _u2 = u; }

        if ( tn.i == 0 )
        {
            _ut.resize(M+1, N+1, 0.0);
            _u1.resize(M+1, N+1, 0.0);

            for (unsigned int m=0; m<=M; m++)
            {
                for (unsigned int n=0; n<=N; n++)
                {
                    SpaceNodePDE sn;
                    sn.i = static_cast<int>(n); sn.x = n*hx;
                    sn.j = static_cast<int>(m); sn.y = m*hy;
                    _u1[m][n] = initial(sn, InitialCondition::InitialValue);
                    _ut[m][n] = initial(sn, InitialCondition::InitialFirstDerivative);
                }
            }
            //printf("%5d %.10f %.10f %.10f\n", ln, integralUP(u), integralUK(u), integralUP(u)+integralUK(u));
        }

        if ( tn.i == 1 )
        {
            for (unsigned int m=0; m<=M; m++)
            {
                for (unsigned int n=0; n<=N; n++)
                {
                    _ut[m][n] = (_u1[m][n] - _u0[m][n])/(ht);
                }
            }
            //printf("%5d %.10f %.10f %.10f\n", ln, integralUP(u), integralUK(u), integralUP(u)+integralUK(u));
        }

        if ( tn.i >= 2)
        {
            for (unsigned int m=0; m<=M; m++)
            {
                for (unsigned int n=0; n<=N; n++)
                {
                    _ut[m][n] = (_u2[m][n] - _u0[m][n])/(2.0*ht);
                }
            }
        }

        //        if (ln==0)
        {
            //        //if (ln%20000==0) IPrinter::printVector(u, nullptr, u.length());
            double k = 1.0;
            //            double k = 9.8681840001001186485325310599663*1.9998766945708295048259473561531;
            //double k = 19.735111566520971090368615144868;
            //        double k = 691.74162213081364457511246502072;
            //        double k=13.790947517451899863330960038627;
            printf("%5d %.10f %.10f %.10f\n", tn.i, integralUP(u), integralUK(u), k*integralUP(u)+integralUK(u));
        }
        return;
    }

    //if (ln == 0) const_cast<CcIHyperbolicIBVP1*>(this)->mu1 = u;
    //if (ln == 1) const_cast<CcIHyperbolicIBVP1*>(this)->mu1 = u;
    //if (ln >= 2) const_cast<CcIHyperbolicIBVP1*>(this)->mu1 = u;

    //if (ln == 0 || ln == 1 || ln == 2 || ln%100==0)
    //{
    //    //IPrinter::printMatrix(u);
    //    //IPrinter::printSeperatorLine();

    //    char filename1[40];
    //    int size1 = sprintf(filename1, "e:/data/img1/image%d.png", ln);
    //    filename1[size1] = 0;

    //    double min = u10.min();
    //    double max = u10.max();

    //    QPixmap pxm;
    //    visualizeMatrixHeat(u10, min, max, pxm, 0, 0);
    //    pxm.save(QString(filename1), "PNG");
    //}
}

double WaveEquationIBVP::U(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
#if defined( x1_t1 )
    return sn.x + tn.t;
#elif defined( x1_t2 )
    return sn.x + tn.t*tn.t;
#elif defined( x1_t3 )
    return sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x2_t1 )
    return sn.x*sn.x + tn.t;
#elif defined( x2_t2 )
    return sn.x*sn.x + tn.t*tn.t;
#elif defined( x2_t3 )
    return sn.x*sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x3_t1 )
    return sn.x*sn.x*sn.x + tn.t;
#elif defined( x3_t2 )
    return sn.x*sn.x*sn.x + tn.t*tn.t;
#elif defined( x3_t3 )
    return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t;
#endif

    //--------- 2D ----------//
#ifdef EXAMPLE_2D_SIMPLE_1
    return sn.x + sn.y + tn.t;
#endif
#ifdef EXAMPLE_2D_1
    return 0.0 - (a*a)*(2.0+2.0);
#endif
#ifdef EXAMPLE_2D_2
    return 2.0 - 4.0*a*a + 2.0*tn.t*alpha;
#endif
    //return 2.0 - 6.0*a*a*(sn.x+sn.y) + 2.0*alpha*tn.t;
    //return 6.0*tn.t - 6.0*a*a*(sn.x+sn.y) + 3.0*alpha*tn.t*tn.t;
    //return 12.0*tn.t*tn.t - 6.0*a*a*(sn.x+sn.y) + 4.0*alpha*tn.t*tn.t*tn.t;
    //return 0.0;
    //return 0.0;
    //return 0.0;
}

double WaveEquationIBVP::Udt1(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
#if defined( x1_t1 ) ||  defined( x2_t1 ) || defined( x3_t1 )
    return 1.0;
#elif defined( x1_t2 ) || defined( x2_t2 ) || defined( x3_t2 )
    return 2.0*tn.t;
#elif defined( x1_t3 ) || defined( x2_t3 ) || defined( x3_t3 )
    return 3.0*tn.t*tn.t;
#endif
}

double WaveEquationIBVP::Udt2(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
#if defined( x1_t1 ) ||  defined( x2_t1 ) || defined( x3_t1 )
    return 0.0;
#elif defined( x1_t2 ) || defined( x2_t2 ) || defined( x3_t2 )
    return 2.0;
#elif defined( x1_t3 ) || defined( x2_t3 ) || defined( x3_t3 )
    return 6.0*tn.t;
#endif
}

double WaveEquationIBVP::Udx1(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
#if defined( x1_t1 ) || defined( x1_t2 ) || defined( x1_t3 )
    return 1.0;
#elif defined( x2_t1 ) || defined( x2_t2 ) || defined( x2_t3 )
    return 2.0*sn.x;
#elif defined( x3_t1 ) || defined( x3_t2 ) ||  defined( x3_t3 )
    return 3.0*sn.x*sn.x;
#endif
}

double WaveEquationIBVP::Udx2(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
#if defined( x1_t1 ) || defined( x1_t2 ) || defined( x1_t3 )
    return 0.0;
#elif defined( x2_t1 ) || defined( x2_t2 ) || defined( x2_t3 )
    return 2.0;
#elif defined( x3_t1 ) || defined( x3_t2 ) ||  defined( x3_t3 )
    return 6.0*sn.x;
#endif
}

void WaveEquationIBVP::saveToImage(const DoubleMatrix &u UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
#ifdef USE_LIB_IMAGING
    //const unsigned int L = static_cast<unsigned int>(_timeDimension.size());
    const unsigned int N = static_cast<unsigned int>(_spaceDimensionX.size());
    const unsigned int M = static_cast<unsigned int>(_spaceDimensionY.size());
    //const double ht = _timeDimension.step();
    //const double hx = _spaceDimensionX.step();
    //const double hy = _spaceDimensionY.step();

    //QDir path("D:/data2");
    //if (!path.exists()) path.mkdir("D:/data2");

    QString filename = QString("data/border/f/png/f_%1.png").arg(tn.i, 4, 10, QChar('0'));
    QPixmap pixmap;
    visualizeMatrixHeat(u, u.min(), u.max(), pixmap, N+1, M+1);
    pixmap.save(filename);
#endif
}

void WaveEquationIBVP::saveToTextF(const DoubleMatrix &u UNUSED_PARAM, const TimeNodePDE &tn) const
{
    std::string txt_number = std::to_string(tn.i);
    std::string filename = std::string("data/border/f/txt/f_") +
            std::string(4 - txt_number.length(), '0') + txt_number + std::string(".txt");
    IPrinter::print(u, filename.c_str());
    printf("Forward: %4d %0.3f %10.8f %10.8f %4d %4d\n", tn.i, tn.t, u.min(), u.max(), 0, 0);
}

double WaveEquationIBVP::integralUP(const DoubleVector &) const
{
    const static double hx = _spaceDimensionX.step();
    const static unsigned int N = static_cast<unsigned int> ( _spaceDimensionX.size() );

    double sum = 0.0;

    sum += 0.50 * vu1[0]*vu1[0];
    for (unsigned int n=1; n<=N-1; n++)
    {
        sum += vu1[n]*vu1[n];
    }
    sum += 0.50 * vu1[N]*vu1[N];

    return sum*hx;
}

double WaveEquationIBVP::integralUK(const DoubleVector &) const
{
    const double hx = _spaceDimensionX.step();
    const unsigned int N = static_cast<unsigned int> ( _spaceDimensionX.size() );

    double sum = 0.0;

    sum += 0.50 * vut[0]*vut[0];
    for (unsigned int n=1; n<=N-1; n++)
    {
        sum += vut[n]*vut[n];
    }
    sum += 0.50 * vut[N]*vut[N];

    return sum*hx;
}

double WaveEquationIBVP::integralUP(const DoubleMatrix &) const
{
    const double hx = _spaceDimensionX.step();
    const double hy = _spaceDimensionY.step();
    const unsigned int N = static_cast<unsigned int> ( _spaceDimensionX.size() );
    const unsigned int M = static_cast<unsigned int> ( _spaceDimensionY.size() );

    double sum = 0.0;

    sum += 0.25 * fabs(mu1[0][0]);// * mu1[0][0];
    sum += 0.25 * fabs(mu1[0][N]);// * mu1[0][N];
    sum += 0.25 * fabs(mu1[M][0]);// * mu1[M][0];
    sum += 0.25 * fabs(mu1[M][N]);// * mu1[M][N];

    for (unsigned int n=1; n<=N-1; n++)
    {
        sum += 0.5 * fabs(mu1[0][n]);// * mu1[0][n];
        sum += 0.5 * fabs(mu1[M][n]);// * mu1[M][n];
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        sum += 0.5 * fabs(mu1[m][0]);// * mu1[m][0];
        sum += 0.5 * fabs(mu1[m][N]);// * mu1[m][N];
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        for (unsigned int n=1; n<=N-1; n++)
        {
            sum += fabs(mu1[m][n]);// * mu1[m][n];
        }
    }

    return sum*(hx*hy);
}

double WaveEquationIBVP::integralUK(const DoubleMatrix &) const
{
    const double hx = _spaceDimensionX.step();
    const double hy = _spaceDimensionY.step();
    const unsigned int N = static_cast<unsigned int> ( _spaceDimensionX.size() );
    const unsigned int M = static_cast<unsigned int> ( _spaceDimensionY.size() );

    double sum = 0.0;

    sum += 0.25 * fabs(mut[0][0]);// * mut[0][0];
    sum += 0.25 * fabs(mut[0][N]);// * mut[0][N];
    sum += 0.25 * fabs(mut[M][0]);// * mut[M][0];
    sum += 0.25 * fabs(mut[M][N]);// * mut[M][N];

    for (unsigned int n=1; n<=N-1; n++)
    {
        sum += 0.5 * fabs(mut[0][n]);// * mut[0][n];
        sum += 0.5 * fabs(mut[M][n]);// * mut[M][n];
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        sum += 0.5 * fabs(mut[m][0]);// * mut[m][0];
        sum += 0.5 * fabs(mut[m][N]);// * mut[m][N];
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        for (unsigned int n=1; n<=N-1; n++)
        {
            sum += fabs(mut[m][n]);// * mut[m][n];
        }
    }

    return sum*(hx*hy);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void WaveEquationFBVP::Main(int argc UNUSED_PARAM, char **argv UNUSED_PARAM)
{
#ifdef USE_LIB_IMAGING
    QGuiApplication app(argc, argv);
#endif

    WaveEquationFBVP w;
    w.setTimeDimension(Dimension(0.01, 0, 100));
    w.setSpaceDimensionX(Dimension(0.01, 100, 200));
    //w.setSpaceDimensionY(Dimension(0.01, 0, 100));

    w.setWaveSpeed(a);//a
    w.setWaveDissipation(d);//d
    w.setUnknownB(b);//b
    w.setRestoration(c);//c


    Benchmark bm;
    bm.tick();
    w.implicit_calculate_D1V1();
    bm.tock();
    bm.printDuration();
}

double WaveEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition condition) const
{
    TimeNodePDE tn; tn.t = timeDimension().max()*timeDimension().step();
    if (condition == FinalCondition::FinalValue) return U(sn, tn);
    if (condition == FinalCondition::FinalFirstDerivative) return Udt1(sn, tn);
    throw std::exception();
}

double WaveEquationFBVP::boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM, BoundaryConditionPDE &condition UNUSED_PARAM) const
{
    if (sn.i == spaceDimensionX().min())
    {
#if defined(__DIRICHLET_LEFT_)
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
        return U(sn, tn)*(condition.alpha()/condition.gamma());
#endif
#if defined(__ROBIN_LEFT__)
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, +4.0, +2.0, +1.0);
        return (condition.alpha()*U(sn, tn)+condition.beta()*Udx1(sn, tn))/condition.gamma();
#endif
    }

    if (sn.i == spaceDimensionX().max())
    {
#if defined(__DIRICHLET_RIGHT_)
        condition = BoundaryConditionPDE(BoundaryCondition::Dirichlet, +1.0, +0.0, +1.0);
        return U(sn, tn)*(condition.alpha()/condition.gamma());
#endif
#if defined(__ROBIN_RIGHT__)
        condition = BoundaryConditionPDE(BoundaryCondition::Robin, +4.0, +2.0, +1.0);
        return (condition.alpha()*U(sn, tn)+condition.beta()*Udx1(sn, tn))/condition.gamma();
#endif
    }

    throw std::exception();
}

double WaveEquationFBVP::f(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    const double a = waveSpeed();//ws
    const double b = unknownB();//wd
    const double c = restoration();//rs
    const double d = waveDissipation();//wd

#if defined (x1_t1) || defined (x2_t1) || defined (x3_t1) || defined (x1_t2) || defined (x2_t2) || defined (x3_t2) || defined (x1_t3) || defined (x2_t3) || defined (x3_t3)
    return Udt2(sn,tn)+d*Udt1(sn,tn)-a*a*Udx2(sn,tn)-b*Udx1(sn,tn)-c*U(sn,tn);
#endif
}

void WaveEquationFBVP::layerInfo(const DoubleVector& u, const TimeNodePDE& tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);

    IPrinter::printVector(u);
}

void WaveEquationFBVP::layerInfo(const DoubleMatrix& p, const TimeNodePDE& tn) const
{
    if (tn.i==200 || tn.i==199 || tn.i==198 || /*tn.i==397 || tn.i==396 ||*/ tn.i==2 || tn.i==1 || tn.i==0)
    {
        IPrinter::printMatrix(p);
        IPrinter::printSeperatorLine();
    }
}

double WaveEquationFBVP::U(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
#if defined( x1_t1 )
    return sn.x + tn.t;
#elif defined( x1_t2 )
    return sn.x + tn.t*tn.t;
#elif defined( x1_t3 )
    return sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x2_t1 )
    return sn.x*sn.x + tn.t;
#elif defined( x2_t2 )
    return sn.x*sn.x + tn.t*tn.t;
#elif defined( x2_t3 )
    return sn.x*sn.x + tn.t*tn.t*tn.t;
#endif

#if defined( x3_t1 )
    return sn.x*sn.x*sn.x + tn.t;
#elif defined( x3_t2 )
    return sn.x*sn.x*sn.x + tn.t*tn.t;
#elif defined( x3_t3 )
    return sn.x*sn.x*sn.x + tn.t*tn.t*tn.t;
#endif

    //--------- 2D ----------//
#ifdef EXAMPLE_2D_SIMPLE_1
    return sn.x + sn.y + tn.t;
#endif
#ifdef EXAMPLE_2D_1
    return 0.0 - (a*a)*(2.0+2.0);
#endif
#ifdef EXAMPLE_2D_2
    return 2.0 - 4.0*a*a + 2.0*tn.t*alpha;
#endif
}

double WaveEquationFBVP::Udt1(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
#if defined( x1_t1 ) ||  defined( x2_t1 ) || defined( x3_t1 )
    return 1.0;
#elif defined( x1_t2 ) || defined( x2_t2 ) || defined( x3_t2 )
    return 2.0*tn.t;
#elif defined( x1_t3 ) || defined( x2_t3 ) || defined( x3_t3 )
    return 3.0*tn.t*tn.t;
#endif
}

double WaveEquationFBVP::Udt2(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
#if defined( x1_t1 ) ||  defined( x2_t1 ) || defined( x3_t1 )
    return 0.0;
#elif defined( x1_t2 ) || defined( x2_t2 ) || defined( x3_t2 )
    return 2.0;
#elif defined( x1_t3 ) || defined( x2_t3 ) || defined( x3_t3 )
    return 6.0*tn.t;
#endif
}

double WaveEquationFBVP::Udx1(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
#if defined( x1_t1 ) || defined( x1_t2 ) || defined( x1_t3 )
    return 1.0;
#elif defined( x2_t1 ) || defined( x2_t2 ) || defined( x2_t3 )
    return 2.0*sn.x;
#elif defined( x3_t1 ) || defined( x3_t2 ) ||  defined( x3_t3 )
    return 3.0*sn.x*sn.x;
#endif
}

double WaveEquationFBVP::Udx2(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
#if defined( x1_t1 ) || defined( x1_t2 ) || defined( x1_t3 )
    return 0.0;
#elif defined( x2_t1 ) || defined( x2_t2 ) || defined( x2_t3 )
    return 2.0;
#elif defined( x3_t1 ) || defined( x3_t2 ) ||  defined( x3_t3 )
    return 6.0*sn.x;
#endif
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
