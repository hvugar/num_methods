#include "waveequationibvp.h"

void WaveEquationIBVP::Main(int argc UNUSED_PARAM, char** argv UNUSED_PARAM)
{
    WaveEquationIBVP w;
    w.setWaveSpeed(1.0);
    w.setWaveDissipation(0.01);
    w.setTimeDimension(Dimension(0.01, 0, 100));
    w.addSpaceDimension(Dimension(0.01, 0, 100));
    w.addSpaceDimension(Dimension(0.01, 0, 100));

    Benchmark bm;
    bm.tick();
    //w.implicit_calculate_D2V1();
    w.implicit_calculate_D2V1X();
    bm.tock();
    bm.printDuration();
}

double WaveEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition condition) const
{
    const Dimension &dimX = spaceDimension(Dimension::DimensionX); const double hx = dimX.step(); const double sigmaX = 5.0 * hx;
    const Dimension &dimY = spaceDimension(Dimension::DimensionY); const double hy = dimY.step(); const double sigmaY = 5.0 * hy;

    if (condition == InitialCondition::InitialValue)
    {
        //--------- 1D ----------//
        //return sn.x*sn.x;
        //return sn.x*sn.x*sn.x;
        //return sn.x*sn.x*sn.x;
        //return 0.0;
        //return 0.0;
        //return (1.0/sqrt(2.0*M_PI*sigmaX*sigmaX)) * exp(-(((sn.x-0.5)*(sn.x-0.5))/(2.0*sigmaX*sigmaX)));

        //--------- 2D ----------//
        //return sn.x*sn.x + sn.y*sn.y;
        //return sn.x*sn.x + sn.y*sn.y;
        return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y;
        //return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y;
        //return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y;
        //return 0.0;
        //return 0.0;
        //return (1.0/(2.0*M_PI*sigmaX*sigmaX)) * exp(-(((sn.x-0.5)*(sn.x-0.5))/(2.0*sigmaX*sigmaX)+((sn.y-0.5)*(sn.y-0.5))/(2.0*sigmaY*sigmaY)));

    }

    if (condition == InitialCondition::FirstDerivative)
    {
        //--------- 1D ----------//
        //return 0.0;
        //return 0.0;
        //return 0.0;
        //return sin(M_PI*sn.x);
        //return (1.0/sqrt(2.0*M_PI*sigmaX*sigmaX)) * exp(-(((sn.x-0.5)*(sn.x-0.5))/(2.0*sigmaX*sigmaX)));
        //return 0.0;

        //--------- 2D ----------//
        //return 1.0;
        //return 0.0;
        return 0.0;
        //return 0.0;
        //return 0.0;
        //return sin(M_PI*sn.x)*sin(M_PI*sn.y);
        //return (1.0/(2.0*M_PI*sigmaX*sigmaX)) * exp(-(((sn.x-0.5)*(sn.x-0.5))/(2.0*sigmaX*sigmaX)+((sn.y-0.5)*(sn.y-0.5))/(2.0*sigmaY*sigmaY)));
        //return 0.0;
    }

    return 0.0;
}

double WaveEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    //--------- 1D ----------//
    //return sn.x*sn.x+tn.t*tn.t;
    //return sn.x*sn.x*sn.x+tn.t*tn.t;
    //return sn.x*sn.x*sn.x+tn.t*tn.t*tn.t;
    //return 0.0;
    //return 0.0;
    //return 0.0;

    //--------- 2D ----------//
    //return sn.x*sn.x + sn.y*sn.y + tn.t;
    //return sn.x*sn.x + sn.y*sn.y + tn.t*tn.t;
    return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y + tn.t*tn.t;
    //return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y + tn.t*tn.t*tn.t;
    //return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y + tn.t*tn.t*tn.t*tn.t;
    //return 0.0;
    //return 0.0;
    //return 0.0;
}

double WaveEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    const double a = waveSpeed();
    const double alpha = waveDissipation();

    //--------- 1D ----------//
    //return sn.x*sn.x+tn.t*tn.t;
    //return 2.0 - 6.0*a*a*sn.x;
    //return 6.0*tn.t - 6.0*a*a*sn.x;
    //return 0.0;
    //return 0.0;
    //return 0.0;

    //--------- 2D ----------//
    //return 0.0 - (a*a)*(2.0+2.0);
    //return 2.0 - (a*a)*(2.0+2.0);
    return 2.0 - 6.0*a*a*(sn.x+sn.y) + 2.0*alpha*tn.t;
    //return 6.0*tn.t - 6.0*a*a*(sn.x+sn.y) + 3.0*alpha*tn.t*tn.t;
    //return 12.0*tn.t*tn.t - 6.0*a*a*(sn.x+sn.y) + 4.0*alpha*tn.t*tn.t*tn.t;
    //return 0.0;
    //return 0.0;
    //return 0.0;
}

void WaveEquationIBVP::layerInfo(const DoubleVector& u, const TimeNodePDE& tn) const
{
    C_UNUSED(u);
    C_UNUSED(tn);
#ifdef USE_LIB_IMAGING
    QPixmap pxm;
    if (method == 1)
    {
        //visualString(u, -0.32, +0.32, 1000, 200, pxm, Qt::transparent, Qt::red, QString("d:/img/1/%1.png").arg(ln));
        visualString(u, -0.60, +0.60, 1000, 200, pxm, Qt::transparent, Qt::red, QString("d:/img/1/%1.png").arg(tn.i));
    }
    if (method == 2)
    {
        //visualString(u, -0.32, +0.32, 1000, 200, pxm, Qt::transparent, Qt::green, QString("d:/img/2/%1.png").arg(ln));
        visualString(u, -0.60, +0.60, 1000, 200, pxm, Qt::transparent, Qt::green, QString("d:/img/2/%1.png").arg(tn.i));
    }
    if (method == 3)
    {
        //visualString(u, -0.32, +0.32, 1000, 200, pxm, Qt::transparent, Qt::blue, QString("d:/img/3/%1.png").arg(ln));
        visualString(u, -0.60, +0.60, 1000, 200, pxm, Qt::transparent, Qt::blue, QString("d:/img/3/%1.png").arg(tn.i));

        QPixmap pxm(1000, 600);
        pxm.fill(Qt::yellow);
        QPainter painter(&pxm);
        QPixmap img1(QString("d:/img/1/%1.png").arg(tn.i), "PNG");
        QPixmap img2(QString("d:/img/2/%1.png").arg(tn.i), "PNG");
        QPixmap img3(QString("d:/img/3/%1.png").arg(tn.i), "PNG");
        painter.drawPixmap(5, 200, img1);
        painter.drawPixmap(5, 200, img2);
        painter.drawPixmap(2, 200, img3);
        pxm.save(QString("d:/img/0/%1.png").arg(tn.i));
    }

    //    QPixmap pxm(1000, 600);
    //    pxm.fill(Qt::yellow);
    //    QPainter painter(&pxm);
    //    QPixmap img1(QString("d:/img/1/%1.png").arg(ln), "PNG");
    //    QPixmap latex(QString("d:/img/latex_1_200.png"), "PNG");
    //    painter.drawPixmap(5, 200, img1);
    //    painter.drawPixmap(1000-540, 10, latex);
    //    pxm.save(QString("d:/img/4/%1.png").arg(ln));

    //printf("\rSaved: %d\n", ln);

    //IPrinter::printVector(u);
    //double w = u.min(); if (w < 0.0) {} else { w = u.max(); } printf("%14.10f\n", w);
    return;

    const static double ht = timeDimension().step();
    const static unsigned int N = static_cast<unsigned int> ( spaceDimension(Dimension::DimensionX).size() );

    DoubleVector &_u0 = const_cast<WaveEquationIBVP*>(this)->vu0;
    DoubleVector &_u1 = const_cast<WaveEquationIBVP*>(this)->vu1;
    DoubleVector &_u2 = const_cast<WaveEquationIBVP*>(this)->vu2;
    DoubleVector &_ut = const_cast<WaveEquationIBVP*>(this)->vut;

    _ut.resize(N+1, 0.0);

    if (tn.i == 0) { const_cast<WaveEquationIBVP*>(this)->vu0 = u; }
    if (tn.i == 1) { const_cast<WaveEquationIBVP*>(this)->vu1 = u; }
    if (tn.i == 2) { const_cast<WaveEquationIBVP*>(this)->vu2 = u; }

    if (tn.i >= 3)
    {
        _u0 = _u1; _u1 = _u2; _u2 = u;
    }

    if (tn.i >= 2)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            _ut[n] = (_u2[n] - _u0[n])/(2.0*ht);
        }
    }

    //if (ln%20000==0) IPrinter::printVector(u, nullptr, u.length());
    //double k = 9.8681840001001186485325310599663;
    //double k = 1.0/523.3314251378401924542377492593;
    double k = 23.788030395009524606979129597455;
    if (tn.i > 1) printf("%5d %.10f %.10f %.10f\n", tn.i, k*integralUP(u), integralUK(u), k*integralUP(u)+integralUK(u));
#endif
}

void WaveEquationIBVP::layerInfo(const DoubleMatrix& u, const TimeNodePDE& tn) const
{
    if (tn.i==1 || tn.i==2 || tn.i==3 || tn.i==4 || tn.i==200)
    {
        IPrinter::printMatrix(u);
        IPrinter::printSeperatorLine();
    }
    return;

    //if (ln%1==0) { double w = u.min(); if (w < 0.0) {} else { w = u.max(); } printf("%14.10f\n", w); }

    //    if (ln%2==0)
    {
        const static double ht = timeDimension().step();
        const Dimension &dimX = spaceDimension(Dimension::DimensionX);
        const Dimension &dimY = spaceDimension(Dimension::DimensionY);

        const static unsigned int N = static_cast<unsigned int> ( dimX.size() );
        const static unsigned int M = static_cast<unsigned int> ( dimY.size() );
        const double hx = dimX.step();
        const double hy = dimY.step();

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
                    _ut[m][n] = initial(sn, InitialCondition::FirstDerivative);
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

void WaveEquationIBVP::saveToImage(const DoubleMatrix &u UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
#ifdef USE_LIB_IMAGING
    const Dimension &time = timeDimension();
    const Dimension &dimX = spaceDimension(Dimension::DimensionX);
    const Dimension &dimY = spaceDimension(Dimension::DimensionY);
    //const unsigned int L = static_cast<unsigned int>(time.size());
    const unsigned int N = static_cast<unsigned int>(dimX.size());
    const unsigned int M = static_cast<unsigned int>(dimY.size());
    //const double ht = time.step();
    //const double hx = dimX.step();
    //const double hy = dimY.step();

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
    const static double hx = spaceDimension(Dimension::DimensionX).step();
    const static unsigned int N = static_cast<unsigned int> ( spaceDimension(Dimension::DimensionX).size() );

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
    const double hx = spaceDimension(Dimension::DimensionX).step();
    const unsigned int N = static_cast<unsigned int> ( spaceDimension(Dimension::DimensionX).size() );

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
    const double hx = spaceDimension(Dimension::DimensionX).step();
    const double hy = spaceDimension(Dimension::DimensionY).step();
    const unsigned int N = static_cast<unsigned int> ( spaceDimension(Dimension::DimensionX).size() );
    const unsigned int M = static_cast<unsigned int> ( spaceDimension(Dimension::DimensionY).size() );

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
    const double hx = spaceDimension(Dimension::DimensionX).step();
    const double hy = spaceDimension(Dimension::DimensionY).step();
    const unsigned int N = static_cast<unsigned int> ( spaceDimension(Dimension::DimensionX).size() );
    const unsigned int M = static_cast<unsigned int> ( spaceDimension(Dimension::DimensionY).size() );

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

void ConjugateCdIHyperbolicIBVP1::Main(int argc UNUSED_PARAM, char **argv UNUSED_PARAM)
{
    ConjugateCdIHyperbolicIBVP1 chibvp;
    chibvp.a = 1.0;
    chibvp.alpha = 0.1;
    chibvp.setTimeDimension(Dimension(0.005, 0, 200));
    chibvp.addSpaceDimension(Dimension(0.01, 0, 100));
    chibvp.addSpaceDimension(Dimension(0.01, 0, 100));

    DoubleMatrix p;
    //chibvp.explicit_calculate_D2V1(p, chibvp.a, chibvp.alpha);
    //chibvp.implicit_calculate_D2V1(p, chibvp.a, chibvp.alpha);
    //chibvp.implicit_calculate_D2V3(p, chibvp.a, chibvp.alpha);
    IPrinter::printSeperatorLine("******************");
    IPrinter::printMatrix(p);
    IPrinter::printSeperatorLine();
}

void ConjugateCdIHyperbolicIBVP1::layerInfo(const DoubleMatrix& p, const TimeNodePDE& tn) const
{
    if (tn.i==200 || tn.i==199 || tn.i==198 || tn.i==2 || tn.i==1 || tn.i==0)
    {
        IPrinter::printMatrix(p);
        IPrinter::printSeperatorLine();
    }
}

double ConjugateCdIHyperbolicIBVP1::initial(const SpaceNodePDE &sn UNUSED_PARAM, InitialCondition condition) const
{
    if (condition == InitialCondition::InitialValue)
    {
        return sn.x*sn.x + sn.y*sn.y + 1.0;
    }
    else
    {
        return 2.0;
    }
}

double ConjugateCdIHyperbolicIBVP1::boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    return sn.x*sn.x + sn.y*sn.y + tn.t*tn.t;
}

double ConjugateCdIHyperbolicIBVP1::f(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    return 2.0 - 4.0*a*a - 2.0*tn.t*alpha;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
