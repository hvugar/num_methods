#include "hyperbolicibvp1.h"

double distribute_blows(const Dimension &dimX, const Dimension &dimY, const SpaceNodePDE &sn)
{
    std::vector<SpacePoint> theta(2);
    DoubleVector q(2);
    q[0] = +0.145; theta[0].x = 0.2500; theta[0].y = 0.7200;
    q[1] = +0.147; theta[1].x = 0.6400; theta[1].y = 0.2700;

    double hx = dimX.step();
    double hy = dimY.step();

    const double sigmaX = 3.0*hx;
    const double sigmaY = 3.0*hy;
    double sum = 0.0;
    for (unsigned int s=0; s<2; s++)
    {
        double factor = 1.0/(2.0*M_PI*sigmaX*sigmaY);
        sum += q[s]*factor*exp(-0.5*(((sn.x-theta[s].x)*(sn.x-theta[s].x))/(sigmaX*sigmaX)+((sn.y-theta[s].y)*(sn.y-theta[s].y))/(sigmaY*sigmaY)));
    }
    return sum;
}

double distribute_control(const Dimension &dimX, const Dimension &dimY, const SpaceNodePDE &sn, const DoubleMatrix &u10)
{
    const static unsigned int N = static_cast<unsigned int> ( dimX.size() );
    const static unsigned int M = static_cast<unsigned int> ( dimY.size() );

    DoubleMatrix k(2, 2);
    DoubleMatrix z(2, 2);
    std::vector<SpacePoint> xi(2);
    std::vector<SpacePoint> eta(2);

    k[0][0] = -0.1820; k[0][1] = -0.1250; k[1][0] = -0.1550; k[1][1] = -0.1310;
    z[0][0] = -0.0262; z[0][1] = -0.0773; z[1][0] = -0.0570; z[1][1] = +0.0653;
    xi[0].x  = +0.3849; xi[0].y  = +0.5442; xi[1].x  = +0.7661; xi[1].y  = +0.6785;
    eta[0].x = +0.6656; eta[0].y = +0.7909; eta[1].x = +0.4856; eta[1].y = +0.3810;

    double hx = dimX.step();
    double hy = dimY.step();

    const double sigmaX = hx;
    const double sigmaY = hy;

    static double u[2] = {0.0,0.0};
    if (sn.i == 1 && sn.j == 1)
    {
        for (unsigned int j=0; j<2; j++)
        {
            SpaceNodePDE sn0;
            for (unsigned int m=0; m<=M; m++)
            {
                sn0.j = static_cast<int>(m); sn0.y = sn0.j*hy;
                for (unsigned int n=0; n<=N; n++)
                {
                    sn0.i = static_cast<int>(n); sn0.x = sn0.i*hx;
                    double factor = 1.0/(2.0*M_PI*sigmaX*sigmaY);
                    double w = factor*exp(-0.5*(((sn0.x-xi[j].x)*(sn0.x-xi[j].x))/(sigmaX*sigmaX)+
                                                ((sn0.y-xi[j].y)*(sn0.y-xi[j].y))/(sigmaY*sigmaY)));
                    u[j] += u10[m][n] * hx * hy * w;
                }
            }
        }
    }

    double v[2] = {0.0, 0.0};
    double sum = 0.0;
    for (unsigned int i=0; i<2; i++)
    {
        for (unsigned int j=0; j<2; j++)
        {
            v[i] += k[i][j]*(u[j]-z[i][j]);
        }
        double factor = 1.0/(2.0*M_PI*sigmaX*sigmaY);
        double w = factor*exp(-0.5*(((sn.x-eta[i].x)*(sn.x-eta[i].x))/(sigmaX*sigmaX)+
                                    ((sn.y-eta[i].y)*(sn.y-eta[i].y))/(sigmaY*sigmaY)));
        sum += v[i] * w;
    }
    return sum;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CcIHyperbolicIBVP1::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    QGuiApplication app(argc, argv);

    CcIHyperbolicIBVP1 hibvp;
    hibvp.a = 1.0;
    hibvp.setTimeDimension(Dimension(0.01, 0, 1000));
    hibvp.addSpaceDimension(Dimension(0.01, 0, 100));
    hibvp.addSpaceDimension(Dimension(0.01, 0, 100));

    //    DoubleVector u;
    //    hibvp.method = 1;
    //    hibvp.explicit_calculate_D1V1(u, hibvp.a);
    //    IPrinter::printVector(u);
    //    hibvp.method = 2;
    //    hibvp.implicit_calculate_D1V1(u, hibvp.a);
    //    IPrinter::printSeperatorLine();
    //    IPrinter::printVector(u);

    DoubleMatrix u;
    //hibvp.explicit_calculate_D2V1(u, hibvp.a);
    //IPrinter::printMatrix(u);
    //IPrinter::printSeperatorLine();

    hibvp.implicit_calculate_D2V1(u, hibvp.a);
    //IPrinter::printMatrix(u);
    //IPrinter::printSeperatorLine();
}

void CcIHyperbolicIBVP1::layerInfo(const DoubleMatrix &u, unsigned int ln) const
{
    //    if (ln==1 || ln==2000)
    //    {
    //        IPrinter::printMatrix(u);
    //        IPrinter::printSeperatorLine();
    //    }

    //    if (ln%1==0) { double w = u.min(); if (w < 0.0) {} else { w = u.max(); } printf("%14.10f\n", w); }
    //    return;

    //    if (ln%2==0)
    {
        const static double ht = timeDimension().step();
        const Dimension &dimX = spaceDimension(Dimension::DimensionX);
        const Dimension &dimY = spaceDimension(Dimension::DimensionY);

        const static unsigned int N = static_cast<unsigned int> ( dimX.size() );
        const static unsigned int M = static_cast<unsigned int> ( dimY.size() );
        const double hx = dimX.step();
        const double hy = dimY.step();

        DoubleMatrix &_u0 = const_cast<CcIHyperbolicIBVP1*>(this)->mu0;
        DoubleMatrix &_u1 = const_cast<CcIHyperbolicIBVP1*>(this)->mu1;
        DoubleMatrix &_u2 = const_cast<CcIHyperbolicIBVP1*>(this)->mu2;

        DoubleMatrix &_ut = const_cast<CcIHyperbolicIBVP1*>(this)->mut;

        if (ln == 0) { _u0 = u; }
        if (ln == 1) { _u1 = u; }
        if (ln == 2) { _u2 = u; }
        if (ln >= 3) { _u0 = _u1; _u1 = _u2; _u2 = u; }

        if ( ln == 0 )
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

        if ( ln == 1 )
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

        if ( ln >= 2)
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
            printf("%5d %.10f %.10f %.10f\n", ln, integralUP(u), integralUK(u), k*integralUP(u)+integralUK(u));
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

void CcIHyperbolicIBVP1::layerInfo(const DoubleVector &u, unsigned int ln) const
{
    QPixmap pxm;
    if (method == 1)
    {
        //visualString(u, -0.32, +0.32, 1000, 200, pxm, Qt::transparent, Qt::red, QString("d:/img/1/%1.png").arg(ln));
        visualString(u, -0.60, +0.60, 1000, 200, pxm, Qt::transparent, Qt::red, QString("d:/img/1/%1.png").arg(ln));
    }
    if (method == 2)
    {
        //visualString(u, -0.32, +0.32, 1000, 200, pxm, Qt::transparent, Qt::green, QString("d:/img/2/%1.png").arg(ln));
        visualString(u, -0.60, +0.60, 1000, 200, pxm, Qt::transparent, Qt::green, QString("d:/img/2/%1.png").arg(ln));
    }
    if (method == 3)
    {
        //visualString(u, -0.32, +0.32, 1000, 200, pxm, Qt::transparent, Qt::blue, QString("d:/img/3/%1.png").arg(ln));
        visualString(u, -0.60, +0.60, 1000, 200, pxm, Qt::transparent, Qt::blue, QString("d:/img/3/%1.png").arg(ln));

        QPixmap pxm(1000, 600);
        pxm.fill(Qt::yellow);
        QPainter painter(&pxm);
        QPixmap img1(QString("d:/img/1/%1.png").arg(ln), "PNG");
        QPixmap img2(QString("d:/img/2/%1.png").arg(ln), "PNG");
        QPixmap img3(QString("d:/img/3/%1.png").arg(ln), "PNG");
        painter.drawPixmap(5, 200, img1);
        painter.drawPixmap(5, 200, img2);
        painter.drawPixmap(2, 200, img3);
        pxm.save(QString("d:/img/0/%1.png").arg(ln));
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

    DoubleVector &_u0 = const_cast<CcIHyperbolicIBVP1*>(this)->vu0;
    DoubleVector &_u1 = const_cast<CcIHyperbolicIBVP1*>(this)->vu1;
    DoubleVector &_u2 = const_cast<CcIHyperbolicIBVP1*>(this)->vu2;
    DoubleVector &_ut = const_cast<CcIHyperbolicIBVP1*>(this)->vut;

    _ut.resize(N+1, 0.0);

    if (ln == 0) { const_cast<CcIHyperbolicIBVP1*>(this)->vu0 = u; }
    if (ln == 1) { const_cast<CcIHyperbolicIBVP1*>(this)->vu1 = u; }
    if (ln == 2) { const_cast<CcIHyperbolicIBVP1*>(this)->vu2 = u; }

    if (ln >= 3)
    {
        _u0 = _u1; _u1 = _u2; _u2 = u;
    }

    if (ln >= 2)
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
    if (ln > 1) printf("%5d %.10f %.10f %.10f\n", ln, k*integralUP(u), integralUK(u), k*integralUP(u)+integralUK(u));
}

double CcIHyperbolicIBVP1::initial(const SpaceNodePDE &sn UNUSED_PARAM, InitialCondition condition) const
{
    const Dimension &dimX = spaceDimension(Dimension::DimensionX); const double hx = dimX.step(); const double sigmaX = 5.0 * hx;
    const Dimension &dimY = spaceDimension(Dimension::DimensionY); const double hy = dimY.step(); const double sigmaY = 5.0 * hy;

    if (condition == InitialCondition::InitialValue)
    {
        // 1D
        //return sn.x*sn.x;
        //return sn.x*sn.x*sn.x;
        //return sn.x*sn.x*sn.x;
        //return 0.0;
        //return 0.0;
        //return (1.0/sqrt(2.0*M_PI*sigmaX*sigmaX)) * exp(-(((sn.x-0.5)*(sn.x-0.5))/(2.0*sigmaX*sigmaX)));

        // 2D
        //return sn.x*sn.x + sn.y*sn.y;
        //return sn.x*sn.x + sn.y*sn.y;
        //return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y;
        //return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y;
        //return 0.0;
        return 0.0;
        //return (1.0/(2.0*M_PI*sigmaX*sigmaX)) * exp(-(((sn.x-0.5)*(sn.x-0.5))/(2.0*sigmaX*sigmaX)+((sn.y-0.5)*(sn.y-0.5))/(2.0*sigmaY*sigmaY)));

        //return 0.0;
    }
    else
    {
        // 1D
        //return 0.0;
        //return 0.0;
        //return 0.0;
        //return sin(M_PI*sn.x);
        //return (1.0/sqrt(2.0*M_PI*sigmaX*sigmaX)) * exp(-(((sn.x-0.5)*(sn.x-0.5))/(2.0*sigmaX*sigmaX)));
        //return 0.0;

        // 2D
        //return 1.0;
        //return 0.0;
        //return 0.0;
        //return 0.0;
        //return sin(M_PI*sn.x)*sin(M_PI*sn.y);
        return (1.0/(2.0*M_PI*sigmaX*sigmaY)) * exp(-(((sn.x-0.5)*(sn.x-0.5))/(2.0*sigmaX*sigmaX) + ((sn.y-0.5)*(sn.y-0.5))/(2.0*sigmaY*sigmaY)));
        //return 0.0;

        //return distribute_blows(dimX, dimY, sn);
    }
}

double CcIHyperbolicIBVP1::boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    //const static Dimension &dimX = spaceDimension(Dimension::DimensionX);
    //const static Dimension &dimY = spaceDimension(Dimension::DimensionY);

    // 1D
    //return sn.x*sn.x+tn.t*tn.t;
    //return sn.x*sn.x*sn.x+tn.t*tn.t;
    //return sn.x*sn.x*sn.x+tn.t*tn.t*tn.t;
    //return 0.0;
    //return 0.0;
    //return 0.0;

    // 2D
    //return sn.x*sn.x + sn.y*sn.y + tn.t;
    //return sn.x*sn.x + sn.y*sn.y + tn.t*tn.t;
    //return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y + tn.t*tn.t;
    //return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y + tn.t*tn.t*tn.t;
    //return 0.0;
    return 0.0;
    //return 0.0;

    //return 0.0;
}

double CcIHyperbolicIBVP1::f(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    //const static Dimension &dimX = spaceDimension(Dimension::DimensionX);
    //const static Dimension &dimY = spaceDimension(Dimension::DimensionY);

    // 1D
    //return 2.0 - 2.0*a*a;
    //return 2.0 - 6.0*a*a*sn.x;
    //return 6.0*tn.t - 6.0*a*a*sn.x;
    //return 0.0;
    //return 0.0;
    //return 0.0;

    // 2D
    //return 0.0 - (a*a)*(2.0+2.0);
    //return 2.0 - (a*a)*(2.0+2.0);
    //return 2.0 - 6.0*a*a*(sn.x+sn.y);
    //return 6.0*tn.t - 6.0*a*a*(sn.x+sn.y);
    //return 0.0;
    return 0.0;
    //return 0.0;

    //return distribute_control(dimX, dimY, sn, this->u10);
}

double CcIHyperbolicIBVP1::integralUP(const DoubleVector &) const
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

double CcIHyperbolicIBVP1::integralUK(const DoubleVector &) const
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

double CcIHyperbolicIBVP1::integralUP(const DoubleMatrix &) const
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

double CcIHyperbolicIBVP1::integralUK(const DoubleMatrix &) const
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

void CdIHyperbolicIBVP1::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    //QGuiApplication app(argc, argv);

    CdIHyperbolicIBVP1 hibvp;
    hibvp.a = 1.0;
    hibvp.alpha = 0.0;
    hibvp.setTimeDimension(Dimension(0.01, 0, 100));
    hibvp.addSpaceDimension(Dimension(0.01, 0, 100));
    hibvp.addSpaceDimension(Dimension(0.01, 0, 100));

    DoubleMatrix u;
    //hibvp.explicit_calculate_D2V1(u, hibvp.a, hibvp.alpha);
    hibvp.implicit_calculate_D2V1(u, hibvp.a, hibvp.alpha);
}

void CdIHyperbolicIBVP1::layerInfo(const DoubleMatrix &u, unsigned int ln) const
{
    //    if (ln==200)
    //    {
    //        IPrinter::printMatrix(u);
    //        IPrinter::printSeperatorLine();

    //        SpaceNodePDE sp; TimeNodePDE tn; tn.i = 200; tn.t = 1.0;
    //        double sum = 0.0;
    //        for (unsigned int m=0; m<=100; m++)
    //        {
    //            sp.y = m*0.01;
    //            for (unsigned int n=0; n<=100; n++)
    //            {
    //                sp.x = n*0.01;
    //                sum += (u[m][n]-boundary(sp, tn))*(u[m][n]-boundary(sp, tn));
    //            }
    //        }
    //        sum = sqrt(0.01*0.01*sum);
    //        printf("Norm: %.8f\n", sum);
    //    }

    CdIHyperbolicIBVP1 *forward = const_cast<CdIHyperbolicIBVP1*>(this);

    if (ln == 0 || ln == 1 || ln == 2 || ln == 3 || ln == 4)
    {
        //        IPrinter::printMatrix(u);
        //        IPrinter::printSeperatorLine();
    }

    if (ln == 194) { forward->u2  = -2.0*u; }
    //    if (ln == 196) { forward->u2  = +1.0*u; }
    if (ln == 196) { forward->u2 += +9.0*u; }
    //    if (ln == 198) { forward->u2 += -4.0*u; }
    if (ln == 198) { forward->u2 += -18.0*u; }
    if (ln == 200)
    {
        forward->u1 = u;
        //        forward->u2 += +3.0*u;
        //        forward->u2 *= +(1.0/(2.0*0.01));
        forward->u2 += +11.0*u;
        forward->u2 *= +(1.0/(6.0*0.01));

        IPrinter::printMatrix(u1);
        IPrinter::printSeperatorLine();
        IPrinter::printMatrix(u2);
        IPrinter::printSeperatorLine();
    }

    //if (ln%1==0) { double w = u.min(); if (w < 0.0) {} else { w = u.max(); } printf("%14.10f\n", w); }
    return;
}

double CdIHyperbolicIBVP1::initial(const SpaceNodePDE &sn, InitialCondition condition) const
{
    C_UNUSED(sn);

    //const Dimension &dimX = spaceDimension(Dimension::DimensionX); const double hx = dimX.step(); const double sigmaX = 5.0 * hx;
    //const Dimension &dimY = spaceDimension(Dimension::DimensionY); const double hy = dimY.step(); const double sigmaY = 5.0 * hy;

    if (condition == InitialCondition::InitialValue)
    {
        // 1D
        //return sn.x*sn.x;
        //return sn.x*sn.x*sn.x;
        //return sn.x*sn.x*sn.x;
        //return 0.0;
        //return 0.0;
        //return (1.0/sqrt(2.0*M_PI*sigmaX*sigmaX)) * exp(-(((sn.x-0.5)*(sn.x-0.5))/(2.0*sigmaX*sigmaX)));

        // 2D
        //return sn.x*sn.x + sn.y*sn.y;
        //return sn.x*sn.x + sn.y*sn.y;
        //return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y;
        return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y;
        //return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y;
        //return 0.0;
        //return 0.0;
        //return (1.0/(2.0*M_PI*sigmaX*sigmaX)) * exp(-(((sn.x-0.5)*(sn.x-0.5))/(2.0*sigmaX*sigmaX)+((sn.y-0.5)*(sn.y-0.5))/(2.0*sigmaY*sigmaY)));

        //return 0.0;
    }
    else
    {
        // 1D
        //return 0.0;
        //return 0.0;
        //return 0.0;
        //return sin(M_PI*sn.x);
        //return (1.0/sqrt(2.0*M_PI*sigmaX*sigmaX)) * exp(-(((sn.x-0.5)*(sn.x-0.5))/(2.0*sigmaX*sigmaX)));
        //return 0.0;

        // 2D
        //return 1.0;
        //return 0.0;
        //return 0.0;
        return 0.0;
        //return 0.0;
        //return sin(M_PI*sn.x)*sin(M_PI*sn.y);
        //return (1.0/(2.0*M_PI*sigmaX*sigmaX)) * exp(-(((sn.x-0.5)*(sn.x-0.5))/(2.0*sigmaX*sigmaX)+((sn.y-0.5)*(sn.y-0.5))/(2.0*sigmaY*sigmaY)));
        //return 0.0;

        //return distribute_blows(dimX, dimY, sn);
    }
}

double CdIHyperbolicIBVP1::boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
    //const static Dimension &dimX = spaceDimension(Dimension::DimensionX);
    //const static Dimension &dimY = spaceDimension(Dimension::DimensionY);

    // 1D
    //return sn.x*sn.x+tn.t*tn.t;
    //return sn.x*sn.x*sn.x+tn.t*tn.t;
    //return sn.x*sn.x*sn.x+tn.t*tn.t*tn.t;
    //return 0.0;
    //return 0.0;
    //return 0.0;

    // 2D
    //return sn.x*sn.x + sn.y*sn.y + tn.t;
    //return sn.x*sn.x + sn.y*sn.y + tn.t*tn.t;
    //return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y + tn.t*tn.t;
    return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y + tn.t*tn.t*tn.t;
    //return sn.x*sn.x*sn.x + sn.y*sn.y*sn.y + tn.t*tn.t*tn.t*tn.t;
    //return 0.0;
    //return 0.0;
    //return 0.0;

    //return 0.0;
}

double CdIHyperbolicIBVP1::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    //const static Dimension &dimX = spaceDimension(Dimension::DimensionX);
    //const static Dimension &dimY = spaceDimension(Dimension::DimensionY);

    // 1D
    //return 2.0 - 2.0*a*a;
    //return 2.0 - 6.0*a*a*sn.x;
    //return 6.0*tn.t - 6.0*a*a*sn.x;
    //return 0.0;
    //return 0.0;
    //return 0.0;

    // 2D
    //return 0.0 - (a*a)*(2.0+2.0);
    //return 2.0 - (a*a)*(2.0+2.0);
    //return 2.0 - 6.0*a*a*(sn.x+sn.y) + alpha*2.0*tn.t;
    return 6.0*tn.t - 6.0*a*a*(sn.x+sn.y) + alpha*3.0*tn.t*tn.t;
    //return 12.0*tn.t*tn.t - 6.0*a*a*(sn.x+sn.y) + alpha*4.0*tn.t*tn.t*tn.t;
    //return 0.0;
    //return 0.0;
    //return 0.0;

    //return distribute_control(dimX, dimY, sn, this->u10);
}

CdIHyperbolicIBVP1::~CdIHyperbolicIBVP1() {}

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

void ConjugateCdIHyperbolicIBVP1::layerInfo(const DoubleMatrix& p, unsigned int ln) const
{
    if (ln==200 or ln==199 or ln==198 or ln==2 or ln==1 or ln==0)
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

ConjugateCdIHyperbolicIBVP1::~ConjugateCdIHyperbolicIBVP1() {}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void WaveEquationIBVP::Main(int argc, char* argv[])
{
    WaveEquationIBVP w;
    w.setTimeDimension(Dimension(0.001, 0, 1000));
    w.addSpaceDimension(Dimension(0.01, 0, 100));
    w.addSpaceDimension(Dimension(0.01, 0, 100));
    DoubleMatrix u;
    Benchmark bm;
    bm.tick();
    w.implicit_calculate_D2V1(u, w.a, w.alpha);
    bm.tock();
    IPrinter::printMatrix(12, 6, u);
    bm.printDuration();
}

double WaveEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition condition) const
{
    const auto x = sn.x;
    const auto y = sn.y;

    if (condition == InitialCondition::InitialValue) return x*x*x + y*y*y;
    if (condition == InitialCondition::FirstDerivative) return 0.0;
    return NAN;
}

double WaveEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    const auto x = sn.x;
    const auto y = sn.y;
    const auto t = tn.t;

    return x*x*x + y*y*y + t*t;
}

double WaveEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    const auto x = sn.x;
    const auto y = sn.y;
    const auto t = tn.t;

    return 2.0 - 6.0*a*a*(x+y) + 2.0*alpha*t;
}

void WaveEquationIBVP::layerInfo(const DoubleVector &, unsigned int) const
{}

void WaveEquationIBVP::layerInfo(const DoubleMatrix &, unsigned int) const
{}
