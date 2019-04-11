#include "testwaveequation.h"

void TestWaveEquation::Main(int argc, char **argv)
{
    QGuiApplication app(argc, argv);

    TestWaveEquation teq;
    teq.p1 = SpacePoint(0.5, 0.4);
    teq._waveSpeed = 1.0;
    teq._dissipation = 0.0;
    teq.setTimeDimension(Dimension(0.0025, 0, 400));
    teq.addSpaceDimension(Dimension(0.005, 0, 200));
    teq.addSpaceDimension(Dimension(0.005, 0, 200));

    double sigmaX = 8.0*teq.spaceDimension(Dimension::DimensionX).step();
    double sigmaY = 8.0*teq.spaceDimension(Dimension::DimensionY).step();
    teq.alpha1 = 1.0/(2.0*M_PI*sigmaX*sigmaY);
    teq.alpha2 = 1.0/(2.0*sigmaX*sigmaY);
    teq.alpha3 = 50.0;

    DoubleMatrix u;
    teq.implicit_calculate_D2V1(u, teq._waveSpeed, teq._dissipation);
    //teq.layerInfo(u, 0);
}

double TestWaveEquation::U(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    return alpha1 * exp( -alpha2*((p1.x-sn.x)*(p1.x-sn.x)+(p1.y-sn.y)*(p1.y-sn.y)) - alpha3*tn.t );
}

double TestWaveEquation::initial(const SpaceNodePDE &sn, InitialCondition condition) const
{
    TimeNodePDE tn; tn.t = 0.0; tn.i = 0;
    if (condition == InitialCondition::InitialValue)    { return U(sn, tn); }
    if (condition == InitialCondition::FirstDerivative) { return -alpha3*U(sn, tn); }
    throw std::exception();
}

double TestWaveEquation::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    return U(sn, tn);
}

double TestWaveEquation::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double _u = U(sn, tn);
    return ( alpha3*alpha3 - _dissipation*alpha3 - _waveSpeed*_waveSpeed *
             (
                 ( 4.0*alpha2*alpha2*(sn.x-p1.x)*(sn.x-p1.x) - 2.0*alpha2) +
                 ( 4.0*alpha2*alpha2*(sn.y-p1.y)*(sn.y-p1.y) - 2.0*alpha2)
                 )
             ) * _u;
}

void TestWaveEquation::layerInfo(const DoubleMatrix &u, unsigned int ln) const
{
    saveToImage(u, ln);

    const Dimension &time = timeDimension();
    const unsigned int L = static_cast<unsigned int>(time.size());
    const double ht = time.step();

    const Dimension &dimX = spaceDimension(Dimension::DimensionX);
    const Dimension &dimY = spaceDimension(Dimension::DimensionY);
    const unsigned int N = static_cast<unsigned int>(dimX.size());
    const unsigned int M = static_cast<unsigned int>(dimY.size());
    const double hx = dimX.step();
    const double hy = dimY.step();

    TestWaveEquation *forward = const_cast<TestWaveEquation*>(this);
    if (ln == 2*(L-2)) { forward->u2  = +1.0*u; }
    if (ln == 2*(L-1)) { forward->u2 += -4.0*u; }
    if (ln == 2*(L-0)) { forward->u2 += +3.0*u; forward->u2 *= 1.0/(2.0*ht); }

    if (ln == 2*(L-0))
    {
        IPrinter::printMatrix(u2);

        DoubleMatrix ut(M+1, N+1);
        SpaceNodePDE sn;
        TimeNodePDE tn; tn.t = 1.0; tn.i = 2*L;
        for (unsigned int m=0; m<=M; m++)
        {
            sn.y = m*hy;
            for (unsigned int n=0; n<=N; n++)
            {
                sn.x = n*hx;
                ut[m][n] = -alpha3*U(sn, tn);
            }
        }

        IPrinter::printSeperatorLine();
        IPrinter::printMatrix(ut);
    }

    //    DoubleMatrix p(M+1, N+1);
    //    {
    //        SpaceNodePDE sn;
    //        TimeNodePDE tn; tn.i = ln; tn.t = ln*ht*0.5;
    //        for (unsigned int m=0; m<=M; m++)
    //        {
    //            sn.y = m*hy;
    //            for (unsigned int n=0; n<=N; n++)
    //            {
    //                sn.x = n*hx;
    //                p[m][n] = U(sn, tn);
    //            }
    //        }
    //    }

    //printf("%10.6f %10.6f\n", p.min(), p.max());

    //std::cout << p.min() << " " << p.max() << "      " << u.min() << " " << u.max() << std::endl;
}

auto TestWaveEquation::saveToImage(const DoubleMatrix &u, unsigned int ln) const -> void
{
    //    const Dimension &time = Problem0HForward::timeDimension();
    const Dimension &dimX = spaceDimension(Dimension::DimensionX);
    const Dimension &dimY = spaceDimension(Dimension::DimensionY);
    //    const unsigned int L = static_cast<unsigned int>(time.size());
    const unsigned int N = static_cast<unsigned int>(dimX.size());
    const unsigned int M = static_cast<unsigned int>(dimY.size());
    //    const double ht = time.step();
    //    const double hx = dimX.step();
    //    const double hy = dimY.step();

    //QString filename1 = QString("data/problem0H/b/txt/b_%1.txt").arg(ln, 4, 10, QChar('0'));
    //IPrinter::print(p,filename1.toLatin1().data());
    //IPrinter::printSeperatorLine();
    //printf("Forward: %4d %0.3f %10.8f %10.8f %4d %4d\n", ln, ln*0.005, u.min(), u.max(), 0, 0);

    QString filename2 = QString("data/problem0H/c/png/c_%1.png").arg(ln, 4, 10, QChar('0'));
    QPixmap pixmap;
    visualizeMatrixHeat(u, u.min(), u.max(), pixmap, N+1, M+1);
    pixmap.save(filename2);
    //IPrinter::printSeperatorLine();
}

/**************************************************************************************/

void TestWaveEquationEx1::Main(int argc, char **argv)
{
    QGuiApplication app(argc, argv);

    TestWaveEquationEx1 teq;
    teq._waveSpeed = 5.0;
    teq._waveDissipation = 0.0;
    teq.setTimeDimension(Dimension(0.001, 0, 2000));
    teq.addSpaceDimension(Dimension(0.001, 0, 1000));

    DoubleMatrix u;
    teq.implicit_calculate_D1V1();
}

double TestWaveEquationEx1::initial(const SpaceNodePDE &sn, InitialCondition condition) const
{
    if (condition == InitialCondition::InitialValue)
    {
        //return 0.1*sin(2.0*M_PI*sn.x);
        return 0.0;
    }

    if (condition == InitialCondition::FirstDerivative)
    {
         return 0.0;
    }

    throw std::exception();
}

double TestWaveEquationEx1::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    return 0.0;
}

double TestWaveEquationEx1::boundary1(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const
{
    if (sn.i == 0)
    {
        //std::cout << sn.i << std::endl;
        condition.type = BoundaryCondition::Neumann;
        condition.coefficientValue = 0.0;
        condition.coefficientDerivative = +1.0;
        return +5.00;
    }
    else
    {
        //std::cout << sn.i << std::endl;
        condition.type = BoundaryCondition::Neumann;
        condition.coefficientValue = 0.0;
        condition.coefficientDerivative = +1.0;
        return -5.0;
    }
}

double TestWaveEquationEx1::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    return 0.0;
}

void TestWaveEquationEx1::layerInfo(const DoubleVector &v, unsigned int ln) const
{
    saveToImage(v, ln);
}

auto TestWaveEquationEx1::saveToImage(const DoubleVector &v, unsigned int ln) const -> void
{
    //const Dimension &time = Problem0HForward::timeDimension();
    //const Dimension &dimX = spaceDimension(Dimension::DimensionX);
    //const unsigned int L = static_cast<unsigned int>(time.size());
    //const unsigned int N = static_cast<unsigned int>(dimX.size());
    //const double ht = time.step();
    //const double hx = dimX.step();
    //const double hy = dimY.step();

    //QString filename1 = QString("data/problem0H/b/txt/b_%1.txt").arg(ln, 4, 10, QChar('0'));
    //IPrinter::print(p,filename1.toLatin1().data());
    //IPrinter::printSeperatorLine();
    //printf("Forward: %4d %0.3f %10.8f %10.8f %4d %4d\n", ln, ln*0.005, u.min(), u.max(), 0, 0);

    QString filename2 = QString("data/problem0H/s/png/c_%1.png").arg(ln, 4, 10, QChar('0'));
    QPixmap pixmap;
    visualString1(v, -0.105, +0.105, 1001, 200, Qt::white, Qt::red, filename2);
    pixmap.save(filename2);
    //IPrinter::printSeperatorLine();
}

