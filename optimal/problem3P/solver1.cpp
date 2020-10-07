#include "solver1.h"
#include "heatequationibvp1.h"

using namespace p3p1;

void Solver1::optimize(int argc, char **argv)
{
    QGuiApplication app(argc, argv);

    HeatEquationIBVP1 he;
    he.setThermalDiffusivity(1.0);
    he.setThermalConvection(-he.lambda0);
    he.implicit_calculate_D2V1();
    return;

    const size_t time_size = 100;
    const double time_step = 0.01;
    const size_t heatSourceNumber = 2;
    const size_t measrPointNumber = 4;



    Solver1 s(Dimension(0.01, 0, 100), Dimension(0.01, 0, 100), Dimension(0.01, 0, 100));
    s._factor = 1.0;
    s.setPointNumber(heatSourceNumber, measrPointNumber);
    s.forward.solver = &s;
    s.backward.solver = &s;

    double a = 0.01;

    s._initialTemperature = s._initialTemperatureList[0];
    s._environmentTemperature = s._environmentTemperatureList[0];
    s.forward.setThermalDiffusivity(a);
    s.forward.setThermalConvection(-s._lambda0);
    s.forward.setThermalConductivity(0.0);
    s.drawImage = true;
    s.forward.findUMode = true;
    s.forward.implicit_calculate_D2V1();
    s.drawImage = false;
    s.V = s.U;
    s.forward.findUMode = false;
    //return;

    IPrinter::printSeperatorLine("matrix");
    IPrinter::printMatrix(10, 4, s.U);
    IPrinter::printSeperatorLine();

    //s.V.resize(101, 101, 15.0);

    s.backward.setThermalDiffusivity(-a);
    s.backward.setThermalConvection(s._lambda0);
    s.backward.setThermalConductivity(0.0);

#ifdef ENABLE_PROGRAM_CONTROL

#else
    const size_t size = s.heatSourceNumber*s.measrPointNumber;
    const size_t indexes[12] = { 0, size, 2*size, 3*size, 4*size, 5*size, 6*size, 7*size, 8*size, 9*size, 10*size, 11*size };
    std::string labels[11] = { "alfa1:  ", "alfa2:  ", "alfa3:  ", "betaX1: ", "betaY1: ", "betaX2: ", "betaY2: ", "betaX3: ", "betaY3: ", "nomU1:  ",  "nomU2:  " };
    {
        puts("------------------------------------------------------------------------------");

        DoubleVector x;
        s.parameterToVector(x);
        for (size_t i=0; i<11; i++)
        {
            if (i==3 || i==9 || i==12) std::cout << "---" << std::endl;
            std::cout << labels[i]; IPrinter::print(x.mid(indexes[i], indexes[i+1]-1), x.mid(indexes[i], indexes[i+1]-1).length(), 8, 4);
        }

#ifdef ENABLE_CHANGING_VALUES
        puts("------------------------------- CHANCING VECTOR ------------------------------");

        for (size_t i=0; i<s.heatSourceNumber; i++)
        {
#ifdef ENABLE_PROGRAM_CONTROL
#else
            for (size_t j=0; j<s.measrPointNumber; j++)
            {
#ifdef ENABLE_ALPHA1_OPTIMIZATION
                s.alpha1[i][j] *= 1.0 + (cos((i+1)*0.1)+sin((j+1)*0.5)) * 0.10;
#endif
#ifdef ENABLE_ALPHA2_OPTIMIZATION
                s.alpha2[i][j] *= 1.0 + (cos((i+1)*0.1)+sin((j+1)*0.5)) * 0.10;
#endif
#ifdef ENABLE_ALPHA3_OPTIMIZATION
                s.alpha3[i][j] *= 1.0 + (cos((i+1)*0.1)+sin((j+1)*0.5)) * 0.10;
#endif
#ifdef ENABLE_BETTA1_OPTIMIZATION
                s.bettaX1[i][j] *= 1.0 + (cos((i+1)*0.1)+sin((j+1)*0.5)) * 0.10;
                s.bettaY1[i][j] *= 1.0 + (cos((i+1)*0.1)+sin((j+1)*0.5)) * 0.10;
#endif
#ifdef ENABLE_BETTA2_OPTIMIZATION
                s.bettaX2[i][j] *= 1.0 + (cos((i+1)*0.1)+sin((j+1)*0.5)) * 0.10;
                s.bettaY2[i][j] *= 1.0 + (cos((i+1)*0.1)+sin((j+1)*0.5)) * 0.10;
#endif
#ifdef ENABLE_BETTA3_OPTIMIZATION
                s.bettaX3[i][j] *= 1.0 + (cos((i+1)*0.1)+sin((j+1)*0.5)) * 0.10;
                s.bettaY3[i][j] *= 1.0 + (cos((i+1)*0.1)+sin((j+1)*0.5)) * 0.10;
#endif
#ifdef ENABLE_NOMIN1_OPTIMIZATION
                s.nomnU1[i][j] *= 1.0 + (cos((i+1)*0.1)+sin((j+1)*0.5)) * 0.10;
#endif
#ifdef ENABLE_NOMIN2_OPTIMIZATION
                s.nomnU2[i][j] *= 1.0 + (cos((i+1)*0.1)+sin((j+1)*0.5)) * 0.10;
#endif
            }
#endif
        }

        s.parameterToVector(x);
        for (size_t i=0; i<11; i++)
        {
            if (i==3 || i==9 || i==12) std::cout << "---" << std::endl;
            std::cout << labels[i]; IPrinter::print(x.mid(indexes[i], indexes[i+1]-1), x.mid(indexes[i], indexes[i+1]-1).length(), 8, 4);
        }

        puts("------------------------------- CHANCING VECTOR ------------------------------");
#endif
    }

    {
#ifdef ENABLE_CHECKING_GRADIENTS
        puts("\n-------------------------- ENABLE CHECKING GRADIENTS -------------------------");
        DoubleVector x0;
        s.parameterToVector(x0);
        DoubleVector g0;
        s.gradient(x0, g0);

        unsigned int p = 8;
        unsigned int d = 4;

#ifdef ENABLE_ALPHA1_OPTIMIZATION
        std::cout << labels[0]; IPrinter::print(g0.mid(indexes[0], indexes[1]-1).L2Normalize(), g0.mid(indexes[0], indexes[1]-1).length(), p, d);
#endif
#ifdef ENABLE_ALPHA2_OPTIMIZATION
        std::cout << labels[1]; IPrinter::print(g0.mid(indexes[1], indexes[2]-1).L2Normalize(), g0.mid(indexes[1], indexes[2]-1).length(), p, d);
#endif
#ifdef ENABLE_ALPHA3_OPTIMIZATION
        std::cout << labels[2]; IPrinter::print(g0.mid(indexes[2], indexes[3]-1).L2Normalize(), g0.mid(indexes[2], indexes[3]-1).length(), p, d);
#endif
        std::cout << "---" << std::endl;
#ifdef ENABLE_BETTA1_OPTIMIZATION
        std::cout << labels[3]; IPrinter::print(g0.mid(indexes[3], indexes[4]-1).L2Normalize(), g0.mid(indexes[3], indexes[4]-1).length(), p, d);
        std::cout << labels[4]; IPrinter::print(g0.mid(indexes[4], indexes[5]-1).L2Normalize(), g0.mid(indexes[4], indexes[5]-1).length(), p, d);
#endif
#ifdef ENABLE_BETTA2_OPTIMIZATION
        std::cout << labels[5]; IPrinter::print(g0.mid(indexes[5], indexes[6]-1).L2Normalize(), g0.mid(indexes[5], indexes[6]-1).length(), p, d);
        std::cout << labels[6]; IPrinter::print(g0.mid(indexes[6], indexes[7]-1).L2Normalize(), g0.mid(indexes[6], indexes[7]-1).length(), p, d);
#endif
#ifdef ENABLE_BETTA3_OPTIMIZATION
        std::cout << labels[7]; IPrinter::print(g0.mid(indexes[7], indexes[8]-1).L2Normalize(), g0.mid(indexes[7], indexes[8]-1).length(), p, d);
        std::cout << labels[8]; IPrinter::print(g0.mid(indexes[8], indexes[9]-1).L2Normalize(), g0.mid(indexes[8], indexes[9]-1).length(), p, d);
#endif
        std::cout << "---" << std::endl;
#ifdef ENABLE_NOMIN1_OPTIMIZATION
        std::cout << labels[9]; IPrinter::print(g0.mid(indexes[9], indexes[10]-1).L2Normalize(), g0.mid(indexes[9], indexes[10]-1).length(), p, d);
#endif
#ifdef ENABLE_NOMIN2_OPTIMIZATION
        std::cout << labels[10]; IPrinter::print(g0.mid(indexes[10], indexes[11]-1).L2Normalize(), g0.mid(indexes[10], indexes[11]-1).length(), p, d);
#endif

        x0.clear();

        puts("------------------------------------------------------------------------------");
        DoubleVector x1;
        s.parameterToVector(x1);
        DoubleVector g1(x1.length(), 0.0);
#ifdef ENABLE_ALPHA1_OPTIMIZATION
        IGradient::Gradient(&s, 0.01, x1, g1, indexes[0], indexes[1]-1); std::cout << labels[0]; IPrinter::print(g1.mid(indexes[0], indexes[1]-1).L2Normalize(), g1.mid(indexes[0], indexes[1]-1).length(), p, d);
#endif
#ifdef ENABLE_ALPHA2_OPTIMIZATION
        IGradient::Gradient(&s, 0.01, x1, g1, indexes[1], indexes[2]-1); std::cout << labels[1]; IPrinter::print(g1.mid(indexes[1], indexes[2]-1).L2Normalize(), g1.mid(indexes[1], indexes[2]-1).length(), p, d);
#endif
#ifdef ENABLE_ALPHA3_OPTIMIZATION
        IGradient::Gradient(&s, 0.01, x1, g1, indexes[2], indexes[3]-1); std::cout << labels[2]; IPrinter::print(g1.mid(indexes[2], indexes[3]-1).L2Normalize(), g1.mid(indexes[2], indexes[3]-1).length(), p, d);
#endif
        std::cout << "---" << std::endl;
#ifdef ENABLE_BETTA1_OPTIMIZATION
        IGradient::Gradient(&s, 0.01, x1, g1, indexes[3], indexes[4]-1); std::cout << labels[3]; IPrinter::print(g1.mid(indexes[3], indexes[4]-1).L2Normalize(), g1.mid(indexes[3], indexes[4]-1).length(), p, d);
        IGradient::Gradient(&s, 0.01, x1, g1, indexes[4], indexes[5]-1); std::cout << labels[4]; IPrinter::print(g1.mid(indexes[4], indexes[5]-1).L2Normalize(), g1.mid(indexes[4], indexes[5]-1).length(), p, d);
#endif
#ifdef ENABLE_BETTA2_OPTIMIZATION
        IGradient::Gradient(&s, 0.01, x1, g1, indexes[5], indexes[6]-1); std::cout << labels[5]; IPrinter::print(g1.mid(indexes[5], indexes[6]-1).L2Normalize(), g1.mid(5*size, 6*size-1).length(), p, d);
        IGradient::Gradient(&s, 0.01, x1, g1, indexes[6], indexes[7]-1); std::cout << labels[6]; IPrinter::print(g1.mid(indexes[6], indexes[7]-1).L2Normalize(), g1.mid(6*size, 7*size-1).length(), p, d);
#endif
#ifdef ENABLE_BETTA3_OPTIMIZATION
        IGradient::Gradient(&s, 0.01, x1, g1, indexes[7], indexes[8]-1); std::cout << labels[7]; IPrinter::print(g1.mid(indexes[7], indexes[8]-1).L2Normalize(), g1.mid(indexes[7], indexes[8]-1).length(), p, d);
        IGradient::Gradient(&s, 0.01, x1, g1, indexes[8], indexes[9]-1); std::cout << labels[8]; IPrinter::print(g1.mid(indexes[8], indexes[9]-1).L2Normalize(), g1.mid(indexes[8], indexes[9]-1).length(), p, d);
#endif
        std::cout << "---" << std::endl;
#ifdef ENABLE_NOMIN1_OPTIMIZATION
        IGradient::Gradient(&s, 0.01, x1, g1, indexes[9], indexes[10]-1); std::cout << labels[9]; IPrinter::print(g1.mid(indexes[9], indexes[10]-1).L2Normalize(), g1.mid(indexes[9], indexes[10]-1).length(), p, d);
#endif
#ifdef ENABLE_NOMIN2_OPTIMIZATION
        IGradient::Gradient(&s, 0.01, x1, g1, indexes[10], indexes[11]-1); std::cout << labels[10]; IPrinter::print(g1.mid(indexes[10], indexes[11]-1).L2Normalize(), g1.mid(indexes[10], indexes[11]-1).length(), p, d);
#endif

        x1.clear();

        puts("-------------------------- ENABLE CHECKING GRADIENTS -------------------------\n");
#endif
    }
#endif


#ifdef ENABLE_PROGRAM_CONTROL

    DoubleVector x2(606, 2.0);
    for (size_t i=0; i<s.heatSourceNumber; i++)
    {
        for (size_t ln=0; ln<101; ln++)
        {
            s.sourceParams[i].q[ln] = 1.0;
            s.sourceParams[i].v[ln].x = 2.0;
            s.sourceParams[i].v[ln].y = 2.0;
        }
    }

#else
    DoubleVector x2;
    s.parameterToVector(x2);
    //    IPrinter::printSeperatorLine();
    //    for (size_t i=0; i<11; i++)
    //    {
    //        if (i==3 || i==9 || i==12) std::cout << "---" << std::endl;
    //        std::cout << labels[i]; IPrinter::print(x2.mid(indexes[i], indexes[i+1]-1), x2.mid(indexes[i], indexes[i+1]-1).length(), 8, 4);
    //    }
    //    IPrinter::printSeperatorLine();
#endif
    //exit(-1);

    //ConjugateGradient g;
    SteepestDescentGradient g;
    g.setFunction(&s);
    g.setGradient(&s);
    g.setPrinter(&s);
    g.setProjection(&s);
    //g.setGradientNormalizer(&prob);
    g.setOptimalityTolerance(0.0000001);
    g.setStepTolerance(0.0000001);
    g.setFunctionTolerance(0.0000001);
    g.setR1MinimizeEpsilon(0.01, 0.001);
    g.setNormalize(true);
    g.showExitMessage(true);
    g.setMaxIterationCount(10);
    s.gradMethod = &g;

    g.calculate(x2);

    IPrinter::printSeperatorLine(nullptr, '=');

    //    IPrinter::printSeperatorLine();
    //    for (size_t i=0; i<11; i++)
    //    {
    //        if (i==3 || i==9 || i==12) std::cout << "---" << std::endl;
    //        std::cout << labels[i]; IPrinter::print(x2.mid(indexes[i], indexes[i+1]-1), x2.mid(indexes[i], indexes[i+1]-1).length(), 8, 4);
    //    }
    //    IPrinter::printSeperatorLine();

    s.drawImage = true;
    s.forward.findUMode = false;
    s.vectorToParameter(x2);
    s.forward.implicit_calculate_D2V1();
    s.drawImage = false;
}

void Solver1::Main(int argc, char **argv)
{
    optimize(argc, argv);
}

Solver1::Solver1(const Dimension &timeDimension, const Dimension &spaceDimensionX, const Dimension &spaceDimensionY)
{
    this->_timeDimension = timeDimension;
    this->_spaceDimensionX = spaceDimensionX;
    this->_spaceDimensionY = spaceDimensionY;

    V.resize(spaceDimensionX.size(), spaceDimensionY.size());
}

Solver1::~Solver1() {}

void Solver1::setPointNumber(size_t heatSourceNumber, size_t measrPointNumber)
{
    this->heatSourceNumber = heatSourceNumber;
    // q
#ifdef ENABLE_PROGRAM_CONTROL

#else
    this->measrPointNumber = measrPointNumber;

    alpha1.resize(heatSourceNumber, measrPointNumber, 0.197);
    alpha2.resize(heatSourceNumber, measrPointNumber, 0.173);
    alpha3.resize(heatSourceNumber, measrPointNumber, 0.165);

    // v
    bettaX1.resize(heatSourceNumber, measrPointNumber, 0.013);
    bettaX2.resize(heatSourceNumber, measrPointNumber, 0.015);
    bettaX3.resize(heatSourceNumber, measrPointNumber, 0.019);

    bettaY1.resize(heatSourceNumber, measrPointNumber, 0.014);
    bettaY2.resize(heatSourceNumber, measrPointNumber, 0.011);
    bettaY3.resize(heatSourceNumber, measrPointNumber, 0.018);

    nomnU1.resize(heatSourceNumber, measrPointNumber, 0.0500);
    nomnU2.resize(heatSourceNumber, measrPointNumber, 0.0400);

    for (size_t i=0; i<heatSourceNumber; i++)
    {
        for (size_t j=0; j<measrPointNumber; j++)
        {
            alpha1[i][j] *= 1.0 + (sin(5.0*(i+1)) + cos(4.0*(j+1))) * 0.10;
            alpha2[i][j] *= 1.0 + (sin(5.0*(i+1)) + cos(4.0*(j+1))) * 0.10;
            alpha3[i][j] *= 1.0 + (sin(5.0*(i+1)) + cos(4.0*(j+1))) * 0.10;

            bettaX1[i][j] *= 1.0 + (sin(5.0*(i+1)) + cos(4.0*(j+1))) * 0.10;
            bettaY1[i][j] *= 1.0 + (sin(5.0*(i+1)) + cos(4.0*(j+1))) * 0.10;
            bettaX2[i][j] *= 1.0 + (sin(5.0*(i+1)) + cos(4.0*(j+1))) * 0.10;
            bettaY2[i][j] *= 1.0 + (sin(5.0*(i+1)) + cos(4.0*(j+1))) * 0.10;
            bettaX3[i][j] *= 1.0 + (sin(5.0*(i+1)) + cos(4.0*(j+1))) * 0.10;
            bettaY3[i][j] *= 1.0 + (sin(5.0*(i+1)) + cos(4.0*(j+1))) * 0.10;

            nomnU1[i][j] *= 1.0 + (sin(5.0*(i+1)) + cos(4.0*(j+1))) * 0.10;
            nomnU2[i][j] *= 1.0 + (sin(5.0*(i+1)) + cos(4.0*(j+1))) * 0.10;
        }
    }

    this->measurePoints = new SpacePoint[measrPointNumber];
    this->measurePoints[0] = SpacePoint(0.50, 0.20);
    this->measurePoints[1] = SpacePoint(0.20, 0.50);
    this->measurePoints[2] = SpacePoint(0.50, 0.80);
    this->measurePoints[3] = SpacePoint(0.80, 0.50);
#endif

    sourceParams = new ProblemParams[timeDimension().size()];
    const size_t time_size = timeDimension().size();
    for (size_t ln=0; ln<time_size; ln++)
    {
        ProblemParams & pp = sourceParams[ln];
        pp.q = new double[heatSourceNumber];
        pp.v = new SpacePoint[heatSourceNumber];

        pp.z  = new SpacePointX[heatSourceNumber];
        pp.f  = new SpacePointX[heatSourceNumber];
        pp.p  = new SpacePoint[heatSourceNumber];
        pp.u  = new SpacePoint[measrPointNumber];
    }
}

double Solver1::frw_initial(const SpaceNodePDE &, InitialCondition) const { return _initialTemperature; }

double Solver1::bcw_final(const SpaceNodePDE &sn, FinalCondition) const
{
    size_t i = static_cast<size_t>(sn.i);
    size_t j = static_cast<size_t>(sn.j);
    return -2.0*(U[j][i]-V[j][i]);
}

double Solver1::frw_boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryConditionPDE &cn) const
{
    //cn = BoundaryConditionPDE::Dirichlet(1.0, 1.0); return 0.0;
    //    cn = BoundaryConditionPDE::Neumann(1.0, 0.0); return 0.0;
    cn = BoundaryConditionPDE::Robin(_lambda1, -1.0, _lambda1); return _environmentTemperature;
}

double Solver1::bcw_boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryConditionPDE &cn) const
{
    //    cn = BoundaryConditionPDE::Neumann(1.0, 0.0); return 0.0;
    cn = BoundaryConditionPDE::Robin(_lambda1, -1.0, 0.0); return 0.0;
}

void Solver1::frw_saveToImage(const DoubleMatrix &u UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
#ifdef USE_LIB_IMAGING
    static double MIN = +10000.0;
    static double MAX = -10000.0;
    if (u.max()>MAX) MAX = u.max();
    if (u.min()<MIN) MIN = u.min();
    printf(">>> %6d %14.6f %14.6f %14.6f %14.6f\n", tn.i, u.min(), u.max(), MIN, MAX);

    //if (tn.i%10==0)
    {
        QString filename = QString("data/problem3P/f/png/%1.png").arg(tn.i, 8, 10, QChar('0'));
        QPixmap pixmap;
        //visualizeMatrixHeat(u, 0.0, 46.052, pixmap, 101, 101);
        //visualizeMatrixHeat(u, 1.50, 2.1807, pixmap, spaceDimensionX().size(), spaceDimensionY().size());
        //visualizeMatrixHeat(u, 0.50, 2.32, pixmap, spaceDimensionX().size(), spaceDimensionY().size());
        visualizeMatrixHeat(u, u.min(), u.max(), pixmap, spaceDimensionX().size(), spaceDimensionY().size());

        //        QPainter painter(&pixmap);
        //        painter.setFont(QFont("Consolas", 12));
        //        painter.setPen(Qt::white);

        //        painter.drawText(2, 30, QString("time: ")+QString::number(tn.t, 'f', 3));
        //        painter.drawText(2, 60, QString("temp: ")+QString::number(u.max(), 'f', 3)+QString(" max temp: ")+QString::number(MAX, 'f', 3));


        pixmap.save(filename);
    }
#endif
}

void Solver1::bcw_saveToImage(const DoubleMatrix &p UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
{
#ifdef USE_LIB_IMAGING
    static double MIN = +10000.0;
    static double MAX = -10000.0;
    if (p.max()>MAX) MAX = p.max();
    if (p.min()<MIN) MIN = p.min();
    //printf(">>> %6d %14.6f %14.6f %14.6f %14.6f\n", tn.i, p.min(), p.max(), MIN, MAX);

    //if (tn.i%10==0)
    {
        QString filename = QString("data/problem3P/b/png/%1.png").arg(tn.i, 8, 10, QChar('0'));
        QPixmap pixmap;
        //visualizeMatrixHeat(u, 0.0, 46.052, pixmap, 101, 101);
        visualizeMatrixHeat(p, p.min(), p.max(), pixmap, spaceDimensionX().size(), spaceDimensionY().size());
        pixmap.save(filename);
    }
#endif
}

auto Solver1::bcw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    double fx = 0.0;

#ifdef ENABLE_PROGRAM_CONTROL
#else
    int ln = static_cast<int>(tn.i);
    PointNodeODE node; node.i = ln; node.x = tn.t;
    const ProblemParams &pp = sourceParams[ln];

    for (size_t j=0; j<measrPointNumber; j++)
    {
        const SpacePoint &mp = measurePoints[j];

        double sum = 0.0;
        for (size_t i=0; i<heatSourceNumber; i++)
        {
            //const SpacePoint &zi = pp.zv[i];
            const SpacePointX &zi = pp.z[i];
            const SpacePoint &pi = pp.p[i];
            const SpacePoint &fi = pp.f[i];
            if (isPointOnPlate(zi))
            {
                double dist = sqrt((zi.x - mp.x)*(zi.x - mp.x) + (zi.y - mp.y)*(zi.y - mp.y));
                sum += (alpha1[i][j]*dist*dist + alpha2[i][j]*dist + alpha3[i][j]) * ( pi.z );

                double betaXij = bettaX1[i][j]*dist*dist + bettaX2[i][j]*dist + bettaX3[i][j];
                double betaYij = bettaY1[i][j]*dist*dist + bettaY2[i][j]*dist + bettaY3[i][j];

                sum += (A4(node, 1, 1, i+1)*fi.x + A4(node, 2, 1, i+1)*fi.y) * betaXij
                        + (A4(node, 1, 2, i+1)*fi.x + A4(node, 2, 2, i+1)*fi.y) * betaYij;
            }
            else
            {
#ifdef ENABLE_ERROR_PRINT
                //throw std::runtime_error("bcw_f: unknown error...");
                std::cout << "bcw_f" << tn.t << ", " << zi.x << ", " << zi.y << std::endl;
#endif
            }
        }
        fx -= sum * DeltaFunction::gaussian(sn, mp, SpacePoint(spaceDimensionX().step()*_factor, spaceDimensionY().step()*_factor));
    }
#endif

    return fx;
}

bool Solver1::isPointOnPlate(const SpacePoint &z) const
{
    return (0.05 <= z.x && z.x <= 0.95) && (0.05 <= z.y && z.y <= 0.95);
}

void Solver1::bcw_layerInfo(const DoubleMatrix &p, const TimeNodePDE &tn) const
{
    HeatEquationFBVP &const_backward = const_cast<HeatEquationFBVP&>(this->backward);
    const int ln = static_cast<int>(tn.i);
    const double ht = timeDimension().step();

    if (ln == timeDimension().max())
    {
        DoubleVector f0(4);
        PointNodeODE n0;
        n0.i = static_cast<int>(ln); n0.x = n0.i*ht;
        ProblemParams &pp0 = sourceParams[ln];

        for (size_t i=0; i<heatSourceNumber; i++)
        {
            const_backward.i = i+1;

            /**********************************************************************************************************/

            const_backward.start( f0, n0 );
            pp0.f[i] = SpacePointX(f0);

            /**********************************************************************************************************/

            const SpacePointX &zi = pp0.z[i];
            SpacePoint &pi = pp0.p[i];
            if (isPointOnPlate(zi))
            {
                SpacePoint d = { 0.0, 0.0, 0.0 };
                pi.z = DeltaFunction::lumpedPoint4(p, zi, spaceDimensionX(), spaceDimensionY(), d);
                pi.x = d.x; pi.y = d.y;
            }
            else
            {
                pi.x = pi.y = pi.z = 0.0;
#ifdef ENABLE_ERROR_PRINT
                std::cout << "bcw_layerInfo" << tn.t << ", " << zi.x << ", " << zi.y << std::endl;
#endif
            }

            /**********************************************************************************************************/
        }
    }

    if (ln != timeDimension().min())
    {
        DoubleVector f0(4), f1(4);
        PointNodeODE n0, n1;
        n0.i = static_cast<int>(ln+0); n0.x = n0.i*ht;
        n1.i = static_cast<int>(ln-1); n1.x = n1.i*ht;
        ProblemParams &pp0 = sourceParams[ln];
        ProblemParams &pp1 = sourceParams[ln-1];

        for (size_t i=0; i<heatSourceNumber; i++)
        {
            const_backward.i = i+1;

            pp0.f[i].toDoubleVector(f0);
            const_backward.next( f0, n0, f1, n1, method );
            pp1.f[i] = SpacePointX(f1);

            /**********************************************************************************************************/

            const SpacePointX &zi = pp1.z[i];
            SpacePoint &pi = pp1.p[i];
            if (isPointOnPlate(zi))
            {
                SpacePoint d = { 0.0, 0.0, 0.0 };
                pi.z = DeltaFunction::lumpedPoint4(p, zi, spaceDimensionX(), spaceDimensionY(), d);
                pi.x = d.x; pi.y = d.y;
            }
            else
            {
                pi.x = pi.y = pi.z = 0.0;
#ifdef ENABLE_ERROR_PRINT
                std::cout << "bcw_layerInfo" << tn.t << ", " << zi.x << ", " << zi.y << std::endl;
#endif
            }

            /**********************************************************************************************************/
        }
    }

    //bcw_saveToImage(p, tn);
}

double Solver1::A1(const PointNodeODE &, size_t r, size_t c, size_t i) const
{
#if defined (TEST_PROBLEM_1) || defined (TEST_PROBLEM_2)
    return 0.0;
#else
    if (i==1) { double data[2][2] = { { +0.24, +0.00 }, { +0.00, +0.35 } }; return data[r-1][c-1]; }
    if (i==2) { double data[2][2] = { { -0.12, +0.00 }, { +0.00, -0.19 } }; return data[r-1][c-1]; }
#endif
}

double Solver1::A2(const PointNodeODE &, size_t r, size_t c, size_t i) const
{
#if defined (TEST_PROBLEM_1) || defined (TEST_PROBLEM_2)
    return 0.0;
#else
    if (i==1) { double data[2][2] = { { +0.00, +0.00 }, { +0.00, +0.00 } }; return data[r-1][c-1]; }
    if (i==2) { double data[2][2] = { { +0.00, +0.00 }, { +0.00, +0.00 } }; return data[r-1][c-1]; }
#endif
}

double Solver1::A3(const PointNodeODE &, size_t r, size_t i) const
{
#if defined (TEST_PROBLEM_1) || defined (TEST_PROBLEM_2)
    return 0.0;
#else
    if (i==1) { double data[2] = { +0.00, +0.00 }; return data[r-1]; }
    if (i==2) { double data[2] = { +0.00, +0.00 }; return data[r-1]; }
#endif
}

double Solver1::A4(const PointNodeODE &, size_t r, size_t c, size_t i) const
{
#if defined (TEST_PROBLEM_1) || defined (TEST_PROBLEM_2)
    if (i==1) { double data[2][2] = { { +1.00, +0.00 }, { +0.00, +1.00 } }; return data[r-1][c-1]; }
    if (i==2) { double data[2][2] = { { +1.00, +0.00 }, { +0.00, +1.00 } }; return data[r-1][c-1]; }
#else
    if (i==1) { double data[2][2] = { { +0.79, +0.00 }, { +0.00, +0.60 } }; return data[r-1][c-1]; }
    if (i==2) { double data[2][2] = { { -0.84, +0.00 }, { +0.00, +0.83 } }; return data[r-1][c-1]; }
#endif
}

auto Solver1::frw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    size_t ln = static_cast<size_t>(tn.i);

    const ProblemParams &pp = sourceParams[ln];
    double fx = _lambda0 * _environmentTemperature;

    double sum = 0.0;
    for (size_t i=0; i<heatSourceNumber; i++)
    {
        const SpacePointX &zi = pp.z[i];
        if ( isPointOnPlate(zi) )
        {
            sum += pp.q[i] * DeltaFunction::gaussian(sn, zi, SpacePoint(spaceDimensionX().step()*_factor, spaceDimensionY().step()*_factor));
        }
        else
        {
#ifdef ENABLE_ERROR_PRINT
            printf("Solver1::frw_f:  %.4f %14.10f, %14.10f\n", tn.t, zi.x, zi.y);
#endif
        }
    }

    return fx + sum;
}

void Solver1::frw_layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const
{
    HeatEquationIBVP &const_forward = const_cast<HeatEquationIBVP&>(this->forward);
    const int ln = static_cast<int>(tn.i);
    const double ht = timeDimension().step();

    if (ln == timeDimension().max()) { const_cast<Solver1*>(this)->U = u; }

    if (ln == timeDimension().min())
    {
        PointNodeODE n0;
        DoubleVector z0(4);
        n0.i = static_cast<int>(ln); n0.x = n0.i*ht;
        ProblemParams &pp0 = sourceParams[ln];

        for (size_t i=0; i<heatSourceNumber; i++)
        {
            const_forward.i = i+1;

            /**********************************************************************************************************/

            forward.start( z0, n0 );
            pp0.z[i] = SpacePointX(z0);

            /**********************************************************************************************************/

            const SpacePointX &z = pp0.z[i];
#ifdef ENABLE_PROGRAM_CONTROL
#else
            double qi  = 0.0;
            double vXi = 0.0;
            double vYi = 0.0;
            for (size_t j=0; j<measrPointNumber; j++)
            {
                double nominal1 = nomnU1[i][j];
                double nominal2 = nomnU2[i][j];
                SpacePoint d;
                pp0.u[j].z = DeltaFunction::lumpedPoint4(u, measurePoints[j], spaceDimensionX(), spaceDimensionY(), d);
                pp0.u[j].x = d.x; pp0.u[j].y = d.y;

                const SpacePoint &mp = measurePoints[j];
                double dist = sqrt((z.x - mp.x)*(z.x - mp.x) + (z.y - mp.y)*(z.y - mp.y));
                qi  += (alpha1[i][j]*dist*dist + alpha2[i][j]*dist + alpha3[i][j]) * (pp0.u[j].z - nominal1);
                vXi += (bettaX1[i][j]*dist*dist + bettaX2[i][j]*dist + bettaX3[i][j]) * (pp0.u[j].z - nominal2);
                vYi += (bettaY1[i][j]*dist*dist + bettaY2[i][j]*dist + bettaY3[i][j]) * (pp0.u[j].z - nominal2);
            }
            pp0.q[i]  = qi;
            pp0.vX[i] = vXi;
            pp0.vY[i] = vYi;
#endif

            if (!isPointOnPlate(z))
            {
                pp0.q[i]  = 0.0;
                //pp0.vX[i] = 0.0;
                //pp0.vY[i] = 0.0;
#ifdef ENABLE_ERROR_PRINT
                std::cout << std::setprecision (15) << "frw_layerInfo" << tn.t << ", " << z.x << ", " << z.y << std::endl;
#endif
            }

            if (forward.findUMode)
            {
                const double t = tn.t;
#if defined(TEST_PROBLEM_2)
                pp0.q[i] = forward.q(i,tn);
                PointNodeODE node; node.i = ln; node.x = node.i*ht;
                forward.v(i, node, pp0.v[i]);
#endif
            }

            /**********************************************************************************************************/
        }
    }

    if (ln != timeDimension().max())
    {
        DoubleVector z0(4), z1(4);
        PointNodeODE n0, n1;
        n0.i = static_cast<int>(ln+0); n0.x = n0.i*ht;
        n1.i = static_cast<int>(ln+1); n1.x = n1.i*ht;
        ProblemParams &pp0 = sourceParams[ln];
        ProblemParams &pp1 = sourceParams[ln+1];

        for (size_t i=0; i<heatSourceNumber; i++)
        {
            const_forward.i = i+1;

            pp0.z[i].toDoubleVector(z0);

            forward.next( z0, n0, z1, n1, method );
            pp1.z[i] = SpacePointX(z1);

            /**********************************************************************************************************/

            const SpacePointX &zi = pp1.z[i];
#ifdef ENABLE_PROGRAM_CONTROL
#else
            double qi  = 0.0;
            double vXi = 0.0;
            double vYi = 0.0;

            for (size_t j=0; j<measrPointNumber; j++)
            {
                double nominal1 = nomnU1[i][j];
                double nominal2 = nomnU2[i][j];

                const SpacePoint &mp = measurePoints[j];
                double dist = sqrt((zi.x - mp.x)*(zi.x - mp.x) + (zi.y - mp.y)*(zi.y - mp.y));
                qi  += (alpha1[i][j]*dist*dist + alpha2[i][j]*dist + alpha3[i][j]) * (pp0.u[j].z - nominal1);
                vXi += (bettaX1[i][j]*dist*dist + bettaX2[i][j]*dist + bettaX3[i][j]) * (pp0.u[j].z - nominal2);
                vYi += (bettaY1[i][j]*dist*dist + bettaY2[i][j]*dist + bettaY3[i][j]) * (pp0.u[j].z - nominal2);

                SpacePoint d;
                pp1.u[j].z = DeltaFunction::lumpedPoint4(u, measurePoints[j], spaceDimensionX(), spaceDimensionY(), d);
                pp1.u[j].x = d.x; pp0.u[j].y = d.y;
            }
            pp1.q[i]  = qi;
            pp1.vX[i] = vXi;
            pp1.vY[i] = vYi;
#endif

            if (!isPointOnPlate(zi))
            {
                pp1.q[i]  = 0.0;
                //pp1.vX[i] = 0.0;
                //pp1.vY[i] = 0.0;
#ifdef ENABLE_ERROR_PRINT
                std::cout << std::setprecision (15) << "frw_layerInfo" << tn.t << ", " << zi.x << ", " << zi.y << std::endl;
#endif
            }

            if (forward.findUMode)
            {
                const double t = tn.t;
#if defined(TEST_PROBLEM_2)
                pp1.q[i] = forward.q(i,tn);
                PointNodeODE node; node.i = ln; node.x = node.i*ht;
                forward.v(i, node, pp1.v[i]);
#endif
            }

            /**********************************************************************************************************/
        }
    }

    if (drawImage) frw_saveToImage(u, tn);

    if (ln == timeDimension().max())
    {
        const_cast<Solver1*>(this)->U = u;

        if (drawImage) { printf("q1: "); for (unsigned int i=0; i<=ln; i++) printf("%10.4f", sourceParams[i].q[0]); puts(""); }
        if (drawImage) { printf("q2: "); for (unsigned int i=0; i<=ln; i++) printf("%10.4f", sourceParams[i].q[1]); puts(""); }
        if (drawImage) { printf("v1: "); for (unsigned int i=0; i<=ln; i++) printf("%10.4f", sourceParams[i].v[0].x); puts(""); }
        if (drawImage) { printf("v2: "); for (unsigned int i=0; i<=ln; i++) printf("%10.4f", sourceParams[i].v[1].x); puts(""); }
        if (drawImage) { printf("z1: "); for (unsigned int i=0; i<=ln; i++) printf("%10.4f", sourceParams[i].z[0].x); puts(""); }
        if (drawImage) { printf("z2: "); for (unsigned int i=0; i<=ln; i++) printf("%10.4f", sourceParams[i].z[1].x); puts(""); }
    }
}

void Solver1::gradient(const DoubleVector &x, DoubleVector &g) const
{
    const int min = timeDimension().min();
    const int max = timeDimension().max();
    const double ht = timeDimension().step();
    const size_t time_size = timeDimension().size();
    const_cast<Solver1*>(this)->vectorToParameter(x);

    g.clear();
    g.resize(x.length(), 0.0);

    for (unsigned int i1=0; i1<_size1; i1++)
    {
        const_cast<Solver1*>(this)->_initialTemperature = _initialTemperatureList[i1];

        for (unsigned int i2=0; i2<_size2; i2++)
        {
            const_cast<Solver1*>(this)->_environmentTemperature = _environmentTemperatureList[i2];

            forward.implicit_calculate_D2V1();
            backward.implicit_calculate_D2V1();

            DoubleVector g0(x.length(), 0.0);


#ifdef ENABLE_PROGRAM_CONTROL
            for (size_t i=0; i<heatSourceNumber; i++)
            {
                size_t offset = i*time_size;
                for (unsigned ln=0; ln<time_size; ln++)
                {
                    PointNodeODE node; node.i = static_cast<int>(ln); node.x = node.i*ht;

                    const SpacePoint  &pi = sourceParams[ln].p[i];
                    const SpacePointX &fi = sourceParams[ln].f[i];
                    g0[offset+0*time_size+ln] = -pi.z;
                    g0[offset+1*time_size+ln] = -(fi.x*A4(node, 1, 1, i+1) + fi.y*A4(node, 1, 2, i+1));
                    g0[offset+2*time_size+ln] = -(fi.x*A4(node, 2, 1, i+1) + fi.y*A4(node, 2, 2, i+1));

                    g[offset+0*time_size+ln] = g0[offset+0*time_size+ln] * _ratio1 * _ratio2;
                    g[offset+1*time_size+ln] = g0[offset+1*time_size+ln] * _ratio1 * _ratio2;
                    g[offset+2*time_size+ln] = g0[offset+2*time_size+ln] * _ratio1 * _ratio2;
                }
            }
#else
            const size_t size = heatSourceNumber*measrPointNumber;

            for (size_t i=0; i<heatSourceNumber; i++)
            {
                for (size_t j=0; j<measrPointNumber; j++)
                {
                    const SpacePoint &mp = measurePoints[j];

                    // alpha
                    g0[0*size + i*measrPointNumber  + j] = 0.0;
                    g0[1*size + i*measrPointNumber  + j] = 0.0;
                    g0[2*size + i*measrPointNumber  + j] = 0.0;

                    // betta
                    g0[3*size + i*measrPointNumber  + j] = 0.0;
                    g0[4*size + i*measrPointNumber  + j] = 0.0;
                    g0[5*size + i*measrPointNumber  + j] = 0.0;
                    g0[6*size + i*measrPointNumber  + j] = 0.0;
                    g0[7*size + i*measrPointNumber  + j] = 0.0;
                    g0[8*size + i*measrPointNumber  + j] = 0.0;

                    // nomnU
                    g0[ 9*size + i*measrPointNumber  + j] = 0.0;
                    g0[10*size + i*measrPointNumber  + j] = 0.0;

                    PointNodeODE node; node.i = min; node.x = node.i*ht;

                    const SpacePointX &zi0 = sourceParams[min].z[i];
                    const SpacePointX &fi0 = sourceParams[min].f[i];
                    const SpacePoint  &uj0 = sourceParams[min].u[j];
                    const SpacePoint  &pi0 = sourceParams[min].p[i];

                    if (isPointOnPlate(zi0))
                    {
                        double dist = sqrt((zi0.x - mp.x)*(zi0.x - mp.x) + (zi0.y - mp.y)*(zi0.y - mp.y));
                        double dif1 = uj0.z-nomnU1[i][j];
                        double dif2 = uj0.z-nomnU2[i][j];
                        double val0 = pi0.z;
                        double valX1 = A4(node, 1, 1, i+1)*fi0.x + A4(node, 2, 1, i+1)*fi0.y;
                        double valY1 = A4(node, 1, 2, i+1)*fi0.x + A4(node, 2, 2, i+1)*fi0.y;

#ifdef ENABLE_ALPHA1_OPTIMIZATION
                        g0[0*size + i*measrPointNumber + j] += 0.5*val0*dif1*dist*dist;
#endif
#ifdef ENABLE_ALPHA2_OPTIMIZATION
                        g0[1*size + i*measrPointNumber + j] += 0.5*val0*dif1*dist;
#endif
#ifdef ENABLE_ALPHA3_OPTIMIZATION
                        g0[2*size + i*measrPointNumber + j] += 0.5*val0*dif1;
#endif
#ifdef ENABLE_BETTA1_OPTIMIZATION
                        g0[3*size + i*measrPointNumber + j] += 0.5*valX1*dif2*dist*dist;
                        g0[4*size + i*measrPointNumber + j] += 0.5*valY1*dif2*dist*dist;
#endif
#ifdef ENABLE_BETTA2_OPTIMIZATION
                        g0[5*size + i*measrPointNumber + j] += 0.5*valX1*dif2*dist;
                        g0[6*size + i*measrPointNumber + j] += 0.5*valY1*dif2*dist;
#endif
#ifdef ENABLE_BETTA3_OPTIMIZATION
                        g0[7*size + i*measrPointNumber + j] += 0.5*valX1*dif2;
                        g0[8*size + i*measrPointNumber + j] += 0.5*valY1*dif2;
#endif
#ifdef ENABLE_NOMIN1_OPTIMIZATION
                        g0[9*size + i*measrPointNumber + j] += 0.5 * val0 * ( alpha1[i][j]*dist*dist + alpha2[i][j]*dist + alpha3[i][j] );
#endif
#ifdef ENABLE_NOMIN2_OPTIMIZATION
                        double betaXij = bettaX1[i][j]*dist*dist + bettaX2[i][j]*dist + bettaX3[i][j];
                        double betaYij = bettaY1[i][j]*dist*dist + bettaY2[i][j]*dist + bettaY3[i][j];
                        g0[10*size + i*measrPointNumber + j] += 0.5 * (valX1 * betaXij + valY1 * betaYij);
#endif
                    }
                    else
                    {
#ifdef ENABLE_ERROR_PRINT
                        std::cout << std::setprecision (15) << "gradient" << node.x << ", " << zi0.x << ", " << zi0.y << std::endl;
#endif
                    }

                    for (int ln=min+1; ln<max-1; ln++)
                    {
                        node.i = ln; node.x = node.i*ht;

                        const SpacePointX &zi = sourceParams[ln].z[i];
                        const SpacePointX &fi = sourceParams[ln].f[i];
                        const SpacePoint  &uj = sourceParams[ln].u[j];
                        const SpacePoint  &pi = sourceParams[ln].p[i];

                        if (isPointOnPlate(zi))
                        {
                            double dist = sqrt((zi.x - mp.x)*(zi.x - mp.x) + (zi.y - mp.y)*(zi.y - mp.y));
                            double dif1 = uj.z-nomnU1[i][j];
                            double dif2 = uj.z-nomnU2[i][j];
                            double val0 = pi.z;
                            double valX1 = A4(node, 1, 1, i+1)*fi.x + A4(node, 2, 1, i+1)*fi.y;
                            double valY1 = A4(node, 1, 2, i+1)*fi.x + A4(node, 2, 2, i+1)*fi.y;

#ifdef ENABLE_ALPHA1_OPTIMIZATION
                            g0[0*size + i*measrPointNumber + j] += val0*dif1*dist*dist;
#endif
#ifdef ENABLE_ALPHA2_OPTIMIZATION
                            g0[1*size + i*measrPointNumber + j] += val0*dif1*dist;
#endif
#ifdef ENABLE_ALPHA3_OPTIMIZATION
                            g0[2*size + i*measrPointNumber + j] += val0*dif1;
#endif
#ifdef ENABLE_BETTA1_OPTIMIZATION
                            g0[3*size + i*measrPointNumber + j] += valX1*dif2*dist*dist;
                            g0[4*size + i*measrPointNumber + j] += valY1*dif2*dist*dist;
#endif
#ifdef ENABLE_BETTA2_OPTIMIZATION
                            g0[5*size + i*measrPointNumber + j] += valX1*dif2*dist;
                            g0[6*size + i*measrPointNumber + j] += valY1*dif2*dist;
#endif
#ifdef ENABLE_BETTA3_OPTIMIZATION
                            g0[7*size + i*measrPointNumber + j] += valX1*dif2;
                            g0[8*size + i*measrPointNumber + j] += valY1*dif2;
#endif
#ifdef ENABLE_NOMIN1_OPTIMIZATION
                            g0[9*size + i*measrPointNumber + j] += val0 * ( alpha1[i][j]*dist*dist + alpha2[i][j]*dist + alpha3[i][j] );
#endif
#ifdef ENABLE_NOMIN2_OPTIMIZATION
                            double betaXij = bettaX1[i][j]*dist*dist + bettaX2[i][j]*dist + bettaX3[i][j];
                            double betaYij = bettaY1[i][j]*dist*dist + bettaY2[i][j]*dist + bettaY3[i][j];
                            g0[10*size + i*measrPointNumber + j] += valX1 * betaXij + valY1 * betaYij;
#endif
                        }
                        else
                        {
#ifdef ENABLE_ERROR_PRINT
                            std::cout << "gradient" << std::setprecision (15) << node.x << ", " << zi0.x << ", " << zi0.y << std::endl;
#endif
                        }
                    }

                    node.i = max; node.x = node.i*ht;

                    const SpacePointX &zi1 = sourceParams[max].z[i];
                    const SpacePointX &fi1 = sourceParams[max].f[i];
                    const SpacePoint  &uj  = sourceParams[max].u[j];
                    const SpacePoint  &pi  = sourceParams[max].p[i];

                    if (isPointOnPlate(zi1))
                    {
                        double dist = sqrt((zi1.x - mp.x)*(zi1.x - mp.x) + (zi1.y - mp.y)*(zi1.y - mp.y));
                        double dif1 = uj.z-nomnU1[i][j];
                        double dif2 = uj.z-nomnU2[i][j];
                        double val0 = pi.z;
                        double valX1 = A4(node, 1, 1, i+1)*fi1.x + A4(node, 2, 1, i+1)*fi1.y;
                        double valY1 = A4(node, 1, 2, i+1)*fi1.x + A4(node, 2, 2, i+1)*fi1.y;

#ifdef ENABLE_ALPHA1_OPTIMIZATION
                        g0[0*size + i*measrPointNumber + j] += 0.5*val0*dif1*dist*dist;
#endif
#ifdef ENABLE_ALPHA2_OPTIMIZATION
                        g0[1*size + i*measrPointNumber + j] += 0.5*val0*dif1*dist;
#endif
#ifdef ENABLE_ALPHA3_OPTIMIZATION
                        g0[2*size + i*measrPointNumber + j] += 0.5*val0*dif1;
#endif
#ifdef ENABLE_BETTA1_OPTIMIZATION
                        g0[3*size + i*measrPointNumber + j] += 0.5*valX1*dif2*dist*dist;
                        g0[4*size + i*measrPointNumber + j] += 0.5*valY1*dif2*dist*dist;
#endif
#ifdef ENABLE_BETTA2_OPTIMIZATION
                        g0[5*size + i*measrPointNumber + j] += 0.5*valX1*dif2*dist;
                        g0[6*size + i*measrPointNumber + j] += 0.5*valY1*dif2*dist;
#endif
#ifdef ENABLE_BETTA3_OPTIMIZATION
                        g0[7*size + i*measrPointNumber + j] += 0.5*valX1*dif2;
                        g0[8*size + i*measrPointNumber + j] += 0.5*valY1*dif2;
#endif
#ifdef ENABLE_NOMIN1_OPTIMIZATION
                        g0[9*size + i*measrPointNumber + j] += 0.5 * val0 * ( alpha1[i][j]*(dist*dist) + alpha2[i][j]*dist + alpha3[i][j] );
#endif
#ifdef ENABLE_NOMIN2_OPTIMIZATION
                        double betaXij = bettaX1[i][j]*dist*dist + bettaX2[i][j]*dist + bettaX3[i][j];
                        double betaYij = bettaY1[i][j]*dist*dist + bettaY2[i][j]*dist + bettaY3[i][j];
                        g0[10*size + i*measrPointNumber + j] += 0.5 * (valX1 * betaXij + valY1 * betaYij);
#endif
                    }
                    else
                    {
#ifdef ENABLE_ERROR_PRINT
                        std::cout << "gradient" << std::setprecision (15) << node.x << ", " << zi0.x << ", " << zi0.y << std::endl;
#endif
                    }

#ifdef ENABLE_ALPHA1_OPTIMIZATION
                    g0[0*size + i*measrPointNumber + j] *= -ht;
#endif
#ifdef ENABLE_ALPHA2_OPTIMIZATION
                    g0[1*size + i*measrPointNumber + j] *= -ht;
#endif
#ifdef ENABLE_ALPHA3_OPTIMIZATION
                    g0[2*size + i*measrPointNumber + j] *= -ht;
#endif
#ifdef ENABLE_BETTA1_OPTIMIZATION
                    g0[3*size + i*measrPointNumber + j] *= +ht;
                    g0[4*size + i*measrPointNumber + j] *= +ht;
#endif
#ifdef ENABLE_BETTA2_OPTIMIZATION
                    g0[5*size + i*measrPointNumber + j] *= +ht;
                    g0[6*size + i*measrPointNumber + j] *= +ht;
#endif
#ifdef ENABLE_BETTA3_OPTIMIZATION
                    g0[7*size + i*measrPointNumber + j] *= +ht;
                    g0[8*size + i*measrPointNumber + j] *= +ht;
#endif
#ifdef ENABLE_NOMIN1_OPTIMIZATION
                    g0[9*size + i*measrPointNumber + j] *= +ht;
#endif
#ifdef ENABLE_NOMIN2_OPTIMIZATION
                    g0[10*size + i*measrPointNumber + j] *= -ht;
#endif
                }
            }

            for (size_t i=0; i<heatSourceNumber; i++)
            {

                for (size_t j=0; j<measrPointNumber; j++)
                {
                    // alpha
                    g[0*size + i*measrPointNumber  + j] += g0[0*size + i*measrPointNumber  + j] * _ratio1 * _ratio2;
                    g[1*size + i*measrPointNumber  + j] += g0[1*size + i*measrPointNumber  + j] * _ratio1 * _ratio2;
                    g[2*size + i*measrPointNumber  + j] += g0[2*size + i*measrPointNumber  + j] * _ratio1 * _ratio2;

                    // betta
                    g[3*size + i*measrPointNumber  + j] += g0[3*size + i*measrPointNumber  + j] * _ratio1 * _ratio2;
                    g[4*size + i*measrPointNumber  + j] += g0[4*size + i*measrPointNumber  + j] * _ratio1 * _ratio2;
                    g[5*size + i*measrPointNumber  + j] += g0[5*size + i*measrPointNumber  + j] * _ratio1 * _ratio2;
                    g[6*size + i*measrPointNumber  + j] += g0[6*size + i*measrPointNumber  + j] * _ratio1 * _ratio2;
                    g[7*size + i*measrPointNumber  + j] += g0[7*size + i*measrPointNumber  + j] * _ratio1 * _ratio2;
                    g[8*size + i*measrPointNumber  + j] += g0[8*size + i*measrPointNumber  + j] * _ratio1 * _ratio2;

                    // nomnU
                    g[9*size + i*measrPointNumber  + j] += g0[9*size + i*measrPointNumber  + j] * _ratio1 * _ratio2;
                    g[10*size + i*measrPointNumber  + j] += g0[10*size + i*measrPointNumber  + j] * _ratio1 * _ratio2;
                }
            }
#endif
        }
    }
}

double Solver1::fx(const DoubleVector &x) const
{
    const_cast<Solver1*>(this)->vectorToParameter(x);
    double sum = 0.0;

    for (unsigned int i1=0; i1<_size1; i1++)
    {
        const_cast<Solver1*>(this)->_initialTemperature = _initialTemperatureList[i1];

        for (unsigned int i2=0; i2<_size2; i2++)
        {
            const_cast<Solver1*>(this)->_environmentTemperature = _environmentTemperatureList[i2];

            forward.implicit_calculate_D2V1();

            sum += integral(U) * _ratio1 * _ratio2;


            double penalty = 0.0;
            unsigned int time_size = timeDimension().size()-1;
            double time_step = timeDimension().step();
            for (size_t i=0; i<heatSourceNumber; i++)
            {
                double penalty_i = 0.0;

                double g0x = fabs(sourceParams[0].z[i].x - 0.5) - 0.5; if (g0x <= 0) g0x = 0.0;
                double g0y = fabs(sourceParams[0].z[i].y - 0.5) - 0.5; if (g0y <= 0) g0y = 0.0;
                penalty_i += 0.5*(g0x*g0x+g0y*g0y);
                for (unsigned int ln=1; ln<time_size; ln++)
                {
                    double gx = fabs(sourceParams[ln].z[i].x - 0.5) - 0.5; if (gx <= 0) gx = 0.0;
                    double gy = fabs(sourceParams[ln].z[i].y - 0.5) - 0.5; if (gy <= 0) gy = 0.0;
                    penalty_i += (gx*gx+gy*gy);
                }
                double g1x = fabs(sourceParams[time_size].z[i].x - 0.5) - 0.5; if (g1x <= 0) g1x = 0.0;
                double g1y = fabs(sourceParams[time_size].z[i].y - 0.5) - 0.5; if (g1y <= 0) g1y = 0.0;
                penalty_i += 0.5*(g1x*g1x+g1y*g1y);

                penalty_i *= time_step;

                penalty += penalty_i;
            }

            sum += R*penalty;
        }
    }

    return sum;
}

auto Solver1::integral(const DoubleMatrix &) const -> double
{
    const Dimension &dimensionX = spaceDimensionX();
    const Dimension &dimensionY = spaceDimensionY();
    const unsigned int N = static_cast<unsigned int> ( dimensionX.size()-1 );
    const unsigned int M = static_cast<unsigned int> ( dimensionY.size()-1 );
    const double hx = dimensionX.step();
    const double hy = dimensionY.step();

    double udiff = 0.0;
    double usum = 0.0;

    udiff = (U[0][0]-V[0][0]); usum += 0.25 * udiff * udiff;// * mu(0, 0);
    udiff = (U[0][N]-V[0][N]); usum += 0.25 * udiff * udiff;// * mu(N, 0);
    udiff = (U[M][0]-V[M][0]); usum += 0.25 * udiff * udiff;// * mu(0, M);
    udiff = (U[M][N]-V[M][N]); usum += 0.25 * udiff * udiff;// * mu(N, M);

    for (unsigned int n=1; n<=N-1; n++)
    {
        udiff = U[0][n]-V[0][n]; usum += 0.5 * udiff * udiff;// * mu(n, 0);
        udiff = U[M][n]-V[M][n]; usum += 0.5 * udiff * udiff;// * mu(n, M);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        udiff = U[m][0]-V[m][0]; usum += 0.5 * udiff * udiff;// * mu(0, m);
        udiff = U[m][N]-V[m][N]; usum += 0.5 * udiff * udiff;// * mu(N, m);
    }

    for (unsigned int m=1; m<=M-1; m++)
    {
        for (unsigned int n=1; n<=N-1; n++)
        {
            udiff = U[m][n]-V[m][n]; usum += udiff * udiff;// * mu(n, m);
        }
    }

    return usum*(hx*hy);
}

auto Solver1::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientBasedMethod::MethodResult result) const -> void
{
    printf_s("%4d %16.8f %16.8f\n", i, f, alpha);
    //if (alpha > 0.0)
    //const_cast<Solver1*>(this)->gradMethod->setR1MinimizeEpsilon(alpha, alpha/100.0);
    //if (fabs(alpha) < 0.0001)
    //const_cast<Solver1*>(this)->gradMethod->setR1MinimizeEpsilon(0.1, 0.001);
}

auto Solver1::project(DoubleVector &x) const -> void
{

}

auto Solver1::project(DoubleVector &x, size_t /*index*/) -> void {}

//--------------------------------------------------------------------------------------------------------------//

HeatEquationIBVP::HeatEquationIBVP(Solver1 *solver) { this->solver = solver; }

HeatEquationIBVP::HeatEquationIBVP(const HeatEquationIBVP &) {}

HeatEquationIBVP & HeatEquationIBVP::operator =(const HeatEquationIBVP &) { return *this; }

HeatEquationIBVP::~HeatEquationIBVP() {}

auto HeatEquationIBVP::initial(const SpaceNodePDE &sn, InitialCondition cn) const -> double { return solver->frw_initial(sn, cn); }

auto HeatEquationIBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &cn) const -> double { return solver->frw_boundary(sn, tn, cn); }

auto HeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const-> double { return solver->frw_f(sn, tn); }

auto HeatEquationIBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const -> void { solver->frw_layerInfo(u, tn); }

auto HeatEquationIBVP::timeDimension() const -> Dimension { return solver->timeDimension(); }

auto HeatEquationIBVP::spaceDimensionX() const -> Dimension { return solver->spaceDimensionX(); }

auto HeatEquationIBVP::spaceDimensionY() const -> Dimension { return solver->spaceDimensionY(); }

auto HeatEquationIBVP::spaceDimensionZ() const -> Dimension { throw std::runtime_error("Not supportted..."); }

auto HeatEquationIBVP::q(size_t i, const TimeNodePDE &tn) const -> double
{
    //if (i==0) return pow(2.7182818284590452353602874713527, -0.2*tn.t);
    //if (i==1) return pow(2.7182818284590452353602874713527, -0.3*tn.t);
    //    return pow(2.7182818284590452353602874713527, -tn.t);
    return 0.25;
}

auto HeatEquationIBVP::v(size_t i, const PointNodeODE &tn, SpacePoint &vl) const -> void
{

    if (i==0) { vl.x = +1.44*(M_PI*M_PI)*sin(2.0*M_PI*tn.x); vl.y = 0.0; }
    if (i==1) { vl.x = -1.44*(M_PI*M_PI)*sin(2.0*M_PI*tn.x); vl.y = 0.0; }
}

//--------------------------------------------------------------------------------------------------------------//

auto HeatEquationIBVP::A(const PointNodeODE &node, size_t r, size_t c) const -> double
{
    return solver->A1(node, r, c, i);
}

auto HeatEquationIBVP::B(const PointNodeODE &node, size_t r, size_t c) const -> double
{
    return solver->A2(node, r, c, i);
}

auto HeatEquationIBVP::C(const PointNodeODE &node, size_t r) const -> double
{
    const size_t ln = static_cast<size_t>(node.i);
    return solver->A3(node, r, i) + ( solver->A4(node, r, 1, i) * solver->sourceParams[ln].v[i-1].x
                                    + solver->A4(node, r, 2, i) * solver->sourceParams[ln].v[i-1].y );
}

auto HeatEquationIBVP::initial(InitialCondition c, size_t r) const -> double
{
    //if (findUMode)
    {
        if (c == InitialCondition::InitialValue)
        {
#if defined (TEST_PROBLEM_1)
            if (i==1) { const double data[2] = { +0.50, +0.15 }; return data[r-1]; }
            if (i==2) { const double data[2] = { +0.92, +0.40 }; return data[r-1]; }
#endif
#if defined (TEST_PROBLEM_2)
            if (i==1) { const double data[2] = { +0.50, +0.25 }; return data[r-1]; }
            if (i==2) { const double data[2] = { +0.50, +0.75 }; return data[r-1]; }
#endif
        }
        if (c == InitialCondition::InitialFirstDerivative)
        {
#if defined (TEST_PROBLEM_1)
            if (i==1) { const double data[2] = { -1.68*M_PI, +0.00*M_PI }; return data[r-1]; }
            if (i==2) { const double data[2] = { +0.00*M_PI, +1.00*M_PI }; return data[r-1]; }
#endif
#if defined (TEST_PROBLEM_2)
            if (i==1) { const double data[2] = { -0.72*M_PI, +0.00 }; return data[r-1]; }
            if (i==2) { const double data[2] = { +0.72*M_PI, +0.00 }; return data[r-1]; }
#endif
        }
    }
    //    else
    //    {
    //        if (c == InitialCondition::InitialValue)
    //        {
    //            //if (i==1) { const double data[2] = { 0.04, 0.04 }; return data[r-1]; }
    //            //if (i==2) { const double data[2] = { 0.96, 0.04 }; return data[r-1]; }
    //            if (i==1) { const double data[2] = { 0.10, 0.10 }; return data[r-1]; }
    //            if (i==2) { const double data[2] = { 0.90, 0.10 }; return data[r-1]; }
    //        }
    //        if (c == InitialCondition::InitialFirstDerivative)
    //        {
    //            //        if (i==1) { const double data[2] = { +0.35, +0.56 }; return data[r-1]; }
    //            //        if (i==2) { const double data[2] = { -0.54, +0.42 }; return data[r-1]; }
    //            if (i==1) { const double data[2] = { +0.00, +0.00 }; return data[r-1]; }
    //            if (i==2) { const double data[2] = { -0.00, +0.00 }; return data[r-1]; }
    //        }
    //    }
    throw std::runtime_error(std::string("auto HeatEquationIBVP::initial(InitialCondition, size_t row) const -> double"));
}

auto HeatEquationIBVP::count() const -> size_t { return 2; }

auto HeatEquationIBVP::dimension() const -> Dimension { return solver->timeDimension(); }

auto HeatEquationIBVP::iterationInfo(const DoubleVector &, const PointNodeODE &) const -> void {}

//--------------------------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------------------------//

HeatEquationFBVP::HeatEquationFBVP(Solver1 *solver) { this->solver = solver; }

HeatEquationFBVP::HeatEquationFBVP(const HeatEquationFBVP &) {}

HeatEquationFBVP & HeatEquationFBVP::operator =(const HeatEquationFBVP &) { return *this; }

HeatEquationFBVP::~HeatEquationFBVP() {}

auto HeatEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition cn) const -> double { return solver->bcw_final(sn, cn); }

auto HeatEquationFBVP::boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &cn) const -> double { return solver->bcw_boundary(sn, tn, cn); }

auto HeatEquationFBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const-> double { return solver->bcw_f(sn, tn); }

auto HeatEquationFBVP::layerInfo(const DoubleMatrix &p, const TimeNodePDE &tn) const -> void { solver->bcw_layerInfo(p, tn); }

auto HeatEquationFBVP::timeDimension() const -> Dimension { return solver->timeDimension(); }

auto HeatEquationFBVP::spaceDimensionX() const -> Dimension { return solver->spaceDimensionX(); }

auto HeatEquationFBVP::spaceDimensionY() const -> Dimension { return solver->spaceDimensionY(); }

auto HeatEquationFBVP::spaceDimensionZ() const -> Dimension { throw std::runtime_error(""); }

//--------------------------------------------------------------------------------------------------------------//

auto HeatEquationFBVP::A(const PointNodeODE &node, size_t r, size_t c) const -> double { return -solver->A1(node, r, c, i); }

auto HeatEquationFBVP::B(const PointNodeODE &node, size_t r, size_t c) const -> double { return +solver->A2(node, r, c, i); }

auto HeatEquationFBVP::C(const PointNodeODE &node, size_t r) const -> double
{
    const size_t ln = static_cast<size_t>(node.i);
    const ProblemParams &pp = solver->sourceParams[ln];
    double sum = 0.0;

    const SpacePointX &zi = pp.z[i-1];
    if (solver->isPointOnPlate(zi))
    {
#ifdef ENABLE_PROGRAM_CONTROL
        // grad
        if (r==1) sum = pp.p[i-1].x * pp.q[i-1];
        if (r==2) sum = pp.p[i-1].y * pp.q[i-1];
#else
        //double val0 = pp.p[i-1].z;
        double valX1 = solver->A4(node, 1, 1, i)*pp.f[i-1].x + solver->A4(node, 1, 2, i)*pp.f[i-1].y;
        double valY1 = solver->A4(node, 2, 1, i)*pp.f[i-1].x + solver->A4(node, 2, 2, i)*pp.f[i-1].y;

        for (size_t j=0; j<solver->measrPointNumber; j++)
        {
            const SpacePoint &mp = solver->measurePoints[j];

            double sum2 = 0.0, sum1 = 0.0, sum3 = 0.0;

            // grad
            if (r==1) sum1 = pp.p[i-1].x * pp.q[i-1];
            if (r==2) sum1 = pp.p[i-1].y * pp.q[i-1];
            //double dist1 = sqrt((zi.x - mp.x)*(zi.x - mp.x) + (zi.y - mp.y)*(zi.y - mp.y));
            //if (r==1) sum1 = pp.p[i-1].x * (solver->alpha1[i-1][j]*dist1*dist1 + solver->alpha2[i-1][j]*dist1 + solver->alpha3[i-1][j]) * diff1;
            //if (r==2) sum1 = pp.p[i-1].y * (solver->alpha1[i-1][j]*dist1*dist1 + solver->alpha2[i-1][j]*dist1 + solver->alpha3[i-1][j]) * diff1;

            //
            if (r==1) sum2 = valX1 * (2.0*solver->bettaX1[i-1][j]*(zi.x - mp.x) + solver->bettaX2[i-1][j]) * (pp.u[j].z-solver->nomnU2[i-1][j]);
            if (r==2) sum2 = valY1 * (2.0*solver->bettaY1[i-1][j]*(zi.y - mp.y) + solver->bettaY2[i-1][j]) * (pp.u[j].z-solver->nomnU2[i-1][j]);

            //
            if (r==1) sum3 = pp.p[i-1].z * (2.0*solver->alpha1[i-1][j]*(zi.x - mp.x) + solver->alpha2[i-1][j]) * (pp.u[j].z-solver->nomnU1[i-1][j]);
            if (r==2) sum3 = pp.p[i-1].z * (2.0*solver->alpha1[i-1][j]*(zi.y - mp.y) + solver->alpha2[i-1][j]) * (pp.u[j].z-solver->nomnU2[i-1][j]);

            sum -= (sum2 + sum1 + sum3);
        }
#endif
    }
    else
    {
#ifdef ENABLE_ERROR_PRINT
        std::cout << "HeatEquationFBVP::C" << std::setprecision (15) << node.x << ", " << zi.x << ", " << zi.y << std::endl;
#endif
    }

    double gx = fabs(zi.x-0.5)-0.5; if (gx <= 0) gx = 0.0;
    double gy = fabs(zi.y-0.5)-0.5; if (gy <= 0) gy = 0.0;
    sum -= solver->R * 2.0*(sgn(zi.x-0.5)*gx + sgn(zi.y-0.5)*gy);

    return sum;
}

auto HeatEquationFBVP::final(FinalCondition c, size_t r) const -> double
{
    if (c == FinalCondition::FinalValue)
    {
        if (i==1) { const double data[2] = { 0.00, 0.00 }; return data[r-1]; }
        if (i==2) { const double data[2] = { 0.00, 0.00 }; return data[r-1]; }
    }

    if (c == FinalCondition::FinalFirstDerivative)
    {
        const int max = timeDimension().max();
        const double ht = timeDimension().size();
        PointNodeODE node; node.i = max; node.x = max*ht;
        if (i==1) { return -( solver->A1(node, r, 1, i)*final(FinalCondition::FinalValue, 1) + solver->A1(node, r, 2, i)*final(FinalCondition::FinalValue, 2) ); }
        if (i==2) { return -( solver->A1(node, r, 1, i)*final(FinalCondition::FinalValue, 1) + solver->A1(node, r, 2, i)*final(FinalCondition::FinalValue, 2) ); }
    }

    return 0.0;
}

auto HeatEquationFBVP::count() const -> size_t { return 2; }

auto HeatEquationFBVP::dimension() const -> Dimension { return solver->timeDimension(); }

void Solver1::vectorToParameter(const DoubleVector &x)
{
#ifdef ENABLE_PROGRAM_CONTROL
    const size_t time_size = timeDimension().size();

    for (size_t i=0; i<heatSourceNumber; i++)
    {
        const size_t offset = i*time_size*3;

        for (size_t ln=0; ln<time_size; ln++)
        {
            sourceParams[i].q[ln]   = x[offset + 0*time_size + ln];
            sourceParams[i].v[ln].x = x[offset + 1*time_size + ln];
            sourceParams[i].v[ln].y = x[offset + 2*time_size + ln];
        }
    }
#else
    for (size_t i=0; i<heatSourceNumber; i++)
    {
        const size_t size = heatSourceNumber*measrPointNumber;

        for (size_t j=0; j<measrPointNumber; j++)
        {
            alpha1[i][j]  = x[0*size + i*measrPointNumber + j];
            alpha2[i][j]  = x[1*size + i*measrPointNumber + j];
            alpha3[i][j]  = x[2*size + i*measrPointNumber + j];
            bettaX1[i][j] = x[3*size + i*measrPointNumber + j];
            bettaY1[i][j] = x[4*size + i*measrPointNumber + j];
            bettaX2[i][j] = x[5*size + i*measrPointNumber + j];
            bettaY2[i][j] = x[6*size + i*measrPointNumber + j];
            bettaX3[i][j] = x[7*size + i*measrPointNumber + j];
            bettaY3[i][j] = x[8*size + i*measrPointNumber + j];
            nomnU1[i][j]  = x[9*size + i*measrPointNumber + j];
            nomnU2[i][j]  = x[10*size + i*measrPointNumber + j];
        }
    }
#endif
}

void Solver1::parameterToVector(DoubleVector &x)
{
#ifdef ENABLE_PROGRAM_CONTROL
    const auto time_size = timeDimension().size();
    const size_t length = heatSourceNumber*time_size*3;

    x.resize(length);

    for (size_t i=0; i<heatSourceNumber; i++)
    {
        const size_t offset = i*time_size*3;
        for (size_t ln=0; ln<time_size; ln++)
        {
            x[offset + 0*time_size + ln] = sourceParams[i].q[ln];
            x[offset + 1*time_size + ln] = sourceParams[i].v[ln].x;
            x[offset + 2*time_size + ln] = sourceParams[i].v[ln].y;
        }
    }
#else
    const size_t size = heatSourceNumber*measrPointNumber;
    x.clear();
    x.resize(size*11);

    for (size_t i=0; i<heatSourceNumber; i++)
    {
        for (size_t j=0; j<measrPointNumber; j++)
        {
            x[0*size + i*measrPointNumber + j] = alpha1[i][j];
            x[1*size + i*measrPointNumber + j] = alpha2[i][j];
            x[2*size + i*measrPointNumber + j] = alpha3[i][j];
            x[3*size + i*measrPointNumber + j] = bettaX1[i][j];
            x[4*size + i*measrPointNumber + j] = bettaY1[i][j];
            x[5*size + i*measrPointNumber + j] = bettaX2[i][j];
            x[6*size + i*measrPointNumber + j] = bettaY2[i][j];
            x[7*size + i*measrPointNumber + j] = bettaX3[i][j];
            x[8*size + i*measrPointNumber + j] = bettaY3[i][j];
            x[9*size + i*measrPointNumber + j] = nomnU1[i][j];
            x[10*size + i*measrPointNumber + j] = nomnU2[i][j];
        }
    }
#endif
}


