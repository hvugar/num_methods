#include "solver5.h"

using namespace p3p5;

void Functional::Main(int /*argc*/, char** /*argv*/)
{
    Functional f(0.01, 0.0, 0.01, 0.0001);

    //    std::cout << std::setprecision(6) << f.LoadedHeatEquationIBVP::thermalDiffusivity() << " "
    //              << f.LoadedHeatEquationIBVP::thermalConductivity() << " " << f.LoadedHeatEquationIBVP::thermalConvection() << " " << f.lambda1 << std::endl;
    //    std::cout << std::setprecision(6) << f.LoadedHeatEquationFBVP::thermalDiffusivity() << " "
    //              << f.LoadedHeatEquationFBVP::thermalConductivity() << " " << f.LoadedHeatEquationFBVP::thermalConvection() << " " << f.lambda1 << std::endl;

    //    f.LoadedHeatEquationIBVP::implicit_calculate_D2V1();
    //    f.LoadedHeatEquationFBVP::implicit_calculate_D2V1();

    size_t time_size = f._timeDimension.size();
    size_t heating_point_number = f.heating_point_number;


    //double x1[] = {4.281320,4.015598,4.094010,4.487715,5.165369,6.065671,7.126982,8.267653,9.400826,10.436499,11.292881,11.895477,12.192561,12.151197,11.771056,11.084692,10.143746,9.023558,7.806453,6.581709,5.432876,4.435138,3.640652,3.098134,2.827804,2.847890,3.158956,3.758098,4.625613,5.742473,7.057111,8.517783,10.044295,11.545338,12.919116,14.068094,14.898064,15.341176,15.350635,14.920771,14.088449,12.916095,11.497830,9.937732,8.349308,6.838864,5.502558,4.406904,3.615646,3.156957,3.058021,3.326407,3.965251,4.957713,6.285717,7.889860,9.706647,11.634315,13.554807,15.334914,16.845764,17.961839,18.592386,18.672837,18.190862,17.188110,15.739137,13.959683,11.979232,9.940281,7.977120,6.212277,4.730436,3.611580,2.889952,2.594968,2.732480,3.299355,4.270740,5.618373,7.271003,9.151123,11.139858,13.099967,14.877105,16.321468,17.290284,17.683081,17.436737,16.553218,15.101373,13.189094,10.974168,8.622295,6.301977,4.156853,2.300565,0.792480,-0.334286,-1.091902,-1.508038,-1.618225,-1.462646,-1.076040,-0.488491,0.261422,1.132275,2.060794,2.971168,3.773892,4.379822,4.706178,4.705529,4.362403,3.709629,2.826567,1.817596,0.814210,-0.063084,-0.717807,-1.097588,-1.198663,-1.070461,-0.793532,-0.468378,-0.183869,-0.001721,0.056426,0.008469,-0.091928,-0.189081,-0.228202,-0.186857,-0.073392,0.072186,0.192901,0.225130,0.132871,-0.098378,-0.447421,-0.860947,-1.271062,-1.594623,-1.765061,-1.745233,-1.541470,-1.203301,-0.824181,-0.508780,-0.364731,-0.460826,-0.815055,-1.388230,-2.094890,-2.818320,-3.451729,-3.909174,-4.154431,-4.193381,-4.070107,-3.847800,-3.601367,-3.375115,-3.205211,-3.103266,-3.059040,-3.049975,-3.033423,-2.969122,-2.827788,-2.601918,-2.304797,-1.981321,-1.678693,-1.455804,-1.347259,-1.360482,-1.474490,-1.643368,-1.803473,-1.912552,-1.945172,-1.921696,-1.886483,-1.889900,-1.956254,-2.072045,-2.168048,-2.170378,-2.027666,-1.740088,-1.370868,-0.999354,-0.693550,-0.463723,-0.255163,0.014587,0.346803,0.605990,0.482308,-0.258083};
    //double x2[] = {3.158009,3.184737,3.233357,3.302677,3.387367,3.482846,3.581328,3.679663,3.768907,3.845587,3.901365,3.932497,3.937601,3.914573,3.863813,3.788365,3.690446,3.574667,3.450023,3.320842,3.196287,3.080140,2.980047,2.900626,2.846679,2.818525,2.816534,2.839257,2.885022,2.948258,3.023320,3.102096,3.180397,3.248726,3.302995,3.334617,3.339575,3.316043,3.262359,3.179218,3.069188,2.935191,2.781939,2.618157,2.447508,2.278792,2.115349,1.962862,1.826106,1.707721,1.608317,1.527951,1.465341,1.419639,1.388185,1.366480,1.351376,1.341406,1.332127,1.323747,1.312707,1.298839,1.281927,1.264286,1.247407,1.231604,1.219140,1.211147,1.210983,1.216489,1.228873,1.245157,1.262592,1.280468,1.295013,1.302047,1.298945,1.283192,1.254290,1.212668,1.156715,1.088962,1.013840,0.933207,0.854155,0.778812,0.712310,0.657749,0.621251,0.605551,0.609400,0.633992,0.677127,0.738125,0.809102,0.885778,0.959606,1.021730,1.068001,1.090251,1.082005,1.039741,0.961162,0.847738,0.703545,0.529702,0.334411,0.126715,-0.086015,-0.291611,-0.483410,-0.652437,-0.793747,-0.898019,-0.961500,-0.987186,-0.973682,-0.925022,-0.844709,-0.743866,-0.630273,-0.515067,-0.411991,-0.326785,-0.271832,-0.255568,-0.283298,-0.358595,-0.481286,-0.645358,-0.850081,-1.083412,-1.333371,-1.588268,-1.829641,-2.044465,-2.216782,-2.336144,-2.384123,-2.351057,-2.237013,-2.036326,-1.752663,-1.390966,-0.967690,-0.499067,-0.006668,0.481381,0.946955,1.362130,1.700628,1.943587,2.075154,2.087622,1.983687,1.761311,1.437981,1.037585,0.583361,0.113844,-0.344427,-0.757655,-1.102206,-1.344889,-1.468641,-1.477577,-1.367288,-1.150847,-0.842006,-0.476174,-0.080762,0.308629,0.654673,0.939638,1.137432,1.233002,1.224681,1.118357,0.930891,0.691085,0.415098,0.136679,-0.113603,-0.320798,-0.458328,-0.529295,-0.534665,-0.493377,-0.405221,-0.286442,-0.167386,-0.049582,0.058734,0.163545,0.239074,0.263544,0.177669,-0.088068,-0.510898,-0.695372,1.048845};

    double x1[] = {
        6.982833  ,4.577156,3.224840,2.832113,3.317968,4.563161,6.440510,8.769984,11.338750,13.898802,16.196012,17.973891,19.039924,19.251916,18.565636,17.039717,14.801229,12.056757,9.033209,5.964676,3.052944,0.462788,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,0.757563,3.870492,7.415416,11.153133,14.771154,17.932641,20.290553,21.591786,21.660065,20.464494,18.125139,14.859536,10.997945,6.886950,2.864585,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,0.651492,4.791701,9.209678,13.524614,17.316072,20.147186,21.715827,21.827217,20.484657,17.889998,14.366513,10.366032,6.331520,2.653958,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,0.170150,3.305231,7.085932,11.228951,15.285245,18.725299,20.986488,21.693148,20.645729,17.951486,14.024282,9.458181,5.012131,1.365083,-0.000000,-0.000000,0.578824,2.477937,4.000119,4.680164,4.072860,1.970710,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,0.449955,2.643053,2.897678,0.961346,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,3.547128,5.296285,3.721917,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,1.730139,3.306135,3.087334,2.536876,1.197331,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,0.191389,0.010250,0.147556,0.102606
    };
    double x2[] = {
        2.648017,2.505403,2.385275,2.286378,2.206906,2.134341,2.065679,1.999303,1.924437,1.844439,1.749854,1.640449,1.513638,1.377635,1.236273,1.087955,0.941280,0.800618,0.677136,0.566082,0.474525,0.398638,0.331141,0.277158,0.225669,0.168990,0.101107,0.015165,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,0.210449,0.649465,1.049184,1.381743,1.629234,1.777462,1.822246,1.776312,1.634469,1.423643,1.175140,0.907776,0.668697,0.475093,0.356439,0.324481,0.412618,0.627599,0.949634,1.376897,1.885948,2.462944,3.054808,3.631668,4.143529,4.538760,4.800168,4.885951,4.769434,4.441963,3.903096,3.175123,2.298490,1.290438,0.214635,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,0.063525,1.129638,2.174282,3.101880,3.866244,4.369444,4.488198,4.210208,3.591814,2.575161,1.190783,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,0.261323,2.416994,4.173365,5.597915,6.284769,6.108468,5.025880,3.234277,0.971175,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,0.061726,1.529046,2.761664,3.857612,3.795101,2.592491,0.368039,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,0.016863,0.192386,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000,-0.000000
    };

    DoubleVector x(heating_point_number*time_size);
//    for (size_t i=0; i<heating_point_number; i++)
//    {
//        for (size_t ln=0; ln<time_size; ln++)
//        {
//            x[i*time_size+ln] = 0.0;//0.05*0.001*ln;
//        }
//    }

    for (size_t ln=0; ln<time_size; ln++)
    {
        x[0*time_size+ln] = x1[ln];
        x[1*time_size+ln] = x2[ln];
    }


    //    IPrinter::printVector(x.mid(   0, 1000));
    //    IPrinter::printVector(x.mid(1001, 2001));
    //    IPrinter::printVector(f.pv.mid(   0, 1000));
    //    IPrinter::printVector(f.pv.mid(1001, 2001));

    //    DoubleVector g1(x.length());
    //    DoubleVector g2(x.length(), 0.0);

    //    f.gradient(x, g1);
    //    IPrinter::printVector(g1.mid(   0, 1000).L2Normalize());
    //    IPrinter::printVector(g1.mid(1001, 2001).L2Normalize());
    //    IPrinter::printSeperatorLine();

    //    std::vector<size_t> i1 = {   0,  100,  200,  300,  400,  500,  600,  700,  800,  900, 1000};
    //    std::vector<size_t> i2 = {1001, 1101, 1201, 1301, 1401, 1501, 1601, 1701, 1801, 1901, 2001};
    //    IGradient::Gradient(__NUM_GRAD_STEP__, &f, x, g2, i1); g2[   0] *= 2; g2[1000] *= 2;
    //    IPrinter::printVector(g2.mid(   0, 1000).L2Normalize());
    //    IGradient::Gradient(__NUM_GRAD_STEP__, &f, x, g2, i2); g2[1001] *= 2; g2[2001] *= 2;
    //    IPrinter::printVector(g2.mid(1001, 2001).L2Normalize());


    //f.gm = new ConjugateGradient;       f.gm->setNormalize(false);
    f.gm = new SteepestDescentGradient; f.gm->setNormalize(true);
    f.gm->setFunction(&f);
    f.gm->setGradient(&f);
    f.gm->setPrinter(&f);
    f.gm->setProjection(&f);
    f.gm->setOptimalityTolerance(0.0);
    f.gm->setStepTolerance(0.0);
    f.gm->setFunctionTolerance(0.0);
    f.gm->setR1MinimizeEpsilon(1.0, 0.1);
    //f.gm->setMaxIterationCount(100);
    f.gm->showExitMessage(true);
    f.gm->calculate(x);
    delete f.gm;
}

Functional::Functional(double diffusivity, double conductivity, double convection, double lambda):
    LoadedHeatEquationIBVP(+diffusivity, conductivity, -convection),
    LoadedHeatEquationFBVP(-diffusivity, conductivity, +convection)
{
    this->lambda1 = lambda;
}

auto Functional::fx(const DoubleVector &x) const -> double
{
    const_cast<Functional*>(this)->qv = x;
    LoadedHeatEquationIBVP::implicit_calculate_D2V1();
    double integralU = integral(U);
    return integralU;
}

auto Functional::integral(const DoubleMatrix &) const -> double
{
    const size_t N = static_cast<size_t> ( _spaceDimensionX.size()-1 );
    const size_t M = static_cast<size_t> ( _spaceDimensionY.size()-1 );
    const double hx = _spaceDimensionX.step();
    const double hy = _spaceDimensionY.step();

    double udiff = 0.0;
    double usum = 0.0;

    udiff = (U[0][0]-V[0][0]); usum += 0.25 * udiff * udiff;// * mu(0, 0);
    udiff = (U[0][N]-V[0][N]); usum += 0.25 * udiff * udiff;// * mu(N, 0);
    udiff = (U[M][0]-V[M][0]); usum += 0.25 * udiff * udiff;// * mu(0, M);
    udiff = (U[M][N]-V[M][N]); usum += 0.25 * udiff * udiff;// * mu(N, M);

    for (size_t n=1; n<=N-1; n++)
    {
        udiff = U[0][n]-V[0][n]; usum += 0.5 * udiff * udiff;// * mu(n, 0);
        udiff = U[M][n]-V[M][n]; usum += 0.5 * udiff * udiff;// * mu(n, M);
    }

    for (size_t m=1; m<=M-1; m++)
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

auto Functional::gradient(const DoubleVector &x, DoubleVector &g) const -> void
{
    const_cast<Functional*>(this)->qv = x;
    g.resize(x.length());
    LoadedHeatEquationIBVP::implicit_calculate_D2V1();
    LoadedHeatEquationFBVP::implicit_calculate_D2V1();

    size_t time_size = _timeDimension.size();
    for (size_t i=0; i<heating_point_number; i++)
    {
        for (size_t ln=0; ln<time_size; ln++)
        {
            /*if (ln%100==0)*/ g[i*time_size+ln] = -pv[i*time_size+ln];
        }
    }
}

auto Functional::project(DoubleVector &x, size_t index) -> void
{
    if (x[index] < -0.0) x[index] = -0.0;
}

auto Functional::print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientBasedMethod::MethodResult /*result*/) const -> void
{
    const_cast<Functional*>(this)->drawImages = true;
    const_cast<Functional*>(this)->qv = x;
    LoadedHeatEquationIBVP::implicit_calculate_D2V1();
    const_cast<Functional*>(this)->drawImages = false;

    if (f <= 0.1) { const_cast<Functional*>(this)->gm->setR1MinimizeEpsilon(0.100, 0.0100); }
    //if (fabs(f - lastFx) <= 0.10000) { const_cast<Functional*>(this)->gm->setR1MinimizeEpsilon(0.100, 0.0100); }
    //if (fabs(f - lastFx) <= 0.00100) { const_cast<Functional*>(this)->gm->setR1MinimizeEpsilon(0.010, 0.0010); }
    //if (fabs(f - lastFx) <= 0.00001) { const_cast<Functional*>(this)->gm->setR1MinimizeEpsilon(0.001, 0.0001); }


    const_cast<Functional*>(this)->lastFx = f;

    size_t time_size = static_cast<size_t>(_timeDimension.max());
    printf("I[%d]: %10.6f | %10.6f %10.6f | %10.6f %10.6f | %10.6f\n", iteration, f, U.min(), U.max(), gm->min_step, gm->min_epsilon, alpha);
    IPrinter::printVector(x.mid(   0, time_size));
    IPrinter::printVector(g.mid(   0, time_size));
    IPrinter::printVector(x.mid(time_size+1, 2*time_size+1));
    IPrinter::printVector(g.mid(time_size+1, 2*time_size+1));

    if (f<=0.1)
    {
        std::cout << "---" << std::endl;
        IPrinter::printVector(x.mid(   0, time_size), nullptr, x.mid(   0, time_size).length());
        IPrinter::printVector(x.mid(time_size+1, 2*time_size+1), nullptr, x.mid(time_size+1, 2*time_size+1).length());
        std::cout << "---" << std::endl;
    }


    std::string path = "E:\\project\\hvugar\\num_methods\\trunk\\optimal\\bin\\data\\problem3P\\f\\png\\";
    system(std::string("mkdir ").append(path).append(std::to_string(iteration)).append(" >nul").data());
    system(std::string("move  ").append(path).append("*.png ").append(path).append(std::to_string(iteration)).append("\\ >nul").data());
}

/*********************************************************************************************************************************************/

LoadedHeatEquationIBVP::LoadedHeatEquationIBVP(double diffusivity, double conductivity, double convection)
{
    setThermalDiffusivity(diffusivity);
    setThermalConductivity(conductivity);
    setThermalConvection(convection);
}

LoadedHeatEquationIBVP::LoadedHeatEquationIBVP(const LoadedHeatEquationIBVP &lhe) : IHeatEquationIBVP(lhe) {}

LoadedHeatEquationIBVP::~LoadedHeatEquationIBVP() {}

LoadedHeatEquationIBVP& LoadedHeatEquationIBVP::operator =(const LoadedHeatEquationIBVP &) { throw std::exception(); }

auto LoadedHeatEquationIBVP::initial(const SpaceNodePDE &/*sn*/, InitialCondition /*ic*/) const -> double { return initial_temperature; }

auto LoadedHeatEquationIBVP::boundary(const SpaceNodePDE &/*sn*/, const TimeNodePDE &/*tn*/, BoundaryConditionPDE &bc) const -> double
{
    //bc = BoundaryConditionPDE::Dirichlet(); return 0.0;
    //bc = BoundaryConditionPDE::Neumann(); return 0.0;
    bc = BoundaryConditionPDE::Robin(lambda1, +1.0); return lambda1 * envrmnt_temperature;
}

auto LoadedHeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    //const double envrmnt_temperature = _functional->envrmnt_temperature;
    //const double lambda0 = _functional->lambda0;
    const double dimX_step = spaceDimensionX().step();
    const double dimY_step = spaceDimensionY().step();
    //const size_t heating_point_number = _functional->heating_point_number;

    double fx = -thermalConvection() * envrmnt_temperature;

    double sum = 0.0;
    for (size_t i=0; i<heating_point_number; i++)
    {
        const SpacePoint &zi = z(tn, i);
        const double w = DeltaFunction::gaussian(sn, zi, SpacePoint(dimX_step, dimY_step));
        if (w > 0.0)
        {
            const double qi = q(tn, i);
            sum += qi * w;
        }
        //sum += qi * DeltaFunction::nearest(sn, zi, dimX_step, dimY_step, dimX_size, dimY_size);
    }

    return fx + sum;
}

auto LoadedHeatEquationIBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const -> void
{
    LoadedHeatEquationIBVP* lhe = const_cast<LoadedHeatEquationIBVP*>(this);
    if (static_cast<int>(tn.i) == timeDimension().max()) { lhe->Shared::U = u; }

    if (drawImages)
    {
        //printf(">>> %6d %14.6f %14.6f\n", tn.i, u.min(), u.max());

        QString filename = QString("data/problem3P/f/png/%1.png").arg(tn.i, 8, 10, QChar('0'));
        QPixmap pixmap;
        visualizeMatrixHeat(u, u.min(), u.max(), pixmap, 101, 101);
        pixmap.save(filename);
    }

    //static double MIN = +10000.0;
    //static double MAX = -10000.0;
    //if (u.max()>MAX) MAX = u.max();
    //if (u.min()<MIN) MIN = u.min();
    //printf(">>> %6d %14.6f %14.6f %14.6f %14.6f\n", tn.i, u.min(), u.max(), MIN, MAX);

    //if (tn.i%10==0)
    //    {
    //QString filename = QString("data/problem3P/f/png/%1.png").arg(tn.i, 8, 10, QChar('0'));
    //QPixmap pixmap;
    //visualizeMatrixHeat(u, 0.0, 3.864, pixmap, 101, 101);
    //visualizeMatrixHeat(u, 1.50, 2.1807, pixmap, spaceDimensionX().size(), spaceDimensionY().size());
    //visualizeMatrixHeat(u, 0.50, 2.32, pixmap, spaceDimensionX().size(), spaceDimensionY().size());
    //visualizeMatrixHeat(u, u.min(), u.max(), pixmap, spaceDimensionX().size(), spaceDimensionY().size());

    //QPainter painter(&pixmap);
    //painter.setFont(QFont("Consolas", 12));
    //painter.setPen(Qt::white);

    //painter.drawText(2, 30, QString("time: ")+QString::number(tn.t, 'f', 3));
    //painter.drawText(2, 60, QString("temp: ")+QString::number(u.max(), 'f', 3)+QString(" max temp: ")+QString::number(MAX, 'f', 3));
    //pixmap.save(filename);
    //    }
}

auto LoadedHeatEquationIBVP::timeDimension() const -> Dimension { return Shared::_timeDimension; }

auto LoadedHeatEquationIBVP::spaceDimensionX() const -> Dimension { return Shared::_spaceDimensionX; }

auto LoadedHeatEquationIBVP::spaceDimensionY() const -> Dimension { return Shared::_spaceDimensionY; }

auto LoadedHeatEquationIBVP::spaceDimensionZ() const -> Dimension { return Shared::_spaceDimensionZ; }

/*********************************************************************************************************************************************/

LoadedHeatEquationFBVP::LoadedHeatEquationFBVP(double diffusivity, double conductivity, double convection)
{
    setThermalDiffusivity(diffusivity);
    setThermalConductivity(conductivity);
    setThermalConvection(convection);
}

LoadedHeatEquationFBVP::LoadedHeatEquationFBVP(const LoadedHeatEquationFBVP &) : IHeatEquationFBVP() { }

LoadedHeatEquationFBVP::~LoadedHeatEquationFBVP() {}

LoadedHeatEquationFBVP& LoadedHeatEquationFBVP::operator =(const LoadedHeatEquationFBVP &) { return *this; }

auto LoadedHeatEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition /*ic*/) const -> double
{
    size_t j = static_cast<size_t>(sn.j), i = static_cast<size_t>(sn.i);
    return -2.0 * ( U[j][i] - V[j][i] );
}

auto LoadedHeatEquationFBVP::boundary(const SpaceNodePDE &/*sn*/, const TimeNodePDE &/*tn*/, BoundaryConditionPDE &bc) const -> double
{
    //bc = BoundaryConditionPDE::Dirichlet(); return 0.0;
    //bc = BoundaryConditionPDE::Neumann(); return 0.0;
    bc = BoundaryConditionPDE::Robin(lambda1, +1.0); return 0.0;
}

auto LoadedHeatEquationFBVP::f(const SpaceNodePDE &/*sn*/, const TimeNodePDE &/*tn*/) const -> double
{
    return 0.0;
}

auto LoadedHeatEquationFBVP::layerInfo(const DoubleMatrix &p, const TimeNodePDE &tn) const -> void
{
    LoadedHeatEquationFBVP* lhe = const_cast<LoadedHeatEquationFBVP*>(this);
    size_t time_size = timeDimension().size();
    for (size_t i=0; i<heating_point_number; i++)
    {
        const SpacePoint &zi = z(tn, i);
        lhe->Shared::pv[i*time_size+tn.i] = DeltaFunction::lumpedPointG(p, zi, _spaceDimensionX, _spaceDimensionY, 1, 4);
    }

    //static double MIN = +10000.0;
    //static double MAX = -10000.0;
    //if (p.max()>MAX) MAX = p.max();
    //if (p.min()<MIN) MIN = p.min();
    //printf(">>> %6d %14.6f %14.6f %14.6f %14.6f\n", tn.i, p.min(), p.max(), MIN, MAX);

    //if (tn.i%10==0)
    //    {
    //QString filename = QString("data/problem3P/b/png/%1.png").arg(tn.i, 8, 10, QChar('0'));
    //QPixmap pixmap;
    //visualizeMatrixHeat(p, p.min(), p.max(), pixmap, 101, 101);
    //visualizeMatrixHeat(u, 1.50, 2.1807, pixmap, spaceDimensionX().size(), spaceDimensionY().size());
    //visualizeMatrixHeat(u, 0.50, 2.32, pixmap, spaceDimensionX().size(), spaceDimensionY().size());
    //visualizeMatrixHeat(u, u.min(), u.max(), pixmap, spaceDimensionX().size(), spaceDimensionY().size());

    //                QPainter painter(&pixmap);
    //                painter.setFont(QFont("Consolas", 12));
    //                painter.setPen(Qt::white);

    //                painter.drawText(2, 30, QString("time: ")+QString::number(tn.t, 'f', 3));
    //                painter.drawText(2, 60, QString("temp: ")+QString::number(u.max(), 'f', 3)+QString(" max temp: ")+QString::number(MAX, 'f', 3));
    //pixmap.save(filename);
    //    }
}

auto LoadedHeatEquationFBVP::timeDimension() const -> Dimension { return Shared::_timeDimension; }

auto LoadedHeatEquationFBVP::spaceDimensionX() const -> Dimension { return Shared::_spaceDimensionX; }

auto LoadedHeatEquationFBVP::spaceDimensionY() const -> Dimension { return Shared::_spaceDimensionY; }

auto LoadedHeatEquationFBVP::spaceDimensionZ() const -> Dimension { return Shared::_spaceDimensionZ; }

/*********************************************************************************************************************************************/

auto Shared::q(const TimeNodePDE &tn, size_t i) const -> double
{
    size_t time_size = _timeDimension.size();
    return qv[i*time_size + tn.i];
    //    return 0.05*tn.t;
}

auto Shared::z(const TimeNodePDE &tn, size_t i) const -> SpacePoint
{
    SpacePoint sp;
    const double t = tn.t;
    const double v = 2.0*M_PI/**t*/;
    switch (i) {
    case 0: {
        sp.x = /*0.2*/0.4*sin(v*t) + 0.5;
        sp.y = /*0.2*/0.4*cos(v*t) + 0.5;
    } break;
    case 1: {
        sp.x = -0.2*sin(v*t) + 0.5;
        sp.y = -0.2*cos(v*t) + 0.5;
    } break;
    }

    return sp;
}
