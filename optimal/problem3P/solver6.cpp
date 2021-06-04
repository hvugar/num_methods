#include "solver6.h"

using namespace p3p6;

void Functional::Main(int /*argc*/, char** /*argv*/)
{
    Functional f(0.01, 0.0, 0.001, 0.0001);

    const size_t heating_source_number = f.heating_source_number;

#ifdef P3P6_OPTIMIZE_Q
    size_t time_size = f._timeDimension.size();

    double x1[] = {4.281320,4.015598,4.094010,4.487715,5.165369,6.065671,7.126982,8.267653,9.400826,10.436499,11.292881,11.895477,12.192561,12.151197,11.771056,11.084692,10.143746,9.023558,7.806453,6.581709,5.432876,4.435138,3.640652,3.098134,2.827804,2.847890,3.158956,3.758098,4.625613,5.742473,7.057111,8.517783,10.044295,11.545338,12.919116,14.068094,14.898064,15.341176,15.350635,14.920771,14.088449,12.916095,11.497830,9.937732,8.349308,6.838864,5.502558,4.406904,3.615646,3.156957,3.058021,3.326407,3.965251,4.957713,6.285717,7.889860,9.706647,11.634315,13.554807,15.334914,16.845764,17.961839,18.592386,18.672837,18.190862,17.188110,15.739137,13.959683,11.979232,9.940281,7.977120,6.212277,4.730436,3.611580,2.889952,2.594968,2.732480,3.299355,4.270740,5.618373,7.271003,9.151123,11.139858,13.099967,14.877105,16.321468,17.290284,17.683081,17.436737,16.553218,15.101373,13.189094,10.974168,8.622295,6.301977,4.156853,2.300565,0.792480,-0.334286,-1.091902,-1.508038,-1.618225,-1.462646,-1.076040,-0.488491,0.261422,1.132275,2.060794,2.971168,3.773892,4.379822,4.706178,4.705529,4.362403,3.709629,2.826567,1.817596,0.814210,-0.063084,-0.717807,-1.097588,-1.198663,-1.070461,-0.793532,-0.468378,-0.183869,-0.001721,0.056426,0.008469,-0.091928,-0.189081,-0.228202,-0.186857,-0.073392,0.072186,0.192901,0.225130,0.132871,-0.098378,-0.447421,-0.860947,-1.271062,-1.594623,-1.765061,-1.745233,-1.541470,-1.203301,-0.824181,-0.508780,-0.364731,-0.460826,-0.815055,-1.388230,-2.094890,-2.818320,-3.451729,-3.909174,-4.154431,-4.193381,-4.070107,-3.847800,-3.601367,-3.375115,-3.205211,-3.103266,-3.059040,-3.049975,-3.033423,-2.969122,-2.827788,-2.601918,-2.304797,-1.981321,-1.678693,-1.455804,-1.347259,-1.360482,-1.474490,-1.643368,-1.803473,-1.912552,-1.945172,-1.921696,-1.886483,-1.889900,-1.956254,-2.072045,-2.168048,-2.170378,-2.027666,-1.740088,-1.370868,-0.999354,-0.693550,-0.463723,-0.255163,0.014587,0.346803,0.605990,0.482308,-0.258083};
    double x2[] = {3.158009,3.184737,3.233357,3.302677,3.387367,3.482846,3.581328,3.679663,3.768907,3.845587,3.901365,3.932497,3.937601,3.914573,3.863813,3.788365,3.690446,3.574667,3.450023,3.320842,3.196287,3.080140,2.980047,2.900626,2.846679,2.818525,2.816534,2.839257,2.885022,2.948258,3.023320,3.102096,3.180397,3.248726,3.302995,3.334617,3.339575,3.316043,3.262359,3.179218,3.069188,2.935191,2.781939,2.618157,2.447508,2.278792,2.115349,1.962862,1.826106,1.707721,1.608317,1.527951,1.465341,1.419639,1.388185,1.366480,1.351376,1.341406,1.332127,1.323747,1.312707,1.298839,1.281927,1.264286,1.247407,1.231604,1.219140,1.211147,1.210983,1.216489,1.228873,1.245157,1.262592,1.280468,1.295013,1.302047,1.298945,1.283192,1.254290,1.212668,1.156715,1.088962,1.013840,0.933207,0.854155,0.778812,0.712310,0.657749,0.621251,0.605551,0.609400,0.633992,0.677127,0.738125,0.809102,0.885778,0.959606,1.021730,1.068001,1.090251,1.082005,1.039741,0.961162,0.847738,0.703545,0.529702,0.334411,0.126715,-0.086015,-0.291611,-0.483410,-0.652437,-0.793747,-0.898019,-0.961500,-0.987186,-0.973682,-0.925022,-0.844709,-0.743866,-0.630273,-0.515067,-0.411991,-0.326785,-0.271832,-0.255568,-0.283298,-0.358595,-0.481286,-0.645358,-0.850081,-1.083412,-1.333371,-1.588268,-1.829641,-2.044465,-2.216782,-2.336144,-2.384123,-2.351057,-2.237013,-2.036326,-1.752663,-1.390966,-0.967690,-0.499067,-0.006668,0.481381,0.946955,1.362130,1.700628,1.943587,2.075154,2.087622,1.983687,1.761311,1.437981,1.037585,0.583361,0.113844,-0.344427,-0.757655,-1.102206,-1.344889,-1.468641,-1.477577,-1.367288,-1.150847,-0.842006,-0.476174,-0.080762,0.308629,0.654673,0.939638,1.137432,1.233002,1.224681,1.118357,0.930891,0.691085,0.415098,0.136679,-0.113603,-0.320798,-0.458328,-0.529295,-0.534665,-0.493377,-0.405221,-0.286442,-0.167386,-0.049582,0.058734,0.163545,0.239074,0.263544,0.177669,-0.088068,-0.510898,-0.695372,1.048845};

    DoubleVector x(heating_source_number*time_size);

    //for (size_t ln=0; ln<time_size; ln++) { x[0*time_size+ln] = x1[ln]; x[1*time_size+ln] = x2[ln]; }

    for (size_t i=0; i<heating_source_number; i++) { for (size_t ln=0; ln<time_size; ln++) { x[i*time_size+ln] = 0.1*(i+1)*ln*P3P6_TIME_STEP; } }

#endif

#ifdef P3P6_OPTIMIZE_Y
    const size_t meausere_point_number = f.meausere_point_number;

    //srand(time(nullptr));
    for (size_t i=0; i<heating_source_number; i++)
    {
        for (size_t j=0; j<meausere_point_number; j++)
        {
            f.k[i][j] = 0.0;//0.05 - (rand() % 1000) * 0.0001;
            f.z[i][j] = 10.0 + f.k[i][j];
        }
    }
#endif

#ifdef P3P6_OPTIMIZE_Y
    DoubleVector x(2 * heating_source_number * meausere_point_number);
    for (size_t i=0; i<heating_source_number; i++)
    {
        for (size_t j=0; j<meausere_point_number; j++)
        {
            x[0*heating_source_number*meausere_point_number + i*meausere_point_number + j] = f.k[i][j];
            x[1*heating_source_number*meausere_point_number + i*meausere_point_number + j] = f.z[i][j];
        }
    }

    // T=500
    //const double x1[] = { -0.517938,-0.554214,-0.561210,-0.573928,-0.498175,-0.549009,-0.537552,-0.519455,10.043862,10.003153,10.015863,10.000000,10.031721,9.978533,10.002108,10.013608 };
    //const double x1[] = { -0.495664,-0.572195,-0.519468,-0.586235,-0.504107,-0.530131,-0.570985,-0.501467,10.047378,10.006916,10.019673,10.003897,10.027930,9.974356,9.998018,10.009655 };
    //const double x1[] = { -0.279171,-0.613798,-0.462846,-0.471367,-0.540759,-0.499876,-0.565659,-0.570900,10.044100,10.002614,10.016456,9.999561,10.032629,9.978964,10.003696,10.014012 };
    //const double x1[] = { -0.261460,-0.591819,-0.450201,-0.442891,-0.517557,-0.510802,-0.544676,-0.589113,10.040061,9.994068,10.010016,9.992880,10.039718,9.985547,10.011200,10.021426 };
    //const double x1[] = { -0.218397,-0.570487,-0.347546,-0.405004,-0.464217,-0.514702,-0.539920,-0.608108,10.034023,9.981388,10.000469,9.982918,10.050805,9.995808,10.022939,10.032965 };
    //const double x1[] = { -0.061823,-0.408228,-0.083657,-0.140866,-0.184070,-0.072747,0.056958,-0.070874,10.040246,9.989278,10.017226,10.000269,10.035124,9.977792,10.001686,10.015058 };
    //const double x1[] = { -0.123125,-0.485238,-0.186630,-0.150516,-0.220890,-0.016917,0.051749,-0.017874,10.058294,10.060698,10.032328,10.019139,10.029638,9.975233,10.003558,10.012117 };
    //const double x1[] = { -0.098624,-0.522315,-0.144864,-0.109226,-0.249009,-0.031725,0.110425,0.081891,10.089100,10.169224,10.078014,10.051277,10.072794,9.980609,9.992992,10.017170 };
    //const double x1[] = { -0.142684,-0.477681,-0.148600,-0.086369,-0.229681,0.140437,0.279019,0.173836,10.156707,10.437467,10.168664,10.123407,10.152389,9.971530,9.957007,9.991577 };

    // T=600
    //const double x1[] = { -0.139751,-0.474833,-0.150734,-0.085685,-0.230342,0.140775,0.279316,0.173955,10.149486,10.412895,10.160875,10.119001,10.151724,9.971919,9.957810,9.992079 };
    const double x1[] = { -0.138460,-0.473968,-0.151473,-0.084952,-0.230539,0.141092,0.279258,0.173815,10.138741,10.376056,10.149122,10.112448,10.150277,9.972811,9.959564,9.993170 };

    // T=800
    // const double x1[] = { -0.138864,-0.471275,-0.152512,-0.088256,-0.231484,0.140026,0.279987,0.175223,10.056367,10.094603,10.059245,10.061847,10.218046,9.930965,9.877373,9.942132 };
    // const double x1[] = { -0.134553,-0.470363,-0.153744,-0.083357,-0.230604,0.142347,0.278852,0.174166,10.056316,10.094430,10.059189,10.061815,10.217860,9.931077,9.877598,9.942273 };

    // T=10000
    //const double x1[] = { -0.134009,-0.472945,-0.152794,-0.079919,-0.229945,0.143364,0.277946,0.172524,10.138670,10.375812,10.149044,10.112404,10.150022,9.972967,9.959873,9.99336 };

    for (size_t i=0; i<2*heating_source_number * meausere_point_number; i++) { x[i] = x1[i]; }



#endif

#if defined(P3P6_OPTIMIZE_Q) && defined(P3P6_CALCULATE_GRAD)
    IPrinter::printVector(x.mid(0*time_size, 1*time_size-1));
    IPrinter::printVector(x.mid(1*time_size, 2*time_size-1));
    IPrinter::printSeperatorLine();

    DoubleVector g1(x.length());
    DoubleVector g2(x.length(), 0.0);

    f.gradient(x, g1);
    IPrinter::printVector(g1.mid(0*time_size, 1*time_size-1).L2Normalize());
    IPrinter::printVector(g1.mid(1*time_size, 2*time_size-1).L2Normalize());
    IPrinter::printSeperatorLine();

    //    std::vector<size_t> i1 = {   0,  100,  200,  300,  400,  500,  600,  700,  800,  900, 1000};
    //    std::vector<size_t> i2 = {1001, 1101, 1201, 1301, 1401, 1501, 1601, 1701, 1801, 1901, 2001};
    //    IGradient::Gradient(P3P6_NUM_GRAD_STEP, &f, x, g2, i1); g2[0*time_size] *= 2; g2[1*time_size-1] *= 2;
    //    IPrinter::printVector(g2.mid(0*time_size, 1*time_size-1).L2Normalize());
    //    IGradient::Gradient(P3P6_NUM_GRAD_STEP, &f, x, g2, i2); g2[1*time_size] *= 2; g2[2*time_size-1] *= 2;
    //    IPrinter::printVector(g2.mid(1*time_size, 2*time_size-1).L2Normalize());
    //    IPrinter::printSeperatorLine();


    IGradient::Gradient(&f, P3P6_NUM_GRAD_STEP, x, g2, 0*time_size, 1*time_size-1); g2[0*time_size] *= 2; g2[1*time_size-1] *= 2;
    IPrinter::printVector(g2.mid(0*time_size, 1*time_size-1).L2Normalize());
    IGradient::Gradient(&f, P3P6_NUM_GRAD_STEP, x, g2, 1*time_size, 2*time_size-1); g2[1*time_size] *= 2; g2[2*time_size-1] *= 2;
    IPrinter::printVector(g2.mid(1*time_size, 2*time_size-1).L2Normalize());
    IPrinter::printSeperatorLine();

#endif

#if defined(P3P6_OPTIMIZE_Y) && defined(CALCULATE_GRAD)
    IPrinter::printVector(x.mid(0, 7), "k:\t", 8);
    IPrinter::printVector(x.mid(8, 15), "z:\t", 8);
    //IPrinter::printVector(f.pi_v.mid(   0, 1000));
    //IPrinter::printVector(f.pi_v.mid(1001, 2001));
    IPrinter::printSeperatorLine();

    DoubleVector g1(x.length(), 0.0);
    DoubleVector g2(x.length(), 0.0);

    f.gradient(x, g1);
    IPrinter::printVector(g1.mid(0,  7).L2Normalize(), "k:\t", 8);
    IPrinter::printVector(g1.mid(8, 15).L2Normalize(), "z:\t", 8);
    IPrinter::printSeperatorLine();

    IGradient::Gradient(&f, 0.01, x, g2, 0, 7);
    IPrinter::printVector(g2.mid(0, 7).L2Normalize(), "k:\t", 8);
    IGradient::Gradient(&f, 0.01, x, g2, 8, 15);
    IPrinter::printVector(g2.mid(8, 15).L2Normalize(), "z:\t", 8);
    IPrinter::printSeperatorLine();
#endif


    f.gm = new ConjugateGradient;       f.gm->setNormalize(false);
//    f.gm = new SteepestDescentGradient; f.gm->setNormalize(false);
    f.gm->setFunction(&f);
    f.gm->setGradient(&f);
    f.gm->setPrinter(&f);
    f.gm->setProjection(&f);
    f.gm->setOptimalityTolerance(0.0);
    f.gm->setStepTolerance(0.0);
    f.gm->setFunctionTolerance(0.0);
    f.gm->setR1MinimizeEpsilon(0.01, 0.001);
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

#ifdef P3P6_OPTIMIZE_Y
    measurePoint = new SpacePoint[meausere_point_number];
    measurePoint[0] = SpacePoint(0.20);
    measurePoint[1] = SpacePoint(0.40);
    measurePoint[2] = SpacePoint(0.60);
    measurePoint[3] = SpacePoint(0.80);

    k = DoubleMatrix(heating_source_number, meausere_point_number,  -0.1);
    z = DoubleMatrix(heating_source_number, meausere_point_number, +10.0);
#endif

    const size_t time_size = _timeDimension.size();

    uj_v = DoubleVector(meausere_point_number*time_size, 0.0);
    pi_v = DoubleVector(heating_source_number*time_size, 0.0);

    qi_v = DoubleVector(heating_source_number*time_size, 0.0);
    gj_v = DoubleVector(meausere_point_number*time_size, 0.0);
}

auto Functional::fx(const DoubleVector &x) const -> double
{
    const_cast<Functional*>(this)->functionCount++;

#ifdef P3P6_OPTIMIZE_Q
    const_cast<Functional*>(this)->qi_v = x;
#endif

#ifdef P3P6_OPTIMIZE_Y
    for (size_t i=0; i<heating_source_number; i++)
    {
        for (size_t j=0; j<meausere_point_number; j++)
        {
            const_cast<Functional*>(this)->k[i][j] = x[0*heating_source_number*meausere_point_number + i*meausere_point_number + j];
            const_cast<Functional*>(this)->z[i][j] = x[1*heating_source_number*meausere_point_number + i*meausere_point_number + j];
        }
    }
#endif

    LoadedHeatEquationIBVP::implicit_calculate_D1V1();
    double integralU = integral(U);
    return integralU;
}

auto Functional::integral(const DoubleMatrix &) const -> double
{
    const size_t N = static_cast<size_t> ( _spaceDimensionX.size()-1 );
    const double hx = _spaceDimensionX.step();

    double udiff = 0.0;
    double usum = 0.0;

    udiff = (U[0]-V[0]); usum += 0.5 * udiff * udiff;// * mu(0, 0);
    udiff = (U[N]-V[N]); usum += 0.5 * udiff * udiff;// * mu(N, 0);

    for (size_t n=1; n<=N-1; n++) { udiff = U[n]-V[n]; usum += udiff * udiff; }

    return usum*hx;
}

auto Functional::gradient(const DoubleVector &x, DoubleVector &g) const -> void
{
    const size_t time_size = _timeDimension.size();

    g.resize(x.length());
    for (size_t n=0; n<g.length(); n++) g[n] = 0.0;

#ifdef P3P6_OPTIMIZE_Q
    const_cast<Functional*>(this)->qi_v = x;
#endif

#ifdef P3P6_OPTIMIZE_Y
    for (size_t i=0; i<heating_source_number; i++)
    {
        for (size_t j=0; j<meausere_point_number; j++)
        {
            const_cast<Functional*>(this)->k[i][j] = x[0*heating_source_number*meausere_point_number + i*meausere_point_number + j];
            const_cast<Functional*>(this)->z[i][j] = x[1*heating_source_number*meausere_point_number + i*meausere_point_number + j];
        }
    }
#endif

    LoadedHeatEquationIBVP::implicit_calculate_D1V1();
    LoadedHeatEquationFBVP::implicit_calculate_D1V1();

#ifdef P3P6_OPTIMIZE_Q
    for (size_t i=0; i<heating_source_number; i++)
    {
        for (size_t ln=0; ln<time_size; ln++)
        {
            /*if (ln%100==0)*/ g[i*time_size+ln] = -pi_v[i*time_size+ln];
        }
    }
#endif

#ifdef P3P6_OPTIMIZE_Y
    const double time_step = _timeDimension.step();
    for (size_t i=0; i<heating_source_number; i++)
    {
        const size_t offset1 = i*time_size;
        for (size_t j=0; j<meausere_point_number; j++)
        {
            const size_t offset2 = j*time_size;

            for (size_t ln=0; ln<time_size; ln++)
            {
                double w = (ln==0 || ln == time_size-1) ? 0.5*time_step : time_step;
                if (optimizeK) g[0*heating_source_number*meausere_point_number + i*meausere_point_number + j] += -pi_v[offset1+ln] * (uj_v[offset2+ln]-z[i][j]) * w;
                if (optimizeZ) g[1*heating_source_number*meausere_point_number + i*meausere_point_number + j] += +pi_v[offset1+ln] * k[i][j] * w;
            }
        }
    }
#endif
}

auto Functional::project(DoubleVector &x, size_t index) -> void
{
#ifdef P3P6_OPTIMIZE_Q
    if (x[index] <=  0.0) x[index] =  0.0;
    //if (x[index] >= 12.0) x[index] = 12.0;
#endif
}

auto Functional::print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientBasedMethod::MethodResult /*result*/) const -> void
{
    //if (f <= 0.1) { const_cast<Functional*>(this)->gm->setR1MinimizeEpsilon(0.100, 0.0100); }
    //if (fabs(f - lastFx) <= 0.10000) { const_cast<Functional*>(this)->gm->setR1MinimizeEpsilon(0.100, 0.0100); }
    //if (fabs(f - lastFx) <= 0.00100) { const_cast<Functional*>(this)->gm->setR1MinimizeEpsilon(0.010, 0.0010); }
    //if (fabs(f - lastFx) <= 0.00001) { const_cast<Functional*>(this)->gm->setR1MinimizeEpsilon(0.001, 0.0001); }

    const_cast<Functional*>(this)->lastFx = f;

    const size_t time_size = _timeDimension.size();
    printf("I[%3d]: %10.6f | %d | %10.6f %10.6f | %10.6f %10.6f | %10.6f | %4zu | %4zu | %4zu | %4zu | %10.6f %10.6f\n", iteration,f,  time_size, U.min(), U.max(),
           gm->min_step, gm->min_epsilon, alpha, functionCount,
           gm->search_function_count, gm->golden_function_count, gm->total__function_count, x.min(), x.max());

#ifdef P3P6_OPTIMIZE_Q
    //    IPrinter::printVector(x.mid(0*time_size, 1*time_size-1));
    //    IPrinter::printVector(x.mid(1*time_size, 2*time_size-1));
    //    IPrinter::printVector(g.mid(0*time_size, 1*time_size-1));
    //    IPrinter::printVector(g.mid(1*time_size, 2*time_size-1));
#endif
#ifdef OPTIMIZE_Y
    IPrinter::printVector(x, "x:\t", x.length());
    IPrinter::printVector(g, "g:\t", g.length());
#endif
    const_cast<Functional*>(this)->functionCount = 0;
    //    const_cast<Functional*>(this)->optimizeK = iteration%2==0;
    //    const_cast<Functional*>(this)->optimizeZ = iteration%2==1;

    //if (iteration!=0) const_cast<Functional*>(this)->gm->setR1MinimizeEpsilon(alpha*0.5, alpha*0.1);

    if (f<=0.0001)
    {
        IPrinter::printSeperatorLine();
#ifdef P3P6_OPTIMIZE_Q
        IPrinter::printVector(x.mid(0*time_size, 1*time_size-1), nullptr, x.mid(0*time_size, 1*time_size-1).length());
        IPrinter::printVector(x.mid(1*time_size, 2*time_size-1), nullptr, x.mid(1*time_size, 2*time_size-1).length());
#endif
#ifdef P3P6_OPTIMIZE_Y
        //        const size_t offset = heating_source_number*meausere_point_number;
        //        const size_t time_size1 = static_cast<size_t>(_timeDimension.size());
        //        IPrinter::printVector(qi_v.mid(0*time_size1, 1*time_size1-1), nullptr, time_size1);
        //        IPrinter::printVector(qi_v.mid(1*time_size1, 2*time_size1-1), nullptr, time_size1);
#endif
        IPrinter::printSeperatorLine();
    }


    const_cast<Functional*>(this)->drawImages = 1;
#ifdef P3P6_OPTIMIZE_Q
    const_cast<Functional*>(this)->qi_v = x;
#endif
#ifdef P3P6_OPTIMIZE_Y
    for (size_t i=0; i<heating_source_number; i++)
    {
        for (size_t j=0; j<meausere_point_number; j++)
        {
            const_cast<Functional*>(this)->k[i][j] = x[0*heating_source_number*meausere_point_number + i*meausere_point_number + j];
            const_cast<Functional*>(this)->z[i][j] = x[1*heating_source_number*meausere_point_number + i*meausere_point_number + j];
        }
    }
#endif
    //    LoadedHeatEquationIBVP::implicit_calculate_D1V1();
    //    const_cast<Functional*>(this)->drawImages = 0;
    //    std::string path = "E:\\project\\hvugar\\num_methods\\trunk\\optimal\\bin\\data\\problem3P\\f\\png\\";
    //    system(std::string("mkdir ").append(path).append(std::to_string(iteration)).append(" >nul").data());
    //    system(std::string("move  ").append(path).append("*.png ").append(path).append(std::to_string(iteration)).append("\\ >nul").data());
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
    bc = BoundaryConditionPDE::Robin(lambda1, +1.0); return lambda1 * envrmnt_temperature;
}

auto LoadedHeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    const double dimX_step = spaceDimensionX().step();

    double sum = 0.0;

#ifdef P3P6_OPTIMIZE_Q
    for (size_t i=0; i<heating_source_number; i++)
    {
        const double si = s(tn, i);
        const double w = DeltaFunction::gaussian(sn.x, si, dimX_step);
        if (w > 0.0) { sum += q(tn, i) * w; }
    }
#endif

#ifdef P3P6_OPTIMIZE_Y
    size_t time_size = timeDimension().size();
    for (size_t i=0; i<heating_source_number; i++)
    {
        const SpacePoint &si = s(tn, i);
        const double w = DeltaFunction::gaussian(sn, si, SpacePoint(dimX_step, dimY_step));
        if (w > 0.0) { sum += qi_v[i*time_size + tn.i] * w; }
    }
#endif

    return sum - thermalConvection() * envrmnt_temperature;
}

auto LoadedHeatEquationIBVP::layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const -> void
{
    LoadedHeatEquationIBVP* lhe = const_cast<LoadedHeatEquationIBVP*>(this);

#ifdef P3P6_OPTIMIZE_Y
    size_t time__min = static_cast<size_t>(timeDimension().min());
    size_t time__max = static_cast<size_t>(timeDimension().max());
    size_t time_size = timeDimension().size();
    for (size_t j=0; j<meausere_point_number; j++)
    {
        const SpacePoint &mp = measurePoint[j];
        lhe->Shared::uj_v[j*time_size+tn.i] = DeltaFunction::lumpedPointG(u, mp, _spaceDimensionX, _spaceDimensionY, 1, 4);
    }

    for (size_t i=0; i<heating_source_number; i++)
    {
        double qi = 0.0;
        size_t index = i*time_size + tn.i;
        for (size_t j=0; j<meausere_point_number; j++) { qi += k[i][j]*(uj_v[j*time_size+tn.i] - z[i][j]); }
        if (tn.i==time__min) lhe->qi_v[index] = qi;
        if (tn.i!=time__max) lhe->qi_v[index+1] = qi;
        if (tn.i==time__max) lhe->qi_v[index] = qi;
    }
#endif

    //saveImage(u, tn);
    if (static_cast<int>(tn.i) == timeDimension().max()) { lhe->Shared::U = u; }
}

auto LoadedHeatEquationIBVP::saveImage(const DoubleVector &u, const TimeNodePDE &tn) const -> void
{
    if (drawImages != 0)
    {
        //printf(">>> %6d %14.6f %14.6f\n", tn.i, u.min(), u.max());

        QString filename = QString("data/problem3P/f/png/%1.png").arg(tn.i, 8, 10, QChar('0'));
        QPixmap pixmap;
        visualizeMatrixHeat(u, u.min(), u.max(), pixmap, 101, 101);
        pixmap.save(filename);
    }
}

auto LoadedHeatEquationIBVP::timeDimension() const -> Dimension { return Shared::_timeDimension; }

auto LoadedHeatEquationIBVP::spaceDimensionX() const -> Dimension { return Shared::_spaceDimensionX; }

auto LoadedHeatEquationIBVP::spaceDimensionY() const -> Dimension { return spaceDimensionX(); }

auto LoadedHeatEquationIBVP::spaceDimensionZ() const -> Dimension { return spaceDimensionX(); }

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
    size_t n = static_cast<size_t>(sn.i);
    return -2.0 * ( U[n] - V[n] );
}

auto LoadedHeatEquationFBVP::boundary(const SpaceNodePDE &/*sn*/, const TimeNodePDE &/*tn*/, BoundaryConditionPDE &bc) const -> double
{
    bc = BoundaryConditionPDE::Robin(lambda1, +1.0); return 0.0;
}

auto LoadedHeatEquationFBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    double sum = 0.0;
#ifdef P3P6_OPTIMIZE_Y
    const double dimX_step = spaceDimensionX().step();
    size_t time_size = timeDimension().size();
    for (size_t j=0; j<meausere_point_number; j++)
    {
        const SpacePoint &mp = measurePoint[j];
        const double w = DeltaFunction::gaussian(sn, mp, SpacePoint(dimX_step, dimY_step));
        if (w > 0.0) { sum -= gj_v[j*time_size+tn.i] * w; }
    }
#endif

    return sum;
}

auto LoadedHeatEquationFBVP::layerInfo(const DoubleVector &p, const TimeNodePDE &tn) const -> void
{
    LoadedHeatEquationFBVP* lhe = const_cast<LoadedHeatEquationFBVP*>(this);

#ifdef P3P6_OPTIMIZE_Q
    size_t time_size = timeDimension().size();
    for (size_t i=0; i<heating_source_number; i++)
    {
        const double si = s(tn, i);
        const double w = DeltaFunction::lumpedPointG(p, si, _spaceDimensionX, 1, 4);
        lhe->Shared::pi_v[i*time_size+tn.i] = w;
    }
#endif

#ifdef P3P6_OPTIMIZE_Y
    size_t time__min = static_cast<size_t>(timeDimension().min());
    size_t time__max = static_cast<size_t>(timeDimension().max());
    size_t time_size = timeDimension().size();
    for (size_t i=0; i<heating_source_number; i++)
    {
        const SpacePoint &si = s(tn, i);
        lhe->Shared::pi_v[i*time_size+tn.i] = DeltaFunction::lumpedPointG(p, si, _spaceDimensionX, _spaceDimensionY, 1, 4);
    }

    for (size_t j=0; j<meausere_point_number; j++)
    {
        double gj = 0.0;
        size_t index = j*time_size + tn.i;
        for (size_t i=0; i<heating_source_number; i++) { gj += k[i][j]*pi_v[i*time_size+tn.i]; }
        if (tn.i==time__max) lhe->gj_v[index] = gj;
        if (tn.i!=time__min) lhe->gj_v[index-1] = gj;
        if (tn.i==time__min) lhe->gj_v[index] = gj;
    }
#endif

    //saveImage(p, tn);
}

auto LoadedHeatEquationFBVP::saveImage(const DoubleVector &p, const TimeNodePDE &tn) const -> void
{
    if (drawImages != 0)
    {
        QString filename = QString("data/problem3P/b/png/%1.png").arg(tn.i, 8, 10, QChar('0'));
        QPixmap pixmap;
        visualizeMatrixHeat(p, p.min(), p.max(), pixmap, 101, 101);
        pixmap.save(filename);
    }

}

auto LoadedHeatEquationFBVP::timeDimension() const -> Dimension { return Shared::_timeDimension; }

auto LoadedHeatEquationFBVP::spaceDimensionX() const -> Dimension { return Shared::_spaceDimensionX; }

auto LoadedHeatEquationFBVP::spaceDimensionY() const -> Dimension { return spaceDimensionX(); }

auto LoadedHeatEquationFBVP::spaceDimensionZ() const -> Dimension { return spaceDimensionX(); }

/*********************************************************************************************************************************************/

auto Shared::q(const TimeNodePDE &tn, size_t i) const -> double
{
    size_t time_size = _timeDimension.size();
    return qi_v[i*time_size + tn.i];
}

auto Shared::s(const TimeNodePDE &tn, size_t i) const -> double
{
    double x;
    const double t = tn.t;
    const double v = 2.0*M_PI/**t*/;
    switch (i) {
    case 0: {
        x = 0.2*cos(v*t) + 0.5;
    } break;
    case 1: {
        x = -0.2*sin(v*t) + 0.5;
    } break;
    }

    return x;
}
