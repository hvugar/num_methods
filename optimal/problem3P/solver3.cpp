#include "solver3.h"

#define SYNTEZ
//#define OMTIMIZE_V

using namespace p3p3;

void Functional::Main(int /*argc*/, char** /*argv*/)
{
    Functional *fn = new Functional();

    const size_t time_size = fn->timeDimension().size();
    const size_t dimX_size = fn->spaceDimensionX().size();
    const size_t dimY_size = fn->spaceDimensionY().size();
    const double time_step = fn->timeDimension().step();
    //const double dimX_step = fn->spaceDimensionX().step();
    //const double dimY_step = fn->spaceDimensionX().step();

    fn->V.resize(dimY_size, dimX_size, 5.0);
    fn->U.resize(dimY_size, dimX_size, 0.0);

#ifdef SYNTEZ
    fn->k.resize(fn->heating_point_number, fn->measure_point_number, 0.0);
    fn->z.resize(fn->heating_point_number, fn->measure_point_number, 0.0);
    fn->measure_point.push_back(SpacePoint(0.25, 0.25));
    fn->measure_point.push_back(SpacePoint(0.25, 0.75));
    fn->measure_point.push_back(SpacePoint(0.75, 0.75));
    fn->measure_point.push_back(SpacePoint(0.75, 0.25));

    DoubleVector x(4*time_size + 2*fn->heating_point_number*fn->measure_point_number + 2*fn->measure_point_number);
#else
#ifdef OMTIMIZE_V
    DoubleVector x(4*time_size);
#else
    DoubleVector x(6*time_size);
#endif
#endif

    unsigned int w = 10, p = 4;

    for (unsigned int ln=0; ln<time_size; ln++)
    {
        const double t = ln*time_step;
        const TimeNodePDE tn(ln, t);

#ifdef SYNTEZ
        x[0*time_size + ln] = +fn->R[0]*sin(2.0*t) + 0.50;/*+fn->R[0]*sin(fn->v(tn, 0)) + 0.50;*/
        x[1*time_size + ln] = +fn->R[0]*cos(2.0*t) + 0.50;/*+fn->R[0]*cos(fn->v(tn, 0)) + 0.50;*/
        x[2*time_size + ln] = -fn->R[1]*sin(2.0*t) + 0.50;/*-fn->R[1]*sin(fn->v(tn, 1)) + 0.50;*/
        x[3*time_size + ln] = +fn->R[1]*cos(2.0*t) + 0.50;/*+fn->R[1]*cos(fn->v(tn, 1)) + 0.50;*/

#else
        x[0*time_size + ln] = 50.00;
        x[1*time_size + ln] = 50.00;
#ifdef OMTIMIZE_V
        x[2*time_size + ln] = 2.0*t*t;
        x[3*time_size + ln] = 2.0*t*t;
#else
        x[2*time_size + ln] = +fn->R[0]*sin(fn->v(tn, 0)) + 0.50; /*+0.8*t+0.1;*/
        x[3*time_size + ln] = +fn->R[0]*cos(fn->v(tn, 0)) + 0.50; /*+0.8*t+0.1;*/
        x[4*time_size + ln] = -fn->R[1]*sin(fn->v(tn, 1)) + 0.50; /*+0.8*t+0.1;*/
        x[5*time_size + ln] = +fn->R[1]*cos(fn->v(tn, 1)) + 0.50; /*-0.8*t+0.9;*/
#endif
#endif
    }

#ifdef SYNTEZ
    x[4*time_size +  0] = -20.4*sin(rand());
    x[4*time_size +  1] = -20.4*sin(rand());
    x[4*time_size +  2] = -20.4*sin(rand());
    x[4*time_size +  3] = -20.4*sin(rand());
    x[4*time_size +  4] = -20.4*sin(rand());
    x[4*time_size +  5] = -20.4*sin(rand());
    x[4*time_size +  6] = -20.4*sin(rand());
    x[4*time_size +  7] = -20.4*sin(rand());

    x[4*time_size +  8] = +4.5*sin(rand());
    x[4*time_size +  9] = +4.5*sin(rand());
    x[4*time_size + 10] = +4.5*sin(rand());
    x[4*time_size + 11] = +4.5*sin(rand());
    x[4*time_size + 12] = +4.5*sin(rand());
    x[4*time_size + 13] = +4.5*sin(rand());
    x[4*time_size + 14] = +4.5*sin(rand());
    x[4*time_size + 15] = +4.5*sin(rand());

    x[4*time_size + 16] = 0.25;
    x[4*time_size + 17] = 0.25;
    x[4*time_size + 18] = 0.25;
    x[4*time_size + 19] = 0.75;
    x[4*time_size + 20] = 0.75;
    x[4*time_size + 21] = 0.75;
    x[4*time_size + 22] = 0.75;
    x[4*time_size + 23] = 0.25;

    IPrinter::printVector(w, p, x.mid(0,   100));
    IPrinter::printVector(w, p, x.mid(101, 201));
    IPrinter::printVector(w, p, x.mid(202, 302));
    IPrinter::printVector(w, p, x.mid(303, 403));
    IPrinter::print(x.mid(404, 411), 8, w, p);
    IPrinter::print(x.mid(412, 419), 8, w, p);
    IPrinter::print(x.mid(420, 427), 8, w, p);
    IPrinter::printSeperatorLine();
#endif

    {
        DoubleVector g0;
        fn->gradient(x, g0);

#ifdef SYNTEZ

        g0[100] = g0[201] = g0[302] = g0[403] = 0.0;
        IPrinter::printVector(w, p, g0.mid(0,   100).L2Normalize());
        IPrinter::printVector(w, p, g0.mid(101, 201).L2Normalize());
        IPrinter::printVector(w, p, g0.mid(202, 302).L2Normalize());
        IPrinter::printVector(w, p, g0.mid(303, 403).L2Normalize());
        IPrinter::print(g0.mid(404, 411).L2Normalize(), 8, w, p);
        IPrinter::print(g0.mid(412, 419).L2Normalize(), 8, w, p);
        IPrinter::print(g0.mid(420, 427).L2Normalize(), 8, w, p);
        IPrinter::printSeperatorLine();

#else
        IPrinter::printVector(w, p, g0.mid(0,   100).L2Normalize());
        IPrinter::printVector(w, p, g0.mid(101, 201).L2Normalize());
        IPrinter::printVector(w, p, g0.mid(202, 302).L2Normalize());
        IPrinter::printVector(w, p, g0.mid(303, 403).L2Normalize());
#ifndef OMTIMIZE_V
        IPrinter::printVector(w, p, g0.mid(404, 504).L2Normalize());
        IPrinter::printVector(w, p, g0.mid(505, 605).L2Normalize());
#endif
        IPrinter::printSeperatorLine();

        //#ifndef OMTIMZIE_V
        //        IPrinter::printVector(w, p, fn->pp[0]/*.L2Normalize()*/);
        //        IPrinter::printVector(w, p, fn->px[0]/*.L2Normalize()*/);
        //        IPrinter::printVector(w, p, fn->py[0]/*.L2Normalize()*/);

        //        IPrinter::printVector(w, p, fn->pp[1]/*.L2Normalize()*/);
        //        IPrinter::printVector(w, p, fn->px[1]/*.L2Normalize()*/);
        //        IPrinter::printVector(w, p, fn->py[1]/*.L2Normalize()*/);
        //#endif


        //        IPrinter::printVector(10, 6, g0.mid(0,   100));
        //        IPrinter::printVector(10, 6, g0.mid(101, 201));
        //        IPrinter::printVector(10, 6, g0.mid(202, 302));
        //        IPrinter::printVector(10, 6, g0.mid(303, 403));
        //        IPrinter::printSeperatorLine();
#endif
    }

    {
        DoubleVector g1(x.length());

#ifdef SYNTEZ
        const double step = 0.01;
        //IGradient::Gradient(fn, step, x, g1,   0, 100); /*g1[  0] *= 2.0; g1[100] *= 2.0;*/ IPrinter::printVector(w, p, g1.mid(  0, 100).L2Normalize());
        //IGradient::Gradient(fn, step, x, g1, 101, 201); /*g1[101] *= 2.0; g1[201] *= 2.0;*/ IPrinter::printVector(w, p, g1.mid(101, 201).L2Normalize());
        //IGradient::Gradient(fn, step, x, g1, 202, 302); /*g1[202] *= 2.0; g1[302] *= 2.0;*/ IPrinter::printVector(w, p, g1.mid(202, 302).L2Normalize());
        //IGradient::Gradient(fn, step, x, g1, 303, 403); /*g1[303] *= 2.0; g1[403] *= 2.0;*/ IPrinter::printVector(w, p, g1.mid(303, 403).L2Normalize());
        IGradient::Gradient(fn, step, x, g1, 404, 411); /*g1[303] *= 2.0; g1[403] *= 2.0;*/ IPrinter::print(g1.mid(404, 411).L2Normalize(), 8, w, p);
        IGradient::Gradient(fn, step, x, g1, 412, 419); /*g1[303] *= 2.0; g1[403] *= 2.0;*/ IPrinter::print(g1.mid(412, 419).L2Normalize(), 8, w, p);
        IGradient::Gradient(fn, step, x, g1, 420, 427); /*g1[303] *= 2.0; g1[403] *= 2.0;*/ IPrinter::print(g1.mid(420, 427).L2Normalize(), 8, w, p);
#else
        const double step = 0.01;
        IGradient::Gradient(fn, step, x, g1,   0, 100); /*g1[  0] *= 2.0; g1[100] *= 2.0;*/ IPrinter::printVector(w, p, g1.mid(  0, 100).L2Normalize());
        IGradient::Gradient(fn, step, x, g1, 101, 201); /*g1[101] *= 2.0; g1[201] *= 2.0;*/ IPrinter::printVector(w, p, g1.mid(101, 201).L2Normalize());
        IGradient::Gradient(fn, step, x, g1, 202, 302); /*g1[202] *= 2.0; g1[302] *= 2.0;*/ IPrinter::printVector(w, p, g1.mid(202, 302).L2Normalize());
        IGradient::Gradient(fn, step, x, g1, 303, 403); /*g1[303] *= 2.0; g1[403] *= 2.0;*/ IPrinter::printVector(w, p, g1.mid(303, 403).L2Normalize());
#ifndef OMTIMIZE_V
        IGradient::Gradient(fn, step, x, g1, 404, 504); /*g1[404] *= 2.0; g1[504] *= 2.0;*/ IPrinter::printVector(w, p, g1.mid(404, 504).L2Normalize());
        IGradient::Gradient(fn, step, x, g1, 505, 605); /*g1[505] *= 2.0; g1[605] *= 2.0;*/ IPrinter::printVector(w, p, g1.mid(505, 605).L2Normalize());
#endif
        //        IPrinter::printSeperatorLine();

        //        IPrinter::printVector(10, 6, g1.mid(0, 100));
        //        IPrinter::printVector(10, 6, g1.mid(101, 201));
        //        IPrinter::printVector(10, 6, g1.mid(202, 302));
        //        IPrinter::printVector(10, 6, g1.mid(303, 403));
#endif
    }
}

Functional::Functional()
{
    ih = new LoadedHeatEquationIBVP(this);
    fh = new LoadedHeatEquationFBVP(this);

    const double a = 0.0001;
    const double lambda0 = 0.0;

    ih->setThermalDiffusivity(a);
    ih->setThermalConductivity(0.0);
    ih->setThermalConvection(-lambda0);

    fh->setThermalDiffusivity(-a);
    fh->setThermalConductivity(0.0);
    fh->setThermalConvection(lambda0);
}

SpacePoint Functional::tr(const TimeNodePDE &tn, size_t i) const
{
    const double t = tn.t;
    const size_t time_size = timeDimension().size();
    SpacePoint sp;

#ifdef SYNTEZ
    if (i==0) { sp.x = s[0*time_size + tn.i]; sp.y = s[1*time_size + tn.i]; }
    if (i==1) { sp.x = s[2*time_size + tn.i]; sp.y = s[3*time_size + tn.i]; }
#else
#ifdef OMTIMIZE_V
    switch (i)
    {

    case 0:
    {
        sp.x = +R[0]*sin(v(tn, i)) + 0.50;
        sp.y = +R[0]*cos(v(tn, i)) + 0.50;
    } break;

    case 1:
    {
        sp.x = -R[1]*sin(v(tn,i)) + 0.50;
        sp.y = +R[1]*cos(v(tn,i)) + 0.50;
    } break;

    default:
    {
        //throw std::exception("i is not valid...");
    }

    }
#else
    if (i==0) { sp.x = x[2*time_size + tn.i]; sp.y = x[3*time_size + tn.i]; }
    if (i==1) { sp.x = x[4*time_size + tn.i]; sp.y = x[5*time_size + tn.i]; }
#endif
#endif

    return sp;
}

double Functional::v(const TimeNodePDE &tn, size_t i) const
{
#ifdef SYNTEZ
    return 2.0*tn.t;
#else
#ifndef OPTIMIZE_V
    return 2.0*tn.t;
#else
    const size_t time_size = timeDimension().size();
    return x[(2+i)*time_size + tn.i];
#endif
#endif
}

double Functional::q(const TimeNodePDE &tn, size_t i) const
{
    const size_t time_size = timeDimension().size();
    return s[i*time_size + tn.i];
}

double Functional::fx(const DoubleVector &x) const
{
    setVector(x);

    ih->implicit_calculate_D2V2();

    return integral(U);
}

auto Functional::integral(const DoubleMatrix &U) const -> double
{
    const Dimension &dimensionX = spaceDimensionX();
    const Dimension &dimensionY = spaceDimensionY();
    const size_t N = static_cast<size_t> ( dimensionX.size()-1 );
    const size_t M = static_cast<size_t> ( dimensionY.size()-1 );
    const double hx = dimensionX.step();
    const double hy = dimensionY.step();

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

void Functional::gradient(const DoubleVector &x, DoubleVector &g) const
{
    setVector(x);

    const double time_step = timeDimension().step();
    const size_t time_size = timeDimension().size();

    ih->implicit_calculate_D2V2();
    fh->implicit_calculate_D2V2();

    g.clear();
    g.resize(x.length());

#ifdef SYNTEZ

    g[404] = g[405] = g[406] = g[407] = g[408] = g[409] = g[410] = g[411] = 0.0;
    g[412] = g[413] = g[414] = g[415] = g[416] = g[417] = g[418] = g[419] = 0.0;

    //g[420 + 2*j + 0] = 0.0;
    //g[420 + 2*j + 1] = 0.0;

    for (unsigned int ln=0; ln<time_size; ln++)
    {
        for (size_t i=0; i<heating_point_number; i++)
        {
            double qi = 0.0;
            for (size_t j=0; j<measure_point_number; j++)
            {
                qi += k.at(i,j) * (uu[j][ln]-z.at(i,j));
            }

            if (i==0)
            {
                g[0*time_size+ln] = -qi * px[0][ln];
                g[1*time_size+ln] = -qi * py[0][ln];
            }

            if (i==1)
            {
                g[2*time_size+ln] = -qi * px[1][ln];
                g[3*time_size+ln] = -qi * py[1][ln];
            }
        }

        for (size_t i=0; i<heating_point_number; i++)
        {
            for (size_t j=0; j<measure_point_number; j++)
            {
                g[404 + i*measure_point_number + j] += -pp[i][ln] * (uu[j][ln]-z.at(i,j)) * time_step;
                g[412 + i*measure_point_number + j] -= -pp[i][ln] * k.at(i,j) * time_step;
            }
        }

        for (size_t j=0; j<measure_point_number; j++)
        {
            for (size_t i=0; i<heating_point_number; i++)
            {
                g[420 + 2*j + 0] += -pp[i][ln] * ux[j][ln] * k.at(i,j) * time_step;
                g[420 + 2*j + 1] += -pp[i][ln] * uy[j][ln] * k.at(i,j) * time_step;
            }
        }

    }


#else

    for (unsigned int ln=0; ln<time_size; ln++)
    {
        const double t = ln*time_step;
        const TimeNodePDE tn = TimeNodePDE(ln, t);

        for (size_t i=0; i<heating_point_number; i++)
        {
            double qi = q(tn, i);
            g[i*time_size+ln] = -pp[i][ln];


#ifdef OMTIMIZE_V
            if (i==0) g[(2+i)*time_size+ln] = -qi * R[0] * ( +px[0][ln]*cos(v(tn, 0)) - py[0][ln]*sin(v(tn, 0)) );
            if (i==1) g[(2+i)*time_size+ln] = -qi * R[1] * ( -px[1][ln]*cos(v(tn, 1)) - py[1][ln]*sin(v(tn, 1)) );
#else
            if (i==0)
            {
                g[2*time_size+ln] = -qi * px[0][ln];
                g[3*time_size+ln] = -qi * py[0][ln];
            }

            if (i==1)
            {
                g[4*time_size+ln] = -qi * px[1][ln];
                g[5*time_size+ln] = -qi * py[1][ln];
            }
#endif
        }
    }

#endif
}

void Functional::setVector(const DoubleVector &x) const
{
    const size_t time_size = timeDimension().size();

    const_cast<Functional*>(this)->s = x.mid(0, 4*time_size-1);

    const_cast<Functional*>(this)->k.at(0,0) = x[4*time_size+0];
    const_cast<Functional*>(this)->k.at(0,1) = x[4*time_size+1];
    const_cast<Functional*>(this)->k.at(0,2) = x[4*time_size+2];
    const_cast<Functional*>(this)->k.at(0,3) = x[4*time_size+3];

    const_cast<Functional*>(this)->k.at(1,0) = x[4*time_size+4];
    const_cast<Functional*>(this)->k.at(1,1) = x[4*time_size+5];
    const_cast<Functional*>(this)->k.at(1,2) = x[4*time_size+6];
    const_cast<Functional*>(this)->k.at(1,3) = x[4*time_size+7];

    const_cast<Functional*>(this)->z.at(0,0) = x[4*time_size+8];
    const_cast<Functional*>(this)->z.at(0,1) = x[4*time_size+9];
    const_cast<Functional*>(this)->z.at(0,2) = x[4*time_size+10];
    const_cast<Functional*>(this)->z.at(0,3) = x[4*time_size+11];

    const_cast<Functional*>(this)->z.at(1,0) = x[4*time_size+12];
    const_cast<Functional*>(this)->z.at(1,1) = x[4*time_size+13];
    const_cast<Functional*>(this)->z.at(1,2) = x[4*time_size+14];
    const_cast<Functional*>(this)->z.at(1,3) = x[4*time_size+15];

    const_cast<Functional*>(this)->measure_point[0].x = x[4*time_size+16];
    const_cast<Functional*>(this)->measure_point[0].y = x[4*time_size+17];
    const_cast<Functional*>(this)->measure_point[1].x = x[4*time_size+18];
    const_cast<Functional*>(this)->measure_point[1].y = x[4*time_size+19];
    const_cast<Functional*>(this)->measure_point[2].x = x[4*time_size+20];
    const_cast<Functional*>(this)->measure_point[2].y = x[4*time_size+21];
    const_cast<Functional*>(this)->measure_point[3].x = x[4*time_size+22];
    const_cast<Functional*>(this)->measure_point[3].y = x[4*time_size+23];
}

void Functional::frw_saveToImage(const DoubleMatrix &u UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const
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

/********************************************************************************************************************************************************/

LoadedHeatEquationIBVP::LoadedHeatEquationIBVP(Functional *function) : _functional(function) {}

LoadedHeatEquationIBVP::LoadedHeatEquationIBVP(const LoadedHeatEquationIBVP &other) { *this = other; }

LoadedHeatEquationIBVP& LoadedHeatEquationIBVP::operator =(const LoadedHeatEquationIBVP &) { return *this; }

LoadedHeatEquationIBVP::~LoadedHeatEquationIBVP() {}

auto LoadedHeatEquationIBVP::initial(const SpaceNodePDE &/*sn*/, InitialCondition /*ic*/) const -> double { return _initial_temperature; }

auto LoadedHeatEquationIBVP::boundary(const SpaceNodePDE &/*sn*/, const TimeNodePDE &/*tn*/, BoundaryConditionPDE &bc) const -> double
{
    const double lambda1 = _functional->_lambda1;
    //bc = BoundaryConditionPDE::Dirichlet(1.0, 1.0); return 0.0;
    //bc = BoundaryConditionPDE::Neumann(1.0, 0.0); return 0.0;
    bc = BoundaryConditionPDE::Robin(lambda1, +1.0, lambda1); return _enviroment_temperature;
}

auto LoadedHeatEquationIBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    const size_t ln = static_cast<size_t>(tn.i);
    const double dimX_step = spaceDimensionX().step();
    const double dimY_step = spaceDimensionY().step();
    const size_t heating_point_number = _functional->heating_point_number;

    double fx = -thermalConvection() * _enviroment_temperature;

#ifdef SYNTEZ
    const size_t measure_point_number = _functional->measure_point_number;
    const DoubleMatrix &k = _functional->k;
    const DoubleMatrix &z = _functional->z;
    const DoubleVector *uu = _functional->uu;


    double sum = 0.0;
    double qi[] = { 0.0, 0.0 };
    //printf("t:%3zu x:%3zu y:%3zu ", tn.i, sn.i, sn.j);
    for (size_t i=0; i<heating_point_number; i++)
    {
        const SpacePoint &zi = _functional->tr(tn, i);
        const double w = DeltaFunction::gaussian(sn, zi, SpacePoint(dimX_step, dimY_step));
        if (w>=0.0)
        {
            qi[i] = 0.0;
            for (size_t j=0; j<measure_point_number; j++)
            {
                qi[i] += k[i][j] * ( uu[j][ln] - z[i][j] ) * w;
            }
            sum += qi[i];
        }
        //if (i==0) printf("%.2f %.2f %.10f ", zi.x, zi.y, qi[i]);
    }
    //puts("");
#else
    //const unsigned dimX_size = spaceDimensionX().size();
    //const unsigned dimY_size = spaceDimensionY().size();

    double sum = 0.0;
    for (size_t i=0; i<heating_point_number; i++)
    {
        const SpacePoint &zi = _functional->tr(tn, i);
        const double qi = _functional->q(tn, i);
        sum += qi * DeltaFunction::gaussian(sn, zi, SpacePoint(dimX_step, dimY_step));
        //sum += qi * DeltaFunction::nearest(sn, zi, dimX_step, dimY_step, dimX_size, dimY_size);
    }
#endif
    return fx + sum;
}

auto LoadedHeatEquationIBVP::layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const -> void
{
    if (static_cast<int>(tn.i) == timeDimension().max()) { _functional->U = u; }

    /****************************************************************************************************************************/

#ifdef SYNTEZ
    const size_t measure_point_number = _functional->measure_point_number;

    if (_functional->uu == nullptr)
    {
        const size_t time_size = timeDimension().size();

        _functional->uu = new DoubleVector[measure_point_number];
        _functional->ux = new DoubleVector[measure_point_number];
        _functional->uy = new DoubleVector[measure_point_number];

        for (size_t j=0; j<measure_point_number; j++)
        {
            _functional->uu[j].resize(time_size);
            _functional->ux[j].resize(time_size);
            _functional->uy[j].resize(time_size);
        }
    }

    for (size_t j=0; j<measure_point_number; j++)
    {
        const SpacePoint &mp = _functional->measure_point[j];
        double dx, dy;
        _functional->uu[j][tn.i] = DeltaFunction::lumpedPointG(u, mp, spaceDimensionX(), spaceDimensionY(), 1, 4, dx, dy);
        _functional->ux[j][tn.i] = dx;
        _functional->uy[j][tn.i] = dy;
    }

#endif

    /****************************************************************************************************************************/
}

auto LoadedHeatEquationIBVP::timeDimension() const -> Dimension { return _functional->timeDimension(); }
auto LoadedHeatEquationIBVP::spaceDimensionX() const -> Dimension { return _functional->spaceDimensionX(); }
auto LoadedHeatEquationIBVP::spaceDimensionY() const -> Dimension { return _functional->spaceDimensionY(); }
auto LoadedHeatEquationIBVP::spaceDimensionZ() const -> Dimension { return _functional->spaceDimensionZ(); }

/********************************************************************************************************************************************************/

LoadedHeatEquationFBVP::LoadedHeatEquationFBVP(Functional *function) :_functional(function) {}

LoadedHeatEquationFBVP::LoadedHeatEquationFBVP(const LoadedHeatEquationFBVP &other) { *this = other; }

LoadedHeatEquationFBVP& LoadedHeatEquationFBVP::operator =(const LoadedHeatEquationFBVP &) { return *this; }

LoadedHeatEquationFBVP::~LoadedHeatEquationFBVP() {}

auto LoadedHeatEquationFBVP::final(const SpaceNodePDE &sn, FinalCondition /*condition*/) const -> double
{
    const size_t i = static_cast<size_t>(sn.i), j = static_cast<size_t>(sn.j);
    return -2.0*(_functional->U[j][i] - _functional->V[j][i]);
}

auto LoadedHeatEquationFBVP::boundary(const SpaceNodePDE &/*sn*/, const TimeNodePDE &/*tn*/, BoundaryConditionPDE &bc) const -> double
{
    const double lambda1 = _functional->_lambda1;
    //bc = BoundaryConditionPDE::Dirichlet(1.0, 1.0); return 0.0;
    //bc = BoundaryConditionPDE::Neumann(1.0, 0.0); return 0.0;
    bc = BoundaryConditionPDE::Robin(lambda1, +1.0, lambda1); return 0.0;
}

auto LoadedHeatEquationFBVP::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double
{
    double fx = 0.0;

    const size_t ln = static_cast<size_t>(tn.i);
    const double dimX_step = spaceDimensionX().step();
    const double dimY_step = spaceDimensionY().step();
    const size_t heating_point_number = _functional->heating_point_number;
    const size_t measure_point_number = _functional->measure_point_number;
    const DoubleMatrix &k = _functional->k;
    const DoubleVector *pp = _functional->pp;

    double sum = 0.0;
#ifdef SYNTEZ
    for (size_t j=0; j<measure_point_number; j++)
    {
        const SpacePoint &mp = _functional->measure_point[j];

        for (size_t i=0; i<heating_point_number; i++)
        {
            sum -= k.at(i,j) * pp[i][ln] * DeltaFunction::gaussian(sn, mp, SpacePoint(dimX_step, dimY_step));
        }
    }
#endif

    return fx + sum;
}

auto LoadedHeatEquationFBVP::layerInfo(const DoubleMatrix &p, const TimeNodePDE &tn) const -> void
{
    /****************************************************************************************************************************/

    const size_t heating_point_number = _functional->heating_point_number;

    if (_functional->pp == nullptr)
    {
        const size_t time_size = timeDimension().size();

        _functional->pp = new DoubleVector[heating_point_number];
        _functional->px = new DoubleVector[heating_point_number];
        _functional->py = new DoubleVector[heating_point_number];

        for (size_t i=0; i<_functional->heating_point_number; i++)
        {
            _functional->pp[i].resize(time_size);
            _functional->px[i].resize(time_size);
            _functional->py[i].resize(time_size);
        }
    }

    for (size_t i=0; i<heating_point_number; i++)
    {
        const SpacePoint &spi = _functional->tr(tn, i);
        double dx, dy;
        _functional->pp[i][tn.i] = DeltaFunction::lumpedPointG(p, spi, spaceDimensionX(), spaceDimensionY(), 1, 4, dx, dy);
        _functional->px[i][tn.i] = dx;
        _functional->py[i][tn.i] = dy;
    }

    /****************************************************************************************************************************/
}

auto LoadedHeatEquationFBVP::timeDimension() const -> Dimension { return _functional->timeDimension(); }
auto LoadedHeatEquationFBVP::spaceDimensionX() const -> Dimension { return _functional->spaceDimensionX(); }
auto LoadedHeatEquationFBVP::spaceDimensionY() const -> Dimension { return _functional->spaceDimensionY(); }
auto LoadedHeatEquationFBVP::spaceDimensionZ() const -> Dimension { return _functional->spaceDimensionZ(); }
