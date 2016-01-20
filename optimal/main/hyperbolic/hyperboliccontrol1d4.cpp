#include "hyperboliccontrol1d4.h"
#include <integral.h>

double gause;

void HyperbolicControl1D4::main()
{
    HyperbolicControl1D4 hc;
    hc.calculateSettings();

        hc.e[0] = 0.4;

        double a, b, t;
        a = 0.7;
        b = 1.0;
        goldenSectionSearch(a, b, t, &hc, 0.00001);
        printf("%f\n", t);

    //hc.fx(1.0);

    //for (double e = 0.1; e < 1.0; e+=0.1)
    {
        //        hc.e[0] = 0.4;
        //for (double t=0.4; t<1.21; t+=0.1) hc.fx(t);
        //        hc.fx(1.2);
    }

    hc.count = 0;
    hc.fx(t);
    printf("Function call count: %u\n", hc.count);
}

HyperbolicControl1D4::HyperbolicControl1D4()
{
    t1 = 1.0;
    //calculateSettings();
    count = 0;
}

void HyperbolicControl1D4::calculateSettings()
{
    t0 = 0.0;

    x0 = 0.0;
    x1 = 1.0;

    N = 100;
    hx = (x1-x0)/N;

    ht = 0.01;
    M = (unsigned int) round((t1-t0)/ht);
    D = 10;

    a = 1.0;
    lamda = 0.25;
    R = 1.0;
    U = 0.0;

    L = 1;
    e.resize(L);
}

double HyperbolicControl1D4::fx(const DoubleVector &v)
{
    count++;
    DoubleMatrix u;
    calculateU(v, u);

    //    struct Integral : public R1Function {
    //        virtual double fx(double x)
    //        {
    //            return x - U;
    //        }
    //    } intf;
    double sum = 0.0;

    double integral = 0.0;
    int count = 0;
    for (unsigned int j=M; j<=M+D-1; j++)
    {
        for (unsigned int i=0; i<=N-1; i++)
        {
            double f1 = u[j+0][i+0] - U;
            double f2 = u[j+0][i+1] - U;
            double f3 = u[j+1][i+0] - U;
            double f4 = u[j+1][i+1] - U;
            count += 4;

            integral = integral + (f1*f1 + f2*f2 + f3*f3 + f4*f4);
        }
    }
    integral = integral * 0.25 * ht *hx;

//    printf("%d %.16f\n", count, integral);

//    integral = 0.0;
//    count = 4;
//    double f1 = (u[M][0] - U);
//    double f2 = (u[M][N] - U);
//    double f3 = (u[M+D][0] - U);
//    double f4 = (u[M+D][N] - U);

//    integral += 0.25*hx*ht*(f1*f1 + f2*f2 + f3*f3 + f4*f4);

//    double integral1=0.0;
//    for (unsigned j=M+1; j<=M+D-1; j++)
//    {
//        double f1 = (u[j][0]-U);
//        double f2 = (u[j][N]-U);
//        integral1 += f1*f1;
//        integral1 += f2*f2;
//        count += 2;
//    }
//    integral1 = integral1*0.50*hx*ht;
//    integral += integral1;

//    double integral2=0.0;
//    for (unsigned i=1; i<=N-1; i++)
//    {
//        double f1 = (u[M][i]-U);
//        double f2 = (u[M+D][i]-U);
//        integral2 += f1*f1;
//        integral2 += f2*f2;
//        count += 2;
//    }
//    integral2 = integral2*0.50*hx*ht;
//    integral += integral2;

//    for (unsigned j=M+1; j<=M+D-1; j++)
//    {
//        for (unsigned i=1; i<=N-1; i++)
//        {
//            integral += hx*ht*(u[j][i]-U)*(u[j][i]-U);
//            count += 1;
//        }
//    }
//    //printf("%d %.16f\n", count, integral);

    return sum + integral;
}

void HyperbolicControl1D4::gradient(const DoubleVector &v, DoubleVector &g)
{
    DoubleMatrix u;
    calculateU(v, u);
    pu = &u;
    calculateP(u, g);


//    for (unsigned int i=0; i<v.size(); i++)
//    {
//        DoubleVector v1 = v;
//        DoubleVector v2 = v;
//        double h = 0.01;
//        v1[i] = v[i] + h;
//        v2[i] = v[i] - h;
//        g[i] = (fx(v1)-fx(v2))/(2.0*h);
//    }

    //    double h = 0.01;
    //    DoubleVector v1 = v;
    //    for (unsigned int j=0; j<=M; j++)
    //    {
    //        double i = 2*(M+D+1)+j;
    //        double _v = v1[i];
    //        v1[i] = _v + h;
    //        double f1 = fx(v1);
    //        v1[i] = _v - h;
    //        double f2 = fx(v1);
    //        g[i] = (f1 - f2) / (2.0*h);
    //        //printf("%.16f %.16f %.16f\n", f1, f2, g[i]);
    //    }
}

double HyperbolicControl1D4::fx(double t)
{
    t1 = t;
    calculateSettings();
    //e[0] = 0.25;
    //e[1] = 0.66;
    //e[2] = 0.75;

    //puts("------------------------------------------------------------------------------------------------------");
    //printf("T: %.8f\n", t);

    DoubleVector v((M+D+1)*(L+2));
    for (unsigned int j=0; j<=(M+D); j++)
    {
        //double t = j*ht;
        v[0*(M+D+1)+j] = 1.0;
        v[1*(M+D+1)+j] = 1.0;
        v[2*(M+D+1)+j] = 1.0;
        //v[3*(M+D+1)+j] = U;
        //v[4*(M+D+1)+j] = U;
    }

    //printf("bashlangic: %.16f\n", fx(v));
    //return 0.0;

    //printf("%.8f %.8f %.8f %.8f %.8f %.8f %u %u %u\n", t0, t1, x0, x1, ht, hx, M, N, D);

    double min_step = 1.0;
    double gold_eps = 0.00001;

    ConjugateGradient cg;
    cg.setFunction(this);
    cg.setGradient(this);
    cg.setEpsilon1(0.0000001);
    cg.setEpsilon2(0.0000001);
    //g.setGradientStep(0.000001);
    //    g.setR1MinimizeEpsilon(20.0, 0.01);
    cg.setR1MinimizeEpsilon(min_step, gold_eps);
    cg.setPrinter(this);
    cg.setProjection(this);
    cg.setNormalize(true);
    cg.showEndMessage(false);
    cg.calculate(v);

    double rf = fx(v);

    DoubleVector gr1(v.size());
    for (unsigned int i=0; i<v.size(); i++)
    {
        DoubleVector v1 = v;
        DoubleVector v2 = v;
        double h = 0.01;
        v1[i] = v[i] + h;
        v2[i] = v[i] - h;
        gr1[i] = (fx(v1)-fx(v2))/(2.0*h);
    }

    DoubleMatrix u;
    DoubleVector gr2(v.size());
    calculateU(v, u);
    pu = &u;
    calculateP(u, gr2);

    FILE* file = fopen("20160120_1.txt", "w");
    gr1.L2Normalize();
    gr2.L2Normalize();

    fprintf(file, "t: %f\n", t);
    Printer::printVector(v, "v1: ", (M+D+1)/10, 0*(M+D+1), 0*(M+D+1)+(M+D), file);
    Printer::printVector(v, "v2: ", (M+D+1)/10, 1*(M+D+1), 1*(M+D+1)+(M+D), file);
    Printer::printVector(v, "v3: ", (M+D+1)/10, 2*(M+D+1), 2*(M+D+1)+(M+D), file);
    fputs("\n", file);
    Printer::printVector(gr1, "g11:", (M+D+1)/10, 0*(M+D+1), 0*(M+D+1)+(M+D), file);
    Printer::printVector(gr1, "g12:", (M+D+1)/10, 1*(M+D+1), 1*(M+D+1)+(M+D), file);
    Printer::printVector(gr1, "g13:", (M+D+1)/10, 2*(M+D+1), 2*(M+D+1)+(M+D), file);
    fputs("\n", file);
    Printer::printVector(gr2, "g21:", (M+D+1)/10, 0*(M+D+1), 0*(M+D+1)+(M+D), file);
    Printer::printVector(gr2, "g22:", (M+D+1)/10, 1*(M+D+1), 1*(M+D+1)+(M+D), file);
    Printer::printVector(gr2, "g23:", (M+D+1)/10, 2*(M+D+1), 2*(M+D+1)+(M+D), file);

    fclose(file);

/*

    //    DoubleVector v1(D+1); for (unsigned j=M; j<=M+D; j++) v1[j] = v[0*(M+D+1)+j];
    //    DoubleVector v2(D+1); for (unsigned j=M; j<=M+D; j++) v2[j] = v[1*(M+D+1)+j];
    //    DoubleVector v3(D+1); for (unsigned j=M; j<=M+D; j++) v3[j] = v[2*(M+D+1)+j];
    //    Printer::printVector(v1, 10, "v1:\t");
    //    Printer::printVector(v2, 10, "v2:\t");
    //    Printer::printVector(v3, 10, "v2:\t");

    DoubleMatrix u;
    DoubleVector gr(3*(M+D+1));
    calculateU(v, u);
    pu = &u;
    calculateP(u, gr);

    printf("Norm: %.12f %d\n", gr.L2Norm(), cg.count());

    FILE* file = fopen("20151220.txt", "a");
    fprintf(file, "------------------------------------------------------------\n");
    fprintf(file, "e1: %f T: %.8f Functional: %.16f hx: %f ht: %f step: %f gold_epsilon: %f N: %d M: %d\n", e[0], t, rf, hx, ht, min_step, gold_eps, N, M);

//    Printer::printVector(v, "v1:\t", M+D+1, 0, M+D, file);
//    Printer::printVector(v, "v2:\t", M+D+1, (M+D+1), (M+D+1)+M+D+1, file);
//    Printer::printVector(v, "v3:\t", M+D+1, 2*(M+D+1), 2*(M+D+1)+M+D+1, file);
//    Printer::printVector(g, "g1:\t", M+D+1, 0, M+D+1, file);
//    Printer::printVector(g, "g2:\t", M+D+1, (M+D+1), (M+D+1)+M+D+1, file);
//    Printer::printVector(g, "g3:\t", M+D+1, 2*(M+D+1), 2*(M+D+1)+M+D+1, file);

    fprintf(file, "v1:\t");
    for (unsigned int j=0; j<=M+D; j++)
    {
        double v1 = v[0*(M+D+1)+j];
        if (v1<0)
            fprintf(file, "%14.8f ", v1);
        else
            fprintf(file, "%+14.8f ", v1);
    }
    fprintf(file, "\n");

    fprintf(file, "v2:\t");
    for (unsigned int j=0; j<=M+D; j++)
    {
        double v2 = v[1*(M+D+1)+j];
        if (v2<0)
            fprintf(file, "%14.8f ", v2);
        else
            fprintf(file, "%+14.8f ", v2);
    }
    fprintf(file, "\n");

    fprintf(file, "v3:\t");
    for (unsigned int j=0; j<=M+D; j++)
    {
        double v3 = v[2*(M+D+1)+j];
        if (v3<0)
            fprintf(file, "%14.8f ", v3);
        else
            fprintf(file, "%+14.8f ", v3);
    }
    fprintf(file, "\n");

    fprintf(file, "g1:\t");
    for (unsigned int j=0; j<=M+D; j++)
    {
        double g1 = gr[0*(M+D+1)+j];
        if (g1<0)
            fprintf(file, "%14.8f ", g1);
        else
            fprintf(file, "%+14.8f ", g1);
    }
    fprintf(file, "\n");

    fprintf(file, "g2:\t");
    for (unsigned int j=0; j<=M+D; j++)
    {
        double g2 = gr[1*(M+D+1)+j];
        if (g2<0)
            fprintf(file, "%14.8f ", g2);
        else
            fprintf(file, "%+14.8f ", g2);
    }
    fprintf(file, "\n");

    fprintf(file, "g3:\t");
    for (unsigned int j=0; j<=M+D; j++)
    {
        double g3 = gr[2*(M+D+1)+j];
        if (g3<0)
            fprintf(file, "%14.8f ", g3);
        else
            fprintf(file, "%+14.8f ", g3);
    }
    fprintf(file, "\n");

    for (unsigned int j=M; j<=M+D; j++)
    {
        //printf("u[%d]:\t", j);
        //Printer::printVector(u[j]);

        fprintf(file, "u[%d]:\t", j);
        for (unsigned int i=0; i<=N; i++)
        {
            double uji = u[j][i];
            if (uji<0)
                fprintf(file, "%10.8f ", uji);
            else
                fprintf(file, "+%10.8f ", uji);
        }
        fprintf(file, "\n");
    }
    fclose(file);
    //Printer::printVector(u[M], 10, "u[M]:\t");
    //Printer::printVector(u[M+D], 10, "u[M+D]:");

    printf("e1: %f T: %.8f Integral: %.16f M: %d\n", e[0], t, rf, M);
*/
    return rf;
}

void HyperbolicControl1D4::print(unsigned int iteration, const DoubleVector &x, const DoubleVector &gradient, double alpha, RnFunction *fn) const
{
    //printf("alpha: %.16f\n", alpha);
}

void HyperbolicControl1D4::project(DoubleVector &v, int index)
{}

double HyperbolicControl1D4::fi1(unsigned int i) const
{
    double x = i*hx;
    return x*x;
}

double HyperbolicControl1D4::fi2(unsigned int i) const
{
    return 0.0;
}

double HyperbolicControl1D4::m1(unsigned int j) const
{
    const DoubleVector &v = *pv;
    double v1 = v[j];
    //if (M <= j && j <= M+D) v = U;
    return v1;
}

double HyperbolicControl1D4::m2(unsigned int j) const
{
    const DoubleVector &v = *pv;
    double v2 = v[M+D+1 + j];
    //if (M <= j && j <= M+D) v = U;
    return v2;
}

double HyperbolicControl1D4::f(unsigned int i, unsigned int j) const
{
    double x = i*hx;
    //double t = j*ht;
    double sum = 0.0;

    const DoubleVector &v = *pv;
    double v3 = v[2*(M+D+1)+j];

//    double sgm = 3.0*hx;
//    double a = 1.0/(sgm*sqrt(2.0*M_PI));
//    double b = 2.0*sgm*sgm;
//    double g = a * exp(-((x-e[0])*(x-e[0]))/b);
//    sum += v3 * g;

//    gause += g;

        if (fabs(x-e[0]) < (hx+0.000001))
        {
            //const DoubleVector &v = *pv;
            //double v3 = v[2*(M+D+1)+j];
            //if (M <= j && j <= M+D) v = U;
            sum += (1.0/hx) * v3 * ((hx-fabs(x-e[0]))/hx);
            //printf("x: %f e: %f sum: %f %.20f %.20f\n", x, e[0], sum, fabs(x-e[0]), hx);
        }
//    if (i==40) sum = (1/hx)*v3;

    //if (fabs(x-e[1]) < (hx+0.000001))
    //{
    //    double v = (*pv)[3*(M+D+1)+j];
    //    sum += (1.0/hx) * v * ((hx-fabs(x-e[1]))/hx);
    //}

    //if (fabs(x-e[2]) < (hx+0.000001))
    //{
    //    double v = (*pv)[4*(M+D+1)+j];
    //    sum += (1.0/hx) * v * ((hx-fabs(x-e[2]))/hx);
    //}

    return sum;
}

void HyperbolicControl1D4::calculateU(const DoubleVector &v, DoubleMatrix &u)
{
    pv = &v;

    u.clear();
    u.resize(M+D+1);
    for (unsigned int j=0; j<=M+D; j++) u[j].resize(N+1);

    DoubleVector u0(N+1);
    DoubleVector u1(N+1);

    DoubleVector da;
    DoubleVector db;
    DoubleVector dc;
    DoubleVector rd;
    DoubleVector rx;

    da.resize(N-1);
    db.resize(N-1);
    dc.resize(N-1);
    rd.resize(N-1);
    rx.resize(N-1);

    double alpha1 = -(lamda*a*a)*((ht*ht)/(hx*hx));
    double beta1  = 1.0 + (2.0*lamda*a*a*(ht*ht))/(hx*hx);
    double alpha2 = (1.0-2.0*lamda)*a*a*((ht*ht)/(hx*hx));
    double alpha3 = +(lamda*a*a)*((ht*ht)/(hx*hx));

    for (unsigned int j=0; j<=M+D-1; j++)
    {

        if (j==0)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u0[i] = fi1(i);
                u1[i] = fi1(i) + ht*fi2(i);
                u[0][i] = u0[i];
                u[1][i] = u1[i];
            }
        }
        else
        {
            gause = 0.0;
            for (unsigned int i=1; i<=N-1; i++)
            {
                da[i-1] = alpha1;
                db[i-1] = beta1;
                dc[i-1] = alpha1;
                rd[i-1] = alpha2*(u1[i-1]-2.0*u1[i]+u1[i+1]) + 2.0*u1[i] + alpha3*(u0[i+1] - 2.0*u0[i] + u0[i-1]) - u0[i] + (ht*ht)*f(i, j);
            }
            //printf("gause: %.8f %.8f\n", gause, gause*hx);

            da[0]   = 0.0;
            dc[N-2] = 0.0;
            rd[0]   -= alpha1 * m1(j+1);
            rd[N-2] -= alpha1 * m2(j+1);
            TomasAlgorithm(da, db, dc, rd, rx);

            u[j+1][0] = m1(j+1);
            for (unsigned int i=1; i<=N-1; i++)
            {
                u[j+1][i] = rx[i-1];
            }
            u[j+1][N] = m2(j+1);

            for (unsigned int i=0; i<=N; i++)
            {
                u0[i] = u1[i];
                u1[i] = u[j+1][i];
            }
        }
    }

    da.clear();
    db.clear();
    dc.clear();
    rd.clear();
    rx.clear();

    u1.clear();
    u0.clear();
}

void HyperbolicControl1D4::calculateP(const DoubleMatrix &u, DoubleVector &g)
{
    DoubleVector p(N+1);
    DoubleVector p0(N+1);
    DoubleVector p1(N+1);

    DoubleVector da;
    DoubleVector db;
    DoubleVector dc;
    DoubleVector rd;
    DoubleVector rx;

    da.resize(N-1);
    db.resize(N-1);
    dc.resize(N-1);
    rd.resize(N-1);
    rx.resize(N-1);

    double alpha1 = -(lamda*a*a)*((ht*ht)/(hx*hx));
    double beta1  = 1.0 + (2.0*lamda*a*a*(ht*ht))/(hx*hx);
    double alpha2 = (1.0-2.0*lamda)*a*a*((ht*ht)/(hx*hx));
    double alpha3 = +(lamda*a*a)*((ht*ht)/(hx*hx));

    for (unsigned int j1=0; j1<=M+D-1; j1++)
    {
        unsigned int j = M+D-j1;

        if (j==M+D)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                p0[i] = pfi1(i);
                p1[i] = pfi1(i) - ht*pfi2(i);
            }
            calculateG(p0, g, M+D);
            calculateG(p1, g, M+D-1);
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                da[i-1] = alpha1;
                db[i-1] = beta1;
                dc[i-1] = alpha1;
                rd[i-1] = alpha2*(p1[i-1]-2.0*p1[i]+p1[i+1]) + 2.0*p1[i] - p0[i] + alpha3*(p0[i+1] - 2.0*p0[i] + p0[i-1]);

                if (M<=j-1 && j-1<=M+D)
                {
                    rd[i-1] -= (ht*ht)*(2*R*(u[j][i]-U));
                }
            }

            da[0]   = 0.0;
            dc[N-2] = 0.0;
            rd[0]   -= alpha1 * pmu1(j-1);
            rd[N-2] -= alpha1 * pmu2(j-1);
            TomasAlgorithm(da, db, dc, rd, rx);

            p[0] = pmu1(j-1);
            for (unsigned int i=1; i<=N-1; i++)
            {
                p[i] = rx[i-1];
            }
            p[N] = pmu2(j-1);

            for (unsigned int i=0; i<=N; i++)
            {
                p0[i] = p1[i];
                p1[i] = p[i];
            }
            calculateG(p, g, j-1);
        }
    }

    da.clear();
    db.clear();
    dc.clear();
    rd.clear();
    rx.clear();

    p0.clear();
    p1.clear();
}

void HyperbolicControl1D4::calculateG(const DoubleVector& psi, DoubleVector& g, unsigned int j)
{
    g[0*(M+D+1)+j] = -(psi[1]-psi[0])/hx;
    g[1*(M+D+1)+j] = +(psi[N]-psi[N-1])/hx;

    //unsigned int i1 = (unsigned int)(round(e[0]/hx));
    g[2*(M+D+1)+j] = -psi[40];
    //if (j==M+D) printf("%.8f\n", -psi[i1]);
}

