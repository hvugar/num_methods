#include "problem22dex5.h"

void Problem22DEx5::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    Problem22DEx5 p22d5;
    p22d5.espilon = 0.0;

    {
        Parameter rpm0(1, 3);

        rpm0.k[0][0] = -10.0; //rpm0.k[1][0] = -2.8;
        rpm0.k[0][1] = -11.0; //rpm0.k[1][1] = +3.2;
        rpm0.k[0][2] = -12.0; //rpm.k[1][2] = -10.7;

        rpm0.z[0][0] = 10.5; //rpm0.z[1][0] = 11.7;
        rpm0.z[0][1] = 12.4; //rpm0.z[1][1] = 12.5;
        rpm0.z[0][2] = 13.4; //rpm0.z[1][1] = 12.5;

        rpm0.eta[0].setPoint(0.40, 0.60);
        //rpm0.eta[1].setPoint(0.70, 0.40);

        rpm0.xi[0].setPoint(0.50, 0.80);
        rpm0.xi[1].setPoint(0.80, 0.70);
        rpm0.xi[2].setPoint(0.50, 0.50);

        p22d5.mParameter0 = rpm0;
        //p22d5.calculateU(rpm0);

        //IPrinter::printSeperatorLine("U");
        //IPrinter::printMatrix(p22d5.U);

        //QPixmap px;
        //visualizeMatrixHeat(p22d5.U, 0.1, p22d5.U.max(), px);
        //px.save("U.png", "PNG");
    }

    Parameter rpm(1, 3);

    rpm.k[0][0] = -12.6; //rpm.k[1][0] = -12.8;
    rpm.k[0][1] = -15.3; //rpm.k[1][1] = -13.2;
    rpm.k[0][2] = -11.2; //rpm.k[1][2] = -10.7;

    rpm.z[0][0] = +14.5; //rpm.z[1][0] = +11.7;
    rpm.z[0][1] = +11.4; //rpm.z[1][1] = +12.5;
    rpm.z[0][2] = +10.5; //rpm.z[1][2] = +12.5;

    //rpm.eta[0].setPoint(0.4000, 0.6000);
    rpm.eta[0].setPoint(0.4264, 0.6354);
    //rpm.eta[1].setPoint(0.7547, 0.4124);

    //rpm.xi[0].setPoint(0.5000, 0.3000);
    rpm.xi[0].setPoint(0.5341, 0.8248);
    rpm.xi[1].setPoint(0.8516, 0.7329);
    rpm.xi[2].setPoint(0.5157, 0.5359);

    p22d5.setParameter(rpm);

    //    FILE *file = fopen("file2.txt", "w");
    //    for (unsigned int n=0; n<100; n++)
    //    {
    //        rpm.eta[0].setPoint(0.4124, n*0.01);
    //        DoubleVector pv1;
    //        rpm.toVector(pv1);
    //        double f = p22d5.fx(pv1);
    //        pv1.clear();
    //        fprintf(file, "%d %.10f\n", n, f);
    //        printf("%d %f\n", n, f);
    //    }
    //    fclose(file);

//    DoubleVector x;
//    rpm.toVector(x);
//    p22d5.optimization(x);

//    return;
    IPrinter::printSeperatorLine();
    DoubleVector pv;
    rpm.toVector(pv);
    DoubleVector ag(pv.length());
    IPrinter::print(pv, pv.length(), 10, 6);
    IPrinter::printSeperatorLine();

    puts("Calculating gradients....");
    p22d5.gradient(pv,ag);
    puts("Gradients are calculated.");

    double functional = p22d5.fx(pv);
    printf("Functional: %f\n", functional);

    DoubleVector ng(pv.length());
    ng.resize(pv.length(), 0.0);

    puts("Calculating numerical gradients....");
    IGradient::Gradient(&p22d5, 0.001, pv, ng);
    puts("Numerical gradients are calculated.");

    //k------------------------------------------------------//
    DoubleVector pk = pv.mid(0, rpm.Lc*rpm.Lo-1);
    DoubleVector ak = ag.mid(0, rpm.Lc*rpm.Lo-1);
    DoubleVector nk = ng.mid(0, rpm.Lc*rpm.Lo-1);

    IPrinter::print(pk);
    IPrinter::print(ak); ak.L2Normalize();
    IPrinter::print(nk); nk.L2Normalize();
    IPrinter::print(ak);
    IPrinter::print(nk);
    IPrinter::printSeperatorLine();

    //z------------------------------------------------------//
    DoubleVector pz = pv.mid(rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo-1);
    DoubleVector az = ag.mid(rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo-1);
    DoubleVector nz = ng.mid(rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo-1);

    IPrinter::print(pz);
    IPrinter::print(az); az.L2Normalize();
    IPrinter::print(nz); nz.L2Normalize();
    IPrinter::print(az);
    IPrinter::print(nz);
    IPrinter::printSeperatorLine();

    //eta------------------------------------------------------//
    DoubleVector pe = pv.mid(2*rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo+2*rpm.Lc-1);
    DoubleVector ae = ag.mid(2*rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo+2*rpm.Lc-1);
    DoubleVector ne = ng.mid(2*rpm.Lc*rpm.Lo, 2*rpm.Lc*rpm.Lo+2*rpm.Lc-1);

    IPrinter::print(pe);
    IPrinter::print(ae); ae.L2Normalize();
    IPrinter::print(ne); ne.L2Normalize();
    IPrinter::print(ae);
    IPrinter::print(ne);
    IPrinter::printSeperatorLine();

    //xi------------------------------------------------------//
    DoubleVector px = pv.mid(2*rpm.Lc*rpm.Lo+2*rpm.Lc, 2*rpm.Lc*rpm.Lo+2*rpm.Lc+2*rpm.Lo-1);
    DoubleVector ax = ag.mid(2*rpm.Lc*rpm.Lo+2*rpm.Lc, 2*rpm.Lc*rpm.Lo+2*rpm.Lc+2*rpm.Lo-1);
    DoubleVector nx = ng.mid(2*rpm.Lc*rpm.Lo+2*rpm.Lc, 2*rpm.Lc*rpm.Lo+2*rpm.Lc+2*rpm.Lo-1);

    IPrinter::print(px);
    IPrinter::print(ax); ax.L2Normalize();
    IPrinter::print(nx); nx.L2Normalize();
    IPrinter::print(ax);
    IPrinter::print(nx);
    IPrinter::printSeperatorLine();
}

Problem22DEx5::Problem22DEx5() : AbstactProblem22D()
{
    forward->setEquationParameters(1.0, 0.01, 0.01, 6.3);
    backward->setEquationParameters(1.0, 0.01, 0.01, 6.3);
    forward->fi = 0.1;

    Dimension time = Dimension(0.005, 0, 200);
    Dimension dimX = Dimension(0.01, 0, 100);
    Dimension dimY = Dimension(0.01, 0, 100);

    setGridParameters(time, dimX, dimY);

    U.resize(dimY.sizeN()+1, dimX.sizeN()+1, 10.0);
}
