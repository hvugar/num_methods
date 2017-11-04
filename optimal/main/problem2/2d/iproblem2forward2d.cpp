#include "iproblem2forward2d.h"

IProblem2Forward2D::IProblem2Forward2D(double a, double lambda0, double lambda, double theta, unsigned int Lc1, unsigned int Lo1)
{
    this->a = a;
    this->lambda0 = lambda0;
    this->lambda = lambda;
    this->theta = theta;

    this->Lc = Lc1;
    this->Lo = Lo1;

    //    this->Lo = 5;
    //    this->Lc = 4;

    this->Lo = 1;
    this->Lc = 1;

    k.resize(Lc, Lo);
    z.resize(Lc, Lo);

    for (unsigned int i=0; i<this->Lc; i++)
    {
        for (unsigned int j=0; j<this->Lo; j++)
        {
            k[i][j] = 5.0;
            z[i][j] = 8.0;
        }
    }

    //    xi.resize(this->Lo);
    //    xi[0].x = 0.25;  xi[0].y = 0.25; //xi[0].i = 25; xi[0].j = 25;
    //    xi[1].x = 0.25;  xi[1].y = 0.75; //xi[1].i = 25; xi[1].j = 75;
    //    xi[2].x = 0.75;  xi[2].y = 0.75; //xi[2].i = 75; xi[2].j = 75;
    //    xi[3].x = 0.75;  xi[3].y = 0.25; //xi[3].i = 75; xi[3].j = 25;
    //    xi[4].x = 0.50;  xi[4].y = 0.50; //xi[4].i = 50; xi[4].j = 50;

    //    eta.resize(this->Lc);
    //    eta[0].x = 0.33; eta[0].y = 0.33; //eta[0].i = 33; eta[0].j = 33;
    //    eta[1].x = 0.33; eta[1].y = 0.66; //eta[1].i = 33; eta[1].j = 66;
    //    eta[2].x = 0.66; eta[2].y = 0.66; //eta[2].i = 66; eta[2].j = 66;
    //    eta[3].x = 0.66; eta[3].y = 0.33; //eta[3].i = 66; eta[3].j = 33;

    xi.resize(this->Lo);
    xi[0].x = 0.10;  xi[0].y = 0.80; xi[0].i = 1; xi[0].j = 8;
    //xi[1].x = 0.80;  xi[1].y = 0.50; xi[1].i = 80; xi[1].j = 50;
    //xi[2].x = 0.40;  xi[2].y = 0.70; xi[2].i = 40; xi[2].j = 70;

    eta.resize(this->Lc);
    //    eta[0].x = 0.70; eta[0].y = 0.30; eta[0].i = 70; eta[0].j = 30;
    //    eta[1].x = 0.20; eta[1].y = 0.30; eta[1].i = 20; eta[1].j = 60;
    eta[0].x = 0.30; eta[0].y = 0.20; eta[0].i = 3; eta[0].j = 2;
}

void IProblem2Forward2D::calculateMVD(DoubleMatrix &u)
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);
    Dimension td = timeDimension();

    unsigned int N = xd.sizeN();
    unsigned int M = yd.sizeN();
    unsigned int L = td.sizeN();
    double hx = xd.step();
    double hy = yd.step();
    double ht = td.step();

    u.clear();
    u.resize(M+1, N+1);

    DoubleMatrix uh(M+1, N+1);

    //************************************************************************************************************
    for (unsigned int m=0; m<=M; m++)
    {
        SpaceNodePDE sn;
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            u[m][n] = initial(sn);
        }
    }
    //************************************************************************************************************

    IPrinter::printMatrix(16, 12, u);
    IPrinter::printSeperatorLine();

    TimeNodePDE tn;
    for (unsigned int l=1; l<=1; l++)
    {
        //************************************************************************************************************
        //        tn.i = l;
        //        tn.t = l*ht - 0.5*ht;
        //        for (unsigned int n=0; n<=N; n++)
        //        {
        //            SpaceNodePDE sn1;
        //            sn1.i = n; sn1.x = n*hx;

        //            DoubleMatrix w1(M+1, M+1, 0.0);
        //            DoubleVector d1(M+1, 0.0);

        //            for (unsigned int m=0; m<=M; m++)
        //            {
        //                sn1.j = m; sn1.y = m*hy;

        //                d1[m] = 2.0*u[m][n] + lambda0*theta*ht + ht*f(sn1, tn);

        //                if (n==0)       d1[m] += ((a*a*ht)/(hx*hx))*(u[m][0]   - 2.0*u[m][1]   + u[m][2]);
        //                if (n>0 && n<N) d1[m] += ((a*a*ht)/(hx*hx))*(u[m][n-1] - 2.0*u[m][n]   + u[m][n+1]);
        //                if (n==N)       d1[m] += ((a*a*ht)/(hx*hx))*(u[m][N-2] - 2.0*u[m][N-1] + u[m][N]);

        //                if (m == 0)
        //                {
        //                    w1[0][0] = 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht + (2.0*a*a*lambda*ht)/(hy);
        //                    w1[0][1] = -(2.0*a*a*ht)/(hy*hy);

        //                    d1[0] += (2.0*a*a*lambda*theta*ht)/(hy) + ((2.0*a*a*ht)/(hy))*g3(sn1, tn);
        //                }
        //                else if (m == M)
        //                {
        //                    w1[M][M-1] = -(2.0*a*a*ht)/(hy*hy);
        //                    w1[M][M-0] = 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht + (2.0*a*a*lambda*ht)/(hy);

        //                    d1[M] += (2.0*a*a*lambda*theta*ht)/(hy) + ((2.0*a*a*ht)/(hy))*g4(sn1, tn);
        //                }
        //                else
        //                {
        //                    w1[m][m-1] = -(a*a*ht)/(hy*hy);
        //                    w1[m][m+0] = 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht;
        //                    w1[m][m+1] = -(a*a*ht)/(hy*hy);
        //                }

        //                // Adding delta part ***********************************************************************************

        //                for (unsigned int i=0; i<Lc; i++)
        //                {
        //                    double _delta = delta(sn1, i, 1);

        //                    for (unsigned int j=0; j<Lo; j++)
        //                    {
        //                        d1[m] += -ht*k[i][j]*z[i][j] * _delta;

        //                        //unsigned int jinx = (unsigned int)(xi[j].x*N);
        //                        unsigned int jiny = xi[j].j;//(unsigned int)(xi[j].y*M);

        //                        w1[m][jiny] += -ht*k[i][j]*_delta;
        //                        //if (_delta > 0.0) printf("%d %d %d %f %f %f %f %f %f %.6f %d %.14f\n", i, j, jiny, xi[j].x, xi[j].y, sn1.x, sn1.y, eta[i].x, eta[i].y, _delta,m, w1[m][jiny]);

        //                        /*
        //                        double hx3 = hx*hx*hx;
        //                        double hx32 = (1.0/(2.0*hx3));
        //                        double hx36 = (1.0/(6.0*hx3));

        //                        double hy3 = hy*hy*hy;
        //                        double hy32 = (1.0/(2.0*hy3));
        //                        double hy36 = (1.0/(6.0*hy3));

        //                        for (unsigned int m1=0; m1<=M; m1++)
        //                        {

        //                            double dh = fabs(m1*hy - xi[j].y);

        //                            if (dh <= hy)
        //                            {
        //                                w1[m][m1] += -ht*k[i][j] * ((2.0*hy-dh)*(hy-dh)*(hy+dh)) * hy32 * _delta;
        //                            }

        //                            if (hy < dh && dh <= 2.0*hy)
        //                            {
        //                                w1[m][m1] += -ht*k[i][j] * ((2.0*hy-dh)*(hy-dh)*(3.0*hy-dh)) * hy36 * _delta;
        //                            }
        //                        }
        //                        */
        //                    }
        //                }

        //                //*****************************************************************************************************
        //            }

        //            DoubleVector x1(M+1);
        //            LinearEquation::GaussianElimination(w1,d1,x1);
        //            //IPrinter::printVector(x);

        //            w1.clear();
        //            d1.clear();

        //            for (unsigned int m=0; m<=M; m++) uh[m][n] = x1[m];
        //            x1.clear();
        //        }
        //************************************************************************************************************

        for (unsigned int m=0; m<=M; m++)
        {
            SpaceNodePDE sn;
            sn.j = m; sn.y = m*hy;
            for (unsigned int n=0; n<=N; n++)
            {
                sn.i = n; sn.x = n*hx;
                uh[m][n] = initial(sn) + 0.5*ht;
            }
        }

        IPrinter::printMatrix(16, 12, uh);
        IPrinter::printSeperatorLine();

        //************************************************************************************************************
        tn.i = l;
        tn.t = l*ht;
        for (unsigned int m=0; m<=M; m++)
        {
            SpaceNodePDE sn2;
            sn2.j = m; sn2.y = m*hy;

            DoubleMatrix w2(N+1, N+1, 0.0);
            DoubleVector d2(N+1, 0.0);

            for (unsigned int n=0; n<=N; n++)
            {
                sn2.i = n; sn2.x = n*hx;

                double pp = ht*f(sn2, tn);
                d2[n] = 2.0*uh[m][n] + lambda0*theta*ht + pp;

                if (m==0)       d2[n] += ((a*a*ht)/(hy*hy))*(uh[0][n]   - 2.0*uh[1][n]   + uh[2][n]);
                if (m>0 && m<M) d2[n] += ((a*a*ht)/(hy*hy))*(uh[m-1][n] - 2.0*uh[m][n]   + uh[m+1][n]);
                if (m==M)       d2[n] += ((a*a*ht)/(hy*hy))*(uh[M-2][n] - 2.0*uh[M-1][n] + uh[M][n]);

                //if (m==20 && n==30) printf("--- %d %f %f %f %f\n", n, pp, 2.0*uh[m][n] + lambda0*theta*ht, ((a*a*ht)/(hy*hy))*(uh[m-1][n] - 2.0*uh[m][n]   + uh[m+1][n]),d2[30]);

                if (n == 0)
                {
                    w2[0][0] = 2.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht + (2.0*a*a*lambda*ht)/(hx);
                    w2[0][1] = -(2.0*a*a*ht)/(hx*hx);

                    d2[0] += (2.0*a*a*lambda*theta*ht)/(hx) + ((2.0*a*a*ht)/(hx))*g1(sn2, tn);
                }
                else if (n == N)
                {
                    w2[N][N-1] = -(2.0*a*a*ht)/(hx*hx);
                    w2[N][N-0] = 2.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht + (2.0*a*a*lambda*ht)/(hx);

                    d2[N] += (2.0*a*a*lambda*theta*ht)/(hx) + ((2.0*a*a*ht)/(hx))*g2(sn2, tn);
                }
                else
                {
                    w2[n][n-1] = -(a*a*ht)/(hx*hx);
                    w2[n][n+0] = 2.0 + (2.0*a*a*ht)/(hx*hx) + lambda0*ht;
                    w2[n][n+1] = -(a*a*ht)/(hx*hx);
                }

                //if (m==20 && n==30) printf("%d %f\n", n, d2[30]);

                // Adding delta part ***********************************************************************************

                for (unsigned int i=0; i<Lc; i++)
                {
                    double _delta = delta(sn2, i, 2);

                    for (unsigned int j=0; j<Lo; j++)
                    {
                        d2[n] += -ht*k[i][j]*z[i][j] * _delta;

                        //if (m==20 && n==30) printf("%d %f\n", n, -ht*k[i][j]*z[i][j]*_delta);


                        unsigned int jinx = xi[j].j;//(unsigned int)(xi[j].x*N);
                        //unsigned int jiny = (unsigned int)(xi[j].y*M);

                        double kk = -ht*k[i][j]*_delta;
                        w2[n][jinx] += kk;

                        if (m==20 && n==30) printf("++++++++++++++++++++ %f\n", ht*k[i][j]*_delta);

                        //if (m==20 && n==30) printf("%d %f %f %f\n", n, d2[30], -ht*k[i][j]*z[i][j] * _delta, kk);

                        //if (_delta > 0.0) printf("%d %d %d %f %f %f %f %f %f %.14f\n", i, j, jinx, xi[j].x, xi[j].y, sn2.x, sn2.y, eta[i].x, eta[i].y, _delta);
                        //if (_delta > 0.0) printf("%d %d %d %f %f %f %f %f %f %.6f %d %.6f %.6f %d\n", i, j, jinx, xi[j].x, xi[j].y, sn2.x, sn2.y, eta[i].x, eta[i].y, _delta,n, w2[n][jinx], kk, m);

                        /*
                        double hx3 = hx*hx*hx;
                        double hx32 = (1.0/(2.0*hx3));
                        double hx36 = (1.0/(6.0*hx3));

                        double hy3 = hy*hy*hy;
                        double hy32 = (1.0/(2.0*hy3));
                        double hy36 = (1.0/(6.0*hy3));

                        for (unsigned int n1=0; n1<=N; n1++)
                        {

                            double dh = fabs(n1*hx - xi[j].x);

                            if (dh <= hx)
                            {
                                w2[n][n1] += -ht*k[i][j] * ((2.0*hx-dh)*(hx-dh)*(hx+dh)) * hx32 * _delta;
                            }

                            if (hx < dh && dh <= 2.0*hx)
                            {
                                w2[n][n1] += -ht*k[i][j] * ((2.0*hx-dh)*(hx-dh)*(3.0*hx-dh)) * hx36 * _delta;
                            }
                        }
                        */
                    }
                }



                //*****************************************************************************************************
            }

            //if (m==20)         d2[30] = -5999.717;
            DoubleVector x2(N+1);
            LinearEquation::GaussianElimination(w2,d2,x2);
            //IPrinter::printVector(x);

            if (m==20)
            {

                double sum = 0.0;
                for (unsigned int n=0; n<=N; n++)
                {
                    sum += w2.at(3, n) * U(n*hx, 0.2, 0.1);
                }

                printf("sum: %f\n", sum);

                IPrinter::printMatrix(18, 6, w2);

                //                printf("%f %f %f\n", sn2.x, sn2.y, tn.t);

                //                FILE* file = fopen("D:\\data_m_1.txt", "w");
                //                IPrinter::printMatrix(10, 4, w2, M, N, NULL, file);
                //                fclose(file);

                //                FILE* file1 = fopen("D:\\data_d_1.txt", "w");
                //                IPrinter::printVector(d2, NULL, N, 0,0, file1);
                //                fclose(file1);
            }

            w2.clear();
            d2.clear();

            for (unsigned int n=0; n<=N; n++) u[m][n] = x2[n];
            x2.clear();
        }

        IPrinter::printMatrix(16, 12, u);
        IPrinter::printSeperatorLine();
    }
}

void IProblem2Forward2D::calculateMVD1(DoubleMatrix &u)
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);
    Dimension td = timeDimension();

    unsigned int N = xd.sizeN();
    unsigned int M = yd.sizeN();
    unsigned int L = td.sizeN();
    double hx = xd.step();
    double hy = yd.step();
    double ht = td.step();

    u.clear();
    u.resize(M+1, N+1);

    DoubleMatrix uh(M+1, N+1);

    //************************************************************************************************************
    for (unsigned int m=0; m<=M; m++)
    {
        SpaceNodePDE sn;
        sn.j = m; sn.y = m*hy;
        for (unsigned int n=0; n<=N; n++)
        {
            sn.i = n; sn.x = n*hx;
            u[m][n] = initial(sn);
        }
    }
    //************************************************************************************************************

    IPrinter::printMatrix(16, 12, u);
    IPrinter::printSeperatorLine();

    TimeNodePDE tn;
    for (unsigned int l=1; l<=1; l++)
    {
        //************************************************************************************************************
        //        tn.i = l;
        //        tn.t = l*ht - 0.5*ht;
        //        for (unsigned int n=0; n<=N; n++)
        //        {
        //            SpaceNodePDE sn1;
        //            sn1.i = n; sn1.x = n*hx;

        //            DoubleMatrix w1(M+1, M+1, 0.0);
        //            DoubleVector d1(M+1, 0.0);

        //            for (unsigned int m=0; m<=M; m++)
        //            {
        //                sn1.j = m; sn1.y = m*hy;

        //                d1[m] = 2.0*u[m][n] + lambda0*theta*ht + ht*f(sn1, tn);

        //                if (n==0)       d1[m] += ((a*a*ht)/(hx*hx))*(u[m][0]   - 2.0*u[m][1]   + u[m][2]);
        //                if (n>0 && n<N) d1[m] += ((a*a*ht)/(hx*hx))*(u[m][n-1] - 2.0*u[m][n]   + u[m][n+1]);
        //                if (n==N)       d1[m] += ((a*a*ht)/(hx*hx))*(u[m][N-2] - 2.0*u[m][N-1] + u[m][N]);

        //                if (m == 0)
        //                {
        //                    w1[0][0] = 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht + (2.0*a*a*lambda*ht)/(hy);
        //                    w1[0][1] = -(2.0*a*a*ht)/(hy*hy);

        //                    d1[0] += (2.0*a*a*lambda*theta*ht)/(hy) + ((2.0*a*a*ht)/(hy))*g3(sn1, tn);
        //                }
        //                else if (m == M)
        //                {
        //                    w1[M][M-1] = -(2.0*a*a*ht)/(hy*hy);
        //                    w1[M][M-0] = 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht + (2.0*a*a*lambda*ht)/(hy);

        //                    d1[M] += (2.0*a*a*lambda*theta*ht)/(hy) + ((2.0*a*a*ht)/(hy))*g4(sn1, tn);
        //                }
        //                else
        //                {
        //                    w1[m][m-1] = -(a*a*ht)/(hy*hy);
        //                    w1[m][m+0] = 2.0 + (2.0*a*a*ht)/(hy*hy) + lambda0*ht;
        //                    w1[m][m+1] = -(a*a*ht)/(hy*hy);
        //                }

        //                // Adding delta part ***********************************************************************************

        //                for (unsigned int i=0; i<Lc; i++)
        //                {
        //                    double _delta = delta(sn1, i, 1);

        //                    for (unsigned int j=0; j<Lo; j++)
        //                    {
        //                        d1[m] += -ht*k[i][j]*z[i][j] * _delta;

        //                        //unsigned int jinx = (unsigned int)(xi[j].x*N);
        //                        unsigned int jiny = xi[j].j;//(unsigned int)(xi[j].y*M);

        //                        w1[m][jiny] += -ht*k[i][j]*_delta;
        //                        //if (_delta > 0.0) printf("%d %d %d %f %f %f %f %f %f %.6f %d %.14f\n", i, j, jiny, xi[j].x, xi[j].y, sn1.x, sn1.y, eta[i].x, eta[i].y, _delta,m, w1[m][jiny]);

        //                        /*
        //                        double hx3 = hx*hx*hx;
        //                        double hx32 = (1.0/(2.0*hx3));
        //                        double hx36 = (1.0/(6.0*hx3));

        //                        double hy3 = hy*hy*hy;
        //                        double hy32 = (1.0/(2.0*hy3));
        //                        double hy36 = (1.0/(6.0*hy3));

        //                        for (unsigned int m1=0; m1<=M; m1++)
        //                        {

        //                            double dh = fabs(m1*hy - xi[j].y);

        //                            if (dh <= hy)
        //                            {
        //                                w1[m][m1] += -ht*k[i][j] * ((2.0*hy-dh)*(hy-dh)*(hy+dh)) * hy32 * _delta;
        //                            }

        //                            if (hy < dh && dh <= 2.0*hy)
        //                            {
        //                                w1[m][m1] += -ht*k[i][j] * ((2.0*hy-dh)*(hy-dh)*(3.0*hy-dh)) * hy36 * _delta;
        //                            }
        //                        }
        //                        */
        //                    }
        //                }

        //                //*****************************************************************************************************
        //            }

        //            DoubleVector x1(M+1);
        //            LinearEquation::GaussianElimination(w1,d1,x1);
        //            //IPrinter::printVector(x);

        //            w1.clear();
        //            d1.clear();

        //            for (unsigned int m=0; m<=M; m++) uh[m][n] = x1[m];
        //            x1.clear();
        //        }
        //************************************************************************************************************

        for (unsigned int m=0; m<=M; m++)
        {
            SpaceNodePDE sn;
            sn.j = m; sn.y = m*hy;
            for (unsigned int n=0; n<=N; n++)
            {
                sn.i = n; sn.x = n*hx;
                uh[m][n] = initial(sn) + 0.5*ht;
            }
        }

        IPrinter::printMatrix(16, 12, uh);
        IPrinter::printSeperatorLine();

        //************************************************************************************************************
        tn.i = l;
        tn.t = l*ht;
        for (unsigned int m=0; m<=M; m++)
        {
            SpaceNodePDE sn2;
            sn2.j = m; sn2.y = m*hy;

            DoubleMatrix w2(N+1, N+1, 0.0);
            DoubleVector d2(N+1, 0.0);

            for (unsigned int n=0; n<=N; n++)
            {
                sn2.i = n; sn2.x = n*hx;

                d2[n] = (2.0/ht)*uh[m][n] + lambda0*theta + f(sn2, tn);

                if (m==0)       d2[n] += ((a*a)/(hy*hy))*(uh[0][n]   - 2.0*uh[1][n]   + uh[2][n]);
                if (m>0 && m<M) d2[n] += ((a*a)/(hy*hy))*(uh[m-1][n] - 2.0*uh[m][n]   + uh[m+1][n]);
                if (m==M)       d2[n] += ((a*a)/(hy*hy))*(uh[M-2][n] - 2.0*uh[M-1][n] + uh[M][n]);

                //if (m==20 && n==30) printf("--- %d %f %f %f %f\n", n, pp, 2.0*uh[m][n] + lambda0*theta, ((a*a)/(hy*hy))*(uh[m-1][n] - 2.0*uh[m][n]   + uh[m+1][n]),d2[30]);

                if (n == 0)
                {
                    w2[0][0] = 2.0/ht + (2.0*a*a)/(hx*hx) + lambda0 + (2.0*a*a*lambda)/(hx);
                    w2[0][1] = -(2.0*a*a)/(hx*hx);

                    d2[0] += (2.0*a*a*lambda*theta)/(hx) + ((2.0*a*a)/(hx))*g1(sn2, tn);
                }
                else if (n == N)
                {
                    w2[N][N-1] = -(2.0*a*a)/(hx*hx);
                    w2[N][N-0] = 2.0/ht + (2.0*a*a)/(hx*hx) + lambda0 + (2.0*a*a*lambda)/(hx);

                    d2[N] += (2.0*a*a*lambda*theta)/(hx) + ((2.0*a*a)/(hx))*g2(sn2, tn);
                }
                else
                {
                    w2[n][n-1] = -(a*a)/(hx*hx);
                    w2[n][n+0] = 2.0/ht + (2.0*a*a)/(hx*hx) + lambda0;
                    w2[n][n+1] = -(a*a)/(hx*hx);
                }

                //if (m==20 && n==30) printf("%d %f\n", n, d2[30]);

                // Adding delta part ***********************************************************************************

                for (unsigned int i=0; i<Lc; i++)
                {
                    double _delta = delta(sn2, i, 2);

                    for (unsigned int j=0; j<Lo; j++)
                    {
                        d2[n] += -k[i][j]*z[i][j] * _delta;

                        //if (m==20 && n==30) printf("%d %f\n", n, -ht*k[i][j]*z[i][j]*_delta);


                        unsigned int jinx = xi[j].j;//(unsigned int)(xi[j].x*N);
                        //unsigned int jiny = (unsigned int)(xi[j].y*M);

                        if (m==2 && n==3) printf("%f %d %d %d\n", -k[i][j]*_delta, jinx, m, n);
                        w2[n][jinx] += -k[i][j]*_delta;

                        //if (m==20 && n==30) printf("%d %f %f %f\n", n, d2[30], -k[i][j]*z[i][j] * _delta, kk);

                        //if (_delta > 0.0) printf("%d %d %d %f %f %f %f %f %f %.14f\n", i, j, jinx, xi[j].x, xi[j].y, sn2.x, sn2.y, eta[i].x, eta[i].y, _delta);
                        //if (_delta > 0.0) printf("%d %d %d %f %f %f %f %f %f %.6f %d %.6f %.6f %d\n", i, j, jinx, xi[j].x, xi[j].y, sn2.x, sn2.y, eta[i].x, eta[i].y, _delta,n, w2[n][jinx], kk, m);

                        /*
                        double hx3 = hx*hx*hx;
                        double hx32 = (1.0/(2.0*hx3));
                        double hx36 = (1.0/(6.0*hx3));

                        double hy3 = hy*hy*hy;
                        double hy32 = (1.0/(2.0*hy3));
                        double hy36 = (1.0/(6.0*hy3));

                        for (unsigned int n1=0; n1<=N; n1++)
                        {

                            double dh = fabs(n1*hx - xi[j].x);

                            if (dh <= hx)
                            {
                                w2[n][n1] += -ht*k[i][j] * ((2.0*hx-dh)*(hx-dh)*(hx+dh)) * hx32 * _delta;
                            }

                            if (hx < dh && dh <= 2.0*hx)
                            {
                                w2[n][n1] += -ht*k[i][j] * ((2.0*hx-dh)*(hx-dh)*(3.0*hx-dh)) * hx36 * _delta;
                            }
                        }
                        */
                    }
                }



                //*****************************************************************************************************
            }

//            if (m==2)         d2[3] = -387.17;
            DoubleVector x2(N+1);
            LinearEquation::GaussianElimination(w2,d2,x2);
            //IPrinter::printVector(x);

            if (m==2)
            {

                double sum = 0.0;
                for (unsigned int n=0; n<=N; n++)
                {
                    sum += w2.at(3, n) * U(n*hx, 0.2, 0.1);
                }

                printf("sum: %f %f\n", sum, d2[3]);

                //                printf("%f %f %f\n", sn2.x, sn2.y, tn.t);

                FILE* file = fopen("D:\\data_m_1.txt", "w");
                IPrinter::printMatrix(10, 4, w2, M, N, NULL, file);
                fclose(file);

                FILE* file1 = fopen("D:\\data_d_1.txt", "w");
                IPrinter::printVector(d2, NULL, N, 0,0, file1);
                fclose(file1);
            }

            w2.clear();
            d2.clear();

            for (unsigned int n=0; n<=N; n++) u[m][n] = x2[n];
            x2.clear();
        }

        IPrinter::printMatrix(16, 12, u);
        IPrinter::printSeperatorLine();
    }
}

double IProblem2Forward2D::initial(const SpaceNodePDE &sn) const
{
    double x = sn.x;
    double y = sn.y;
    return x*x + y*y;
}

double IProblem2Forward2D::boundary(const SpaceNodePDE &sn UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM, BoundaryType boundary UNUSED_PARAM) const
{
    return NAN;
}

double IProblem2Forward2D::f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double x = sn.x;
    double y = sn.y;
    double t = tn.t;

    //printf("t: %f\n", t);

    double res = 1.0 - 4.0*a*a + lambda0*(U(x,y,t) - theta);

    for (unsigned int i=0; i<Lc; i++)
    {
        double _delta = delta(sn, i, 2);

        for (unsigned int j=0; j<Lo; j++)
        {
            //printf("---- %f %f %f %f\n", xi[j].x, xi[j].y, t, k[i][j]*U(xi[j].x, xi[j].y, t)*_delta);
            res += -k[i][j]*U(xi[j].x, xi[j].y, t) * _delta;
        }
    }

    for (unsigned int i=0; i<Lc; i++)
    {
        double _delta = delta(sn, i, 2);

        for (unsigned int j=0; j<Lo; j++)
        {
            res += +k[i][j]*z[i][j] * _delta;
        }
    }


//    double W = 0.0;
//    for (unsigned int i=0; i<Lc; i++)
//    {
//        if ( sn.i == eta[i].i && sn.j == eta[i].j )
//        {

//            double _delta = delta(sn, i, 0);

//            double vi = 0.0;
//            for (unsigned int j=0; j<Lo; j++)
//            {
//                //SpaceNodePDE xij = xi[j];

//                if (_delta > 0.0) printf("------------------------------- %d %d %f %f %.14f\n", i, j, xi[j].x, xi[j].y, _delta);

//                double u = U(xi[j].x, xi[j].y, t);
//                vi += k[i][j] * ( u  - z[i][j]) * _delta;

//                printf("-------------- %f\n", k[i][j] * u * _delta);
//            }
//            //printf("---------- %d %f\n", i, vi);
//            W += vi;
//        }
//    }
//    res -= W;

    return res;
}

double IProblem2Forward2D::delta(const SpaceNodePDE &sn UNUSED_PARAM, unsigned int i UNUSED_PARAM, unsigned int source) const
{
    //return 0.0;

    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);
    double hx = xd.step();
    double hy = yd.step();

    // Approximation delta function using normal distribution formula
    static double sigmaX = hx;
    static double sigmaY = hy;
    //double res = (1.0/(2.0*M_PI*sigmaX*sigmaY)) * exp( - (((sn.x-eta[i].x)*(sn.x-eta[i].x))/(2.0*sigmaX*sigmaX) + ((sn.y-eta[i].y)*(sn.y-eta[i].y))/(2.0*sigmaY*sigmaY)) );
    //double res = (1.0/(2.0*M_PI*sigmaX*sigmaY)) * exp(-0.5*(((sn.x-eta[i].x)*(sn.x-eta[i].x))/(sigmaX*sigmaX)+((sn.y-eta[i].y)*(sn.y-eta[i].y))/(sigmaY*sigmaY)));
    //return res;

    // Approximation delta function using L4 Lagrange interpolation
    //    double h3 = hx*hx*hx;
    //    double h32 = (1.0/(2.0*h3));
    //    double h36 = (1.0/(6.0*h3));
    //    double dh = fabs(sn.x - eta[i]);
    //    if (dh <= hx)                return ((2.0*hx-dh)*(hx-dh)*(hx+dh)) * h32;
    //    if (hx < dh && dh <= 2.0*hx) return ((2.0*hx-dh)*(hx-dh)*(3.0*hx-dh)) * h36;
    //    return 0.0;

    double res = 0.0;
    //if ( sn.x == eta[i].x && sn.y == eta[i].y ) res = 1.0/(hx*hy);
    if ( sn.i == eta[i].i && sn.j == eta[i].j ) res = 1.0/(hx*hy);

    return res;
}

double IProblem2Forward2D::g1(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double y = sn.y;
    double t = tn.t;
    //return lambda*(U(0.0, y, t) - theta);
    return lambda*(y*y+t - theta);
}

double IProblem2Forward2D::g2(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double y = sn.y;
    double t = tn.t;
    //return 2.0 + lambda*(U(1.0, y, t) - theta);
    return 2.0 + lambda*(1.0 + y*y + t - theta);
}

double IProblem2Forward2D::g3(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double x = sn.x;
    double t = tn.t;
    //return lambda*(U(x, 0.0, t) - theta);
    return lambda*(x*x + t - theta);
}

double IProblem2Forward2D::g4(const SpaceNodePDE &sn, const TimeNodePDE &tn) const
{
    double x = sn.x;
    double t = tn.t;
    //return 2.0 + lambda*(U(x, 1.0, t) - theta);
    return 2.0 + lambda*(1.0 + x*x + t - theta);
}

double IProblem2Forward2D::U(double x, double y, double t) const
{
    return x*x + y*y + t;
}

