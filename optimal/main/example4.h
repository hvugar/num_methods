#ifndef EXAMPLE4_H
#define EXAMPLE4_H

#include <vector2d.h>
#include <matrix2d.h>
#include <matrix3d.h>
#include <printer.h>
#include <math.h>
#include <vector>

#define SAMPLE_1
//#define SAMPLE_2

typedef std::vector<DoubleMatrix> stdDoubleMatrixVector;
typedef std::vector<DoubleVector> stdDoubleVectorVector;

class Example4
{
public:
    Example4();

    double h;
    unsigned int N = 1000;
    unsigned int K = 4;
    unsigned int n = 3;
    unsigned int F = N/10;
    //unsigned int L = 5;
    unsigned int w = 12;
    unsigned int p = 8;

    void static Main(int argc, char *argv[]);

    double fx(unsigned int n, unsigned int i) const;
//    double fx1(unsigned int) const;
//    double fx2(unsigned int) const;
//    double fx3(unsigned int) const;

    void calculateRX(DoubleMatrix &rx);
    void calculateNX(const DoubleMatrix &rx, DoubleMatrix &nx);

    void calculateRS(std::vector<DoubleVector> &rx);

    void initAMatrices(stdDoubleMatrixVector &A);
    void updateAMatrices(stdDoubleMatrixVector &A, unsigned int k);
    void clearAMatrices(stdDoubleMatrixVector &A);

    void calculateNX(const std::vector<DoubleVector> &rx, DoubleVector &x1, DoubleVector &x2, DoubleVector &x3, std::vector<DoubleVector> &nx);

    void calculateM1(const std::vector<unsigned int> *s, const DoubleMatrix &rx, const DoubleMatrix &nx);
    void calculateM1BE(unsigned int c, const std::vector<unsigned int> s, const DoubleMatrix &nx,
                       DoubleMatrix &M, DoubleVector &B);
    void qovmaM1R2L(const std::vector<DoubleMatrix> &GAMMA, DoubleVector &ETA, const std::vector<unsigned int> &s, std::vector<DoubleMatrix> &BETTA);

    void calculateM2(const std::vector<unsigned int> *s, const DoubleMatrix &rx UNUSED_PARAM, const DoubleMatrix &nx UNUSED_PARAM);
//    void calculateM2BE(unsigned int c, const std::vector<unsigned int> s, const DoubleMatrix &nx, DoubleMatrix &M, DoubleVector &B,
//                       const stdDoubleMatrixVector &P3, const stdDoubleMatrixVector &P2, const stdDoubleMatrixVector &P1, const stdDoubleMatrixVector &P0,
//                       const stdDoubleMatrixVector &Q);
    void calculateM2BE(unsigned int c, const std::vector<unsigned int> s, const DoubleMatrix &nx, DoubleMatrix &M, DoubleVector &B,
                       const std::vector<stdDoubleMatrixVector> &P, const stdDoubleMatrixVector &Q);

    double a(unsigned int i, unsigned int j, unsigned int k) const;
    double b(unsigned int i, unsigned int k) const;

    //void calculatePQ(std::vector<stdDoubleMatrixVector> &P3, stdDoubleMatrixVector &P2, stdDoubleMatrixVector &P1, stdDoubleMatrixVector &P0, stdDoubleMatrixVector &Q);
    void calculatePQ(std::vector<stdDoubleMatrixVector> &P, stdDoubleMatrixVector &Q);

//    void calculateNS(DoubleMatrix &nx, const DoubleMatrix &rx,
//                     const stdDoubleMatrixVector &P3, const stdDoubleMatrixVector &P2, const stdDoubleMatrixVector &P1,
//                     const stdDoubleMatrixVector &P0, const stdDoubleMatrixVector &Q);
    void calculateNS(DoubleMatrix &nx, const DoubleMatrix &rx, const std::vector<stdDoubleMatrixVector> &P, const stdDoubleMatrixVector &Q);

    void fillGamma(stdDoubleMatrixVector &GAMMA, DoubleVector &ETA, unsigned int s, unsigned int k);
};

#endif // EXAMPLE4_H
