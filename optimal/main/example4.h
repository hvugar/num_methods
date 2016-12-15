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

class Example4
{
public:
    Example4();

    double h;
    unsigned int N = 1000;
    unsigned int K = 4;
    unsigned int n = 3;
    unsigned int F = N/10;
    unsigned int L = 5;
    unsigned int w = 12;
    unsigned int p = 8;

    void static Main(int argc, char *argv[]);

    double fx1(unsigned int) const;
    double fx2(unsigned int) const;
    double fx3(unsigned int) const;

    void calculateRX(DoubleMatrix &rx);
    void calculateNX(const DoubleMatrix &rx, DoubleMatrix &nx);

    void calculateRS(std::vector<DoubleVector> &rx);

    void initAMatrices(std::vector<DoubleMatrix> &A);
    void updateAMatrices(std::vector<DoubleMatrix> &A, unsigned int k);
    void clearAMatrices(std::vector<DoubleMatrix> &A);
    void calculateNX(const std::vector<DoubleVector> &rx, DoubleVector &x1, DoubleVector &x2, DoubleVector &x3, std::vector<DoubleVector> &nx);

    void calculateM1(const unsigned int s[][5], const DoubleMatrix &rx, const DoubleMatrix &nx);
    void calculateM1BE(unsigned int c, const unsigned int s[], unsigned int L, const DoubleMatrix &nx,
                       DoubleMatrix &M, DoubleVector &B);

    void calculateM2(const unsigned int s[][5], const DoubleMatrix &rx UNUSED_PARAM, const DoubleMatrix &nx UNUSED_PARAM);
    void calculateM2BE(unsigned int c, const unsigned int s[], unsigned int L, const DoubleMatrix &nx, DoubleMatrix &M, DoubleVector &B,
                       const std::vector<DoubleMatrix> &P3, const std::vector<DoubleMatrix> &P2, const std::vector<DoubleMatrix> &P1, const std::vector<DoubleMatrix> &P0,
                       const std::vector<DoubleVector> &Q);

    double a(unsigned int i, unsigned int j, unsigned int k) const;
    double b(unsigned int i, unsigned int k) const;

    void calculatePQ(std::vector<DoubleMatrix> &P3, std::vector<DoubleMatrix> &P2, std::vector<DoubleMatrix> &P1, std::vector<DoubleMatrix> &P0, std::vector<DoubleVector> &Q);

    void calculateNS(DoubleMatrix &nx, const DoubleMatrix &rx,
                     const std::vector<DoubleMatrix> &P3, const std::vector<DoubleMatrix> &P2, const std::vector<DoubleMatrix> &P1,
                     const std::vector<DoubleMatrix> &P0, const std::vector<DoubleVector> &Q);

    void fillGamma(std::vector<DoubleMatrix> &GAMMA, DoubleVector &ETA, unsigned int s, unsigned int k);
};

#endif // EXAMPLE4_H
