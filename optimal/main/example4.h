#ifndef EXAMPLE4_H
#define EXAMPLE4_H

#include <vector2d.h>
#include <matrix2d.h>
#include <matrix3d.h>
#include <printer.h>
#include <math.h>
#include <vector>

#define SAMPLE_1

typedef std::vector<DoubleMatrix> stdDoubleMatrixVector;
typedef std::vector<DoubleVector> stdDoubleVectorVector;

class Example4
{
public:
    double h;
    unsigned int N;
    unsigned int K;
    unsigned int n;
    unsigned int F;
    unsigned int w;
    unsigned int p;

    void static Main(int argc, char *argv[]);

    /**
     * @brief fx
     * @param i index of function
     * @param k index of time
     * @return
     */
    double fx(unsigned int i, unsigned int k) const;

    void calculateRX(DoubleMatrix &rx);
    void calculateNX(const DoubleMatrix &rx, DoubleMatrix &nx);
    void calculateRS(stdDoubleMatrixVector &rx);

    void initAMatrices(stdDoubleMatrixVector &A);
    void updateAMatrices(std::vector<DoubleMatrix> &A, unsigned int k);
    void clearAMatrices(stdDoubleMatrixVector &A);

    void calculateNX(const stdDoubleMatrixVector &rx, DoubleVector &x1, DoubleVector &x2, DoubleVector &x3, stdDoubleVectorVector &nx);

    void calculateM1(const std::vector<unsigned int> *s, const DoubleMatrix &rx, const DoubleMatrix &nx);
    void calculateM1BE(unsigned int c, const std::vector<unsigned int> s, const DoubleMatrix &rx, const DoubleMatrix &nx, DoubleMatrix &M, DoubleVector &B);
    void qovmaM1R2L(const std::vector<DoubleMatrix> &GAMMA, DoubleVector &ETA, const std::vector<unsigned int> &s, std::vector<DoubleMatrix> &BETTA);

    void calculateM2(const std::vector<unsigned int> *s, const DoubleMatrix &rx UNUSED_PARAM, const DoubleMatrix &nx UNUSED_PARAM);
    void calculateM2BE(unsigned int c, const std::vector<unsigned int> s, const DoubleMatrix &nx, DoubleMatrix &M, DoubleVector &B, const std::vector<stdDoubleMatrixVector> &P, const stdDoubleMatrixVector &Q);

    double a(unsigned int i, unsigned int j, unsigned int k) const;
    double b(unsigned int i, unsigned int k) const;

    void calculatePQ(std::vector<stdDoubleMatrixVector> &P, stdDoubleMatrixVector &Q);
    void calculateNS(DoubleMatrix &nx, const DoubleMatrix &rx, const std::vector<stdDoubleMatrixVector> &P, const stdDoubleMatrixVector &Q);
    void fillGamma(stdDoubleMatrixVector &GAMMA, DoubleVector &ETA, unsigned int s, unsigned int k);
};

#endif // EXAMPLE4_H
