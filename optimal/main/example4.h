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
    unsigned int N;
    unsigned int K;
    unsigned int n;
//    unsigned int F;

    void static Main(int argc, char *argv[]);

    double fx1(unsigned int) const;
    double fx2(unsigned int) const;
    double fx3(unsigned int) const;

    void calculateM1();
    void calculateM1BE(unsigned int c, unsigned int s[], unsigned int L, const std::vector<DoubleVector> &rx, DoubleMatrix &M, DoubleVector &B,
                       const std::vector<DoubleMatrix> &P3, const std::vector<DoubleMatrix> &P2, const std::vector<DoubleMatrix> &P1, const std::vector<DoubleMatrix> &P0,
                       const std::vector<DoubleVector> &Q);

    void calculateM2();
    void calculateM2BE(unsigned int c, unsigned int s[], unsigned int L, const std::vector<DoubleVector> &rx, DoubleMatrix &M, DoubleVector &B);

    double a(unsigned int i, unsigned int j, unsigned int k) const;
    double b(unsigned int i, unsigned int k) const;

    void calculatePQ(std::vector<DoubleMatrix> &P3, std::vector<DoubleMatrix> &P2, std::vector<DoubleMatrix> &P1, std::vector<DoubleMatrix> &P0, std::vector<DoubleVector> &Q);
    void calculateRS(std::vector<DoubleVector> &rx);
    void calculateNS(std::vector<DoubleVector> &nx, const std::vector<DoubleMatrix> &P3, const std::vector<DoubleMatrix> &P2, const std::vector<DoubleMatrix> &P1,
                     const std::vector<DoubleMatrix> &P0, const std::vector<DoubleVector> &Q);

    unsigned int w = 12;
    unsigned int p = 8;
};

#endif // EXAMPLE4_H
