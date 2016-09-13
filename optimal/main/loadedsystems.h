#ifndef LOADEDSYSTEMS_H
#define LOADEDSYSTEMS_H

#include <vector2d.h>
#include <matrix2d.h>
#include <matrix3d.h>
#include <rungekutta.h>
#include <printer.h>
#include <vector>

using namespace std;

class LoadedSystems
{
public:
    LoadedSystems();
    virtual ~LoadedSystems() {}

    void init();
    void calculate(unsigned int j, DoubleMatrix &m, unsigned int N,
                   unsigned int k2, const std::vector<DoubleMatrix> &alpha,
                   unsigned int k1, const std::vector<DoubleMatrix> &betta,
                   unsigned int n,  const DoubleVector &qamma);

    void calculate(DoubleMatrix &m, unsigned int N,
                   const DoubleVector &alpha,
                   const std::vector<DoubleVector> &betta,
                   double qamma, unsigned int n);

    double A(unsigned int row, unsigned int col, double t) const;
    double B(unsigned int row, unsigned int col, unsigned int s, double t) const;
    double C(unsigned int row, double t) const;

    double x1(double t) const;
    double x2(double t) const;
    double x3(double t) const;

    double R(const DoubleVector &args) const;
    double S(const DoubleVector &args, double t) const;

    DoubleVector t1;
    DoubleVector t2;

    // loaded points count
    unsigned int k1;
    // border points count
    unsigned int k2;
    unsigned int n;

    DoubleVector qamma;
    vector<DoubleMatrix> alpha;
    vector<DoubleMatrix> betta;

    double ht;
    unsigned int N;
};

#endif // LOADEDSYSTEMS_H
