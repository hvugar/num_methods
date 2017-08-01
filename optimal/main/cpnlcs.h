#ifndef CAUCHYPROBLEMNONLOCALCONTIONS_H
#define CAUCHYPROBLEMNONLOCALCONTIONS_H

#include <vector2d.h>
#include <matrix2d.h>
#include <vector>
#include <grid/cauchyp.h>

using namespace std;

class CauchyProblemNonLocalContions : public CauchyProblemM1stOrder, public SystemLinearODE1stOrder
{
public:
    static void Main(int agrc, char *argv[]);

    struct NonSeparatedCondition
    {
        double time;
        DoubleMatrix alpha;
        unsigned int rows;
        unsigned int cols;
        unsigned int n;
    };

    struct SeparatedCondition
    {
        double time;
        DoubleMatrix alpha;
        unsigned int rows;
        unsigned int cols;
        unsigned int n;
    };

    DoubleVector M;


    CauchyProblemNonLocalContions(const Dimension &grid);

    virtual void initialize();

    /**
     * @brief nlcs Неразделенные условия.
     */
    std::vector<NonSeparatedCondition> nscs;
    /**
     * @brief nlcs Разделенные условия.
     */
    std::vector<SeparatedCondition> scs;

    unsigned int n = 2;
    unsigned int L = 3;

    //std::vector<DoubleMatrix> alpha;
    DoubleVector betta;

    double x1(unsigned int k) const;
    double x2(unsigned int k) const;

    void calculate1();
    void calculateInterval(unsigned int start, unsigned int end, unsigned int r);

    virtual double S0(unsigned int i, double t, const DoubleVector &x, unsigned int k) const;
    \
protected:
    virtual double A(double t UNUSED_PARAM, unsigned int k, unsigned int row, unsigned int col) const;
    virtual double B(double t UNUSED_PARAM, unsigned int k, unsigned int row) const;

protected:
    virtual double f(double x, const DoubleVector &y, unsigned int k, unsigned int i) const;
};

#endif // CAUCHYPROBLEMNONLOCALCONTIONS_H
