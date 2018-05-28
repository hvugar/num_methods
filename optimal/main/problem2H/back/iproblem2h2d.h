#ifndef IPROBLEM2H2D_H
#define IPROBLEM2H2D_H

#include <grid/hibvp.h>
#include <vector>
#include <algorithm>
#include <gradient_cjt.h>
#include <gradient_sd.h>

#define METHOD_1
#define METHOD_3

#define _INFO_ROWS_ 1
#define _INFO_COLS_ 1

#define _DISTRIBUTION_METHOD_ 4

using namespace std;

namespace IProblem2H
{

struct OptimizeParameter
{
    DoubleMatrix k;
    DoubleMatrix z;
    std::vector<SpacePoint> xi;
    std::vector<SpacePoint> eta;
};

struct EquationParameter
{
    unsigned int No;
    unsigned int Nc;
    unsigned int Ns;
    std::vector<SpacePoint> theta;
    DoubleVector q;
    double a;
    double lambda;
};

struct ExtendedSpacePointNode
{
    SpacePoint pt;
    unsigned int id;
    unsigned int i;
    unsigned int j;
    double x;
    double y;
    double w;
};

struct ExtendedSpacePoint : public SpacePoint
{
public:
    unsigned int id;
private:
    unsigned int rows;
    unsigned int cols;
    DoubleMatrix w;
    unsigned int layerCount;
    DoubleVector uv;
    DoubleVector ux;
    DoubleVector uy;
    DoubleVector ri;
    DoubleVector ci;

    void spreadPoint(const Dimension &dimX, const Dimension &dimY, unsigned int type);
};

struct SpacePointInfo : public SpaceNodePDE
{
    unsigned int id;
    double *u;
    double *ux;
    double *uy;
    unsigned int layerNumber;

    void clearWeights();
    double value(unsigned int layer) const { return u[layer]; }
    double valueDx(unsigned int layer) const { return ux[layer]; }
    double valueDy(unsigned int layer) const { return uy[layer]; }
};

//class ExtendedSpaceNode2DH : public SpaceNodePDE
//{
//public:
//    //struct WISpaceNodePDE : public SpaceNodePDE
//    //{
//    //    double w;
//    //    double *u;
//    //};

//    explicit ExtendedSpaceNode2DH();
//    virtual ~ExtendedSpaceNode2DH();

//    void setSpaceNode(const SpacePoint& sn);

//    void extendWeights(const Dimension &dimX, const Dimension &dimY, unsigned int layerNumber, unsigned int rows, unsigned int cols);
//    void clearWeights();

//    double value(unsigned int layer) const;
//    //double value3(unsigned int ln) const;
//    //double value1(unsigned int layer) const;
//    //double value2(unsigned int layer) const;

//    double valueDx(unsigned int layer) const;
//    double valueDy(unsigned int layer) const;

//    //double valueDx1(unsigned int layer) const;
//    //double valueDy1(unsigned int layer) const;

//    //double valueDx2(unsigned int layer) const;
//    //double valueDy2(unsigned int layer) const;

//    unsigned int id;
//    //unsigned int rows;
//    //unsigned int cols;
//    //WISpaceNodePDE **wi;
//    unsigned int layerNumber;

//    double *u;
//    double *ux;
//    double *uy;
//};

class IProblem2H2D
{
public:
    static void Main(int argc, char* argv[]);

    static void forward();
    static void checkGradient();
    static void optimization1();
    static void forwardS();

    static void distributeDelta(const SpacePoint &pt, std::vector<ExtendedSpacePointNode> &nodes, unsigned int id, const Dimension &xd, const Dimension &yd, unsigned int k);
    static void distributeDeltaP(const SpacePoint &pt, unsigned int id, std::vector<ExtendedSpacePointNode> &nodes, const Dimension &dimX, const Dimension &dimY, unsigned int k=0);
    static void distributeDeltaR(const SpacePoint &pt, unsigned int id, std::vector<ExtendedSpacePointNode> &nodes, const Dimension &dimX, const Dimension &dimY, unsigned int k=0);
    static void distributeDeltaG(const SpacePoint &pt, unsigned int id, std::vector<ExtendedSpacePointNode> &nodes, const Dimension &dimX, const Dimension &dimY, unsigned int k);

    static void distributeDelta3(const SpacePoint &pt, unsigned int id, std::vector<ExtendedSpacePointNode> &nodes, const Dimension &dimX, const Dimension &dimY);
};

}

#endif // IPROBLEM2H2D_H
