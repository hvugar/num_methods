#ifndef IPROBLEM2H2D_H
#define IPROBLEM2H2D_H

#include <grid/hibvp.h>
#include <vector>

using namespace std;

class IProblem2H2D
{
public:
    struct OptimizeParameter
    {
//        unsigned int No;
//        unsigned int Nc;
//        unsigned int Ns;
        DoubleMatrix k;
        DoubleMatrix z;
        std::vector<SpacePoint> xi;
        std::vector<SpacePoint> eta;
//        std::vector<SpacePoint> theta;
//        DoubleVector q;

//        double a;
//        double lambda;
//        double lambda1;
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
        double lambda1;
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

    static void Main(int argc, char* argv[]);

    static void Forward();

    static void distributeDelta(const SpacePoint &pt, std::vector<ExtendedSpacePointNode> &nodes, unsigned int id,
                                const Dimension &xd, const Dimension &yd);
};

class ExtendedSpaceNode2DH : public SpaceNodePDE
{
public:
    struct WISpaceNodePDE : public SpaceNodePDE
    {
        double w;
        double *u;
    };

    explicit ExtendedSpaceNode2DH();
    virtual ~ExtendedSpaceNode2DH();

    void setSpaceNode(const SpacePoint& sn);

    void extendWeights(const Dimension &dimX, const Dimension &dimY, unsigned int rows = 4, unsigned int cols = 4);
    void clearWeights();

    void extendLayers(unsigned int layerNumber);
    void clearLayers();

    double value(unsigned int layer) const;
    double valueDx(unsigned int layer) const;
    double valueDy(unsigned int layer) const;

    unsigned int id;
    unsigned int rows;
    unsigned int cols;
    WISpaceNodePDE **wi;
    unsigned int layerNumber;

private:
    double value(double x, double y, unsigned int layer) const;
};

#endif // IPROBLEM2H2D_H
