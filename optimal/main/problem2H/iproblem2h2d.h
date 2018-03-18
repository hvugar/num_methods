#ifndef IPROBLEM2H2D_H
#define IPROBLEM2H2D_H

#include <grid/hibvp.h>
#include <vector>

using namespace std;

class IProblem2H2D
{
public:
    struct Parameter
    {
        unsigned int No;
        unsigned int Nc;
        unsigned int Ns;
        DoubleMatrix k;
        DoubleMatrix z;
        std::vector<SpacePoint> xi;
        std::vector<SpacePoint> eta;
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

    struct ObservationPointNode : public ExtendedSpacePointNode {};

    struct ObservationDeltaNode : public ExtendedSpacePointNode {};

    struct ControlPointNode : public ExtendedSpacePointNode {};

    struct ControlDeltaNode : public ExtendedSpacePointNode {};

    static void Main(int argc, char* argv[]);

    static void distributeDelta(const SpacePoint &pt, std::vector<ExtendedSpacePointNode> &nodes, unsigned int id,
                                const Dimension &xd, const Dimension &yd);

};

#endif // IPROBLEM2H2D_H
