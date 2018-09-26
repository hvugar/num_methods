#include "common.h"

SpacePointInfo::SpacePointInfo()
    : SpaceNodePDE(), id(0), layerNumber(0)
{
    _vl = NULL;
    _dx = NULL;
    _dy = NULL;
}

SpacePointInfo::SpacePointInfo(unsigned int id, unsigned int layerNumber)
    : SpaceNodePDE(), id(id), layerNumber(layerNumber)
{
    _vl = new double[layerNumber];
    _dx = new double[layerNumber];
    _dy = new double[layerNumber];
}

SpacePointInfo::~SpacePointInfo()
{
    delete [] _vl;  _vl = NULL;
    delete [] _dx;  _dx = NULL;
    delete [] _dy;  _dy = NULL;
}

void SpacePointInfo::createSpacePointInfos(unsigned int layerNumber)
{
    _vl = new double[layerNumber];
    _dx = new double[layerNumber];
    _dy = new double[layerNumber];
    this->layerNumber = layerNumber;
}

void SpacePointInfo::clearWeights()
{
    delete [] _vl;  _vl = NULL;
    delete [] _dx;  _dx = NULL;
    delete [] _dy;  _dy = NULL;
}
