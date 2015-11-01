#include "doublevector.h"
#include <math.h>

DoubleVector::DoubleVector(unsigned int size, double value) : std::vector<double>(size, value)
{
}

DoubleVector::~DoubleVector()
{
    clear();
}

double DoubleVector::EuclideanNorm() const
{
    return L2Norm();
}

double DoubleVector::L2Norm() const
{
    double norm = 0.0;
    for (unsigned int i=0; i<size(); i++)
    {
        norm += at(i)*at(i);
    }
    return sqrt(norm);
}

double DoubleVector::L1Norm() const
{
    double norm = 0.0;
    for (unsigned int i=0; i<size(); i++)
    {
        norm += fabs(at(i));
    }
    return norm;
}

double DoubleVector::LInfNorm() const
{
    double norm = NAN;
    for (unsigned int i=0; i<size(); i++)
    {
        if (norm == NAN)
            norm = fabs(at(i));
        else
            if (norm < fabs(at(i))) norm = fabs(at(i));
    }
    return norm;
}

double DoubleVector::EuclideanDistance(const DoubleVector &p) const
{
    if ( size() != p.size() ) return INFINITY;

    double distance = 0.0;
    for (unsigned int i=0; i<size(); i++)
    {
        distance += (at(i)-p.at(i))*(at(i)-p.at(i));
    }
    return sqrt(distance);
}

void DoubleVector::L2Normalize()
{
    double norm = L2Norm();
    for (unsigned int i=0; i<size(); i++) (*this)[i] /= norm;
}

void DoubleVector::L1Normalize()
{
    double norm = L1Norm();
    for (unsigned int i=0; i<size(); i++) (*this)[i] /= norm;
}

void DoubleVector::EuclideanNormalize()
{
    double norm = L2Norm();
    for (unsigned int i=0; i<size(); i++) (*this)[i] /= norm;
}

double DoubleVector::L2Norm(const DoubleVector& p)
{
    return p.L2Norm();
}

double DoubleVector::L1Norm(const DoubleVector& p)
{
    return p.L1Norm();
}

double DoubleVector::LInfNorm(const DoubleVector& p)
{
    return p.LInfNorm();
}

double DoubleVector::EuclideanNorm(const DoubleVector& p)
{
    return L2Norm(p);
}

double DoubleVector::EuclideanDistance(const DoubleVector& p, const DoubleVector& q)
{
    return p.EuclideanDistance(q);
}

void DoubleVector::L2Normalize(DoubleVector& p)
{
    p.L2Normalize();
}

void DoubleVector::L1Normalize(DoubleVector& p)
{

    p.L1Normalize();
}

void DoubleVector::EuclideanNormalize(DoubleVector& p)
{
    p.EuclideanNormalize();
}

/**
 * @brief DoubleMatrix::DoubleMatrix
 * @param m Number of rows
 * @param n Number of columns in row
 * @param value
 */
DoubleMatrix::DoubleMatrix(unsigned int m, unsigned int n, double value) : std::vector<DoubleVector>()
{
    resize(m);
    for (unsigned int j=0; j<m; j++) this[j].resize(n);
}

DoubleMatrix::~DoubleMatrix()
{}

void DoubleMatrix::Clear()
{
    for (unsigned int j=0; j<size(); j++)
    {
        this[j].clear();
    }
    this->clear();
}

/**
 * @brief DoubleMatrix::Resize
 * @param m Number of rows
 * @param n Number of columns in row
 */
void DoubleMatrix::Resize(unsigned int m, unsigned n)
{
    Clear();

    resize(m);
    for (unsigned int j=0; j<size(); j++)
    {
        this[j].resize(n);
    }
}

DoubleCube::DoubleCube() : std::vector<DoubleMatrix>()
{}

DoubleCube::~DoubleCube()
{}

void DoubleCube::Resize(unsigned int Nz, unsigned int Ny, unsigned Nx)
{
    Clear();

    resize(Nz);
    for (unsigned int k=0; k<size(); k++)
    {
        this[k].resize(Ny);
        for (unsigned int m=0; m<this[k].size(); m++)
        {
            this[k][m].resize(Nx);
        }
    }
}

void DoubleCube::Clear()
{
    for (unsigned int k=0; k<size(); k++)
    {
        for (unsigned int m=0; m<this[k].size(); m++)
        {
            this[k][m].clear();
        }
        this[k].clear();
    }
    this->clear();
}
