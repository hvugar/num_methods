#include "matrix.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>

Vector::Vector(unsigned int size, double value) : sz(size), pdata(NULL)
{
    if (size > 0) pdata = (double*)malloc(sizeof(double) * size);
    for (unsigned int i=0; i<size; i++) pdata[i] = value;
}

Vector::~Vector()
{
    clear();
}

double& Vector::at(unsigned int n)
{
    return pdata[n];
}

const double& Vector::at(unsigned int n) const
{
    return pdata[n];
}


void Vector::clear()
{
    if (pdata != NULL)
    {
        free(pdata);
        pdata = NULL;
        sz = 0;
    }
}

double* Vector::data() noexcept
{
    return pdata;
}

const double* Vector::data() const noexcept
{
    return pdata;
}

bool Vector::empty() const
{
    return sz == 0;
}

void Vector::resize (unsigned int size, double val)
{
    if (size == sz) return;
    if (size > sz)
    {
        pdata = (double*)realloc(pdata, sizeof(double)*size);
        for (unsigned int i=sz; i<size; i++) pdata[i] = val;
        sz = size;
        return;
    }
    if (size < sz)
    {
        pdata = (double*)realloc(pdata, sizeof(double)*size);
        sz = size;
        return;
    }
}

unsigned int Vector::size() const
{
    return sz;
}

Vector& Vector::operator= (const Vector& x)
{
    resize(x.sz, 0.0);
    return *this;
}

double& Vector::operator[] (unsigned int n)
{
    return pdata[n];
}

double Vector::operator[] (unsigned int n) const
{
    return pdata[n];
}

double Vector::L2Norm() const
{
    double norm = 0.0;
    if (!empty())
    {
        for (unsigned int i=0; i<sz; i++)
        {
            double item = pdata[i];
            norm += item*item;
        }
        return sqrt(norm);
    }
    return norm;
}

double Vector::L1Norm() const
{
    double norm = 0.0;
    if (!empty())
    {
        for (unsigned int i=0; i<sz; i++)
        {
            norm += fabs(pdata[i]);
        }
    }
    return norm;
}

double Vector::LInfNorm() const
{
    double norm = 0.0;
    if (!empty())
    {
        for (unsigned int i=0; i<sz; i++)
        {
            if ( i == 0 )
                norm = fabs(pdata[i]);
            else if (norm < fabs(pdata[i]))
                norm = fabs(pdata[i]);
        }
    }
    return norm;
}

double Vector::EuclideanNorm() const
{
    return L2Norm();
}

double Vector::EuclideanDistance(const Vector &p) const
{
    if ( size() != p.size() ) return INFINITY;

    double distance = 0.0;
    for (unsigned int i=0; i<size(); i++)
    {
        distance += (at(i)-p.at(i))*(at(i)-p.at(i));
    }
    return sqrt(distance);
}

void Vector::L2Normalize()
{
    double norm = L2Norm();
    if (norm != 0.0) for (unsigned int i=0; i<sz; i++) pdata[i] /= norm;
}

void Vector::L1Normalize()
{
    double norm = L1Norm();
    for (unsigned int i=0; i<sz; i++) pdata[i] /= norm;
}

void Vector::EuclideanNormalize()
{
    L2Normalize();
}

double Vector::min() const
{
    double _min = NAN;
    if (!empty())
    {
        _min = at(0);
        for (unsigned int i=1; i<sz; i++)
            if (_min > pdata[i]) _min = pdata[i];
    }
    return _min;
}

double Vector::max() const
{
    double _max = 0.0;
    if (!empty())
    {
        _max = pdata[0];
        for (unsigned int i=1; i<sz; i++)
            if (_max < pdata[i]) _max = pdata[i];
    }
    return _max;
}

Vector Vector::mid(unsigned int s, unsigned int e) const
{
    Vector vector(e-s+1);
    for (unsigned int i=s; i<=e; i++) vector[i] = (*this)[i];
    return vector;
}

Matrix::Matrix(unsigned int rows, unsigned int cols, double val) : mrows(rows), mcols(cols), pdata(NULL)
{
    if (mcols > 0 && mrows > 0)
    {
        pdata = new Vector[rows];
        for (unsigned int j=0; j<rows; j++)
        {
            pdata[j].resize(cols, val);
        }
    }
}

Matrix::~Matrix()
{
    for (unsigned int j=0; j<mrows; j++)
    {
        pdata[j].clear();
    }
    delete pdata;
}

Vector& Matrix::row(unsigned int r)
{
    return pdata[r];
}

const Vector& Matrix::row(unsigned int r) const
{
    return pdata[r];
}

unsigned int Matrix::rows() const { return mrows; }
unsigned int Matrix::cols() const { return mcols; }
unsigned int Matrix::size() const { return mrows; }

Vector& Matrix::operator[] (unsigned int j) { return pdata[j]; }
const Vector& Matrix::operator[] (unsigned int j) const { return pdata[j]; }

void Matrix::clear()
{

}

void Matrix::resize(unsigned int rows)
{
    if (rows > mrows)
    {
        Vector *cdata = pdata;
        pdata = new Vector[rows];
        for (unsigned int i=0; i<mrows; i++) pdata[i] = cdata[i];
        for (unsigned int i=mrows; i<rows; i++) pdata[i].resize(mcols);
        delete [] cdata;
    }
    if (rows < mrows)
    {
        Vector *cdata = pdata;
        pdata = new Vector[rows];
        for (unsigned int i=0; i<rows; i++) pdata[i] = cdata[i];
        delete [] cdata;
    }
    mrows = rows;
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

