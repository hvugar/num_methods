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
    if (norm != 0.0)
    {
        for (unsigned int i=0; i<size(); i++) (*this)[i] /= norm;
    }
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

double DoubleVector::min() const
{
    double _min = NAN;
    if (size() != 0)
    {
        _min = at(0);
        for (unsigned int i=1; i<size(); i++)
        {
            if (_min > at(i)) _min = at(i);
        }
    }
    return _min;
}

double DoubleVector::max() const
{
    double _max = NAN;
    if (size()!=0)
    {
        _max = at(0);
        for (unsigned int i=1; i<size(); i++)
        {
            if (_max < at(i)) _max = at(i);
        }
    }
    return _max;
}

DoubleVector* DoubleVector::mid(unsigned int start, unsigned int end) const
{
    unsigned int size = end - start + 1;
    DoubleVector *vector = new DoubleVector(size);
    for (unsigned int i=start; i<=end; i++) (*vector)[i-start] = (*this)[i];
    return vector;
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
    for (unsigned int j=0; j<m; j++) this->at(j).resize(n, value);
}

DoubleMatrix::~DoubleMatrix()
{}

double DoubleMatrix::min() const
{
    double _min = NAN;
    if (size() != 0)
    {
        _min = at(0).min();
        for (unsigned int i=1; i<size(); i++)
        {
            double _vmin = at(i).min();
            if (_min > _vmin) _min = _vmin;
        }
    }
    return _min;
}

double DoubleMatrix::max() const
{
    double _max = NAN;
    if (size() != 0)
    {
        _max = at(0).max();
        for (unsigned int i=1; i<size(); i++)
        {
            double _vmax = at(i).max();
            if (_max < _vmax) _max = _vmax;
        }
    }
    return _max;
}

double DoubleMatrix::determinant() const
{
    double det = 0.0;

    //    if (m->rows != m->columns)
    //        return NAN;

    //    if (m->rows == 1 && m->columns == 1)
    //        return m->items[0][0];
    //    else if (m->rows == 2 && m->columns == 2)
    //        return m->items[0][0] * m->items[1][1] - m->items[0][1] * m->items[1][0];
    //    else
    //    {
    //        int i;
    //        for (i=0; i<m->columns; i++)
    //        {
    //            struct Matrix *mnr = matrix_minor(m, 0, i);
    //            det += (i%2==0 ? +1 : -1) * m->items[0][i] * matrix_det(mnr);
    //            matrix_free(mnr);
    //        }
    //    }

    return det;
}

DoubleMatrix* DoubleMatrix::transpose() const
{
    return NULL;
}

DoubleMatrix* DoubleMatrix::inverse() const
{
    return NULL;
}

DoubleMatrix* DoubleMatrix::minor(size_t row, size_t col) const
{
    C_UNUSED(row);
    C_UNUSED(col);
    return NULL;
}

DoubleMatrix* DoubleMatrix::multiply(const DoubleMatrix &m) const
{
    C_UNUSED(m);
    return NULL;
}

DoubleMatrix* DoubleMatrix::multiply(const DoubleVector &v) const
{
    C_UNUSED(v);
    return NULL;
}

DoubleMatrix DoubleMatrix::operator+(const DoubleMatrix &A) const
{
    C_UNUSED(A);
    return (*this);
}

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

Vector::Vector(unsigned int size, double value) : sz(size), pdata(NULL)
{
    if (size > 0) pdata = (double*)malloc(sizeof(double) * size);
    for (unsigned int i=0; i<size; i++) pdata[i] = value;
}

Vector::~Vector()
{
    free(pdata);
}

double& Vector::at(unsigned int n)
{
    //    if (n < 0) return 0.0;
    //    if (n > sz-1) return 0.0;
    return pdata[n];
}

const double& Vector::at(unsigned int n) const
{
    return pdata[n];
}


void Vector::clear()
{
    if (pdata != NULL) free(pdata);
    sz = 0;
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

double Vector::EuclideanDistance(const DoubleVector &p) const
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

Matrix::Matrix(unsigned int rows, unsigned int cols, double val) : rows(rows), cols(cols), pdata(NULL)
{
    if (cols > 0 && rows > 0)
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
    for (unsigned int j=0; j<rows; j++)
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
