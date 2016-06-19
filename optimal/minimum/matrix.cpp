#include "matrix.h"

CMatrix::CMatrix(size_t rows, size_t cols) : mrows(rows), mcols(cols)
{
    pvector = new CVector[rows];
    for (unsigned int j=0; j<rows; j++) pvector[j].resize(cols);
}

CMatrix::~CMatrix()
{
    for (size_t j=0; j<mrows; j++) pvector[j].clear();
    delete [] pvector;
}

size_t CMatrix::rows() const
{
    return mrows;
}

size_t CMatrix::columns() const
{
    return mcols;
}

void CMatrix::clear()
{
    size_t sz = rows();
    for (size_t j=0; j<sz; j++)
    {
        this[j].clear();
    }
    delete [] pvector;
    pvector = NULL;
}

void CMatrix::resize(size_t rows, size_t cols)
{
    clear();

    pvector = new CVector[rows];
    for (size_t j=0; j<rows; j++)
    {
        pvector[j].resize(cols);
    }
}

double CMatrix::min() const
{
    double _min = NAN;
    //    if (rows() != 0)
    //    {
    //        _min = at(0);
    //        for (unsigned int i=1; i<rows(); i++)
    //        {
    //            if (_min > (*this)(i)) _min = at(i);
    //        }
    //    }
    return _min;
}

double CMatrix::max() const
{
    double _max = NAN;
    //    if (rows()!=0)
    //    {
    //        _max = at(0);
    //        for (unsigned int i=1; i<size(); i++)
    //        {
    //            if (_max < at(i)) _max = at(i);
    //        }
    //    }
    return _max;
}

CVector& CMatrix::operator[](size_t i)
{
    return pvector[i];
}

const CVector& CMatrix::operator[](size_t i) const
{
    return pvector[i];
}

CVector::CVector(size_t size)
{
    msize = size;
    pdata = (double*)malloc(sizeof(double)*size);
}

CVector::~CVector()
{
    clear();
}

size_t CVector::size() const
{
    return msize;
}

void CVector::resize(size_t size)
{
    pdata = (double*)realloc(pdata, size);
}

double& CVector::operator[](unsigned int i)
{
    return pdata[i];
}

const double& CVector::operator[](unsigned int i) const
{
    return pdata[i];
}

void CVector::clear()
{
    free(pdata);
    pdata=NULL;
    msize=0;
}

double CVector::L2Norm() const
{
    double norm = 0.0;
    size_t sz = size();
    for (size_t i=0; i<sz; i++)
    {
        double item = operator [](i);
        norm += item*item;
    }
    return sqrt(norm);
}

double CVector::L1Norm() const
{
    double norm = 0.0;
    size_t sz = size();
    for (size_t i=0; i<sz; i++)
    {
        norm += fabs(operator [](i));
    }
    return norm;
}

double CVector::LInfNorm() const
{
    double norm = 0.0;
    size_t sz = size();
    for (unsigned int i=0; i<sz; i++)
    {
        double item = operator [](i);

        if (norm == NAN)
            norm = fabs(item);
        else
            if (norm < fabs(item)) norm = fabs(item);
    }
    return norm;
}

double CVector::EuclideanDistance(const CVector &v) const
{
    if ( size() != v.size() ) return INFINITY;

    double distance = 0.0;
    for (unsigned int i=0; i<size(); i++)
    {
        double item = operator [](i);
        distance += (item-v[i])*(item-v[i]);
    }
    return sqrt(distance);
}

void CVector::L2Normalize()
{
    double norm = L2Norm();
    if (norm != 0.0)
    {
        for (size_t i=0; i<size(); i++) (*this)[i] /= norm;
    }
}

void CVector::L1Normalize()
{
    double norm = L1Norm();
    for (size_t i=0; i<size(); i++) (*this)[i] /= norm;
}

void CVector::EuclideanNormalize()
{
    double norm = L2Norm();
    for (size_t i=0; i<size(); i++) (*this)[i] /= norm;
}

double CVector::min() const
{
    double _min = NAN;
    if (size() != 0)
    {
        _min = (*this)[0];
        for (size_t i=1; i<size(); i++)
        {
            if (_min > (*this)[i]) _min = (*this)[i];
        }
    }
    return _min;
}

double CVector::max() const
{
    double _max = NAN;
    if (size()!=0)
    {
        _max = (*this)[0];
        for (unsigned int i=1; i<size(); i++)
        {
            if (_max < (*this)[i]) _max = (*this)[i];
        }
    }
    return _max;
}

CVector CVector::mid(size_t start, size_t end) const
{
    size_t size = end - start + 1;
    CVector vector(size);
    for (size_t i=start; i<=end; i++) vector[i-start] = (*this)[i];
    return vector;
}

//double Vector::at(unsigned int i) const
//{
//    //if (i<msize) throw std::
//    return pdata[i];
//}

//double* Vector::data() const
//{
//    return this->pdata;
//}

//void Vector::add(double d)
//{
//    pdata = (double*)realloc(pdata, sizeof(double)*(msize+1));
//    pdata[msize] = d;
//    msize++;
//}

//void Vector::insert(unsigned int i, double d)
//{
//    pdata = (double*)realloc(pdata, sizeof(double)*(msize+1));
//    memcpy(pdata+i+1, pdata+i, sizeof(double)*(msize-i));
//    pdata[i]=d;
//    msize++;
//}

//void Vector::remove(unsigned int i)
//{
//    memcpy(pdata+i, pdata+i+1, sizeof(double)*(msize-i-1));
//    pdata = (double*)realloc(pdata, sizeof(double)*(msize-1));
//    msize--;
//}

