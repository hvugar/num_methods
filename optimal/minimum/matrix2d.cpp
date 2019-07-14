#include "matrix2d.h"
#include "printer.h"
#include "exceptions.h"
#include <time.h>
#include <math.h>
#include <float.h>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

DoubleMatrix::DoubleMatrix(unsigned int rows, unsigned int cols, double value) : mRows(0), mCols(0), mData(nullptr)
{
    if (rows > 0 && cols > 0)
    {
        mRows = rows;
        mCols = cols;
        mData = static_cast<double**>(malloc(sizeof(double*)*rows));
        for (unsigned int i=0; i<rows; i++)
        {
            mData[i] = static_cast<double*>(malloc(sizeof(double)*cols));
            for (unsigned int j=0; j<cols; j++) mData[i][j] = value;
        }
    }
}

DoubleMatrix::DoubleMatrix(const DoubleMatrix &matrix) : mRows(0), mCols(0), mData(nullptr)
{
    //puts("DoubleMatrix::DoubleMatrix(const DoubleMatrix &matrix)");
    if (matrix.mRows > 0 && matrix.mCols > 0)
    {
        mRows = matrix.mRows;
        mCols = matrix.mCols;
        mData = static_cast<double**>(malloc(sizeof(double*)*mRows));
        for (unsigned int i=0; i<mRows; i++)
        {
            mData[i] = static_cast<double*>(malloc(sizeof(double)*mCols));
            memcpy(mData[i], matrix.mData[i], sizeof(double)*mCols);
        }
    }
}

DoubleMatrix::DoubleMatrix(const DoubleVector &vector) : mRows(0), mCols(0), mData(nullptr)
{
    if (vector.length() > 0)
    {
        mRows = vector.length();
        mCols = 1;
        mData = static_cast<double**>(malloc(sizeof(double*)*mRows));
        for (unsigned int i=0; i<mRows; i++)
        {
            mData[i] = static_cast<double*>(malloc(sizeof(double)*mCols));
            mData[i][0] = vector.at(i);
        }
    }
}

DoubleMatrix::~DoubleMatrix()
{
    clear();
}

unsigned int DoubleMatrix::rows() const
{
    return mRows;
}

unsigned int DoubleMatrix::cols() const
{
    return mCols;
}

bool DoubleMatrix::empty() const
{
    return (mRows == 0 || mCols == 0 || mData == nullptr);
}

void DoubleMatrix::clear()
{
    if (mData == nullptr) return;

    for (unsigned int row=0; row<mRows; row++)
    {
        free(mData[row]);
        mData[row] = nullptr;
    }

    free(mData);
    mData = nullptr;
    mRows = 0;
    mCols = 0;
}

void DoubleMatrix::resize(unsigned int rows, unsigned int cols, double value)
{
    if (rows == 0 && cols == 0) clear();

    if (rows == 0 || cols == 0) return;

    if (rows > 0 && cols > 0)
    {
        if (mData == nullptr)
        {
            mData = static_cast<double**>(malloc(sizeof(double*)*rows));
            for (unsigned int row=0; row<rows; row++)
            {
                mData[row] = static_cast<double*>(malloc(sizeof(double)*cols));
                for (unsigned int col=0; col<cols; col++) mData[row][col] = value;
            }
            mRows = rows;
            mCols = cols;
        }
        else
        {
            double **ptr = static_cast<double**>(malloc(sizeof(double*) * rows));
            for (unsigned int row=0; row<rows; row++)
            {
                ptr[row] = static_cast<double*>(malloc(sizeof(double)*cols));

                if (row < mRows)
                {
                    for (unsigned int col=0; col<cols; col++)
                        if (col < mCols) ptr[row][col] = mData[row][col]; else ptr[row][col] = value;
                }
                else
                {
                    for (unsigned int col=0; col<cols; col++)
                        ptr[row][col] = value;
                }
            }

            clear();

            mData = ptr;
            mRows = rows;
            mCols = cols;
        }
    }
}

void DoubleMatrix::reset(double value)
{
    for (unsigned int r=0; r<mRows; r++)
    {
        for (unsigned int c=0; c<mCols; c++)
        {
            mData[r][c] = value;
        }
    }
}

double& DoubleMatrix::at(unsigned int row, unsigned int col)
{
    //if (row>=mRows) { throw std::out_of_range("Error: Row index out of range"); }
    //if (col>=mCols) { throw std::out_of_range("Error: Column index out of range"); }
    if (row>=mRows) { throw double_matrix_exception(5); }
    if (col>=mCols) { throw double_matrix_exception(6); }
    return mData[row][col];
}

const double& DoubleMatrix::at(unsigned int row, unsigned int col) const
{
    //if (row>=mRows) { throw std::out_of_range("Error: Row index out of range"); }
    //if (col>=mCols) { throw std::out_of_range("Error: Column index out of range"); }
    if (row>=mRows) { throw double_matrix_exception(5); }
    if (col>=mCols) { throw double_matrix_exception(6); }
    return mData[row][col];
}

DoubleVector DoubleMatrix::row(unsigned int r) const
{
    return DoubleVector(mData[r], mCols);
}

DoubleVector DoubleMatrix::col(unsigned int c) const
{
    DoubleVector v(mRows);
    for (unsigned int r=0; r<mRows; r++) v[r] = mData[r][c];
    return v;
}

void DoubleMatrix::setColumn(unsigned int c, const DoubleVector &col)
{
    if (mRows == col.length() && c < mCols)
    {
        for (unsigned int r=0; r<mRows; r++) mData[r][c] = col.at(r);
    }
    else
    {
        throw double_matrix_exception(0);
    }
}

void DoubleMatrix::setRow(unsigned int r, const DoubleVector &row)
{
    if (mCols == row.length() && r < mRows)
    {
        for (unsigned int c=0; c<mCols; c++) mData[r][c] = row.at(c);
    }
    else
    {
        throw double_matrix_exception(0);
    }
}

double** DoubleMatrix::data() const
{
    return mData;
}

double** DoubleMatrix::data()
{
    return mData;
}

///////////////////////////////////////////////////////////////////////////////////////////////

DoubleMatrix DoubleMatrix::DiagonalMatrix(const DoubleVector& v)
{
    unsigned int n = v.length();
    DoubleMatrix m(n, n, 0.0);
    for (unsigned int i=0; i<n; i++) m.mData[i][i] = v[i];
    return m;
}

DoubleMatrix DoubleMatrix::IdentityMatrix(unsigned int n)
{
    DoubleMatrix m(n, n, 0.0);
    for (unsigned int i=0; i<n; i++) m.mData[i][i] = 1.0;
    return m;
}

DoubleMatrix DoubleMatrix::ZeroMatrix(unsigned int rows, unsigned int cols)
{
    DoubleMatrix matrix(rows, cols, 0.0);
    return matrix;
}

DoubleMatrix DoubleMatrix::HilbertMatrix(unsigned int rows, unsigned int cols)
{
    DoubleMatrix matrix(rows, cols);
    for (unsigned int row=0; row<rows; row++)
    {
        for (unsigned int col=0; col<cols; col++)
        {
            matrix.mData[row][col] = 1.0/static_cast<double>(row+col+1);
        }
    }
    return matrix;
}

///////////////////////////////////////////////////////////////////////////////////////////////

bool DoubleMatrix::dioqonalMatrix() const
{
    for (unsigned int row=0; row < mRows; row++)
    {
        for (unsigned int col=0; col < mCols; col++)
        {
            if (row != col && fabs(mData[row][col]) >= DBL_EPSILON) return false;
        }
    }
    return true;
}

bool DoubleMatrix::identityMatrix() const
{
    for (unsigned int row=0; row<mRows; row++)
    {
        for (unsigned int col=0; col<mCols; col++)
        {
            if (row != col && fabs(mData[row][col]) >= DBL_EPSILON) return false;
            if (row == col && fabs(mData[row][col]) != 1.0) return false;
        }
    }
    return true;
}

bool DoubleMatrix::zeroMatrix() const
{
    for (unsigned int row=0; row < mRows; row++)
    {
        for (unsigned int col=0; col < mCols; col++)
        {
            if (fabs(mData[row][col]) >= DBL_EPSILON) return false;
        }
    }
    return true;
}

bool DoubleMatrix::squareMatrix() const
{
    return mRows == mCols;
}

///////////////////////////////////////////////////////////////////////////////////////////////

bool DoubleMatrix::dimEquals(const DoubleMatrix &matrix) const
{
    return (mRows == matrix.mRows && mCols == matrix.mCols);
}

bool DoubleMatrix::equals(const DoubleMatrix &matrix) const
{
    if (!dimEquals(matrix)) return false;

    bool equals = true;
    for (unsigned int i=0; i<mRows; i++)
    {
        for (unsigned int j=0; j<mCols; j++)
        {
            if (fabs(mData[i][j] - matrix.mData[i][j]) >= DBL_EPSILON)
            {
                equals = false;
                break;
            }
        }
        if (!equals) break;
    }
    return equals;
}

double DoubleMatrix::min() const
{
    double minimum = DBL_MAX;
    for (unsigned int row=0; row<mRows; row++)
    {
        for (unsigned int col=0; col<mCols; col++)
        {
            if (mData[row][col] < minimum) minimum = mData[row][col];
        }
    }
    return minimum;
}

double DoubleMatrix::max() const
{
    double maximum = DBL_MIN;
    for (unsigned int row=0; row<mRows; row++)
    {
        for (unsigned int col=0; col<mCols; col++)
        {
            if (mData[row][col] > maximum) maximum = mData[row][col];
        }
    }
    return maximum;
}

double DoubleMatrix::determinant2() const
{
    double det = 0.0;

    if (mRows != mCols) { throw double_matrix_exception(2); }

    if (mRows == 1 && mCols == 1)
        return mData[0][0];
    else if (mRows == 2 && mCols == 2)
        return mData[0][0]*mData[1][1] - mData[0][1]*mData[1][0];
    else
    {
        unsigned int i;
        for (i=0; i<mCols; i++)
        {
            DoubleMatrix mnr = minor(0, i);
            det += (i%2==0 ? +1 : -1) * mData[0][i] * mnr.determinant2();
        }
    }

    return det;
}

double DoubleMatrix::determinant() const
{
    if (mRows != mCols) { throw double_matrix_exception(2); }

    // Checking properties of the determinant ---------------------------------------------------------------
    for (unsigned int r=0; r<mRows; r++)
    {
        bool row_items_are_zero = true;
        for (unsigned int c=0; c<mCols; c++) if (fabs(mData[r][c]) >= DBL_EPSILON) row_items_are_zero = false;
        if (row_items_are_zero) return 0.0;
    }
    for (unsigned int c=0; c<mRows; c++)
    {
        bool col_items_are_zero = true;
        for (unsigned int r=0; r<mCols; r++) if (fabs(mData[r][c]) >= DBL_EPSILON) col_items_are_zero = false;
        if (col_items_are_zero) return 0.0;
    }

    if (identityMatrix()) return 1.0;

    // ------------------------------------------------------------------------------------------------------

    double det = 0.0;

    // Gaussian elimination ---------------------------------------------------------------------------------
    DoubleMatrix mx = *this;

    for (unsigned k=0; k < mRows; k++)
    {
        if (fabs(mx.at(k,k)) <= DBL_EPSILON)
        {
            bool swiched = false;
            for (unsigned int p=k+1; p<mRows; p++)
            {
                if (fabs(mx.at(p,k)) > DBL_EPSILON)
                {
                    mx.switchRows(k, p);
                    swiched = true;
                    //break;
                }
                if (swiched==true) break;
            }
        }

        for (unsigned int j=k+1; j<mRows; j++)
        {
            double c = -(mx.at(j,k)/mx.at(k,k));
            for (unsigned int i=k; i<mCols; i++)
            {
                mx.at(j,i) = mx.at(j,i) + mx.at(k,i) * c;
            }
        }
    }
    det = 1.0;
    for (unsigned int i=0; i<mRows; i++) det *= mx.at(i,i);
    mx.clear();
    // Gaussian elimination ---------------------------------------------------------------------------------

    return det;
}

void DoubleMatrix::transpose()
{
    if (empty()) return;

    unsigned int rows = mCols;
    unsigned int cols = mRows;
    double **data = static_cast<double**>(malloc(sizeof(double*)*rows));
    for (unsigned int i=0; i<rows; i++)
    {
        data[i] = static_cast<double*>(malloc(sizeof(double)*cols));
        for (unsigned int j=0; j<cols; j++) data[i][j] = mData[j][i];
    }
    clear();
    mRows = rows;
    mCols = cols;
    mData = data;
}

void DoubleMatrix::inverse()
{
    double det = determinant();
    if (fabs(det) <= DBL_EPSILON) throw double_matrix_exception(4);

    double idet = 1.0/det;

    DoubleMatrix m = *this;

    for (unsigned int r=0; r<rows(); r++)
    {
        for (unsigned int c=0; c<cols(); c++)
        {
            DoubleMatrix minor = m.minor(r,c);
            mData[r][c] = minor.determinant();
            if ((r+c)%2==1) mData[r][c] *= -1.0;
            minor.clear();
        }
    }

    transpose();

    for (unsigned int r=0; r<rows(); r++)
    {
        for (unsigned int c=0; c<cols(); c++)
        {
            mData[r][c] *= idet;
        }
    }

    m.clear();
}

DoubleMatrix DoubleMatrix::minor(unsigned int row, unsigned int col) const
{
    DoubleMatrix m(mRows-1, mCols-1);
    for (unsigned int i=0; i<mRows; i++)
    {
        for (unsigned int j=0; j<mCols; j++)
        {
            if (i == row && j == col) continue;
            if (i < row  && j <  col) { m.mData[i][j]     = mData[i][j]; continue; }
            if (i < row  && j >  col) { m.mData[i][j-1]   = mData[i][j]; continue; }
            if (i > row  && j <  col) { m.mData[i-1][j]   = mData[i][j]; continue; }
            if (i > row  && j >  col) { m.mData[i-1][j-1] = mData[i][j]; continue; }
        }
    }
    return m;
}

///////////////////////////////////////////////////////////////////////////////////////////////

double* DoubleMatrix::operator [](unsigned int row) const
{
    return mData[row];
}

double* DoubleMatrix::operator [](unsigned int row)
{
    return mData[row];
}

double& DoubleMatrix::operator()(unsigned int row, unsigned int col)
{
    return mData[row][col];
}

const double& DoubleMatrix::operator()(unsigned int row, unsigned int col) const
{
    return mData[row][col];
}

///////////////////////////////////////////////////////////////////////////////////////////////

DoubleMatrix& DoubleMatrix::operator =(const DoubleMatrix &other)
{
    // the matrix object holds reusable storage, such as a heap-allocated buffer mData
    //puts("DoubleMatrix& DoubleMatrix::operator =(const DoubleMatrix &m)");

    if (this == &other) return *this; // self-assignment check expected

    clear();
    if (other.mRows > 0 && other.mCols > 0)
    {
        mRows = other.mRows;
        mCols = other.mCols;
        mData = static_cast<double**>(malloc(sizeof(double*)*mRows));
        for (unsigned int i=0; i<mRows; i++)
        {
            mData[i] = static_cast<double*>(malloc(sizeof(double)*mCols));
            memcpy(mData[i], other.mData[i], sizeof(double)*mCols);
        }
    }
    return *this;
}

DoubleMatrix& DoubleMatrix::operator =(const DoubleVector &v)
{
    if (v.length() > 0)
    {
        clear();
        mRows = v.length();
        mCols = 1;
        mData = static_cast<double**>(malloc(sizeof(double*)*mRows));
        for (unsigned int i=0; i<mRows; i++)
        {
            mData[i] = static_cast<double*>(malloc(sizeof(double)*mCols));
            mData[i][0] = v.at(i);
        }
    }
    return *this;
}

DoubleMatrix& DoubleMatrix::operator +=(const DoubleMatrix &m)
{
    if (!dimEquals(m)) { throw double_matrix_exception(1); }

    for (unsigned int rw=0; rw<mRows; rw++)
    {
        for (unsigned int cl=0; cl<mCols; cl++)
        {
            mData[rw][cl] += m.mData[rw][cl];
        }
    }
    return *this;
}

DoubleMatrix& DoubleMatrix::operator -=(const DoubleMatrix &m)
{
    if (!dimEquals(m)) { throw double_matrix_exception(1); }

    for (unsigned int rw=0; rw<mRows; rw++)
    {
        for (unsigned int cl=0; cl<mCols; cl++)
        {
            mData[rw][cl] -= m.mData[rw][cl];
        }
    }
    return *this;
}

DoubleMatrix& DoubleMatrix::operator *=(const DoubleMatrix &m)
{
    if (mCols != m.mRows) { throw double_matrix_exception(3); }

    //    DoubleMatrix m;
    //    m.resize(m1.rows(), m2.cols());

    //    for (unsigned int i=0; i<m.rows(); i++)
    //    {
    //        for (unsigned int j=0; j<m.cols(); j++)
    //        {
    //            double sum = 0.0;
    //            for (unsigned int k=0; k<m1.cols(); k++) sum += m1.at(i,k)*m2.at(k,j);
    //            m.at(i,j) = sum;
    //        }
    //    }
    //    return m;

    //    for (unsigned int rw=0; rw < mRows; rw++)
    //    {
    //        for (unsigned int cl=0; cl < mCols; cl++)
    //        {
    //            mData[rw][cl] -= m.mData[rw][cl];
    //        }
    //    }
    return *this;
}

DoubleMatrix& DoubleMatrix::operator *=(double scalar)
{
    for (unsigned int row=0; row<mRows; row++)
    {
        for (unsigned int col=0; col<mCols; col++)
        {
            mData[row][col] *= scalar;
        }
    }
    return *this;
}

DoubleMatrix& DoubleMatrix::operator *=(const DoubleVector& v)
{
    if (mCols != v.mLength) { throw double_matrix_exception(3); }

    DoubleMatrix m = *this;

    resize(mRows, 1);
    for (unsigned int row=0; row<mRows; row++)
    {
        mData[row][0] = 0.0;
        for (unsigned int col=0; col<m.mCols; col++) mData[row][0] += m.mData[row][col]*v[col];
    }
    m.clear();

    return *this;
}

///////////////////////////////////////////////////////////////////////////////////////////////

DoubleMatrix operator +(DoubleMatrix m1, const DoubleMatrix& m2)
{
    if (!m1.dimEquals(m2)) { throw double_matrix_exception(1); }

    for (unsigned int row=0; row<m1.mRows; row++)
    {
        for (unsigned int col=0; col<m1.mCols; col++)
        {
            m1.mData[row][col] += m2.mData[row][col];
        }
    }
    return m1;
}

DoubleMatrix operator -(DoubleMatrix m1, const DoubleMatrix& m2)
{
    if (!m1.dimEquals(m2)) { throw double_matrix_exception(1); }

    for (unsigned int row=0; row<m1.mRows; row++)
    {
        for (unsigned int col=0; col<m1.mCols; col++)
        {
            m1.mData[row][col] -= m2.mData[row][col];
        }
    }
    return m1;
}

DoubleMatrix operator *(const DoubleMatrix &m1, const DoubleMatrix &m2)
{
    if (m1.cols() != m2.rows()) { throw double_matrix_exception(3); }

    DoubleMatrix m;
    m.resize(m1.rows(), m2.cols());

    for (unsigned int i=0; i<m.rows(); i++)
    {
        for (unsigned int j=0; j<m.cols(); j++)
        {
            double sum = 0.0;
            for (unsigned int k=0; k<m1.cols(); k++) sum += m1.at(i,k)*m2.at(k,j);
            m.at(i,j) = sum;
        }
    }
    return m;
}

DoubleMatrix operator *(double scalar, DoubleMatrix m)
{
    unsigned int m_row = m.rows();
    unsigned int m_col = m.cols();

    for (unsigned int row=0; row < m_row; row++)
    {
        for (unsigned int col=0; col < m_col; col++)
        {
            m.mData[row][col] *= scalar;
        }
    }

    return m;
}

DoubleMatrix operator *(DoubleMatrix m, double scalar)
{
    return scalar*m;
}

DoubleMatrix operator *(DoubleMatrix m, const DoubleVector& v)
{
    if (m.mCols != v.length()) { throw double_matrix_exception(3); }

    DoubleMatrix m1(m.mRows, 1);

    for (unsigned int row=0; row<m.mRows; row++)
    {
        m1.mData[row][0] = 0.0;
        for (unsigned int col=0; col<m.mCols; col++) m1[row][0] += m.mData[row][col]*v[col];
    }

    return m1;
}

///////////////////////////////////////////////////////////////////////////////////////////////

bool operator ==(const DoubleMatrix &matrix1, const DoubleMatrix& matrix2)
{
    if (!matrix1.dimEquals(matrix2)) return false;

    for (unsigned int row = 0; row < matrix1.mRows; row++)
    {
        for (unsigned int col = 0; col < matrix1.mCols; col++)
        {
            if (fabs(matrix1.mData[row][col] - matrix2.mData[row][col]) >= DBL_EPSILON) return false;
        }
    }
    return true;
}

bool operator !=(const DoubleMatrix &matrix1, const DoubleMatrix& matrix2)
{
    return !(matrix1==matrix2);
}

///////////////////////////////////////////////////////////////////////////////////////////////

void DoubleMatrix::switchRows(unsigned int row1, unsigned int row2)
{
    double *row = static_cast<double*>(malloc(sizeof(double) * mCols));
    memcpy(row, mData[row1], sizeof(double)*mCols);
    memcpy(mData[row1], mData[row2], sizeof(double)*mCols);
    memcpy(mData[row2], row, sizeof(double)*mCols);
    free(row);
}

void DoubleMatrix::switchCols(unsigned int, unsigned int)
{
}
