#include "matrix2d.h"
#include "printer.h"
#include "exceptions.h"
#include <time.h>
#include <math.h>
#include <float.h>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

DoubleMatrix::DoubleMatrix(size_t rows, size_t cols, double value) : mRows(0), mCols(0), mData(nullptr)
{
    if (rows > 0 && cols > 0)
    {
        mRows = rows;
        mCols = cols;
        mData = static_cast<double**>(malloc(sizeof(double*)*rows));
        for (size_t i=0; i<rows; i++)
        {
            mData[i] = static_cast<double*>(malloc(sizeof(double)*cols));
            for (size_t j=0; j<cols; j++) mData[i][j] = value;
        }
    }
}

DoubleMatrix::DoubleMatrix(const DoubleMatrix &matrix) : mRows(0), mCols(0), mData(nullptr)
{
    if (matrix.mRows > 0 && matrix.mCols > 0)
    {
        mRows = matrix.mRows;
        mCols = matrix.mCols;
        mData = static_cast<double**>(malloc(sizeof(double*)*mRows));
        for (size_t i=0; i<mRows; i++)
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
        for (size_t i=0; i<mRows; i++)
        {
            mData[i] = static_cast<double*>(malloc(sizeof(double)*mCols));
            mData[i][0] = vector.at(i);
        }
    }
}

DoubleMatrix::DoubleMatrix(const double* const data, size_t length, size_t rows, size_t cols) : mRows(rows), mCols(cols), mData(nullptr)
{
    if (data == nullptr || length) return;


}

DoubleMatrix::DoubleMatrix(double** data, size_t rows, size_t cols)
{
    if (data == nullptr || rows == 0 || cols == 0) throw std::exception();

    mRows = rows;
    mCols = cols;
    mData = static_cast<double**>(malloc(sizeof(double*)*mRows));
    for (size_t r=0; r<rows; r++)
    {
        mData[r] = static_cast<double*>(malloc(sizeof(double)*mCols));
        memcpy(mData[r], data[r], sizeof(double)*mCols);
    }
}

DoubleMatrix::~DoubleMatrix()
{
    clear();
}

size_t DoubleMatrix::rows() const
{
    return mRows;
}

size_t DoubleMatrix::cols() const
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

    for (size_t row=0; row<mRows; row++)
    {
        free(mData[row]);
        mData[row] = nullptr;
    }

    free(mData);
    mData = nullptr;
    mRows = 0;
    mCols = 0;
}

void DoubleMatrix::resize(size_t rows, size_t cols, double value)
{
    if (rows == 0 && cols == 0) clear();

    if (rows == 0 || cols == 0) return;

    if (rows > 0 && cols > 0)
    {
        if (mData == nullptr)
        {
            mData = static_cast<double**>(malloc(sizeof(double*)*rows));
            for (size_t row=0; row<rows; row++)
            {
                mData[row] = static_cast<double*>(malloc(sizeof(double)*cols));
                for (size_t col=0; col<cols; col++) mData[row][col] = value;
            }
            mRows = rows;
            mCols = cols;
        }
        else
        {
            double **ptr = static_cast<double**>(malloc(sizeof(double*) * rows));
            for (size_t row=0; row<rows; row++)
            {
                ptr[row] = static_cast<double*>(malloc(sizeof(double)*cols));

                if (row < mRows)
                {
                    for (size_t col=0; col<cols; col++)
                        if (col < mCols) ptr[row][col] = mData[row][col]; else ptr[row][col] = value;
                }
                else
                {
                    for (size_t col=0; col<cols; col++)
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
    for (size_t r=0; r<mRows; r++)
    {
        for (size_t c=0; c<mCols; c++)
        {
            mData[r][c] = value;
        }
    }
}

double& DoubleMatrix::at(size_t row, size_t col)
{
    //if (row>=mRows) { throw std::out_of_range("Error: Row index out of range"); }
    //if (col>=mCols) { throw std::out_of_range("Error: Column index out of range"); }
    if (row>=mRows) { throw double_matrix_exception(5); }
    if (col>=mCols) { throw double_matrix_exception(6); }
    return mData[row][col];
}

const double& DoubleMatrix::at(size_t row, size_t col) const
{
    //if (row>=mRows) { throw std::out_of_range("Error: Row index out of range"); }
    //if (col>=mCols) { throw std::out_of_range("Error: Column index out of range"); }
    if (row>=mRows) { throw double_matrix_exception(5); }
    if (col>=mCols) { throw double_matrix_exception(6); }
    return mData[row][col];
}

DoubleVector DoubleMatrix::row(size_t r) const
{
    return DoubleVector(mData[r], mCols);
}

DoubleVector DoubleMatrix::col(size_t c) const
{
    DoubleVector v(mRows);
    for (size_t r=0; r<mRows; r++) v[r] = mData[r][c];
    return v;
}

void DoubleMatrix::setColumn(size_t c, const DoubleVector &col)
{
    if (mRows == col.length() && c < mCols)
    {
        for (size_t r=0; r<mRows; r++) mData[r][c] = col.at(r);
    }
    else
    {
        throw double_matrix_exception(0);
    }
}

void DoubleMatrix::setRow(size_t r, const DoubleVector &row)
{
    if (mCols == row.length() && r < mRows)
    {
        for (size_t c=0; c<mCols; c++) mData[r][c] = row.at(c);
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
    size_t n = v.length();
    DoubleMatrix m(n, n, 0.0);
    for (size_t i=0; i<n; i++) m.mData[i][i] = v[i];
    return m;
}

DoubleMatrix DoubleMatrix::IdentityMatrix(size_t n)
{
    DoubleMatrix m(n, n, 0.0);
    for (size_t i=0; i<n; i++) m.mData[i][i] = 1.0;
    return m;
}

DoubleMatrix DoubleMatrix::ZeroMatrix(size_t rows, size_t cols)
{
    DoubleMatrix matrix(rows, cols, 0.0);
    return matrix;
}

DoubleMatrix DoubleMatrix::HilbertMatrix(size_t rows, size_t cols)
{
    DoubleMatrix matrix(rows, cols);
    for (size_t row=0; row<rows; row++)
    {
        for (size_t col=0; col<cols; col++)
        {
            matrix.mData[row][col] = 1.0/static_cast<double>(row+col+1);
        }
    }
    return matrix;
}

///////////////////////////////////////////////////////////////////////////////////////////////

bool DoubleMatrix::dioqonalMatrix() const
{
    for (size_t row=0; row < mRows; row++)
    {
        for (size_t col=0; col < mCols; col++)
        {
            if (row != col && fabs(mData[row][col]) >= DBL_EPSILON) return false;
        }
    }
    return true;
}

bool DoubleMatrix::identityMatrix() const
{
    for (size_t row=0; row<mRows; row++)
    {
        for (size_t col=0; col<mCols; col++)
        {
            if (row != col && fabs(mData[row][col]) >= DBL_EPSILON) return false;
            if (row == col && fabs(mData[row][col]) != 1.0) return false;
        }
    }
    return true;
}

bool DoubleMatrix::zeroMatrix() const
{
    for (size_t row=0; row<mRows; row++)
    {
        for (size_t col=0; col<mCols; col++)
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
    for (size_t i=0; i<mRows; i++)
    {
        for (size_t j=0; j<mCols; j++)
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
    for (size_t row=0; row<mRows; row++)
    {
        for (size_t col=0; col<mCols; col++)
        {
            if (mData[row][col] < minimum) minimum = mData[row][col];
        }
    }
    return minimum;
}

double DoubleMatrix::max() const
{
    double maximum = DBL_MIN;
    for (size_t row=0; row<mRows; row++)
    {
        for (size_t col=0; col<mCols; col++)
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
        size_t i;
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

    if (mRows == 0) { throw double_matrix_exception(0); }

    if (mRows == 1) return mData[0][0];

    // Checking properties of the determinant ---------------------------------------------------------------
    for (size_t r=0; r<mRows; r++)
    {
        bool row_items_are_zero = true;
        for (size_t c=0; c<mCols; c++) if (fabs(mData[r][c]) >= DBL_EPSILON) row_items_are_zero = false;
        if (row_items_are_zero) return 0.0;
    }
    for (size_t c=0; c<mRows; c++)
    {
        bool col_items_are_zero = true;
        for (size_t r=0; r<mCols; r++) if (fabs(mData[r][c]) >= DBL_EPSILON) col_items_are_zero = false;
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
            for (size_t p=k+1; p<mRows; p++)
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

        for (size_t j=k+1; j<mRows; j++)
        {
            double c = -(mx.at(j,k)/mx.at(k,k));
            for (size_t i=k; i<mCols; i++)
            {
                mx.at(j,i) = mx.at(j,i) + mx.at(k,i) * c;
            }
        }
    }
    det = 1.0;
    for (size_t i=0; i<mRows; i++) det *= mx.at(i,i);
    mx.clear();
    // Gaussian elimination ---------------------------------------------------------------------------------

    return det;
}

void DoubleMatrix::transpose()
{
    if (empty()) return;

    size_t rows = mCols;
    size_t cols = mRows;
    double **data = static_cast<double**>(malloc(sizeof(double*)*rows));
    for (size_t i=0; i<rows; i++)
    {
        data[i] = static_cast<double*>(malloc(sizeof(double)*cols));
        for (size_t j=0; j<cols; j++) data[i][j] = mData[j][i];
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

    if (mRows == 1)
    {
        mData[0][0] = idet;
        return;
    }

    DoubleMatrix m = *this;

    for (size_t r=0; r<rows(); r++)
    {
        for (size_t c=0; c<cols(); c++)
        {
            DoubleMatrix minor = m.minor(r,c);
            mData[r][c] = minor.determinant();
            if ((r+c)%2==1) mData[r][c] *= -1.0;
            minor.clear();
        }
    }

    transpose();

    for (size_t r=0; r<rows(); r++)
    {
        for (size_t c=0; c<cols(); c++)
        {
            mData[r][c] *= idet;
        }
    }

    m.clear();
}

DoubleMatrix DoubleMatrix::minor(size_t row, size_t col) const
{
    DoubleMatrix m(mRows-1, mCols-1);
    for (size_t i=0; i<mRows; i++)
    {
        for (size_t j=0; j<mCols; j++)
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

double* DoubleMatrix::operator [](size_t row) const
{
    return mData[row];
}

double* DoubleMatrix::operator [](size_t row)
{
    return mData[row];
}

double& DoubleMatrix::operator()(size_t row, size_t col)
{
    return mData[row][col];
}

const double& DoubleMatrix::operator()(size_t row, size_t col) const
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
        for (size_t i=0; i<mRows; i++)
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
        for (size_t i=0; i<mRows; i++)
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

    for (size_t rw=0; rw<mRows; rw++)
    {
        for (size_t cl=0; cl<mCols; cl++)
        {
            mData[rw][cl] += m.mData[rw][cl];
        }
    }
    return *this;
}

DoubleMatrix& DoubleMatrix::operator -=(const DoubleMatrix &m)
{
    if (!dimEquals(m)) { throw double_matrix_exception(1); }

    for (size_t rw=0; rw<mRows; rw++)
    {
        for (size_t cl=0; cl<mCols; cl++)
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

    //    for (size_t i=0; i<m.rows(); i++)
    //    {
    //        for (size_t j=0; j<m.cols(); j++)
    //        {
    //            double sum = 0.0;
    //            for (size_t k=0; k<m1.cols(); k++) sum += m1.at(i,k)*m2.at(k,j);
    //            m.at(i,j) = sum;
    //        }
    //    }
    //    return m;

    //    for (size_t rw=0; rw < mRows; rw++)
    //    {
    //        for (size_t cl=0; cl < mCols; cl++)
    //        {
    //            mData[rw][cl] -= m.mData[rw][cl];
    //        }
    //    }
    return *this;
}

DoubleMatrix& DoubleMatrix::operator *=(double scalar)
{
    for (size_t row=0; row<mRows; row++)
    {
        for (size_t col=0; col<mCols; col++)
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
    for (size_t row=0; row<mRows; row++)
    {
        mData[row][0] = 0.0;
        for (size_t col=0; col<m.mCols; col++) mData[row][0] += m.mData[row][col]*v[col];
    }
    m.clear();

    return *this;
}

///////////////////////////////////////////////////////////////////////////////////////////////

DoubleMatrix operator +(DoubleMatrix m1, const DoubleMatrix& m2)
{
    if (!m1.dimEquals(m2)) { throw double_matrix_exception(1); }

    for (size_t row=0; row<m1.mRows; row++)
    {
        for (size_t col=0; col<m1.mCols; col++)
        {
            m1.mData[row][col] += m2.mData[row][col];
        }
    }
    return m1;
}

DoubleMatrix operator -(DoubleMatrix m1, const DoubleMatrix& m2)
{
    if (!m1.dimEquals(m2)) { throw double_matrix_exception(1); }

    for (size_t row=0; row<m1.mRows; row++)
    {
        for (size_t col=0; col<m1.mCols; col++)
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

    for (size_t i=0; i<m.rows(); i++)
    {
        for (size_t j=0; j<m.cols(); j++)
        {
            double sum = 0.0;
            for (size_t k=0; k<m1.cols(); k++) sum += m1.at(i,k)*m2.at(k,j);
            m.at(i,j) = sum;
        }
    }
    return m;
}

DoubleMatrix operator *(double scalar, DoubleMatrix m)
{
    size_t m_row = m.rows();
    size_t m_col = m.cols();

    for (size_t row=0; row < m_row; row++)
    {
        for (size_t col=0; col < m_col; col++)
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

    for (size_t row=0; row<m.mRows; row++)
    {
        m1.mData[row][0] = 0.0;
        for (size_t col=0; col<m.mCols; col++) m1[row][0] += m.mData[row][col]*v[col];
    }

    return m1;
}

///////////////////////////////////////////////////////////////////////////////////////////////

bool operator ==(const DoubleMatrix &matrix1, const DoubleMatrix& matrix2)
{
    if (!matrix1.dimEquals(matrix2)) return false;

    for (size_t row = 0; row < matrix1.mRows; row++)
    {
        for (size_t col = 0; col < matrix1.mCols; col++)
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

void DoubleMatrix::switchRows(size_t row1, size_t row2)
{
    double *row = static_cast<double*>(malloc(sizeof(double) * mCols));
    memcpy(row, mData[row1], sizeof(double)*mCols);
    memcpy(mData[row1], mData[row2], sizeof(double)*mCols);
    memcpy(mData[row2], row, sizeof(double)*mCols);
    free(row);
}

void DoubleMatrix::switchCols(size_t, size_t)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////

Matrix2D::Matrix2D(size_t rows, size_t cols, double /*value*/) : m_rows(rows), m_cols(cols), m_data(nullptr)
{
    const size_t length = m_rows * m_cols;
    m_data = static_cast<double*>(std::calloc(length, sizeof(double)));
    if (m_data == nullptr) { throw std::exception(); }
}

Matrix2D::Matrix2D(const Matrix2D &mtx) : m_rows(mtx.m_rows), m_cols(mtx.m_cols)
{
    const size_t length = m_rows * m_cols;
    //m_data = new double[length];
    m_data = static_cast<double*>(std::malloc(length*sizeof(double)));
    std::memcpy(m_data, mtx.m_data, length);
}

Matrix2D::~Matrix2D()
{
    clear();
}

auto Matrix2D::rows() const -> size_t { return m_rows; }

auto Matrix2D::cols() const -> size_t { return m_cols; }

auto Matrix2D::at(size_t row, size_t col) const -> double
{
    if ((row >= m_rows) || (col >= m_cols)) { throw std::exception(); }
    return value(row, col);
}

auto Matrix2D::at(size_t row, size_t col) -> double &
{
    if ((row >= m_rows) || (col >= m_cols)) { throw std::exception(); }
    return value(row, col);
}

auto Matrix2D::value(size_t row, size_t col) const -> double
{
    size_t index = row*m_cols + col;
    return m_data[index];
}

auto Matrix2D::value(size_t row, size_t col) -> double &
{
    size_t index = row*m_cols + col;
    return m_data[index];
}

auto Matrix2D::min() const -> double
{
    const size_t length = m_rows * m_cols;
    double minimum = 0.0;
    if (length == 0) { std::exception(); } else { minimum = m_data[0]; }
    for (size_t i = 1; i < length; ++i)
    {
        if (m_data[i] < minimum) minimum = m_data[i];
    }
    return minimum;
}

auto Matrix2D::max() const -> double
{
    const size_t length = m_rows * m_cols;
    double maximum = 0.0;
    if (length == 0) { std::exception(); } else { maximum = m_data[0]; }
    for (size_t i = 1; i < length; ++i)
    {
        if (m_data[i] > maximum) maximum = m_data[i];
    }
    return maximum;
}

auto Matrix2D::fill(double value) -> void
{
    const size_t length = m_rows * m_cols;
    for (size_t i = 0; i < length; ++i) { m_data[i] = value; }
}

auto Matrix2D::clear() -> void
{
    if (m_data != nullptr) std::free(m_data);
    m_data = nullptr;
    m_rows = m_cols = 0;
}


auto Matrix2D::reset(size_t rows, size_t cols) -> void
{
    if (m_data != nullptr) std::free(m_data);
    m_rows = rows;
    m_cols = cols;
    const size_t length = m_rows * m_cols;
    m_data = static_cast<double*>(std::calloc(length, sizeof(double)));
    if (m_data == nullptr) { throw std::exception(); }
}

