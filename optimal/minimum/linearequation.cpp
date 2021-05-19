#include "linearequation.h"
#include <math.h>
#include <float.h>
#include <cmethods.h>
#include <vector>

/**
 * @brief GaussianElimination1
 * @param m основная матрица системы
 * @param b столбец свободных членов
 * @param x
 */
void GaussianElimination1(const DoubleMatrix& m, const DoubleVector& b, DoubleVector &x)
{
    DoubleMatrix M = m;
    DoubleVector B = b;

    const unsigned int ui = (unsigned)0-1;

    unsigned int n = x.length();

    for (unsigned k=0; k<n; k++)
    {
        if (fabs(M.at(k,k)) <= DBL_EPSILON)
        {
            for (unsigned int p=k+1; p<n; p++)
            {
                if (fabs(M[p][k]) > DBL_EPSILON)
                {
                    M.switchRows(k, p);
                    double bk = B[k];
                    B[k] = B[p];
                    B[p] = bk;
                    break;
                }
            }
        }

        for (unsigned int j=(k+1); j<n; j++)
        {
            double c = M.at(j,k)/M.at(k,k);
            for (unsigned int i=k; i<n; i++)
            {
                M.at(j,i) = M.at(j,i) - M.at(k,i) * c;
            }
            B[j] = B[j] - B[k] *c;
        }
    }

    for (unsigned int i=(n-1); i!=ui; i--)
    {
        for (unsigned int j=(n-1); j>i; j--) B[i] -= (M.at(i,j) * x[j]);
        x[i] = B[i] / M.at(i,i);
    }
}

/**
 * @brief Метод Гаусса с выбором главного элемента
 * @param A - основная матрица системы
 * @param b - вектор свободных членов
 * @param x - вектор искомых значений
 */
void GaussianElimination2(const DoubleMatrix &A, const DoubleVector &b, DoubleVector &x)
{
    if (A.rows() == 0 || A.cols() == 0) throw double_matrix_exception(0);
    if (A.rows() != A.cols()) throw double_matrix_exception(0);
    if (b.length() != A.rows()) throw double_matrix_exception(0);

    DoubleMatrix A1 = A;
    DoubleVector b1 = b;
    unsigned int size = b.length();

    // Поставленная задача будет решаться методом Гаусса с выбором главного элемента по столбцу

    //    IPrinter::printSeperatorLine(nullptr, '=');
    //    IPrinter::print(A1, A1.rows(), A1.cols());
    //    IPrinter::printSeperatorLine(nullptr, '=');
    //

    //    IPrinter::print(b1, b1.length());

    for (unsigned int c=0; c<size; c++)
    {
        unsigned int i = c;
        double max = fabs(A1[i][c]);
        for (unsigned int r=c; r<size; r++)
        {
            if (fabs(A1[r][c]) > max) { max = fabs(A1[r][c]); i = r; }
        }
        if (i != c)
        {
            A1.switchRows(i, c);
            const double sw = b1[i]; b1[i] = b1[c]; b1[c] = sw;
            //printf(">>> %d > %d\n", i, c);
        }

        b1[c] /= A1[c][c];
        for (unsigned int s=c+1; s<size; s++)
        {
            A1[c][s] /= A1[c][c];
        }
        A1[c][c] = 1.0;

        for (unsigned int r=c+1; r<size; r++)
        {
            b1[r] -= A1[r][c]*b1[c];
            for (unsigned int s=c+1; s<size; s++)
            {
                A1[r][s] -= A1[r][c]*A1[c][s];
            }
            A1[r][c] = 0.0;
        }
        //        IPrinter::printSeperatorLine(nullptr, '=');
        //        IPrinter::print(A1, A1.rows(), A1.cols());
        //        IPrinter::printSeperatorLine();
        //        IPrinter::print(b1, b1.length());
    }

    for (unsigned int i=0, r=size-1; i<size; i++, r--)
    {
        x[r] = b1[r];
        for (unsigned int j=0, c=size-1; j<i; j++, c--) x[r] -= A1[r][c]*x[c];
    }
}

void GaussianElimination3(const double** A, const double* b, double* x, size_t N)
{
    double** A1 = new double*[N];
    for (size_t i=0; i<N; i++)
    {
        A1[i] = new double[N];
        memcpy(A1[i], A[i], sizeof (double)*N);
    }
    double* b1 = new double[N];
    memcpy(b1, b, sizeof (double)*N);

    // Поставленная задача будет решаться методом Гаусса с выбором главного элемента по столбцу
    for (size_t c=0; c<N; c++)
    {
        size_t i = c;
        double max = fabs(A1[i][c]);
        for (size_t r=c; r<N; r++)
        {
            if (fabs(A1[r][c]) > max) { max = fabs(A1[r][c]); i = r; }
        }
        if (i != c)
        {
            double *tmp = A1[i]; A1[i] = A1[c]; A1[c] = tmp;
            double  sw  = b1[i]; b1[i] = b1[c]; b1[c] = sw;
        }

        b1[c] /= A1[c][c];
        for (size_t s=c+1; s<N; s++)
        {
            A1[c][s] /= A1[c][c];
        }
        A1[c][c] = 1.0;

        for (size_t r=c+1; r<N; r++)
        {
            b1[r] -= A1[r][c]*b1[c];
            for (size_t s=c+1; s<N; s++)
            {
                A1[r][s] -= A1[r][c]*A1[c][s];
            }
            A1[r][c] = 0.0;
        }
    }

    for (size_t i=0, r=N-1; i<N; i++, r--)
    {
        x[r] = b1[r];
        for (size_t j=0, c=N-1; j<i; j++, c--) x[r] -= A1[r][c]*x[c];
    }

    for (size_t i=0; i<N; i++)
    {
        delete [] A1[i];
    }
    delete [] A1;
    delete [] b1;
}

void LinearEquation::GaussianElimination(const DoubleMatrix& m, const DoubleVector& b, DoubleVector& x)
{
    GaussianElimination2(m,b,x);
    //GaussianElimination3(m.data(), b.data(), x.data(), x.length());
}

void LinearEquation::FirstRowLoaded(const double *e, double f, const double *a, const double *b, const double *c, const double *d, unsigned int N)
{
    double *p = (double*) malloc(sizeof(double)*N);
    double *q = (double*) malloc(sizeof(double)*N);
    double *r = (double*) malloc(sizeof(double)*N);
    double *x = (double*) malloc(sizeof(double)*N);

    p[0] = e[0];
    q[0] = e[1];
    p[0] = f;

    for (unsigned int n=1; n<N-1; n++)
    {
        p[n] = -p[n-1]*(b[n]/a[n]) + q[n-1];
        q[n] = -p[n-1]*(c[n]/a[n]) + e[n+1];
        r[n] = -p[n-1]*(d[n]/a[n]) + r[n-1];
    }

    p[N-1] = -p[N-2]*(b[N-1]/a[N-1]) + q[N-2];
    q[N-1] = 0.0;
    r[N-1] = -p[N-2]*(d[N-1]/a[N-1]) + r[N-2];

    x[N-1] = r[N-1]/p[N-1];

    printf("%f\n", x[N-1]);

    free(r);
    free(q);
    free(p);
}

void LinearEquation::func1(const double *a, const double *b, const double *c, const double *d, double **e, double *x, unsigned int N)
{
    printf(">>> %d\n", N);
    std::vector<unsigned int> selectedCols;
    for (unsigned int col=0; col<N; col++)
    {
        for (unsigned int row=0; row<N; row++)
        {
            if (e[row][col] != 0.0)
            {
                selectedCols.push_back(col);
                break;
            }
        }
    }
    size_t selectedColsSize = selectedCols.size();

    if (selectedColsSize == 0)
    {
        printf(">>> %d %d\n", N, selectedColsSize);
        tomasAlgorithm(a, b, c, d, x, N);
    }
    else
    {

        double *v = new double[N];
        //double *v = (double*) malloc(sizeof(double)*N);
        tomasAlgorithm(a, b, c, d, v, N);

        //double **w = (double**) malloc(sizeof(double*) * selectedColsSize);
        double **w = new double*[selectedColsSize];

        for (unsigned int i=0; i<selectedColsSize; i++)
        {
            //w[i] = (double*) malloc(sizeof(double)*N);
            w[i] = new double[N];
        }

        for (unsigned int sc=0; sc<selectedColsSize; sc++)
        {
            for (unsigned int row=0; row<N; row++)
            {
                w[sc][row] = -e[row][selectedCols[sc]];
            }
            tomasAlgorithm(a, b, c, w[sc], w[sc], N);
        }

        DoubleMatrix M(selectedColsSize, selectedColsSize, 0.0);
        DoubleVector A(selectedColsSize);
        DoubleVector u(selectedColsSize, 0.0);
        for (unsigned int scr=0; scr<selectedColsSize; scr++)
        {
            A[scr] = v[selectedCols[scr]];
            for (unsigned int scc=0; scc<selectedColsSize; scc++)
            {
                M[scr][scc] = -w[scc][selectedCols[scr]];
                if (scr==scc) M[scr][scc] += 1.0;
            }
        }

        LinearEquation::GaussianElimination(M, A, u);

        for (unsigned int i=0; i<N; i++)
        {
            x[i] = v[i];
            for (unsigned int sc=0; sc<selectedColsSize; sc++) x[i] += w[sc][i]*u[sc];
        }

        M.clear();
        A.clear();
        u.clear();
        //    for (unsigned int sc=0; sc<selectedColsSize; sc++) free(w[sc]);
        for (unsigned int sc=0; sc<selectedColsSize; sc++) delete [] w[sc];
        delete [] w;
        delete [] v;
        selectedCols.clear();
    }
}
