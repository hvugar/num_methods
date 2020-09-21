#include "printer.h"
#include <stdio.h>
#include <time.h>

#ifdef _WIN32
#include <windows.h>
#endif

IPrinter::~IPrinter() {}

void IPrinter::printMatrix(const DoubleMatrix &x, unsigned int m, unsigned int n, const char* s, FILE* f)
{
    C_UNUSED(s);

    size_t rows = x.rows();
    size_t cols = x.cols();
    size_t M = rows / m;

    for (size_t j=0; j<rows; j++)
    {
        size_t N = cols / n;
        if (j%M==0)
        {
            for (size_t i=0; i<cols; i++)
            {
                if (i%N==0) fprintf(f, "%14.10f ", x.at(j,i));
            }
            fputs("\n", f);
        }
    }
    fflush(f);
}

void IPrinter::printMatrix(unsigned int width, unsigned int presicion, const DoubleMatrix &x, unsigned int m, unsigned int n, const char* s, FILE* f)
{
    C_UNUSED(s);

    char format[10] = {0};
    int sz = sprintf(format, "%%%d.%df ", width, presicion);
    format[sz] = '\0';

    size_t rows = x.rows();
    size_t cols = x.cols();
    size_t M = rows / m;

    for (size_t j=0; j<rows; j++)
    {
        size_t N = cols / n;
        if (j%M==0)
        {
            for (size_t i=0; i<cols; i++)
            {
                if (i%N==0) fprintf(f, format, x.at(j,i));
            }
            fputs("\n", f);
        }
    }
    fflush(f);
}

void IPrinter::printAsMatrix(const DoubleVector &x, unsigned int M, unsigned int N, unsigned int m, unsigned int n, const char* s, FILE* f)
{
    C_UNUSED(s);

    for (size_t j=0; j<=M; j++)
    {
        if (j%(M/m)==0)
        {
            for (size_t i=0; i<=N; i++)
            {
                if (i%(N/n)==0) fprintf(f, "%12.8f ", x[j*(N+1)+i]);
            }
            fputs("\n", f);
        }
    }
    fflush(f);
}

void IPrinter::printVector(const DoubleVector &x, const char *s, unsigned int n, unsigned int start, unsigned int end, FILE *file)
{
    if (s!=nullptr) fprintf(file, "%s", s);
    if (start != 0 || end != 0)
    {
        size_t N = (end-start+1) / n;
        for (size_t i=start; i<=end; i++)
        {
            if ((i-start)%N==0) fprintf(file, "%10.6f ", x[i]);
        }
    }
    else
    {
        size_t N = x.length() / n;
        for (size_t i=0; i<x.length(); i++)
        {
            if (i%N==0) fprintf(file, "%10.6f ", x[i]);
        }
    }
    fputs("\n", file);
    fflush(file);
}

void IPrinter::printVector(const std::vector<double> &x, const char *s, unsigned int n, unsigned int start, unsigned int end, FILE *file)
{
    if (s!=nullptr) fprintf(file, "%s", s);
    if (start != 0 || end != 0)
    {
        unsigned int N = (end-start+1) / n;
        for (unsigned int i=start; i<=end; i++)
        {
            if ((i-start)%N==0) fprintf(file, "%14.10f ", x[i]);
        }
    }
    else
    {
        unsigned int size = static_cast<unsigned int>(x.size());
        unsigned int N = size / n;
        for (unsigned int i=0; i<size; i++)
        {
            if (i%N==0) fprintf(file, "%14.10f ", x[i]);
        }
    }
    fputs("\n", file);
    fflush(file);
}

void IPrinter::printVector(double *x, unsigned int size, const char *s, unsigned int n, unsigned int start, unsigned int end, FILE *file)
{
    if (s!=nullptr) fprintf(file, "%s", s);
    if (start != 0 || end != 0)
    {
        unsigned int N = (end-start+1) / n;
        for (unsigned int i=start; i<=end; i++)
        {
            if ((i-start)%N==0) fprintf(file, "%14.10f ", x[i]);
        }
    }
    else
    {
        unsigned int N = size / n;
        for (unsigned int i=0; i<size; i++)
        {
            if (i%N==0) fprintf(file, "%14.10f ", x[i]);
        }
    }
    fputs("\n", file);
    fflush(file);
}

void IPrinter::printVector(double *x, unsigned int size, const char *s, unsigned int n, unsigned int start, unsigned int end, const char *filename)
{
    FILE *file = fopen(filename, "w");
    if (s!=nullptr) fprintf(file, "%s", s);
    if (start != 0 || end != 0)
    {
        unsigned int N = (end-start+1) / n;
        for (unsigned int i=start; i<=end; i++)
        {
            if ((i-start)%N==0) fprintf(file, "%6.2f ", x[i]);
        }
    }
    else
    {
        unsigned int N = size / n;
        for (unsigned int i=0; i<size; i++)
        {
            if (i%N==0) fprintf(file, "%6.2f ", x[i]);
        }
    }
    fputs("\n", file);
    fflush(file);
}

void IPrinter::printVector(unsigned int width, unsigned int presicion, const DoubleVector &x, const char *s, unsigned int n, unsigned int start, unsigned int end, FILE *file)
{
    char format[10] = {0};
    int sz = sprintf(format, "%%%d.%df ", width, presicion);
    format[sz] = '\0';

    if (s!=nullptr) fprintf(file, "%s", s);
    if (start != 0 || end != 0)
    {
        size_t N = static_cast<size_t>(end-start+1) / n;
        for (size_t i=start; i<=end; i++)
        {
            if ((i-start)%N==0) fprintf(file, format, x[i]);
        }
    }
    else
    {
        size_t N = x.length() / n;
        for (size_t i=0; i<x.length(); i++)
        {
            if (i%N==0) fprintf(file, format, x[i]);
        }
    }
    fputs("\n", file);
    fflush(file);
}

void IPrinter::printCube(const DoubleVector& x, unsigned int M, unsigned int N2, unsigned int N1, FILE *file)
{
    for (unsigned int k=0; k<=M; k++)
    {
        if (k%(M/10)==0)
        {
            fprintf(file, "k: %d\n", k);
            for (unsigned int j=0; j<=N2; j++)
            {
                if (j%(N2/10)==0)
                {
                    for (unsigned int i=0; i<=N1; i++)
                    {
                        if (i%(N1/10)==0)
                        {
                            fprintf(file, "%14.10f ", x[k*(N2+1)*(N1+1) + j*(N1+1) + i]);
                        }
                    }
                    fprintf(file, "\n");
                }
            }
        }
    }
}

void IPrinter::printDateTime(FILE *file)
{
    time_t t = time(nullptr);
    struct tm * now = localtime( & t );
    char buf[80];
    strftime(buf, sizeof(buf), "%Y-%m-%d %X", now);
    fprintf(file, "%s\n", buf);
}

void IPrinter::print(const DoubleMatrix &m, size_t /*M*/, size_t /*N*/, size_t width, size_t presicion, FILE *file)
{
    char format[10] = {0};
    int sz = sprintf(format, "%%%d.%df ", width, presicion);
    format[sz] = '\0';

    size_t rows = m.rows();
    size_t cols = m.cols();

    for (unsigned int i=0; i<rows; i++)
    {
        for (unsigned int j=0; j<cols; j++)
        {
            fprintf(file, format, m.at(i,j));
        }
        fputs("\n", file);
    }
    fflush(file);
}

void IPrinter::print(const DoubleMatrix &m, const char *filename, size_t M, size_t N, size_t width, size_t presicion)
{
    C_UNUSED(M);
    C_UNUSED(N);

    char format[10] = {0};
    int sz = sprintf(format, "%%%d.%df ", width, presicion);
    format[sz] = '\0';

    size_t rows = m.rows();
    size_t cols = m.cols();

    FILE* file = fopen(filename, "w");
    for (unsigned int i=0; i<rows; i++)
    {
        for (unsigned int j=0; j<cols; j++)
        {
            fprintf(file, format, m.at(i,j));
        }
        fputs("\n", file);
        fflush(file);
    }
    fclose(file);
}

void IPrinter::print(const DoubleVector &v, size_t N, size_t width, size_t presicion, FILE *file)
{
    C_UNUSED(N);

    char format[10] = {0};
    int sz = sprintf(format, "%%%d.%df ", width, presicion);
    format[sz] = '\0';

    //unsigned int size = v.length();

    for (unsigned int i=0; i<N; i++)
    {
        fprintf(file, format, v.at(i));
    }
    fputs("\n", file);
    fflush(file);
}

void IPrinter::print(const double *v, unsigned int N, unsigned int width, unsigned int presicion, FILE *file)
{
    C_UNUSED(N);

    char format[10] = {0};
    int sz = sprintf(format, "%%%d.%df ", width, presicion);
    format[sz] = '\0';

    //    unsigned int size = v.length();

    for (unsigned int i=0; i<N; i++)
    {
        fprintf(file, format, v[i]);
    }
    fputs("\n", file);
    fflush(file);
}

void IPrinter::print(std::vector<DoubleVector> &rv, size_t k, size_t N, size_t width, size_t presicion, FILE *file)
{
    char format[10] = {0};
    int sz = sprintf(format, "%%%d.%df ", width, presicion);
    format[sz] = '\0';

    const auto size = rv.size();
    const auto modN = size/N;

    for (unsigned int row=0; row<k; row++)
    {
        for (unsigned int col=0; col<size; col++)
        {
            if (col%modN==0) fprintf_s(file, format, rv[col][row]);
        }
        fputs("\n", file);
    }
    fflush(file);
}


void IPrinter::printSeperatorLine(const char* msg, char c, FILE* file)
{
    size_t columns=10;

#ifdef _INC_WINDOWS
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
    columns = static_cast<size_t>(csbi.srWindow.Right - csbi.srWindow.Left) + 1;
#endif

    size_t start = 0;
    if (msg != nullptr)
    {
        fprintf(file, "%s ", msg);
        start = strlen(msg)+1;
    }

    for (size_t i=start; i<columns-1; i++) fprintf(file,"%c", c);

#ifdef _INC_WINDOWS
    fputs("\n", file);
#else
    fprintf(file, "\n");
#endif
    fflush(file);
}
