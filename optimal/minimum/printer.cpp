#include "printer.h"
#include <stdio.h>
#include <time.h>

void IPrinter::print(GradientIterationInfo &info) const
{
    C_UNUSED(info);
}

void IPrinter::printMatrix(const DoubleMatrix &x, unsigned int m, unsigned int n, const char* s, FILE* f)
{
    C_UNUSED(s);
    unsigned int size = x.size();
    unsigned int M = x.size() / m;

    for (unsigned int j=0; j<size; j++)
    {
        unsigned int xj_size = x[j].size();
        unsigned int N = xj_size / n;
        if (j%M==0)
        {
            for (unsigned int i=0; i<xj_size; i++)
            {
                if (i%N==0) fprintf(f, "%14.10f ", x[j][i]);
            }
            fputs("\n", f);
        }
    }
    fflush(f);
}

void IPrinter::printAsMatrix(const DoubleVector &x, unsigned int M, unsigned int N, unsigned int m, unsigned int n, const char* s, FILE* f)
{
    C_UNUSED(s);
    for (unsigned int j=0; j<=M; j++)
    {
        if (j%(M/m)==0)
        {
            for (unsigned int i=0; i<=N; i++)
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
    if (s!='\0') fprintf(file, "%s", s);
    if (start != 0 || end != 0)
    {
        unsigned int N = (end-start+1) / n;
        for (unsigned int i=start; i<=end; i++)
        {
            if ((i-start)%N==0) fprintf(file, "%12.8f ", x[i]);
        }
    }
    else
    {
        unsigned int N = x.size() / n;
        for (unsigned int i=0; i<x.size(); i++)
        {
            if (i%N==0) fprintf(file, "%12.8f ", x[i]);
        }
    }
    fputs("\n", file);
    fflush(file);
}

void IPrinter::printVector(double *x, unsigned int size, const char *s, unsigned int n, unsigned int start, unsigned int end, FILE *file)
{
    if (s!='\0') fprintf(file, "%s", s);
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
    if (s!='\0') fprintf(file, "%s", s);
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
    time_t t = time(0);
    struct tm * now = localtime( & t );
    char buf[80];
    strftime(buf, sizeof(buf), "%Y-%m-%d %X", now);
    fprintf(file, "%s\n", buf);
}
