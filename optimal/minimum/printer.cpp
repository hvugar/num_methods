#include "printer.h"

void Printer::printMatrix(const DoubleMatrix &x, unsigned int m, unsigned int n, const char* s, FILE* f)
{
    unsigned int M = x.size() / m;

    for (unsigned int j=0; j<x.size(); j++)
    {
        unsigned int N = x[j].size() / n;
        if (j%M==0)
        {
            for (unsigned int i=0; i<x[j].size(); i++)
            {
                if (i%N==0) fprintf(f, "%14.10f ", x[j][i]);
            }
            fputs("\n", f);
        }
    }
    fflush(f);
}

void Printer::printVector(const DoubleVector& x, unsigned int n, const char *s, FILE *f)
{
    unsigned int N = x.size() / n;
    if (s!=NULL)
    {
        fprintf(f, "%s", s);
    }
    for (unsigned int i=0; i<x.size(); i++)
    {
        if (i%N==0) fprintf(f, "%14.10f ", x[i]);
    }
    fputs("\n", f);
    fflush(f);
}

void Printer::printAsMatrix(const DoubleVector &x, unsigned int M, unsigned int N, unsigned int m, unsigned int n, const char* s, FILE* f)
{
    for (unsigned int j=0; j<=M; j++)
    {
        if (j%(M/m)==0)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                if (i%(N/n)==0) fprintf(f, "%14.10f ", x[j*(N+1)+i]);
            }
            fputs("\n", f);
        }
    }
    fflush(f);
}

void Printer::printVector(const DoubleVector &x, const char *s, unsigned int n, unsigned int start, unsigned int end, FILE *file)
{
    if (s!='\0') fprintf(file, "%s", s);
    if (start != 0 || end != 0)
    {
        unsigned int N = (end-start+1) / n;
        for (unsigned int i=start; i<=end; i++)
        {
            if ((i-start)%N==0) fprintf(file, "%14.8f ", x[i]);
        }
    }
    else
    {
        unsigned int N = x.size() / n;
        for (unsigned int i=0; i<x.size(); i++)
        {
            if (i%N==0) fprintf(file, "%14.8f ", x[i]);
        }
    }
    fputs("\n", file);
    fflush(file);
}
