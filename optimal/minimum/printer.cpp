#include "printer.h"

void Printer::printMatrix(const DoubleMatrix& x, unsigned int m, unsigned int n, const char* s, FILE* f)
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

void Printer::printVector(const DoubleVector& x, unsigned int n, const char* s, FILE* f)
{
    unsigned int N = x.size() / n;
    fprintf(f, "%s", s);
    for (unsigned int i=0; i<x.size(); i++)
    {
        if (i%N==0) fprintf(f, "%14.10f ", x[i]);
    }
    fputs("\n", f);
    fflush(f);
}
