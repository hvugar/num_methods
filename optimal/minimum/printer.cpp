#include "printer.h"

void Printer::printMatrix(const DoubleMatrix& x, int m, int n, FILE* f, const char* s)
{
    for (unsigned int j=0; j<x.size(); j++)
    {
        if (j%m==0)
        {
            for (unsigned int i=0; i<x[j].size(); i++)
            {
                if (i%n==0) fprintf(f, "%14.10f ", x[j][i]);
            }
            fputs("\n", f);
        }
    }
    fflush(f);
}

void Printer::printVector(const DoubleVector& x, int n, FILE* f, const char* s)
{
    for (unsigned int i=0; i<x.size(); i++)
    {
        if (i%n==0) fprintf(f, "%14.10f ", x[i]);
    }
    fputs("\n", f);
    fflush(f);
}
