#include "printer.h"

void printMatrix1(const DoubleMatrix& x)
{
    int m = x.size()/10;
    for (unsigned int j=0; j<x.size(); j++)
    {
        if (j%m==0)
        {
            for (unsigned int i=0; i<x[j].size(); i++)
            {
                if (i%m==0) printf("%14.10f ", x[j][i]);
            }
            puts("");
        }
    }
}
