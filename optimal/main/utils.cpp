#include "utils.h"

void printX(const char* s, const std::vector<double>& x)
{
    unsigned int i;
    unsigned int n = x.size();
    printf("%s =\t{", s);
    for (i=0; i<n; i++)
    {
        if ( i%((n-1)/10) == 0 )
        {
            if (x[i] < 0)
            {
                printf("%12.8f", x[i]);
            }
            else
            {
                printf("%+12.8f", x[i]);
            }
        }
        if ( i%((n-1)/10) == 0 && i != n-1 )
        {
            printf(", ");
        }
    }
    printf("};");
    printf("\n");
}
