#include "printer.h"
#include <stdio.h>
#include <time.h>
#include <windows.h>

void IPrinter::printMatrix(const DoubleMatrix &x, unsigned int m, unsigned int n, const char* s, FILE* f)
{
    C_UNUSED(s);

    unsigned int rows = x.rows();
    unsigned int cols = x.cols();
    unsigned int M = rows / m;

    for (unsigned int j=0; j<rows; j++)
    {
        unsigned int N = cols / n;
        if (j%M==0)
        {
            for (unsigned int i=0; i<cols; i++)
            {
                //if (i%N==0) fprintf(f, "%14.10f ", x.at(j,i));
                if (i%N==0) fprintf(f, "%18.14f ", x.at(j,i));
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

    unsigned int rows = x.rows();
    unsigned int cols = x.cols();
    unsigned int M = rows / m;

    for (unsigned int j=0; j<rows; j++)
    {
        unsigned int N = cols / n;
        if (j%M==0)
        {
            for (unsigned int i=0; i<cols; i++)
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
            if ((i-start)%N==0) fprintf(file, "%14.10f ", x[i]);
        }
    }
    else
    {
        unsigned int N = x.size() / n;
        for (unsigned int i=0; i<x.size(); i++)
        {
            if (i%N==0) fprintf(file, "%14.10f ", x[i]);
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

void IPrinter::printVector(unsigned int width, unsigned int presicion, const DoubleVector &x, const char *s, unsigned int n, unsigned int start, unsigned int end, FILE *file)
{
    char format[10] = {0};
    int sz = sprintf(format, "%%%d.%df ", width, presicion);
    format[sz] = '\0';

    if (s!='\0') fprintf(file, "%s", s);
    if (start != 0 || end != 0)
    {
        unsigned int N = (end-start+1) / n;
        for (unsigned int i=start; i<=end; i++)
        {
            if ((i-start)%N==0) fprintf(file, format, x[i]);
        }
    }
    else
    {
        unsigned int N = x.size() / n;
        for (unsigned int i=0; i<x.size(); i++)
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
    time_t t = time(0);
    struct tm * now = localtime( & t );
    char buf[80];
    strftime(buf, sizeof(buf), "%Y-%m-%d %X", now);
    fprintf(file, "%s\n", buf);
}

void IPrinter::print(const DoubleMatrix &m, unsigned int M, unsigned int N, unsigned int width, unsigned int presicion, FILE *file)
{
    C_UNUSED(M);
    C_UNUSED(N);

    char format[10] = {0};
    int sz = sprintf(format, "%%%d.%df ", width, presicion);
    format[sz] = '\0';

    unsigned int rows = m.rows();
    unsigned int cols = m.cols();

    for (unsigned int i=0; i<rows; i++)
    {
        for (unsigned int j=0; j<cols; j++)
        {
            fprintf(file, format, m.at(i,j));
        }
        fputs("\n", file);
    }
    fflush(file);


//    unsigned int M = rows / m;

//    for (unsigned int j=0; j<rows; j++)
//    {
//        unsigned int N = cols / n;
//        if (j%M==0)
//        {
//            for (unsigned int i=0; i<cols; i++)
//            {
//                if (i%N==0) fprintf(f, format, x.at(j,i));
//            }
//            fputs("\n", f);
//        }
//    }
//    fflush(f);
}

void IPrinter::print(const DoubleVector &v, unsigned int N, unsigned int width, unsigned int presicion, FILE *file)
{
    C_UNUSED(N);

    char format[10] = {0};
    int sz = sprintf(format, "%%%d.%df ", width, presicion);
    format[sz] = '\0';

    unsigned int size = v.size();

    for (unsigned int i=0; i<size; i++)
    {
        fprintf(file, format, v.at(i));
    }
    fputs("\n", file);
    fflush(file);
}

void IPrinter::printSeperatorLine(const char* msg, char c, FILE* file)
{
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    int columns;
    GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
    columns = csbi.srWindow.Right - csbi.srWindow.Left + 1;
    //rows = csbi.srWindow.Bottom - csbi.srWindow.Top + 1;

    int start = 0;
    if (msg != NULL)
    {
        fprintf(file, "%s ", msg);
        start = strlen(msg)+1;
    }

    for (int i=start; i<columns; i++) fprintf(file,"%c", c);
    fputs("", file);
    fflush(file);
}
