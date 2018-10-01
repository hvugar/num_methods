#ifndef PRINTER_H
#define PRINTER_H

#include "global.h"
#include "function.h"
#include "gradient.h"
#include "matrix2d.h"

class MINIMUMSHARED_EXPORT IPrinter
{
public:
    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientMethod::MethodResult result) const = 0;

public:

    static void printVector(const DoubleVector &x, const char *s = NULL, unsigned int n = 10, unsigned int start = 0, unsigned int end = 0, FILE *file = stdout);
    static void printVector(double *x, unsigned int size, const char *s = NULL, unsigned int n = 10, unsigned int start = 0, unsigned int end = 0, FILE *file = stdout);
    static void printVector(double *x, unsigned int size, const char *s, unsigned int n, unsigned int start, unsigned int end, const char *filename);
    static void printVector(unsigned int width, unsigned int presicion, const DoubleVector &x, const char *s = NULL, unsigned int n = 10, unsigned int start = 0, unsigned int end = 0, FILE *file = stdout);
    static void printMatrix(const DoubleMatrix &x, unsigned int m = 10, unsigned int n = 10, const char* s = NULL, FILE* f = stdout);
    static void printMatrix(unsigned int width, unsigned int presicion, const DoubleMatrix &x, unsigned int m = 10, unsigned int n = 10, const char* s = NULL, FILE* f = stdout);
    static void printAsMatrix(const DoubleVector &x, unsigned int M, unsigned int N, unsigned int m = 10, unsigned int n = 10, const char* s = NULL, FILE* f = stdout);
    static void print(const DoubleMatrix &m, unsigned int M=10, unsigned int N=10, unsigned int width=14, unsigned int presicion=10, FILE *file=stdout);
    static void print(const DoubleVector &v, unsigned int N=10, unsigned int width=14, unsigned int presicion=10, FILE *file=stdout);
    static void printCube(const DoubleVector& x, unsigned int M, unsigned int N2, unsigned int N1, FILE *file = stdout);
    static void printDateTime(FILE *file = stdout);
    static void printSeperatorLine(const char* msg = NULL, char c='-', FILE* file=stdout);
};

#endif // PRINTER_H
