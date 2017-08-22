#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

struct Vector
{
    unsigned int size;
    double* items;
} ;

struct Matrix
{
    unsigned int rows;
    unsigned int columns;
    double** items;
};

struct Vector* vector_new(unsigned int size);
void vector_free(struct Vector *vector);
double vector_L2Norm(struct Vector *v);
void vector_Normalize(struct Vector *v);
double vector_max(struct Vector *v);
double vector_min(struct Vector *v);
double vector_maxmin(const struct Vector *v);


struct Matrix* matrix_new(unsigned int rows, unsigned int columns);
void matrix_free(struct Matrix *matrix);
double matrix_det(struct Matrix *m);
struct Matrix* matrix_transpose(struct Matrix *m);
struct Matrix* matrix_inverse(struct Matrix *m);
struct Matrix* matrix_minor(struct Matrix *m, unsigned int row, unsigned int column);
struct Matrix* matrix_mult(struct Matrix *m1, struct Matrix *m2);
struct Matrix* matrix_mult_v1(struct Matrix *m, struct Vector *v);
struct Matrix* matrix_mult_v2(struct Vector *v, struct Matrix *m);

//void matrix_transpose(struct Matrix *m, struct Matrix *transpose);
//void matrix_inverse(struct Matrix *m, struct Matrix *inverse);
//void matrix_minor(struct Matrix *m, unsigned int row, unsigned int column, struct Matrix *minor);
//void matrix_mult(struct Matrix *m1, struct Matrix *m2, struct Matrix *mult);
void matrix_print(struct Matrix *m);
void matrix_rand(struct Matrix *m);

#ifdef __cplusplus
}
#endif

#endif // MATRIX_H
