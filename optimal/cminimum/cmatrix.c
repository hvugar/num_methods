#include "cmatrix.h"
#include <float.h>

struct Vector* vector_new(unsigned int size)
{
    struct Vector* vector = (struct Vector*) malloc(sizeof(struct Vector));
    vector->size = size;
    vector->items = (double*)malloc(sizeof(double)*size);
    return vector;
}

void vector_free(struct Vector *vector)
{
    free(vector->items);
    vector->size = 0;
    vector->items = NULL;
    free(vector);
}

struct Matrix* matrix_new(unsigned int rows, unsigned int columns)
{
    struct Matrix *matrix = (struct Matrix*)malloc(sizeof(struct Matrix));
    matrix->rows = rows;
    matrix->columns = columns;
    matrix->items = (double**) malloc(sizeof(double*)*rows);
    unsigned int i;
    for (i=0; i<rows; i++) matrix->items[i] = (double*) malloc(sizeof(double)*columns);
    return matrix;
}

double vector_L2Norm(struct Vector *v)
{
    double norm = 0.0;
    unsigned int i;
    for (i=0; i<v->size; i++) norm += v->items[i]*v->items[i];
    return sqrt(norm);
}

void vector_Normalize(struct Vector *v)
{
    double norm = vector_L2Norm(v);
    unsigned int i;
    for (i=0; i<v->size; i++) v->items[i] = v->items[i]/norm;
}

double vector_max(struct Vector *v)
{
    double max = DBL_MIN;
    unsigned int i;
    for (i=0; i<v->size; i++) if (max < v->items[i]) max = v->items[i];
    return max;
}

double vector_min(struct Vector *v)\
{
    double min = DBL_MAX;
    unsigned int i;
    for (i=0; i<v->size; i++) if (min > v->items[i]) min = v->items[i];
    return min;
}

double vector_maxmin(const struct Vector *v)
{
    double min = DBL_MAX;
    double max = DBL_MIN;
    unsigned int i;
    for (i=0; i<v->size; i++)
    {
        double item = v->items[i];
        if (min > item) min = item;
        if (max < item) max = item;
        v->items[i] = 0.0;
    }
    return max;
}


void matrix_free(struct Matrix *matrix)
{
    unsigned int i;
    unsigned int rows = matrix->rows;
    for (i=0; i<rows; i++) free(matrix->items[i]);
    free(matrix->items);
    matrix->rows = 0;
    matrix->columns = 0;
    free(matrix);
}

/**
 * @brief Детерминант матрицы
 * @param m
 * @return
 */
double matrix_det(struct Matrix *m)
{
    double det = 0.0;

    if (m->rows != m->columns)
        return NAN;

    if (m->rows == 1 && m->columns == 1)
        return m->items[0][0];
    else if (m->rows == 2 && m->columns == 2)
        return m->items[0][0] * m->items[1][1] - m->items[0][1] * m->items[1][0];
    else
    {
        int i;
        for (i=0; i<m->columns; i++)
        {
            struct Matrix *mnr = matrix_minor(m, 0, i);
            det += (i%2==0 ? +1 : -1) * m->items[0][i] * matrix_det(mnr);
            matrix_free(mnr);
        }
    }

    return det;
}

/**
 * @brief Транспонированная матрица
 * @param m
 * @return
 */
struct Matrix* matrix_transpose(struct Matrix *m)
{
    struct Matrix *transpose = matrix_new(m->columns, m->rows);
    int i,j;
    for (i=0; i<m->rows; i++)
    {
        for (j=0; j<m->columns; j++)
        {
            transpose->items[j][i] = m->items[i][j];
        }
    }
    return transpose;
}

struct Matrix* matrix_inverse(struct Matrix *m)
{
    struct Matrix *inverse = matrix_new(m->columns, m->rows);
    double det = matrix_det(m);
    if (det != 0.0)
    {
        if (m->rows == 1 && m->columns == 1)
        {
            inverse->items[0][0] = 1.0 / det;
        }
        else
        {
            int j,i;
            for (j=0; j<m->rows; j++)
            {
                for (i=0; i<m->columns; i++)
                {
                    struct Matrix *mnr = matrix_minor(m, j, i);
                    double adj = matrix_det(mnr);
                    inverse->items[i][j] = ((i+j)%2==0 ? +1 : -1) * (adj/det);
                    matrix_free(mnr);
                }
            }
        }
    }
    return inverse;
}

/**
 * @brief Дополнительный минор
 * @param m
 * @param row
 * @param column
 * @return
 */
struct Matrix* matrix_minor(struct Matrix *m, unsigned int row, unsigned int column)
{
    struct Matrix *minor = matrix_new(m->rows-1, m->columns-1);
    int i,j;
    for (j=0; j<m->rows; j++)
    {
        for (i=0; i<m->columns; i++)
        {
            if (j == row || i == column) continue;
            if (j < row  && i <  column) { minor->items[j][i] = m->items[j][i];     continue; }
            if (j < row  && i >  column) { minor->items[j][i-1] = m->items[j][i];   continue; }
            if (j > row  && i <  column) { minor->items[j-1][i] = m->items[j][i];   continue; }
            if (j > row  && i >  column) { minor->items[j-1][i-1] = m->items[j][i]; continue; }
        }
    }
    return minor;
}

struct Matrix* matrix_mult(struct Matrix *m1, struct Matrix *m2)
{
    if (m1->columns != m2->rows) return NULL;

    struct Matrix *mult = matrix_new(m1->rows, m2->columns);
    int i,j,k;
    for (j=0; j<mult->rows; j++)
    {
        for (i=0; i<mult->columns; i++)
        {
            double sum = 0.0;
            for (k=0; k<m1->columns; k++) sum += m1->items[j][k]*m2->items[k][i];
            mult->items[j][i] = sum;
        }
    }
    return mult;
}

void matrix_print(struct Matrix *m)
{
    int i,j;
    for (j=0; j<m->rows; j++)
    {
        for (i=0; i<m->columns; i++)
        {
            printf("%14.6f", m->items[j][i]);
        }
        printf("\n");
    }
}

void matrix_rand(struct Matrix *m)
{
    int i,j;
    for (j=0; j<m->rows; j++)
    {
        for (i=0; i<m->columns; i++)
        {
            m->items[j][i] = rand() % 10;
        }
    }
}
