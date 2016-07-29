#ifndef MATRIX_H
#define MATRIX_H

#include "global.h"
#include <vector>

class MINIMUMSHARED_EXPORT Vector
{
public:
    explicit Vector(unsigned int size = 0, double value = 0.0);
    virtual ~Vector();

    double L2Norm() const;
    double L1Norm() const;
    double LInfNorm() const;
    double EuclideanNorm() const;
    double EuclideanDistance(const Vector&) const;
    void L2Normalize();
    void L1Normalize();
    void EuclideanNormalize();
    double min() const;
    double max() const;
    Vector mid(unsigned int s, unsigned int e) const;

    ///////
    void assign (unsigned int n, const double& value);
    // Access element. Returns a reference to the element at position n in the vector.
    double& at (unsigned int n);
    const double& at (unsigned int n) const;
    // Access last element. Returns a reference to the last element in the vector.
    double& back();
    const double& back() const;
    // Return iterator to beginning. Returns an iterator pointing to the first element in the vector.
    double& begin();
    const double& begin() const;
    // Return size of allocated storage capacity. Returns the size of the storage space currently allocated for the vector, expressed in terms of elements.
    unsigned int capacity() const;
    // Return const_iterator to beginning. Returns a const_iterator pointing to the first element in the container.
    double& cbegin() const noexcept;
    // Return const_iterator to end. Returns a const_iterator pointing to the past-the-end element in the container.
    double& cend() const noexcept;
    // Clear content. Removes all elements from the vector (which are destroyed), leaving the container with a size of 0.
    void clear();
    // Access data. Returns a direct pointer to the memory array used internally by the vector to store its owned elements.
    double* data() noexcept;
    const double* data() const noexcept;
    // Test whether vector is empty. Returns whether the vector is empty (i.e. whether its size is 0).
    bool empty() const;
    // Return iterator to end. Returns an iterator referring to the past-the-end element in the vector container.
    double& end();
    const double& end() const;
    // Erase elements. Removes from the vector either a single element (position) or a range of elements ([first,last)).
    double& erase (unsigned int position);
    double& erase (unsigned int first, unsigned int last);
    // Access first element. Returns a reference to the first element in the vector.
    double& front();
    const double& front() const;
    // Insert elements. The vector is extended by inserting new elements before the element at the specified position,
    // effectively increasing the container size by the number of elements inserted.
    // single element
    double& insert (unsigned int position, const double& val);
    // fill
    void insert (unsigned int position, unsigned int n, const double& val);
    // range
    void insert (unsigned int position, unsigned int first, unsigned int last);
    // Return maximum size. Returns the maximum number of elements that the vector can hold.
    unsigned int max_size() const;
    //Assign content. Assigns new contents to the container, replacing its current contents, and modifying its size accordingly.
    Vector& operator= (const Vector& x);
    // Access element. Returns a reference to the element at position n in the vector container.
    double& operator[] (unsigned int n);
    double operator[] (unsigned int n) const;
    // Change size. Resizes the container so that it contains n elements.
    void resize (unsigned int size, double val = 0.0);
    // Return size. Returns the number of elements in the vector.
    unsigned int size() const;
    //Swap content. Exchanges the content of the container by the content of x, which is another vector object of the same type. Sizes may differ.
    void swap(Vector& x);

private:
    unsigned int sz;
    double *pdata;
};

class MINIMUMSHARED_EXPORT Matrix
{
public:
    Matrix(unsigned int rows=0, unsigned int cols=0, double val=0.0);
    virtual ~Matrix();

    Vector& row(unsigned int r);
    const Vector& row(unsigned int r) const;

    unsigned int rows() const;
    unsigned int cols() const;
    unsigned int size() const;

    Vector& operator[] (unsigned int j);
    const Vector& operator[] (unsigned int j) const;

    void clear();
    void resize(unsigned int rows);

    double min() const { return 0.0; }
    double max() const { return 0.0; }

private:
    unsigned int mrows;
    unsigned int mcols;
    Vector *pdata;
};

typedef Vector DoubleVector;
typedef Matrix DoubleMatrix;

class MINIMUMSHARED_EXPORT DoubleCube : public std::vector<DoubleMatrix>
{
public:
    DoubleCube();
    virtual ~DoubleCube();

    void Resize(unsigned int Nz, unsigned int Ny, unsigned Nx);
    void Clear();
};

#endif // MATRIX_H
