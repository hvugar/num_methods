#ifndef R1FUNCTIONVECTOR_H
#define R1FUNCTIONVECTOR_H

#include <function.h>
#include <vector>
#include <new>

using namespace std;

class R1FunctionVector
{
public:
    R1FunctionVector(unsigned int size);
    virtual ~R1FunctionVector();

    R1Function* at(unsigned int) const;
    unsigned int size() const;
private:
    R1Function* *pfx;
    unsigned int sz;
};

//class R1FunctionMatrix
//{
//public:
//    R1FunctionMatrix(unsigned int row, unsigned int col);
//    virtual ~R1FunctionMatrix();

//    R1Function* at(unsigned int row, unsigned int col) const;
//    unsigned int rows() const;
//    unsigned int cols() const;
//private:
//    R1Function* **pfx;
//    unsigned int mRows;
//    unsigned int mCols;
//};

#endif // R1FUNCTIONVECTOR_H
