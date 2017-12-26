#include "expol.h"

ExpOptimalLetters::ExpOptimalLetters() : AbstactProblem22D()
{
    Dimension time = Dimension(0.005, 0, 200);
    Dimension dimX = Dimension(0.01, 0, 100);
    Dimension dimY = Dimension(0.01, 0, 100);
    setGridParameters(time, dimX, dimY);

    U.resize(dimY.sizeN()+1, dimX.sizeN()+1, 10.0);
}

void ExpOptimalLetters::table1()
{

}
