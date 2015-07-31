#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

extern void calculate();
extern void sample_gradient1();
extern void _calculate();
extern void calculate_grid();

int main(int argc, char** argv)
{
    //sample_gradient1();
    //calculate();
    _calculate();
	//calculate_grid();
	return 0;
}
