#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "methods.h"
#include "print.h"
#include "optimal.h"

extern void calculate();
extern void sample_gradient1();

int main(int argc, char** argv)
{
	//sample_gradient1();
	calculate();
	return 0;
}
