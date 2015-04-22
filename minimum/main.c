#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "methods.h"
#include "print.h"
#include "optimal.h"

extern void calculate();
extern void __calculate();

int main(int argc, char** argv)
{
	__calculate();
	//calculate();
	return 0;
}
