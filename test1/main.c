#include <stdio.h>
#include "test_lib.h"

int main(int argc, char** argv)
{
	double x,y,z;
	
	double min = -10.0;
	double max = +10.0;
	for (x=min; x<max; x+=0.001)
	{
		for (y=min; y<max; y+=0.001)
		{
			for (z=min; z<max; z+=0.001)
			{
				double a1 = x*y + x*z + z*y;
				double a2 = x*x*y*y + x*x*z*z + z*z*y*y;
				
				//if ((23.0 < a1 && a1 < 25.0) && (299.0 < a2 && a2 < 301.0))				
				if ((a1 == 24.0) && (a2 < 300.0))						
					printf("%f %f %f %f %f\n", x, y, z, a1, a2);
			}
		}
	}
	puts("ok");
	
	return 0;
}