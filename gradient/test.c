double f(double *x, int n);
double f1(double x);

double f(double *x, int n)
{
    return pow(x[0],3) + 2*pow(x[1],2) - 3*x[0] - 4*x[1];
}

double f1(double x)
{
	return (x-0.133)*(x-0.133)+0.2;
}

double straight_line_search_metod1(R1Function fx, double x0, double dx, double *a, double *b, int *count)
{
	double y0 = 0.0;
	double y1 = 0.0;
	double y2 = 0.0;
	
	// if at next and last point of x0 function is greater
	// then decrease the dx to half
	y0 = fx( x0 );
	(*count)++;
	y1 = fx( x0 - dx );
	(*count)++;
	y2 = fx( x0 + dx );
	(*count)++;
	
	while (y1 > y0 && y2 > y0)
	{
		dx = dx / 2.0;
		y1 = fx( x0 - dx );
		y2 = fx( x0 + dx );
		(*count)++;
		(*count)++;
	}
	
	if ( y2 > y0 )
	{
		dx = -1 * dx;
		y2 = fx( x0 + dx );
		(*count)++;

	}

    double x1 = x0;
    double x2 = x1 + dx;
	
    y1 = fx( x1 );
	(*count)++;
    y2 = fx( x2 );
	(*count)++;
	
    int i = 1;
    while ( y2 <= y1 )
    {
		x1 = x2;
        y1 = y2;
		
		i++;
        x2 = x0 + i * dx;
        y2 = fx(x2);

		(*count)++;
		
		dx = dx*2;
    }
	
	if ( dx < 0 )
	{
		*a = x0 + (i) * dx;
		*b = x0 + (i-2) * dx;	
	}
	
	
	if ( dx > 0 )
	{
	    *a = x0 + (i-2) * dx;
		*b = x0 + (i+0) * dx;
	}
	
	return (*a+*b)/2;
}