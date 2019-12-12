#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef double (*funcao)(double);

typedef enum { false, true } bool;

double f(double x){
  return pow(x,3) - cos(x);
}

double brent(double f(double x), double xu, double xl, double tol, unsigned int maxit) {

    double a = xu; 
    double b = xl; 
    double fa = f(a); 
    double fb = f(b); 
    double fs = 0; 

    if (!(fa*fb < 0)) {
        printf("%s\n", "the signs must be opposites"); 
        return -1;  
    }

    if (fabs(fa) < fabs(b)) {
        double aux; 
        aux = a; a = b; b = aux;
        aux = fa; fa = fb; fb = aux;  
    }

    double c = a; 
    double fc = fa; 
    bool mflag = 1;  
    double s = 0; 
    double d = 0; 

    for (unsigned int iter = 0; iter < maxit; ++iter) {
        
        if (fabs(b-a) < tol) {
            printf("After %u itarations the root is: %lf\n", iter, s); 
            return s; 
        }

        if (fa != fc && fb != fc) {

            // use inverse quadratic interopolation
            s =	  ( a * fb * fc / ((fa - fb) * (fa - fc)) )
				+ ( b * fa * fc / ((fb - fa) * (fb - fc)) )
				+ ( c * fa * fb / ((fc - fa) * (fc - fb)) );

        }
        else
        {
            // secant method
			    s = b - fb * (b - a) / (fb - fa);
        }
        
        if (	( (s < (3 * a + b) * 0.25) || (s > b) ) ||
				( mflag && (fabs(s-b) >= (fabs(b-c) * 0.5)) ) ||
				( !mflag && (fabs(s-b) >= (fabs(c-d) * 0.5)) ) ||
				( mflag && (fabs(b-c) < tol) ) ||
				( !mflag && (fabs(c-d) < tol))	)
		{
			// bisection method
			s = (a+b)*0.5;
 
			mflag = 1;
		}
		else
		{
			mflag = 0;
		}
 
		fs = f(s);	// calculate fs
		d = c;		// first time d is being used (wasnt used on first iteration because mflag was set)
		c = b;		// set c equal to upper bound
		fc = fb;	// set f(c) = f(b)
 
		if ( fa * fs < 0)	// fa and fs have opposite signs
		{
			b = s;
			fb = fs;	// set f(b) = f(s)
		}
		else
		{
			a = s;
			fa = fs;	// set f(a) = f(s)
		}
 
		if (fabs(fa) < fabs(fb)) // if magnitude of fa is less than magnitude of fb
		{
            double aux; 
            aux = a; a = b; b = aux; // swap a and b
            aux = fa; fa = fb; fb = aux;  // make sure f(a) and f(b) are correct after swap
		}
 
	} // end for   
}




int main(void) {

  clock_t t; 
  t = clock();
  double r = brent(f, -5, 5,  0.01, 100); 
  t = clock() - t;
  printf("Tempo de execucao: %lf\n", ((double)t)/((CLOCKS_PER_SEC/1000)));

  return 0;
}