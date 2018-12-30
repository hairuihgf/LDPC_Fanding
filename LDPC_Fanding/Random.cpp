#include "Random.h"
#include <cmath>

////////////////////////
#define TINY                1.e-100

unsigned long state = 57;

int RandomInt(void)
{
	static const unsigned long M = 2147483647L;
	static const int A = 48271;
	static const int Q = M / A;
	static const int R = M % A;

	int tmpState = A * (state % Q) - R * (state / Q);
	if (tmpState >= 0)
		state = tmpState;
	else
		state = tmpState + M;

	return state;
}

double UniformRV(void)
{
	static const int A = 48271;
	static const unsigned long M = 2147483647L;
	static const int Q = M / A;
	static const int R = M % A;

	int tmpState = A * (state % Q) - R * (state / Q);
	if (tmpState >= 0)
		state = tmpState;
	else
		state = tmpState + M;

	return ((double)state / (double)M);
}

double GaussRand(double sigma)
{
	double u1, u0;

	u0 = UniformRV();
	u1 = UniformRV();
	if (u0 < TINY) u0 = TINY;
	return (sigma * sqrt(-2.0*log(u0)) * cos(2 * PI*u1));
}

void GaussComplex(double sigma, double *pReal, double *pImag)
{
	double x, y, z;

	x = UniformRV();
	if (x < TINY) x = TINY;
	y = 2.0 * PI * UniformRV();
	z = sigma * sqrt(-2.0 * log(x));
	*pReal = z * cos(y);
	*pImag = z * sin(y);
}