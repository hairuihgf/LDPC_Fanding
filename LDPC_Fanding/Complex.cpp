#include "Complex.h"
#include <cmath>


double cDistance2(dcomplex& a, dcomplex& b)
{
	double c, d;

	c = a.re - b.re;
	d = a.im - b.im;
	return(c*c + d*d);
}
double cDistance1(dcomplex &a, dcomplex &b)
{
	double c, d;
	c = fabs(a.re - b.re);
	d = fabs(a.im - b.im);
	return (c + d);
}

dcomplex cAdd(dcomplex& a, dcomplex& b)
{
	dcomplex result;

	result.re = a.re + b.re;
	result.im = a.im + b.im;
	return result;
}

dcomplex cSub(dcomplex& a, dcomplex& b)
{
	dcomplex result;

	result.re = a.re - b.re;
	result.im = a.im - b.im;
	return result;
}

/* Form the product of two complex numbers a b- */
dcomplex cMultConj(dcomplex& a, dcomplex& b)
{
	dcomplex result;

	result.re = a.re*b.re + a.im*b.im;
	result.im = a.im*b.re - a.re*b.im;
	return result;
}

/* Form the product of two complex numbers */
dcomplex cMult(dcomplex& a, dcomplex& b)
{
	dcomplex result;

	result.re = a.re*b.re - a.im*b.im;
	result.im = a.re*b.im + a.im*b.re;
	return result;
}

/* Form the quotient of two complex numbers */
dcomplex cDiv(dcomplex& a, dcomplex& b)
{
	dcomplex c;
	double r, den;

	if (fabs(b.re) >= fabs(b.im)) {
		r = b.im / b.re;
		den = b.re + r*b.im;
		c.re = (a.re + r*a.im) / den;
		c.im = (a.im - r*a.re) / den;
	}
	else {
		r = b.re / b.im;
		den = b.im + r*b.re;
		c.re = (a.re*r + a.im) / den;
		c.im = (a.im*r - a.re) / den;
	}

	return c;
}

dcomplex cMakeComplex(double re, double im)
{
	dcomplex result;

	result.re = re;
	result.im = im;
	return result;
}

/* Compute the modulus of a complex number */
double cNorm(dcomplex& z)
{
	return (z.re*z.re + z.im*z.im);
}

dcomplex cConjugate(dcomplex& z)
{
	dcomplex result;

	result.re = z.re;
	result.im = -z.im;
	return result;
}

/* Form the product of a real number and a complex number */
dcomplex cRealMult(double x, dcomplex& a)
{
	dcomplex result;

	result.re = x * a.re;
	result.im = x * a.im;
	return result;
}
