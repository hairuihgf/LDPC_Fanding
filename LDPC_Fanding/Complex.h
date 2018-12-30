#pragma once
#ifndef _COMPLEX_H
#define _COMPLEX_H

#define  CZERO  cMakeComplex(0.0, 0.0)

typedef struct DCOMPLEX {
	double re;
	double im;
} dcomplex;

double cDistance2(dcomplex& a, dcomplex& b);
double cDistance1(dcomplex& a, dcomplex& b);
dcomplex cAdd(dcomplex& a, dcomplex& b);
dcomplex cSub(dcomplex& a, dcomplex& b);
dcomplex cMultConj(dcomplex& a, dcomplex& b);
dcomplex cMult(dcomplex& a, dcomplex& b);
dcomplex cDiv(dcomplex& a, dcomplex& b);
dcomplex cMakeComplex(double re, double im);
double   cNorm(dcomplex& z);
dcomplex cConjugate(dcomplex& z);
dcomplex cRealMult(double x, dcomplex& a);

#endif

