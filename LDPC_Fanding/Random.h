#pragma once
#ifndef _RANDOM_H
#define _RANDOM_H

#define PI                  3.1415926536

int RandomInt(void);
double UniformRV(void);
double GaussRand(double sigma);
void GaussComplex(double sigma, double *pReal, double *pImag);

#endif
