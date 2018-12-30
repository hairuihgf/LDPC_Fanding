#pragma once
#define _CRT_SECURE_NO_DEPRECATE
#include "randomc.h"
#include <math.h>
#include "struct.h"
#include "Random.h"
#define PI 3.14159
#define LLR_CLIP 100
#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))



double gauss_g(struct AWGN *gauss)
{
	double x;
	double x1;

	static double x2;
	static int x2_valid = 0;

	//extern CRandomMersenne rg;


	if (x2_valid) {
		x2_valid = 0;
		x2 *= gauss->sigma;
		return x2;
	}

	do {
		x1 = 2.0*UniformRV() - 1.0;
		x2 = 2.0*UniformRV() - 1.0;
		x = x1*x1 + x2*x2;
	} while (x >= 1.0 || x == 0);

	x1 *= sqrt((-2.0)*log(x) / x);
	x2 *= sqrt((-2.0)*log(x) / x);

	x2_valid = 1;

	x1 *= gauss->sigma;

	return x1;

}


void awgn_channel(int input[], double output[], struct AWGN *gauss, int len)
{
	int i;
	double add;
	int COLUMN;
	COLUMN = len;
	for (i = 1; i <= COLUMN; i++)
	{
		add = gauss_g(gauss);
		output[i] = input[i] + add;
	}
}