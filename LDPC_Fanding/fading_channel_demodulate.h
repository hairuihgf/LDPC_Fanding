#pragma once
#include "Complex.h"
#include <cmath>
#include <stdio.h>
#include <stdlib.h>

void fading_channel_demodulate(dcomplex* modulate, double* demodulated, int mod, double sigma, int moded_len, dcomplex *SideInfo);
void CreateSignalConstellation(int N_BITS, int *SignalMapper, dcomplex *SignalSet, int ConstellationSize);
void NormalizeSignalPower(dcomplex *SignalSet, int ConstellationSize);
void Normalize(double *data, int num);
double clip(double x);