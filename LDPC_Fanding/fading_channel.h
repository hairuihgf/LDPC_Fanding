#pragma once
#include <complex>
#include "Complex.h"
#include "Random.h"
using namespace std;

void fading_channel(complex<float> *modulated_symbols, dcomplex *moded_symbols, dcomplex *received_symbols, dcomplex *SideInfo, int len_symbol, double DOPPLER, double T_SAMPLE, double RAYLEIGH_M, double sigma);
void Rayfadsim(double Fd, double Ts, int Ns, int M, dcomplex *Z);
void Rayfadsim2(double Fd, double Ts, int Ns, int M, dcomplex *Z);