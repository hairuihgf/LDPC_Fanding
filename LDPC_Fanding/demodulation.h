#pragma once
#include <complex>
#include <cmath>
#include <algorithm>

class demodulator {
public:
	demodulator() {};
	~demodulator() {};
	void demodulation(complex<float> *modulated_symbols, int symbol_length, int modulation_mode, float *soft_bits, float *hard_decision);
};