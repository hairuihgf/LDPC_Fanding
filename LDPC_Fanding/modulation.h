#pragma once
#include <cmath>
#include <iostream>
#include <complex>
using namespace std;

class modulator {
public:
	modulator() {};
	~modulator() {};
	void modulation(int *input, int bits_num, int modulate_mode, complex<float> *modulated_symbols);
};