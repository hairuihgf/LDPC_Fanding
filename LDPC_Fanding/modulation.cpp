#include <cmath>
#include <iostream>
#include <complex>
#define M_PI 3.1415926
using namespace std;

class modulator {
public:
	modulator();
	~modulator();
	void modulation(int *input, int bits_num, int modulate_mode, complex<float> *modulated_symbols);
};

void modulator::modulation(int *input, int bits_num, int modulate_mode, complex<float> *modulated_symbols) {
	float GAIN_BPSK2 = 1 / sqrt((float)2);
	float GAIN_QPSK4 = 1 / sqrt((float)2);
	float GAIN_QAM16 = 1 / sqrt((float)10);
	float GAIN_QAM64 = 1 / sqrt((float)42);
	float GAIN_QAM256 = 1 / sqrt((float)170);
	float RealBpsk[] = { 1, -1 };
	float ImagBpsk[] = { 1, -1 };
	float RealQpsk[] = { 1, 1, -1, -1 };
	float ImagQpsk[] = { 1, -1, 1, -1 };
	float Real16Qam[] = { 1, 1, 3, 3, 1, 1, 3, 3, -1, -1, -3, -3, -1, -1, -3, -3 };
	float Imag16Qam[] = { 1, 3, 1, 3, -1, -3, -1, -3, 1, 3, 1, 3, -1, -3, -1, -3 };
	float Real64Qam[] = { 3, 3, 1, 1, 3, 3, 1, 1, 5, 5, 7, 7, 5, 5, 7, 7, 3, 3, 1, 1, 3, 3, 1, 1, 5, 5, 7, 7, 5, 5, 7, 7, -3, -3, -1, -1, -3, -3, -1, -1, -5, -5, -7, -7, -5, -5, -7, -7, -3, -3, -1, -1, -3, -3, -1, -1, -5, -5, -7, -7, -5, -5, -7, -7 };
	float Imag64Qam[] = { 3, 1, 3, 1, 5, 7, 5, 7, 3, 1, 3, 1, 5, 7, 5, 7, -3, -1, -3, -1, -5, -7, -5, -7, -3, -1, -3, -1, -5, -7, -5, -7, 3, 1, 3, 1, 5, 7, 5, 7, 3, 1, 3, 1, 5, 7, 5, 7, -3, -1, -3, -1, -5, -7, -5, -7, -3, -1, -3, -1, -5, -7, -5, -7 };
	float A[] = { 5, 5, 7, 7, 5, 5, 7, 7, 3, 3, 1, 1, 3, 3, 1, 1, 5, 5, 7, 7, 5, 5, 7, 7, 3, 3, 1, 1, 3, 3, 1, 1, 11, 11, 9, 9, 11, 11, 9, 9, 13, 13, 15, 15, 13, 13, 15, 15, 11, 11, 9, 9, 11, 11, 9, 9, 13, 13, 15, 15, 13, 13, 15, 15 };
	float B[] = { 5, 7, 5, 7, 3, 1, 3, 1, 5, 7, 5, 7, 3, 1, 3, 1, 11, 9, 11, 9, 13, 15, 13, 15, 11, 9, 11, 9, 13, 15, 13, 15, 5, 7, 5, 7, 3, 1, 3, 1, 5, 7, 5, 7, 3, 1, 3, 1, 11, 9, 11, 9, 13, 15, 13, 15, 11, 9, 11, 9, 13, 15, 13, 15 };
	float Real256Qam[256];
	float Imag256Qam[256];
	for (int i = 0; i<64; i++) {
		Real256Qam[i] = A[i];
		Real256Qam[64 + i] = A[i];
		Real256Qam[128 + i] = -1 * A[i];
		Real256Qam[192 + i] = -1 * A[i];
		Imag256Qam[i] = B[i];
		Imag256Qam[64 + i] = -1 * B[i];
		Imag256Qam[128 + i] = B[i];
		Imag256Qam[192 + i] = -1 * B[i];
	}
	int symbol_num;
	float Real8psk[] = { cos(M_PI * 11 / 8),cos(M_PI * 9 / 8),cos(M_PI * 5 / 8),cos(M_PI * 7 / 8),cos(M_PI * 13 / 8),cos(M_PI * 15 / 8),cos(M_PI * 3 / 8),cos(M_PI * 1 / 8) };
	float Imag8psk[] = { sin(M_PI * 11 / 8),sin(M_PI * 9 / 8),sin(M_PI * 5 / 8),sin(M_PI * 7 / 8),sin(M_PI * 13 / 8),sin(M_PI * 15 / 8),sin(M_PI * 3 / 8),sin(M_PI * 1 / 8) };

	switch (modulate_mode)
	{
	case 1:
		for (int i = 0; i<bits_num; i++) {
			modulated_symbols[i].real(RealBpsk[(int)input[i]] * GAIN_BPSK2);
			modulated_symbols[i].imag(ImagBpsk[(int)input[i]] * GAIN_BPSK2);
		}
		break;
	case 2:
		symbol_num = bits_num / 2;
		for (int symbol_index = 0; symbol_index<symbol_num; symbol_index++)
		{
			int temp = (int)input[2 * symbol_index] * 2 + input[2 * symbol_index + 1];
			float C_re = RealQpsk[temp] * GAIN_QPSK4;
			float C_im = ImagQpsk[temp] * GAIN_QPSK4;
			modulated_symbols[symbol_index].real(C_re);
			modulated_symbols[symbol_index].imag(C_im);
		}
		break;
	case 3:
		symbol_num = bits_num / 3;
		for (int i = 0; i<symbol_num; i++)
		{
			int index = (int)(input[3 * i] * 4 + input[3 * i + 1] * 2 + input[3 * i + 2] * 1);
			modulated_symbols[i].real(Real8psk[i]);
			modulated_symbols[i].imag(Imag8psk[i]);
		}
		break;
	case 4:
		symbol_num = bits_num / 4;
		for (int symbol_index = 0; symbol_index<symbol_num; symbol_index++)
		{
			int temp = input[4 * symbol_index] * 8 + input[4 * symbol_index + 1] * 4 + input[4 * symbol_index + 2] * 2 + input[4 * symbol_index + 3] * 1;
			modulated_symbols[symbol_index].real(Real16Qam[temp] * GAIN_QAM16);
			modulated_symbols[symbol_index].imag(Imag16Qam[temp] * GAIN_QAM16);
		}
		break;
	case 6:
		symbol_num = bits_num / 6;
		for (int symbol_index = 0; symbol_index<symbol_num; symbol_index++)
		{
			int temp = input[6 * symbol_index] * 32 + input[6 * symbol_index + 1] * 16 + input[6 * symbol_index + 2] * 8 + input[6 * symbol_index + 3] * 4 + input[6 * symbol_index + 4] * 2 + input[6 * symbol_index + 5];
			modulated_symbols[symbol_index].real(Real64Qam[temp] * GAIN_QAM64);
			modulated_symbols[symbol_index].imag(Imag64Qam[temp] * GAIN_QAM64);
		}
		break;
	case 8:
		symbol_num = bits_num / 8;
		for (int symbol_index = 0; symbol_index<symbol_num; symbol_index++)
		{
			int temp = input[8 * symbol_index] * 128 + input[8 * symbol_index + 1] * 64 + input[8 * symbol_index + 2] * 32 + input[8 * symbol_index + 3] * 16 + input[8 * symbol_index + 4] * 8 + input[8 * symbol_index + 5] * 4 + input[8 * symbol_index + 6] * 2 + input[8 * symbol_index + 7];
			modulated_symbols[symbol_index].real(Real256Qam[temp] * GAIN_QAM256);
			modulated_symbols[symbol_index].imag(Imag256Qam[temp] * GAIN_QAM256);
		}
		break;
	}
}