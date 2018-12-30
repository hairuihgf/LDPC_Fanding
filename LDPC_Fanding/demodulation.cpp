#include <complex>
#include <cmath>
#include <algorithm>
using namespace std;
#define M_PI 3.14159

class demodulator {
public:
	demodulator();
	~demodulator();
	void demodulation(complex<float> *modulated_symbols, int symbol_length, int modulation_mode, float *soft_bits_out, float *hard_decision);
};

void demodulator::demodulation(complex<float> *modulated_symbols, int symbol_length, int modulate_mode, float *soft_bits_out, float *hard_decision)
{
	float *C_re = new float[symbol_length];
	float *C_im = new float[symbol_length];
	//float *soft_bits_out = new float[symbol_length*modulation_mode];
	int symbol_row = 1;
	for (int i = 0; i<symbol_length; i++)
	{
		C_re[i] = modulated_symbols[i].real();
		C_im[i] = modulated_symbols[i].imag();
	}
	float D;
	switch (modulate_mode)
	{
	case 1:
	{
		D = 1.0;
		complex<float> pi34(cos((float)M_PI * 3 / 4), sin((float)M_PI * 3 / 4));
		for (int i = 0; i<symbol_length; i++)
		{
			complex<float> temp;
			temp = modulated_symbols[i] * pi34;
			float data_rotate = temp.real();
			soft_bits_out[i] = data_rotate*D * 4;
		}
		break;
	}
	case 2:
	{
		D = 1 / sqrt(2.0);
		for (int i = 0; i<symbol_length; i++)
		{
			soft_bits_out[2 * i] = -4 * D*C_re[i];
			soft_bits_out[2 * i + 1] = -4 * D*C_im[i];
		}
		break;
	}
	case 3:
	{
		float modu_data[3][8];
		for (int i = 0; i<8; i++)
		{
			if (i<4) modu_data[0][i] = 0;
			else modu_data[0][i] = 1;
			if (i / 2 % 2 == 0) modu_data[1][i] = 0;
			else modu_data[1][i] = 1;
			if (i % 2 == 0) modu_data[2][i] = 0;
			else modu_data[2][i] = 1;
		}
		int pos01[4] = { 4,5,6,7 };
		int pos00[4] = { 0,1,2,3 };
		int pos11[4] = { 2,3,6,7 };
		int pos10[4] = { 0,1,4,5 };
		int pos21[4] = { 1,3,5,7 };
		int pos20[4] = { 0,2,4,6 };
		complex<float> j(0, 1);
		complex<float> modu_cons[8];
		modu_cons[0].real(cos(11 / 8 * M_PI));
		modu_cons[0].imag(sin(11 / 8 * M_PI));
		modu_cons[1].real(cos(9 / 8 * M_PI));
		modu_cons[1].imag(sin(9 / 8 * M_PI));
		modu_cons[2].real(cos(5 / 8 * M_PI));
		modu_cons[2].imag(sin(5 / 8 * M_PI));
		modu_cons[3].real(cos(7 / 8 * M_PI));
		modu_cons[3].imag(sin(7 / 8 * M_PI));
		modu_cons[4].real(cos(13 / 8 * M_PI));
		modu_cons[4].imag(sin(13 / 8 * M_PI));
		modu_cons[5].real(cos(15 / 8 * M_PI));
		modu_cons[5].imag(sin(15 / 8 * M_PI));
		modu_cons[6].real(cos(3 / 8 * M_PI));
		modu_cons[6].imag(sin(3 / 8 * M_PI));
		modu_cons[7].real(cos(1 / 8 * M_PI));
		modu_cons[7].imag(sin(1 / 8 * M_PI));
		complex<float> pi18;
		pi18.real(cos(M_PI / 8));
		pi18.imag(sin(M_PI / 8));
		for (int i = 0; i<8; i++) {
			modu_cons[i] *= pi18;
		}
		float *bit0 = new float[symbol_length];
		float *bit1 = new float[symbol_length];
		float *bit2 = new float[symbol_length];
		for (int m = 0; m<symbol_length; m++)
		{
			complex<float> soft_sym(modulated_symbols[m]);
			float d[8];
			for (int l = 0; l<8; l++)
			{
				d[l] = abs(soft_sym - modu_cons[l])*abs(soft_sym - modu_cons[l]);
			}
			float min_pos00 = d[pos00[0]];
			float min_pos01 = d[pos01[0]];
			float min_pos10 = d[pos10[0]];
			float min_pos11 = d[pos11[0]];
			float min_pos20 = d[pos20[0]];
			float min_pos21 = d[pos21[0]];
			for (int l = 1; l<4; l++)
			{
				min_pos00 = min(min_pos00, d[pos00[l]]);
				min_pos01 = min(min_pos01, d[pos01[l]]);
				min_pos10 = min(min_pos10, d[pos10[l]]);
				min_pos11 = min(min_pos11, d[pos11[l]]);
				min_pos20 = min(min_pos20, d[pos20[l]]);
				min_pos21 = min(min_pos21, d[pos21[l]]);
			}
			bit0[m] = min_pos00 - min_pos01;
			bit1[m] = min_pos10 - min_pos11;
			bit2[m] = min_pos20 - min_pos21;
		}
		for (int i = 0; i<symbol_length; i++)
		{
			soft_bits_out[3 * i] = bit0[i];
			soft_bits_out[3 * i + 1] = bit1[i];
			soft_bits_out[3 * i + 2] = bit2[i];
		}
		delete[] bit0;
		delete[] bit1;
		delete[] bit2;
		break;
	}
	case 4:
	{
		D = 1 / sqrt(10.0);
		for (int i = 0; i<symbol_length; i++)
		{
			soft_bits_out[4 * i] = -4 * D*C_re[i];
			soft_bits_out[4 * i + 1] = -4 * D*C_im[i];
			soft_bits_out[4 * i + 2] = 4 * D*(abs(C_re[i]) - 2 * D);
			soft_bits_out[4 * i + 3] = 4 * D*(abs(C_im[i]) - 2 * D);
		}
		break;
	}
	case 6:
	{
		D = 1 / sqrt(42.0);
		for (int i = 0; i<symbol_length; i++)
		{
			soft_bits_out[6 * i] = -4 * D*C_re[i];
			soft_bits_out[6 * i + 1] = -4 * D*C_im[i];
			soft_bits_out[6 * i + 2] = 4 * D*(abs(C_re[i]) - 4 * D);
			soft_bits_out[6 * i + 3] = 4 * D*(abs(C_im[i]) - 4 * D);
			soft_bits_out[6 * i + 4] = -4 * D*(2 * D - abs(abs(C_re[i]) - 4 * D));
			soft_bits_out[6 * i + 5] = -4 * D*(2 * D - abs(abs(C_im[i]) - 4 * D));
		}
		break;
	}
	case 8:
	{
		D = 1 / sqrt(170.0);
		for (int i = 0; i<symbol_length; i++)
		{
			soft_bits_out[8 * i] = -4 * D*C_re[i];
			soft_bits_out[8 * i + 1] = -4 * D*C_im[i];
			soft_bits_out[8 * i + 2] = 4 * D*(abs(C_re[i]) - 8 * D);
			soft_bits_out[8 * i + 3] = 4 * D*(abs(C_im[i]) - 8 * D);
			soft_bits_out[8 * i + 4] = -4 * D*(4 * D - abs(abs(C_re[i]) - 8 * D));
			soft_bits_out[8 * i + 5] = -4 * D*(4 * D - abs(abs(C_im[i]) - 8 * D));
			soft_bits_out[8 * i + 6] = -4 * D*(2 * D - abs(4 * D - abs(abs(C_re[i]) - 8 * D)));
			soft_bits_out[8 * i + 7] = -4 * D*(2 * D - abs(4 * D - abs(abs(C_im[i]) - 8 * D)));
		}
		break;
	}
	}
	for (int i = 0; i<symbol_length*modulate_mode; i++)
	{
		if (soft_bits_out[i] >= 0)
		{
			hard_decision[i] = 1;
		}
		else
		{
			hard_decision[i] = 0;
		}
	}
	delete[] C_re;
	delete[] C_im;
}