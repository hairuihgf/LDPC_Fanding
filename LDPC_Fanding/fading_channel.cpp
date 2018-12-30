#include "fading_channel.h"

void fading_channel(complex<float> *modulated_symbols, dcomplex *moded_symbols, dcomplex *received_symbols, dcomplex *SideInfo, int len_symbol, double DOPPLER, double T_SAMPLE, double RAYLEIGH_M, double sigma)
{
	//这部分是fading channel后加入的内容
	//complex转换成dcomplex
	for (int i = 0; i<len_symbol; i++) {
		moded_symbols[i].re = modulated_symbols[i].real();
		moded_symbols[i].im = modulated_symbols[i].imag();
	}
	//fading channel
	double Fd = DOPPLER;
	double Ts = T_SAMPLE;
	dcomplex temp_complex = cMakeComplex(0.0, 0.0);
	int num_of_symbol = len_symbol;
	Rayfadsim(Fd, Ts, num_of_symbol, RAYLEIGH_M, SideInfo);
	//Rayfadsim(Fd, Ts, num_of_symbol, RAYLEIGH_M, SideInfo);
	//************************************************//
	/*
	dcomplex *SideInfo2 = SideInfo;
	int num_of_symbol2 = 15000;
	dcomplex *SideInfo3 = new dcomplex[num_of_symbol2];
	Rayfadsim2(25, 0.001, num_of_symbol2, 8, SideInfo3);*/
	//Rayfadsim(25, 0.001, 15000, 8, SideInfo3);
	/*
	FILE *fading_re,*fading_im;
	char title_fading_re[300], title_fading_im[300];
	title_fading_re[0] = '\0';
	title_fading_im[0] = '\0';
	strcat_s(title_fading_re, "reylaigh_re.rcd");
	strcat_s(title_fading_im, "reylaigh_im.rcd");
	fopen_s(&fading_re, title_fading_re, "a");
	fopen_s(&fading_im, title_fading_im, "a");
	for (int i = 0; i < num_of_symbol2; i++)
	{
		fprintf_s(fading_re, "%f\t", SideInfo3[i].re);
		fprintf_s(fading_im, "%f\t", SideInfo3[i].im);

	}
	fprintf_s(fading_re, "\n\n");
	fprintf_s(fading_im, "\n\n");
	fclose(fading_re);
	fclose(fading_im);
	fading_re = NULL;
	fading_im = NULL;
	*/
	//***************************************************//
	for (int i = 0; i < num_of_symbol; i++)
	{
		/*GaussComplex(sqrt_, &rvReal, &rvImag);  //独立瑞利
		SideInfo[i].re = rvReal;
		SideInfo[i].im = rvImag;*/

		temp_complex = cMult(moded_symbols[i], SideInfo[i]);
		received_symbols[i].re = temp_complex.re;
		received_symbols[i].im = temp_complex.im;
	}

	double rvReal, rvImag;
	for (int i = 0; i < num_of_symbol; i++)
	{
		GaussComplex(sigma, &rvReal, &rvImag);
		received_symbols[i].re += rvReal;
		received_symbols[i].im += rvImag;
	}

}

void Rayfadsim(double Fd, double Ts, int Ns, int M, dcomplex *Z)
{
	int i, j, NN = 4 * M+2;
	double dop_gain, theta;

	double *doppler = new double[M];
	double *dop_update = new double[M];
	double *phi = new double[M];
	double *varphi = new double[M];
	double *state = new double[M];
	//	dcomplex *Z = new dcomplex[ Ns];

	// set up fading channel parameters
	dop_gain = sqrt(2. / M);          //  This gain differs sqrt(2), it is  from the paper to get Var(Z)=1.
	theta = (2.0 * UniformRV() - 1.0) * PI;
	for (i = 0; i < M; i++)
	{
		doppler[i] = Fd * cos(((2 * PI) / NN) * (i + 1) + (theta / NN) - (PI / NN));
		phi[i] = (2 * UniformRV() - 1) * PI;
		varphi[i] = (2 * UniformRV() - 1) * PI;
		state[i] = 0;
		dop_update[i] = 2 * PI * doppler[i] * Ts;
	}
	for (i = 0; i < Ns; i++)
	{
		Z[i] = cMakeComplex(0.0, 0.0);
	}

	//generate fading channel samples
	for (i = 0; i < Ns; i++)
	{
		for (j = 0; j < M; j++)
		{
			Z[i].re += dop_gain * cos(state[j] + phi[j]);
			Z[i].im += dop_gain * sin(state[j] + varphi[j]);
			state[j] += dop_update[j];
		}
	}


	delete[] doppler; doppler = NULL;
	delete[] dop_update; dop_update = NULL;
	delete[] phi; phi = NULL;
	delete[] varphi; varphi = NULL;
	delete[] state; state = NULL;

	//return Z;

}
//xiao-zheng paper "Simulation Models With Correct Statistical Properties for Rayleigh Fading Channels"
void Rayfadsim2(double Fd, double Ts, int Ns, int M, dcomplex *Z)
{
	int i, j, NN = 4 * M;
	double dop_gain, wd, theta, phi,psi_n, alpha_n,wn;
	//double *psi_n = new double[M];
	//double *alpha_n = new double[M];
	//double *wn = new double[M];

	for (i = 0; i < Ns; i++)
	{
		Z[i] = cMakeComplex(0.0, 0.0);
	}
	dop_gain = 2*sqrt(1./M);
	wd		= 2 * PI*Fd;
	theta	= (2 * UniformRV() - 1) * PI;
	phi		= (2 * UniformRV() - 1) * PI;
	
	for (i = 0; i < M; i++)
	{
		psi_n   = (2 * UniformRV() - 1) * PI;
		alpha_n  = (2 * PI *(i + 1) - PI + theta) / NN;
		wn = wd*cos(alpha_n);
		for (j = 0; j < Ns; j++)
		{
			Z[j].re += dop_gain*cos(psi_n)*cos(wn * j*Ts + phi);
			Z[j].im += dop_gain*sin(psi_n)*cos(wn * j*Ts + phi);
		}
	}
	//delete[] psi_n; psi_n = NULL;
	//delete[] alpha_n; alpha_n = NULL;
	//delete[] wn; wn = NULL;
}