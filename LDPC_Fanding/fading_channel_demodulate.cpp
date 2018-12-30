#include "fading_channel_demodulate.h"

void fading_channel_demodulate(dcomplex* modulate, double* demodulated, int mod, double sigma, int moded_len, dcomplex *SideInfo)
{
	int N_BITS = mod;
	int num_of_symbol = moded_len;
	int M_ARY = 1 << N_BITS;
	int ConstellationSize = M_ARY;
	int i, k, ii, j, index, label;

	//double likelihood[M_ARY];      //likelihood[]长度为星座点数目
	double *likelihood = new double[M_ARY];

	//double TEMP_received_symbol[num_of_symbol][N_BITS];
	//double *TEMP_received_symbol_total = new double[num_of_symbol*N_BITS];
	double *TEMP_received_symbol_total, **TEMP_received_symbol;
	TEMP_received_symbol_total=(double *)malloc(sizeof(double)*(num_of_symbol*N_BITS));
	TEMP_received_symbol = (double **)malloc(sizeof(double*)*(num_of_symbol));
	//double **TEMP_received_symbol = new double*[num_of_symbol];
	for (int i = 0; i<num_of_symbol; i++)
//		TEMP_received_symbol[i] = TEMP_received_symbol_total;
		TEMP_received_symbol[i] = TEMP_received_symbol_total + i*N_BITS;

	int modulation_type = N_BITS;
	int constellation_size = ConstellationSize;
	double pb[2];
	double Lc = 2 / sigma / sigma;
	double n0 = 2 * sigma *sigma;
	double n0_var;

	dcomplex x, y, h;//
	int *SignalMapper = new int[M_ARY];
	dcomplex *SignalSet = new dcomplex[M_ARY];
	CreateSignalConstellation(N_BITS, SignalMapper, SignalSet, ConstellationSize);


	for (int i = 0; i < num_of_symbol; i++)
	{
		y = modulate[i];

		//if(CHANNELTYPE==RAYLEIGH){
		//ChannelBlockDeinterleaver(SideInfo, Interleaver_size);
		h.re = SideInfo[i].re;
		h.im = SideInfo[i].im;
		//}


		//get symbol Probability
		double minDistance = 1.e100;
		for (label = 0; label < ConstellationSize; label++)
		{
			index = SignalMapper[label];
			x = SignalSet[index];
			//-----
			//if(CHANNELTYPE==RAYLEIGH)
			x = cMult(x, h);


			likelihood[label] = cDistance2(x, y);///欧式距离
			if (likelihood[label] < minDistance) minDistance = likelihood[label];

		}

		for (label = 0; label < ConstellationSize; label++) {
			likelihood[label] = likelihood[label] - minDistance;

			likelihood[label] = exp(-likelihood[label] / (n0));
			//if (likelihood[label] < 1.e-20) likelihood[label] = 1.e-20;
		}

		Normalize(likelihood, ConstellationSize);


		for (k = 0; k < N_BITS; k++) {
			pb[0] = pb[1] = 0.0;
			for (label = 0; label < ConstellationSize; label++) {
				j = (label >> k) & 0x01;
				pb[j] += likelihood[label];
			}
			//if (pb[0] < 1.e-30) pb[0] = 1.e-30;
			//if (pb[1] < 1.e-30) pb[1] = 1.e-30; 
			TEMP_received_symbol[i][k] = pb[1] / pb[0]; /////////////////////////////////////////////////////////////////////////
			//if (TEMP_received_symbol[i][k] < 1.e-30) TEMP_received_symbol[i][k] = 1.e-30;
			//if (TEMP_received_symbol[i][k] < 1.e-30) TEMP_received_symbol[i][k] = 1.e-30;
		}
	}

	k = 0;//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for (ii = 0; ii < num_of_symbol; ii++) {///////////////////////////////////////////////
		for (j = N_BITS - 1; j >= 0; j--) {
			//for (j = 0;j < N_BITS;j++){
			demodulated[k] = log(TEMP_received_symbol[ii][j]);
			k++;
		}
	}///////////////////////////////////////////////////////////////////////////////////////////////

	 //	    BCDToBinary(label, N_BITS, &es_bit_tr[i * N_BITS]);

	 /*	demodulated[i * N_BITS] =	( likelihood[2] + likelihood[3] ) / ( likelihood[0] + likelihood[1] );
	 demodulated[i * N_BITS] = log(demodulated[i * N_BITS]);

	 demodulated[i * N_BITS + 1] = ( likelihood[1] + likelihood[3] ) / ( likelihood[0] + likelihood[2] );
	 demodulated[i * N_BITS + 1] = log(demodulated[i * N_BITS + 1] );
	 */
	delete[] likelihood;
	delete[] TEMP_received_symbol_total;
	delete[] TEMP_received_symbol;
	delete[] SignalMapper;
	delete[] SignalSet;



}

void CreateSignalConstellation(int N_BITS, int *SignalMapper, dcomplex *SignalSet, int ConstellationSize)

{
	if (N_BITS == 1) {
		SignalMapper[0] = 0;    SignalMapper[1] = 1;

		SignalSet[0].re = 1;    SignalSet[0].im = 1;
		SignalSet[1].re = -1;    SignalSet[1].im = -1;

	}
	if (N_BITS == 2) {//同 QPSK
		SignalMapper[0] = 0;    SignalMapper[2] = 2;
		SignalMapper[1] = 1;    SignalMapper[3] = 3;

		SignalSet[0].re = 1;    SignalSet[0].im = 1;
		SignalSet[1].re = 1;    SignalSet[1].im = -1;
		SignalSet[2].re = -1;    SignalSet[2].im = 1;
		SignalSet[3].re = -1;    SignalSet[3].im = -1;

	}
	if (N_BITS == 3) {
		SignalMapper[0] = 0;    SignalMapper[1] = 1;
		SignalMapper[2] = 3;    SignalMapper[3] = 2;
		SignalMapper[4] = 6;    SignalMapper[5] = 7;
		SignalMapper[6] = 5;    SignalMapper[7] = 4;

		SignalSet[0].re = (double)sqrt((double)2);    SignalSet[0].im = 0;
		SignalSet[1].re = 1;          SignalSet[1].im = 1;
		SignalSet[2].re = 0;          SignalSet[2].im = (double)sqrt((double)2);
		SignalSet[3].re = -1;         SignalSet[3].im = 1;
		SignalSet[4].re = -(double)sqrt((double)2);   SignalSet[4].im = 0;
		SignalSet[5].re = -1;         SignalSet[5].im = -1;
		SignalSet[6].re = 0;          SignalSet[6].im = -(double)sqrt((double)2);
		SignalSet[7].re = 1;          SignalSet[7].im = -1;
	}
	if (N_BITS == 4) {
		/*SignalMapper[0] =   0;    SignalMapper[2] =   1;
		SignalMapper[10] =  2;    SignalMapper[8] =   3;
		SignalMapper[1] =   4;    SignalMapper[3] =   5;
		SignalMapper[11] =  6;    SignalMapper[9] =   7;
		SignalMapper[5] =   8;    SignalMapper[7] =   9;
		SignalMapper[15] = 10;    SignalMapper[13] = 11;
		SignalMapper[4] =  12;    SignalMapper[6] =  13;
		SignalMapper[14] = 14;    SignalMapper[12] = 15;

		SignalSet[ 0].re = 3;    SignalSet[ 0].im = 3;
		SignalSet[ 1].re = 1;    SignalSet[ 1].im = 3;
		SignalSet[ 2].re =-1;    SignalSet[ 2].im = 3;
		SignalSet[ 3].re =-3;    SignalSet[ 3].im = 3;
		SignalSet[ 4].re = 3;    SignalSet[ 4].im = 1;
		SignalSet[ 5].re = 1;    SignalSet[ 5].im = 1;
		SignalSet[ 6].re =-1;    SignalSet[ 6].im = 1;
		SignalSet[ 7].re =-3;    SignalSet[ 7].im = 1;
		SignalSet[ 8].re = 3;    SignalSet[ 8].im =-1;
		SignalSet[ 9].re = 1;    SignalSet[ 9].im =-1;
		SignalSet[10].re =-1;    SignalSet[10].im =-1;
		SignalSet[11].re =-3;    SignalSet[11].im =-1;
		SignalSet[12].re = 3;    SignalSet[12].im =-3;
		SignalSet[13].re = 1;    SignalSet[13].im =-3;
		SignalSet[14].re =-1;    SignalSet[14].im =-3;
		SignalSet[15].re =-3;    SignalSet[15].im =-3;*/
		for (int i = 0; i<16; i++)
			SignalMapper[i] = i;
		SignalSet[0].re = 1;    SignalSet[0].im = 1;
		SignalSet[1].re = 1;    SignalSet[1].im = 3;
		SignalSet[2].re = 3;    SignalSet[2].im = 1;
		SignalSet[3].re = 3;    SignalSet[3].im = 3;
		SignalSet[4].re = 1;    SignalSet[4].im = -1;
		SignalSet[5].re = 1;    SignalSet[5].im = -3;
		SignalSet[6].re = 3;    SignalSet[6].im = -1;
		SignalSet[7].re = 3;    SignalSet[7].im = -3;
		SignalSet[8].re = -1;    SignalSet[8].im = 1;
		SignalSet[9].re = -1;    SignalSet[9].im = 3;
		SignalSet[10].re = -3;    SignalSet[10].im = 1;
		SignalSet[11].re = -3;    SignalSet[11].im = 3;
		SignalSet[12].re = -1;    SignalSet[12].im = -1;
		SignalSet[13].re = -1;    SignalSet[13].im = -3;
		SignalSet[14].re = -3;    SignalSet[14].im = -1;
		SignalSet[15].re = -3;    SignalSet[15].im = -3;
	}
	if (N_BITS == 5) {
		SignalMapper[31] = 0;    SignalMapper[27] = 1;
		SignalMapper[19] = 2;    SignalMapper[23] = 3;
		SignalMapper[26] = 4;    SignalMapper[15] = 5;
		SignalMapper[11] = 6;    SignalMapper[3] = 7;
		SignalMapper[7] = 8;    SignalMapper[18] = 9;
		SignalMapper[30] = 10;    SignalMapper[14] = 11;
		SignalMapper[10] = 12;    SignalMapper[2] = 13;
		SignalMapper[6] = 14;    SignalMapper[22] = 15;

		SignalMapper[28] = 16;    SignalMapper[12] = 17;
		SignalMapper[8] = 18;    SignalMapper[0] = 19;
		SignalMapper[4] = 20;    SignalMapper[20] = 21;
		SignalMapper[24] = 22;    SignalMapper[13] = 23;
		SignalMapper[9] = 24;    SignalMapper[1] = 25;
		SignalMapper[5] = 26;    SignalMapper[16] = 27;
		SignalMapper[29] = 28;    SignalMapper[25] = 29;
		SignalMapper[17] = 30;    SignalMapper[21] = 31;

		SignalSet[0].re = 3;    SignalSet[0].im = 5;
		SignalSet[1].re = 1;    SignalSet[1].im = 5;
		SignalSet[2].re = -1;    SignalSet[2].im = 5;
		SignalSet[3].re = -3;    SignalSet[3].im = 5;

		SignalSet[4].re = 5;    SignalSet[4].im = 3;
		SignalSet[5].re = 3;    SignalSet[5].im = 3;
		SignalSet[6].re = 1;    SignalSet[6].im = 3;
		SignalSet[7].re = -1;    SignalSet[7].im = 3;
		SignalSet[8].re = -3;    SignalSet[8].im = 3;
		SignalSet[9].re = -5;    SignalSet[9].im = 3;

		SignalSet[10].re = 5;    SignalSet[10].im = 1;
		SignalSet[11].re = 3;    SignalSet[11].im = 1;
		SignalSet[12].re = 1;    SignalSet[12].im = 1;
		SignalSet[13].re = -1;    SignalSet[13].im = 1;
		SignalSet[14].re = -3;    SignalSet[14].im = 1;
		SignalSet[15].re = -5;    SignalSet[15].im = 1;


		SignalSet[16].re = 5;    SignalSet[16].im = -1;
		SignalSet[17].re = 3;    SignalSet[17].im = -1;
		SignalSet[18].re = 1;    SignalSet[18].im = -1;
		SignalSet[19].re = -1;    SignalSet[19].im = -1;
		SignalSet[20].re = -3;    SignalSet[20].im = -1;
		SignalSet[21].re = -5;    SignalSet[21].im = -1;

		SignalSet[22].re = 5;    SignalSet[22].im = -3;
		SignalSet[23].re = 3;    SignalSet[23].im = -3;
		SignalSet[24].re = 1;    SignalSet[24].im = -3;
		SignalSet[25].re = -1;    SignalSet[25].im = -3;
		SignalSet[26].re = -3;    SignalSet[26].im = -3;
		SignalSet[27].re = -5;    SignalSet[27].im = -3;

		SignalSet[28].re = 3;    SignalSet[28].im = -5;
		SignalSet[29].re = 1;    SignalSet[29].im = -5;
		SignalSet[30].re = -1;    SignalSet[30].im = -5;
		SignalSet[31].re = -3;    SignalSet[31].im = -5;
	}
	if (N_BITS == 6) {
		/*SignalMapper[36] = 9;    SignalMapper[37] = 41;
		SignalMapper[39] = 57;    SignalMapper[38] = 25;
		SignalMapper[34] = 17;    SignalMapper[35] = 49;
		SignalMapper[33] = 33;    SignalMapper[32] = 1;

		SignalMapper[44] = 13;    SignalMapper[45] = 45;
		SignalMapper[47] = 61;    SignalMapper[46] = 29;
		SignalMapper[42] = 21;    SignalMapper[43] = 53;
		SignalMapper[41] = 37;    SignalMapper[40] = 5;

		SignalMapper[60] = 15;    SignalMapper[61] = 47;
		SignalMapper[63] = 63;    SignalMapper[62] = 31;
		SignalMapper[58] = 23;    SignalMapper[59] = 55;
		SignalMapper[57] = 39;    SignalMapper[56] = 7;

		SignalMapper[52] = 11;    SignalMapper[53] = 43;
		SignalMapper[54] = 27;    SignalMapper[55] = 59;
		SignalMapper[50] = 19;    SignalMapper[51] = 51;
		SignalMapper[49] = 35;    SignalMapper[48] = 3;
		//-------------------------------------------

		SignalMapper[28] = 14;    SignalMapper[29] = 46;
		SignalMapper[31] = 62;    SignalMapper[30] = 30;
		SignalMapper[26] = 22;    SignalMapper[27] = 54;
		SignalMapper[25] = 38;    SignalMapper[24] = 6;

		SignalMapper[20] = 10;    SignalMapper[21] = 42;
		SignalMapper[23] = 58;    SignalMapper[22] = 26;
		SignalMapper[18] = 18;    SignalMapper[19] = 50;
		SignalMapper[17] = 34;    SignalMapper[16] = 2;

		SignalMapper[12] = 12;    SignalMapper[13] = 44;
		SignalMapper[15] = 60;    SignalMapper[14] = 28;
		SignalMapper[10] = 20;    SignalMapper[11] = 52;
		SignalMapper[9] = 36;    SignalMapper[8] = 4;

		SignalMapper[4] = 8;    SignalMapper[5] = 40;
		SignalMapper[7] = 56;    SignalMapper[6] = 24;
		SignalMapper[2] = 16;    SignalMapper[3] = 48;
		SignalMapper[1] = 32;    SignalMapper[0] = 0;


		SignalSet[0].re = 3;    SignalSet[0].im = 3;
		SignalSet[1].re = 3;    SignalSet[1].im = 1;
		SignalSet[2].re = 1;    SignalSet[2].im = 3;
		SignalSet[3].re = 1;    SignalSet[3].im = 1;
		SignalSet[4].re = 3;    SignalSet[4].im = 5;
		SignalSet[5].re = 3;    SignalSet[5].im = 7;
		SignalSet[6].re = 1;    SignalSet[6].im = 5;
		SignalSet[7].re = 1;    SignalSet[7].im = 7;

		SignalSet[8].re = 5;    SignalSet[8].im = 3;
		SignalSet[9].re = 5;    SignalSet[9].im = 1;
		SignalSet[10].re = 7;    SignalSet[10].im = 3;
		SignalSet[11].re = 7;    SignalSet[11].im = 1;
		SignalSet[12].re = 5;    SignalSet[12].im = 5;
		SignalSet[13].re = 5;    SignalSet[13].im = 7;
		SignalSet[14].re = 7;    SignalSet[14].im = 5;
		SignalSet[15].re = 7;    SignalSet[15].im = 7;


		SignalSet[16].re = 3;    SignalSet[16].im = -3;
		SignalSet[17].re = 3;    SignalSet[17].im = -1;
		SignalSet[18].re = 1;    SignalSet[18].im = -3;
		SignalSet[19].re = 1;    SignalSet[19].im = -1;
		SignalSet[20].re = 3;    SignalSet[20].im = -5;
		SignalSet[21].re = 3;    SignalSet[21].im = -7;
		SignalSet[22].re = 1;    SignalSet[22].im = -5;
		SignalSet[23].re = 1;    SignalSet[23].im = -7;

		SignalSet[24].re = 5;    SignalSet[24].im = -3;
		SignalSet[25].re = 5;    SignalSet[25].im = -1;
		SignalSet[26].re = 7;    SignalSet[26].im = -3;
		SignalSet[27].re = 7;    SignalSet[27].im = -1;
		SignalSet[28].re = 5;    SignalSet[28].im = -5;
		SignalSet[29].re = 5;    SignalSet[29].im = -7;
		SignalSet[30].re = 7;    SignalSet[30].im = -5;
		SignalSet[31].re = 7;    SignalSet[31].im = -7;
		//-----------------------------------------------
		SignalSet[32].re = -3;    SignalSet[32].im = 3;
		SignalSet[33].re = -3;    SignalSet[33].im = 1;
		SignalSet[34].re = -1;    SignalSet[34].im = 3;
		SignalSet[35].re = -1;    SignalSet[35].im = 1;
		SignalSet[36].re = -3;    SignalSet[36].im = 5;
		SignalSet[37].re = -3;    SignalSet[37].im = 7;
		SignalSet[38].re = -1;    SignalSet[38].im = 5;
		SignalSet[39].re = -1;    SignalSet[39].im = 7;

		SignalSet[40].re = -5;    SignalSet[40].im = 3;
		SignalSet[41].re = -5;    SignalSet[41].im = 1;
		SignalSet[42].re = -7;    SignalSet[42].im = 3;
		SignalSet[43].re = -7;    SignalSet[43].im = 1;
		SignalSet[44].re = -5;    SignalSet[44].im = 5;
		SignalSet[45].re = -5;    SignalSet[45].im = 7;
		SignalSet[46].re = -7;    SignalSet[46].im = 5;
		SignalSet[47].re = -7;    SignalSet[47].im = 7;

		SignalSet[48].re = -3;    SignalSet[48].im = -3;
		SignalSet[49].re = -3;   SignalSet[49].im = -1;
		SignalSet[50].re = -1;    SignalSet[50].im = -3;
		SignalSet[51].re = -1;    SignalSet[51].im = -1;
		SignalSet[52].re = -3;    SignalSet[52].im = -5;
		SignalSet[53].re = -3;    SignalSet[53].im = -7;
		SignalSet[54].re = -1;    SignalSet[54].im = -5;
		SignalSet[55].re = -1;   SignalSet[55].im = -7;

		SignalSet[56].re = -5;    SignalSet[56].im = -3;
		SignalSet[57].re = -5;    SignalSet[57].im = -1;
		SignalSet[58].re = -7;    SignalSet[58].im = -3;
		SignalSet[59].re = -7;    SignalSet[59].im = -1;
		SignalSet[60].re = -5;    SignalSet[60].im = -5;
		SignalSet[61].re = -5;   SignalSet[61].im = -7;
		SignalSet[62].re = -7;    SignalSet[62].im = -5;
		SignalSet[63].re = -7;    SignalSet[63].im = -7;*/
		for (int i = 0; i<64; i++)
			SignalMapper[i] = i;
		SignalSet[0].re = 3;    SignalSet[0].im = 3;
		SignalSet[1].re = 3;    SignalSet[1].im = 1;
		SignalSet[2].re = 1;    SignalSet[2].im = 3;
		SignalSet[3].re = 1;    SignalSet[3].im = 1;
		SignalSet[4].re = 3;    SignalSet[4].im = 5;
		SignalSet[5].re = 3;    SignalSet[5].im = 7;
		SignalSet[6].re = 1;    SignalSet[6].im = 5;
		SignalSet[7].re = 1;    SignalSet[7].im = 7;

		SignalSet[8].re = 5;    SignalSet[8].im = 3;
		SignalSet[9].re = 5;    SignalSet[9].im = 1;
		SignalSet[10].re = 7;    SignalSet[10].im = 3;
		SignalSet[11].re = 7;    SignalSet[11].im = 1;
		SignalSet[12].re = 5;    SignalSet[12].im = 5;
		SignalSet[13].re = 5;    SignalSet[13].im = 7;
		SignalSet[14].re = 7;    SignalSet[14].im = 5;
		SignalSet[15].re = 7;    SignalSet[15].im = 7;


		SignalSet[16].re = 3;    SignalSet[16].im = -3;
		SignalSet[17].re = 3;    SignalSet[17].im = -1;
		SignalSet[18].re = 1;    SignalSet[18].im = -3;
		SignalSet[19].re = 1;    SignalSet[19].im = -1;
		SignalSet[20].re = 3;    SignalSet[20].im = -5;
		SignalSet[21].re = 3;    SignalSet[21].im = -7;
		SignalSet[22].re = 1;    SignalSet[22].im = -5;
		SignalSet[23].re = 1;    SignalSet[23].im = -7;

		SignalSet[24].re = 5;    SignalSet[24].im = -3;
		SignalSet[25].re = 5;    SignalSet[25].im = -1;
		SignalSet[26].re = 7;    SignalSet[26].im = -3;
		SignalSet[27].re = 7;    SignalSet[27].im = -1;
		SignalSet[28].re = 5;    SignalSet[28].im = -5;
		SignalSet[29].re = 5;    SignalSet[29].im = -7;
		SignalSet[30].re = 7;    SignalSet[30].im = -5;
		SignalSet[31].re = 7;    SignalSet[31].im = -7;
		//-----------------------------------------------
		SignalSet[32].re = -3;    SignalSet[32].im = 3;
		SignalSet[33].re = -3;    SignalSet[33].im = 1;
		SignalSet[34].re = -1;    SignalSet[34].im = 3;
		SignalSet[35].re = -1;    SignalSet[35].im = 1;
		SignalSet[36].re = -3;    SignalSet[36].im = 5;
		SignalSet[37].re = -3;    SignalSet[37].im = 7;
		SignalSet[38].re = -1;    SignalSet[38].im = 5;
		SignalSet[39].re = -1;    SignalSet[39].im = 7;

		SignalSet[40].re = -5;    SignalSet[40].im = 3;
		SignalSet[41].re = -5;    SignalSet[41].im = 1;
		SignalSet[42].re = -7;    SignalSet[42].im = 3;
		SignalSet[43].re = -7;    SignalSet[43].im = 1;
		SignalSet[44].re = -5;    SignalSet[44].im = 5;
		SignalSet[45].re = -5;    SignalSet[45].im = 7;
		SignalSet[46].re = -7;    SignalSet[46].im = 5;
		SignalSet[47].re = -7;    SignalSet[47].im = 7;

		SignalSet[48].re = -3;    SignalSet[48].im = -3;
		SignalSet[49].re = -3;   SignalSet[49].im = -1;
		SignalSet[50].re = -1;    SignalSet[50].im = -3;
		SignalSet[51].re = -1;    SignalSet[51].im = -1;
		SignalSet[52].re = -3;    SignalSet[52].im = -5;
		SignalSet[53].re = -3;    SignalSet[53].im = -7;
		SignalSet[54].re = -1;    SignalSet[54].im = -5;
		SignalSet[55].re = -1;   SignalSet[55].im = -7;

		SignalSet[56].re = -5;    SignalSet[56].im = -3;
		SignalSet[57].re = -5;    SignalSet[57].im = -1;
		SignalSet[58].re = -7;    SignalSet[58].im = -3;
		SignalSet[59].re = -7;    SignalSet[59].im = -1;
		SignalSet[60].re = -5;    SignalSet[60].im = -5;
		SignalSet[61].re = -5;   SignalSet[61].im = -7;
		SignalSet[62].re = -7;    SignalSet[62].im = -5;
		SignalSet[63].re = -7;    SignalSet[63].im = -7;

	}
	if (N_BITS > 7) {
		printf("The modulation_type %d is not supported currently\n", N_BITS);
		exit(1);
	}
	NormalizeSignalPower(SignalSet, ConstellationSize);
}

void NormalizeSignalPower(dcomplex *SignalSet, int ConstellationSize)
{
	double sum = 0.;
	int i = 0;

	for (i = 0; i < ConstellationSize; i++)
		sum += (SignalSet[i].re * SignalSet[i].re +
			SignalSet[i].im * SignalSet[i].im);

	sum = sqrt(sum / ConstellationSize);
	for (i = 0; i < ConstellationSize; i++)
	{
		SignalSet[i].re /= sum;
		SignalSet[i].im /= sum;
	}
}

void Normalize(double *data, int num)
{
	int m;
	double sum = 0;

	for (m = 0; m < num; m++)
		sum += data[m];
	for (m = 0; m < num; m++) {
		data[m] /= sum;
		data[m] = clip(data[m]);
	}
}

double clip(double x)                //avoid flowimg
{
#define TOOSMALL 1e-50
#define TOOLARGE 1e50
	if (x >= 0) {
		if (x > TOOLARGE)
			x = TOOLARGE;
		else if (x < TOOSMALL)
			x = TOOSMALL;
	}
	if (x < 0) {
		if (x > -TOOSMALL)
			x = -TOOSMALL;
		if (x < -TOOLARGE)
			x = -TOOLARGE;
	}
	return x;
}
