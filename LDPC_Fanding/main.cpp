#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <iostream>
#include <complex>

#include "randomc.h"
#include "struct.h"
#include "awgn_channel.h"
#include "function.h"
#include "pn_seq.h"
#include "Complex.h"
#include "modulation.h"
#include "demodulation.h"
#include "fading_channel.h"
#include "fading_channel_demodulate.h"
//CRandomMersenne rg(179);


using namespace std;

int main()
{
	int i;
	matrix sparse;
	profile profiles;
	char **gauss;
	int *exchange;
	int rank;
	//char matrixfile[200] = "2048_384_mackay.code";
	v_node *v_line;		//Variable node
	c_node *c_line;		//Check node

	//读出参数
	ReadProfile(&profiles);

	//读出校验矩阵并高斯消去
	init_sparse(&sparse, profiles.matrixfile);

	//exchange = (int *)malloc(sizeof(int)*(sparse.col_num + 1)); //Exchanging Column Index
	exchange = new int[sparse.col_num + 1];


	for (i = 1; i <= sparse.col_num; i++)
		exchange[i] = i;

	gauss = (char**)malloc(sizeof(char *)*(sparse.row_num + 1));
	for (i = 0; i <= sparse.row_num; i++)
		gauss[i] = (char*)malloc(sizeof(char)*(sparse.col_num + 1));

	init_gauss(&sparse, gauss, exchange);
	count_gauss(&sparse, gauss, exchange, rank);
	check_gauss(&sparse, gauss, rank);

	//	v_line=new v_node[sparse.col_num+1];
	v_line = (v_node *)malloc(sizeof(v_node)*(sparse.col_num + 1));
	//	c_line=new c_node[sparse.row_num+1];
	c_line = (c_node *)malloc(sizeof(c_node)*(sparse.row_num + 1));
	init_decode(&sparse, v_line, c_line);
	//写文件
	char titile[300];
	char record[300];
	FILE *fout, *rout;
	AWGN p, *pi;
	p.rate = (double)sparse.message_num / sparse.col_num;
	p.reset = 1;
	titile[0] = '\0';

	sprintf_s(titile, "%d_%d_%.4f_%d_%d_%d", sparse.col_num, sparse.message_num, profiles.q_step, profiles.q_bits, profiles.codeword,profiles.channel_style);

	for (i = 0; i < 30; i++)
		record[i] = titile[i];
	strcat_s(record, ".rcd");

	strcat_s(titile, "_PER.txt");

	fopen_s(&fout, titile, "a");

	fprintf(fout, "%s ;", profiles.matrixfile);
	fprintf(fout, "R: %f; EM: %d; MIN ERR: %d; MAX ITR: %d; ER_EA: %d; Br_LEN: %d q_step:%f\n\n", p.rate, profiles.codemethod, profiles.leasterror, profiles.max_iter, profiles.err_era, profiles.burstlen, profiles.q_step);
	fprintf(fout, "Q_bits: %d; Q_step: %f;Algorithm: %d;Channel: %d\n\n", profiles.q_bits, profiles.q_step, profiles.codeword, profiles.channel_style);
	fprintf(fout, "SNR	    iBER			 BER		   FER			ITR#\n");
	fclose(fout);
	fout = NULL;

	fopen_s(&rout, record, "a");
	fprintf(rout, "R: %f; EM: %d; MIN ERR: %d; MAX ITR: %d; ER_EA: %d; Br_LEN: %d\n\n", p.rate, profiles.codemethod, profiles.leasterror, profiles.max_iter, profiles.err_era, profiles.burstlen);
	fprintf(rout, "Q_bits: %d; Q_step: %f;Algorithm: %d;modulate_mod: %d; Channel: %d\n\n", profiles.q_bits, profiles.q_step, profiles.codeword, profiles.modulate_mode, profiles.channel_style);
	fclose(rout);
	rout = NULL;

	printf("\n*Successfully Initiate code; Simulation Start. Written by R.H, 2018*\n");
	printf("Q_bits: %d; Q_step: %f; Algorithm: %d; modulate_mode:%d\n\n", profiles.q_bits, profiles.q_step, profiles.codeword, profiles.modulate_mode);

	printf( "SNR    iBER	     BER	   FER	   ITR#\n");

	//进入主循环
	int step;
	char *message, *code;
	int  *cw_bits_int;
	char reg[24] = { 0 };
	double	ierror_bit_rate, error_bit_rate, error_frame_rate, ave_iter;
	bool	rightframe;
	double  ierror_bit, error_bit, error_frame;
	double	frame;
	int		len_symbol;
	double  snr;
	int		iter;

	double	*lp, *lp0;	//	double **lq, **lr;	int **r_list, **q_list;
	double	*lp1;	// The second largest
	int		*lp_int;		//quantized LLR  
	char	*resultcode;
	double	*transcode;
	int		*min_indx;	//The minimum indx of each check sum
	int		*decodedresult;
	int		*chk_sum, *code_sgn, *weight_1, *weight_2;
	double	q_step;		//quantized step
	int		q_bits;
	float	delta;					  // d of DD-BMP
	bool	decode_valid = true;
	int		max_iter = 50;
	double	noise;
	float	scale;
	double	d_quasi;                         //parameter of (q+1)-bit quasi-uniform quantizer
	int		dsum;
	int *polar_code;

	modulator mod;
	demodulator demod;

	len_symbol = sparse.col_num / profiles.modulate_mode;

	
	reg[0] = 1;
	reg[21] = 1;
	reg[17] = 1;

	message = (char*)malloc(sizeof(char)*(sparse.message_num + 1));
	code = (char*)malloc(sizeof(char)*(sparse.col_num + 1));
	cw_bits_int = (int*)malloc(sizeof(int)*(sparse.col_num));

	polar_code = (int*)malloc(sizeof(int)*(sparse.col_num + 1));
	transcode = (double*)malloc(sizeof(double)*(sparse.col_num + 1));
	resultcode = (char*)malloc(sizeof(char)*(sparse.col_num + 1));
	lp = (double *)malloc(sizeof(double)*(sparse.col_num + 1));
	lp0 = (double *)malloc(sizeof(double)*(2 * sparse.col_num + 1));
	lp_int = (int *)malloc(sizeof(int)*(sparse.col_num + 1));


	lp1 = (double *)malloc(sizeof(double)*(sparse.row_num + 1));
	min_indx = (int *)malloc(sizeof(int)*(sparse.row_num + 1));

	code_sgn = (int *)malloc(sizeof(int)*(sparse.col_num + 1));
	chk_sum = (int *)malloc(sizeof(int)*(sparse.row_num + 1));
	weight_1 = (int *)malloc(sizeof(int)*(sparse.row_num + 1));
	weight_2 = (int *)malloc(sizeof(int)*(sparse.row_num + 1));

	snr		=	profiles.snr;
	q_step	=	profiles.q_step;
	q_bits	=	profiles.q_bits;
	delta	=	profiles.delta;
	max_iter =	profiles.max_iter;
	scale	=	profiles.scale;
	d_quasi =   profiles.d_quasi;                         //parameter of (q+1)-bit quasi-uniform quantizer
	dsum	=	profiles.dsum;

	//double *soft_bits_out = new double[sparse.col_num];

	for (step = 0; step < profiles.test_num; step++)
	{
//		double sigma = sqrt(1.0 / (2 * pow(10.0, snr / 10.0))); //Es/No if Es=1
//		double sigma = sqrt(1.0 / (2 * pow(10.0, snr / 10.0)*p.rate));//Eb/N0

		p.sigma = sqrt(1.0 / (2 * pow(10.0, snr / 10.0)*p.rate));//Eb/N0
		fopen_s(&rout, record, "a");
		fprintf(rout, "\n%s%1.4f\n", "SNR:",snr);
		fclose(rout);
		rout = NULL;
		p.I_sigma = 1.0 / p.sigma;
		noise = 1.0 / p.I_sigma;
		frame = 0;
		error_bit_rate = 0;
		error_bit = 0;
		ierror_bit_rate = 0;
		ierror_bit = 0;
		error_frame = 0;
		error_frame_rate = 0;
		iter = 0;
		double const1 = 2 / pow(p.sigma, 2);

		p.ix = profiles.seed1;
		p.iy = profiles.seed2;
		p.iz = profiles.seed3;
		pi = &p;

		while ((frame < profiles.leastframe) || (error_frame < profiles.leasterror))
		{
			frame++;
			rightframe = true;

			for (i = 0; i <= sparse.message_num; i++)
				message[i] = 0;
			for (i = 1; i <= sparse.col_num; i++)
				code[i] = 0;

			if (profiles.codemethod == 1)
			{
				pn_init(message, reg, (sparse.message_num), false);
				encode(&sparse, gauss, exchange, message, code);

				if (check_code(&sparse, code) == false)
				{
					printf("Code Invalid!!!\n");
					break;
				}
				for (i = 0; i < sparse.col_num; i++)
				{
//					cw_bits_int[i] = (code[i+1] == 0 ? 0 : 1);
					cw_bits_int[i] = (code[i + 1] == 0 ? 1 : 0);
					polar_code[i+1] = (code[i+1] == 0 ? 1 : -1);
				}
			}
			switch (profiles.channel_style)
			{
			case 0:
			{
				complex<float> *modulated_symbols = new complex<float>[len_symbol];
				complex<float> *receiver_frame_signal = new complex<float>[len_symbol];

				dcomplex *moded_symbols = new dcomplex[len_symbol];
				dcomplex *received_symbols = new dcomplex[len_symbol];
				dcomplex *SideInfo = new dcomplex[len_symbol];
				//调制
				mod.modulation(cw_bits_int, sparse.col_num, profiles.modulate_mode, modulated_symbols);
				//衰落信道
				fading_channel(modulated_symbols, moded_symbols, received_symbols, SideInfo, len_symbol, profiles.DOPPLER, profiles.T_SAMPLE, profiles.RAYLEIGH_M, p.sigma);
				//rayleigh demodulation
				double *soft_bits_out_double = new double[sparse.col_num];
				fading_channel_demodulate(received_symbols, soft_bits_out_double, profiles.modulate_mode, p.sigma, len_symbol, SideInfo);
				//for (int i = 0; i<sparse.col_num; i++)
				//	soft_bits_out[i] = (float)soft_bits_out_double[i];
				for (int i = 1; i <= sparse.col_num; i++)
				//	transcode[i] = soft_bits_out_double[i - 1];
					transcode[i] = soft_bits_out_double[i - 1];
				delete[] received_symbols;
				delete[] soft_bits_out_double;
				delete[] moded_symbols;
				delete[] SideInfo;
			} break;
			case 1:
				awgn_channel(polar_code, transcode, pi, sparse.col_num + 1); break;
			default:break;
			}
			/*
			complex<float> *modulated_symbols = new complex<float>[len_symbol];
			complex<float> *receiver_frame_signal = new complex<float>[len_symbol];

			dcomplex *moded_symbols = new dcomplex[len_symbol];
			dcomplex *received_symbols = new dcomplex[len_symbol];
			dcomplex *SideInfo = new dcomplex[len_symbol];
			//调制
			mod.modulation(cw_bits_int, sparse.col_num, profiles.modulate_mode, modulated_symbols);
			//衰落信道
			fading_channel(modulated_symbols, moded_symbols, received_symbols, SideInfo, len_symbol, profiles.DOPPLER, profiles.T_SAMPLE, profiles.RAYLEIGH_M, p.sigma);
			//rayleigh demodulation
			double *soft_bits_out_double = new double[sparse.col_num];
			fading_channel_demodulate(received_symbols, soft_bits_out_double, profiles.modulate_mode, p.sigma, len_symbol, SideInfo);
			//for (int i = 0; i<sparse.col_num; i++)
			//	soft_bits_out[i] = (float)soft_bits_out_double[i];
			for (int i = 1; i<=sparse.col_num; i++)
				transcode[i] = soft_bits_out_double[i-1];
			delete[] received_symbols;
			delete[] soft_bits_out_double;
			delete[] moded_symbols;
			delete[] SideInfo;
			*/
			//高斯信道
			//awgn_channel(polar_code, transcode, pi, sparse.col_num + 1);
			//译码
			switch (profiles.codeword)
			{
			case 0: decode_ms(&sparse, resultcode, v_line, c_line, lp, lp0, iter, transcode, noise, decode_valid, max_iter, delta, q_step, q_bits); break;
			case 1: decode_qua_ms(&sparse, resultcode, v_line, c_line, lp, lp0, iter, transcode, noise, decode_valid, max_iter, delta, q_step, q_bits); break;
				//case 1: decode_ms_q(&sparse,resultcode,v_line,c_line,lp,lp0,iter,transcode,noise,decode_valid,max_iter,delta, v_step, v_max, c_step, c_max); break;
			case 2: decode_rbi_mlgd(&sparse, resultcode, v_line, c_line, lp, lp0, weight_1, weight_2, chk_sum, code_sgn, iter, transcode, noise, decode_valid, max_iter, delta, q_step, q_bits); break;
			case 3: decode_quai_rbi_mlgd(&sparse, resultcode, v_line, c_line, lp, lp0, weight_1, weight_2, chk_sum, code_sgn, iter, transcode, noise, decode_valid, max_iter, delta, q_step, q_bits); break;
			case 4: decode_srbmld(&sparse, resultcode, v_line, c_line, lp, lp0, weight_1, weight_2, chk_sum, code_sgn, iter, transcode, noise, decode_valid, max_iter, delta, scale, q_step, q_bits); break;
			case 5: decode_mrbmld(&sparse, resultcode, v_line, c_line, lp, lp0, weight_1, weight_2, chk_sum, code_sgn, iter, transcode, noise, decode_valid, max_iter, delta, q_step, q_bits); break;
				//case 6: decode_weighted_ddwmld(&sparse,resultcode,v_line,c_line,lp,lp0,chk_sum, code_sgn, iter,transcode, noise,decode_valid,max_iter, delta,q_bits,chk_d); break;
			case 6: decode_improved_srbmld(&sparse, resultcode, v_line, c_line, lp, lp0, weight_1, weight_2, chk_sum, code_sgn, iter, transcode, noise, decode_valid, max_iter, delta, scale, q_step, q_bits); break;
			case 7: decode_iosmld(&sparse, resultcode, v_line, c_line, lp, lp0, chk_sum, code_sgn, iter, transcode, noise, decode_valid, max_iter, delta, q_bits); break;
			case 8: decode_osmld(&sparse, resultcode, v_line, c_line, lp, lp0, chk_sum, code_sgn, iter, transcode, noise, decode_valid, max_iter, delta, q_bits); break;
			case 9: decode_bf(&sparse, resultcode, v_line, c_line, lp, lp0, chk_sum, code_sgn, iter, transcode, noise, decode_valid, max_iter, delta, q_bits); break;
			case 10: decode_ddbmp(&sparse, resultcode, v_line, c_line, lp, lp0, iter, transcode, noise, decode_valid, max_iter, q_step, q_bits); break;
			//case 11: vMSDecoder(&sparse, resultcode, dppCNValue, ippVNSgn, dppVNValue, dpLLRTotal, transcode, lp0, iter, decode_valid, max_iter, delta, code, q_bits, q_step, d_quasi); break;
			case 12: decode_wbf(&sparse, resultcode, v_line, c_line, lp, lp0, chk_sum, code_sgn, iter, transcode, noise, decode_valid, max_iter, delta, q_bits); break;
			case 13:
				if (profiles.channel_style == 0){
					for (i = 1; i <= sparse.col_num; i++) {
						transcode[i] *= const1*0.86;
					}
				
				}
				else {
					for (i = 1; i <= sparse.col_num; i++) {
						transcode[i] *= const1;
					}
				}


				decodeon(&sparse, resultcode, v_line, c_line, lp, lp0, iter, transcode, noise, decode_valid, max_iter);
				break;
			case 14: decode_quasi_ms(&sparse, resultcode, v_line, c_line, lp, lp0, iter, transcode, noise, decode_valid, max_iter, delta, q_step, q_bits, d_quasi); break;
			case 15: // hard-decision min-sum
				for (i = 1; i <= sparse.col_num; i++)
				{
					if (transcode[i] > 0) {
						transcode[i] = 1;
					}
					else {
						transcode[i] = -1;
					}
				}
				decode_ms(&sparse, resultcode, v_line, c_line, lp, lp0, iter, transcode, noise, decode_valid, max_iter, delta, q_step, q_bits);
				break;

			case 16: decode_dsms(&sparse, resultcode, v_line, c_line, lp, lp0, iter, transcode, noise, decode_valid, max_iter, delta, q_step, q_bits, dsum); break;
			case 17:
				if (profiles.channel_style == 0) {
					for (i = 1; i <= sparse.col_num; i++) {
						transcode[i] *= const1*0.86;
					}

				}
				else {
					for (i = 1; i <= sparse.col_num; i++) {
						transcode[i] *= const1;
					}
				}


				spa_qua(&sparse, resultcode, v_line, c_line, lp, lp0, iter, transcode, noise, decode_valid, max_iter,0.125,8);
				break;

			default: break;

			}
			rightframe = true;
			for (i = 1; i <= sparse.col_num; i++)
			{
				if (resultcode[i] != code[i])
				{
					rightframe = false;
					error_bit++;
				}
			}

			for (i = 1; i <= sparse.col_num; i++)
			{
				if (transcode[i] > 0) {
					transcode[i] = 0;
				}
				else {
					transcode[i] = 1;
				}
				if (transcode[i] != code[i])
				{
					ierror_bit++;
				}
			}


			if (!rightframe)
			{
				error_frame++;
			}

			ierror_bit_rate = ierror_bit / (frame*sparse.col_num);
			error_bit_rate = error_bit / (frame*sparse.col_num);
			error_frame_rate = error_frame / frame;
			ave_iter = (double)iter / frame;

			if (!rightframe)
			{
				fopen_s(&rout, record, "a");
				fprintf(rout, "\n%s%1.2e%s%1.2e%s%2.0f%s%8.2e%s%8.2e%s%8.2e%s%5.2f", "F:", frame, " EB:", error_bit, " EF:", error_frame, " FER:", error_frame_rate, " BER:", error_bit_rate, " IER:", ierror_bit_rate, " ITR:", ave_iter);
				fclose(rout);
			}
			
		}


		fopen_s(&fout, titile, "a");

		fprintf(fout, "%.4f	%e %e %e %f\n", snr, ierror_bit_rate, error_bit_rate, error_frame_rate, ave_iter);

		printf("%.4f %e %e %e %f\n", snr, ierror_bit_rate, error_bit_rate, error_frame_rate, ave_iter);
		fclose(fout);
		snr += profiles.snr_step;

		}
		free(v_line);
		free(c_line);
		for (i = 0; i <= sparse.row_num; i++)
			free(gauss[i]);
		free(gauss);

		getchar();
		return 0;
	}