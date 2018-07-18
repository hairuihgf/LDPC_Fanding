#pragma once
#ifndef FUNCTION_H
#define FUNCTION_H
#include <vector>
#include <stdio.h>
#include "struct.h"



using std::vector;


void init_sparse(struct matrix *sparse, char*filename);
void init_gauss(struct matrix *sparse, char **gauss, int *exchange);
void check_gauss(struct matrix *sparse, char **gauss, int& rank);
void count_gauss(struct matrix *sparse, char **gauss, int *exchange, int& rank);
void ReadProfile(AWGN p, AWGN *pi, char ch, char matrixfile[200], int test_num, float readtemp, int	leastframe, int	leasterror, int	codemethod,
	int	codeword, int max_iter, int	err_era, int burstlen, int	burstnum, int de_alg, float delta, float	scale, int	q_bits, double	alpha, double	d_quasi,
	int	dsum, int	max_q_level, double	q_step, double  snr_step);
void encode(struct matrix *sparse, char **gauss, int *exchange, char*message, char *code);
bool check_code(struct matrix *sparse, char*code);
void init_decode(struct matrix *sparse, struct v_node *v_line, struct c_node *c_line);
void check_result2(struct matrix *sparse, bool &decode_valid, double *lp, char *resultcode);
double LLR_add(double in1, double in2);
double LLR_sub(double in1, double in2);
void check_result2(struct matrix *sparse, char*resultcode, double *lp, bool &decode_valid);
void decodeon(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter);
void decode_ms(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, double q_step, int q_bits);
void decode_dsms(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, double q_step, int q_bits, int dsum);

void decode_qua_ms(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, double q_step, int q_bits);


void decode_quasi_ms(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, double q_step, int q_bits, double d_quasi);
void decode_ms_q(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, double scale, int v_step, int v_max, int c_step, int c_max);

void vMSDecoder(struct matrix *S_Sparse, char *szResultCode, double **dppCNValue, int **ippVNSgn, double **dppVNValue, double *dpLLRTotal, double *dpChannelLLR, double *lp0, int &iIter, bool &bDecodeValid, int iMaxIter, double delta, char *szCode, int q_bits, double q_step, double d_quasi);

void decode_bst(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter);
void decode_ddbmp(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, double q_step, int q_bits);
void decode_ddwmld(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, int q_bits);
void decode_weighted_ddwmld(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, float scale, int q_bits);
void decode_srbmld(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *weight_1, int *weight_2, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, float scale, double q_step, int q_bits);
void decode_improved_srbmld(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *weight_1, int *weight_2, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, float scale, double q_step, int q_bits);
void decode_mrbmld(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *weight_1, int *weight_2, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, double q_step, int q_bits);
void decode_rbi_mlgd(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *weight_1, int *weight_2, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, double q_step, int q_bits);
void decode_quai_rbi_mlgd(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *weight_1, int *weight_2, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, double q_step, int q_bits);
void decode_iosmld(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, int q_bits);
void decode_osmld(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, int q_bits);
void decode_bf(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, int q_bits);
void decode_wbf(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, int q_bits);
void decode_mwbf(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, int q_bits, double alpha);
void decode_imwbf(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, int q_bits, double alpha, double *lp1, int *min_indx);


unsigned short bin2gray(unsigned short num);

unsigned short gray2bin(unsigned short num);


void dec2bin(int dec, int d, char* bin);
int bin2dec(int d, char* bin);
void LLR_cal(double voltage, struct CMFM *cm, double *d_tmp, double *llr);


// quantization of LLR vector
int llr_quantize(const double* llr, int* llr_q, int n, double q_step, int max_q_level);


#endif
