#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "randomc.h"
#include "struct.h"
#include "function.h"

//CRandomMersenne rg(179);
int main()
{
	int i;
	matrix sparse;
	char **gauss;
	int *exchange;
	int rank;
	char matrixfile[200] = "2048_384_mackay.code";
	//读出参数
	FILE	*fin;
	AWGN    p;
	AWGN    *pi;
	char	ch;
	char	matrixfile[200];
	int		test_num;
	float	readtemp;
	int		leastframe;
	int		leasterror;
	int		codemethod;
	int		codeword;
	int		max_iter;
	int		err_era;                 // Generate Errors '0' or Erasures '1'
	int		burstlen;
	int		burstnum;
	int		de_alg;					  // Algorithm Type: '0' BP; '1' DD-BMP
	float	delta;					  // d of DD-BMP
	float	scale;
	int		q_bits;					  // Quantizing Bits
	double	alpha;					  // Alpha Factor
	double	d_quasi;                  //parameter of (q+1)-bit quasi-uniform quantizer
	int		dsum;
	int		max_q_level;				  //maximum value of quantization level
	double	q_step;					//quantized step
	double  snr_step;

	ReadProfile(p,*pi,);

	

	//读出校验矩阵并高斯消去
	init_sparse(&sparse, matrixfile);

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

}