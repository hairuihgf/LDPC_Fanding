#include <string.h>
#include <stdlib.h>
#include "function.h"

void ReadProfile(
	AWGN    p, 
	AWGN    *pi,
	char	matrixfile[200], 
	int		test_num,
	int		leastframe,
	int		leasterror,
	int		codemethod,
	int		codeword,
	int		max_iter,
	int		err_era,                  // Generate Errors '0' or Erasures '1'
	int		burstlen, 
	int		burstnum,
	int		de_alg,					  // Algorithm Type: '0' BP; '1' DD-BMP
	float	delta,					  // d of DD-BMP
	float	scale,
	int		q_bits,					  // Quantizing Bits
	double	alpha,					  // Alpha Factor
	double	d_quasi,                  //parameter of (q+1)-bit quasi-uniform quantizer
	int		dsum,
	int		max_q_level,				  //maximum value of quantization level
	double	q_step,					//quantized step
	double  snr_step
	)
{
	FILE	*fin;
	char	ch;
	float	readtemp;

	fopen_s(&fin, "profile.txt", "r");

	/*Read Profile*/
	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%s", matrixfile);


	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &test_num);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%f", &readtemp);
	p.snr = (double)readtemp;

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%f", &readtemp);
	snr_step = (double)readtemp;

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d ", &p.ix);
	//seed1 = p.ix;
	fscanf(fin, "%d ", &p.iy);
	//seed2 = p.iy;
	fscanf(fin, "%d ", &p.iz);
	//seed3 = p.iz;


	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &leastframe);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &leasterror);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &codemethod);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &max_iter);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &err_era);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &burstlen);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &burstnum);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &de_alg);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%f", &readtemp);
	delta = readtemp;

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%f", &scale);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &codeword);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &q_bits);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%f", &readtemp);
	alpha = readtemp;


	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%f", &readtemp);
	q_step = readtemp;

	fscanf(fin, "%f", &readtemp);
	d_quasi = readtemp;

	fscanf(fin, "%d", &max_q_level);

	fscanf(fin, "%d", &dsum);

	fclose(fin);
	fin = NULL;
}