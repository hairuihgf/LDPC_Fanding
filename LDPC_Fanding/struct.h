#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>


struct matrix
{
	int  message_num;
	int  col_num;
	int  row_num;
	int  max_colweigh;
	int  max_rowweigh;
	int *col_weigh;
	int *row_weigh;
	int **col;
	int **row;
};

struct   AWGN
{
	long      reset;
	long      ix;
	long      iy;
	long      iz;
	double     snr;
	double     sigma;
	double     I_sigma;      //(1/(sigma*sigma))
	double     rate;
};

struct profile
{
	char	matrixfile[200];
	int		test_num;
	double  snr;
	double  snr_step;
	int		seed1;
	int		seed2;
	int		seed3;
	int		leastframe;
	int		leasterror;
	int		codemethod;
	int		max_iter;
	int		err_era;                 // Generate Errors '0' or Erasures '1'
	int		burstlen;
	int		burstnum;
	int		de_alg;					  // Algorithm Type: '0' BP; '1' DD-BMP
	float	delta;					  // d of DD-BMP
	float	scale;
	int		codeword;
	int		q_bits;					  // Quantizing Bits
	double	alpha;					  // Alpha Factor
	double	q_step;					//quantized step
	double	d_quasi;                  //parameter of (q+1)-bit quasi-uniform quantizer
	int		max_q_level;				  //maximum value of quantization level
	int		dsum;
	int		modulate_mode;
	double	DOPPLER;
	double	T_SAMPLE;
	int		RAYLEIGH_M;
	int		channel_style;
};
struct v_node  //Variable Node
{
	struct c_node **link; //Connected Check Node
	int sybol;            //Node Symbol
	int step;             //CN Number
	double *r_value;     //CN Values

	int *m_value;		  //M Values for DD BMP
};



struct c_node
{
	struct v_node **link;		//Connected Variable
	int sybol;					//Node Symbol
	int step;					//CV Numbers
	double *q_value;			//CV Values
	int *m_value;               //M Values for DD BMP
};

