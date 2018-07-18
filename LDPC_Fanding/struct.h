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

