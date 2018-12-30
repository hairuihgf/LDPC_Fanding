#include "function.h"
#include <time.h>
#include <math.h>
#include <random>

#define LLR_CLIP 100
#define A_THRSD 15
#define B_THRSD 20


int trap_set[A_THRSD][B_THRSD] = { 0 };

// quantization of LLR vector
// double rounding
//double round_q(double x, double q_step, int max_q_level)
//{
//	int y;
//
//	if(x>0)
//		y=(int)(x/q_step+0.5);
//	else
//		y=(int)(x/q_step-0.5);
//
//	y=std::min(y,max_q_level-1);
//	y=std::max(y,-max_q_level);
//
//	return y;
//}


double min_s(double *x, int y, int z)             //min
{
	int i;
	int j = 1;
	double temp;
	if (y != 1)
	{
		temp = x[1];
		i = 1;
		for (i; i<y - 1; i++)
		{
			if (fabs(temp)>fabs(x[i + 1]))
				temp = fabs(x[i + 1]);
			else
				temp = fabs(temp);

			j++;
		}
		j = j + 1;
		for (j; j<z; j++)
		{
			if (fabs(temp)>fabs(x[j + 1]))
				temp = fabs(x[j + 1]);
			else
				temp = fabs(temp);
		}
		return temp;
	}
	else
	{
		temp = x[2];
		j = 2;
		for (j; j<z; j++)
		{
			if (fabs(temp)>fabs(x[j + 1]))
				temp = fabs(x[j + 1]);
			else
				temp = fabs(temp);
		}
		return temp;

	}

}

int sgn(double x)
{
	int y;

	if (x>0)
		y = 1;
	else if (x == 0)
		y = 1 - ((clock() & 1) << 1);	// Random
	else
		y = -1;

	return y;
}

double LLR_add(double in1, double in2)
{
	double	result;
	double	tp1 = 0, tp2 = 0;
	double	max1, max2;

	if (in1<in2)	max1 = in2;
	else		max1 = in1;

	if ((in1 + in2)<0)	max2 = 0;
	else			max2 = in1 + in2;

	//tp1 = exp(-fabs(in1-in2));	//if 0 modulate -1 and 1modulate 1
	tp1 = exp(-fabs(in1 + in2));		//if 0 modulate 1 and 1modulate -1

	tp1 = log(1 + tp1);
	//if 0 modulate -1 and 1modulate 1
	//tp2 = exp(-fabs(in1+in2));
	tp2 = exp(-fabs(in1 - in2));		//if 0 modulate 1 and 1modulate -1

	tp2 = log(1 + tp2);

	//result = max1 - max2 + tp1 - tp2;//if 0 modulate -1 and 1modulate 1
	result = -max1 + max2 + tp1 - tp2;//if 0 modulate 1 and 1 modulate -1

	return(result);
}

double LLR_sub(double in1, double in2)
{
	double	result;
	double	tp1 = 0, tp2 = 0;
	double	max1, max2;

	if (fabs(in1) == fabs(in2))
	{
		if (in1 >= 0) in1 = in1 + 0.000000000000000000000000000000000000000000000000000000000000000001;
		else in1 = in1 - 0.000000000000000000000000000000000000000000000000000000000000000001;
	}


	if (in1<in2)	max1 = in2;
	else		max1 = in1;

	if ((in1 + in2)<0)	max2 = 0;
	else			max2 = in1 + in2;

	//tp1 = exp(-fabs(in1-in2));	//if 0 modulate -1 and 1modulate 1
	tp1 = exp(-fabs(in1 + in2));		//if 0 modulate 1 and 1modulate -1
	if (tp1>0.9999999999999999)
		tp1 = 0.9999999999999999;
	tp1 = log(1 - tp1);

	//tp2 = exp(-fabs(in1+in2));	//if 0 modulate -1 and 1modulate 1
	tp2 = exp(-fabs(in1 - in2));		//if 0 modulate 1 and 1modulate -1
	if (tp2>0.9999999999999999)
		tp2 = 0.9999999999999999;
	tp2 = log(1 - tp2);

	//result = max1 - max2 + tp1 - tp2;//if 0 modulate -1 and 1modulate 1
	result = -max1 + max2 + tp1 - tp2;//if 0 modulate 1 and 1modulate -1

	return(result);
}


void check_result2(struct matrix *sparse, char *resultcode, double *lp, bool &decode_valid)
{
	int i, j, temp;

	decode_valid = true;

	for (i = 1; i <= sparse->col_num; i++)
	{
		if (lp[i] >= 0.0) resultcode[i] = 0;
		else resultcode[i] = 1;
	}

	for (i = 0; i<sparse->row_num; i++)
	{
		temp = 0;
		for (j = 0; j<sparse->row_weigh[i]; j++)
		{
			temp = (temp + resultcode[sparse->row[i][j]]);
			if (temp == 2)
				temp = 0;
		}

		if (temp == 0) {}
		else
		{
			decode_valid = false;
			break;
		}
	}
}



//SPA
void decodeon(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter)
{
	int i, j, k, k1;
	double temp, temp1;

	for (i = 1; i <= sparse->col_num; i++)			//Initiate
		resultcode[i] = 0;

	for (i = 0; i <= sparse->col_num; i++)
		lp[i] = 0.0;

	for (i = 0; i <= sparse->col_num; i++)
		lp0[i] = 0.0;

	for (i = 1; i <= sparse->col_num; i++)
		lp0[i] = transcode[i];

	for (i = 1; i <= sparse->col_num; i++)
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			v_line[i].link[j] = &c_line[sparse->col[i - 1][j - 1]];

	for (i = 1; i <= sparse->row_num; i++)
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
			c_line[i].link[j] = &v_line[sparse->row[i - 1][j - 1]];

	for (i = 1; i <= sparse->col_num; i++)		//Clean
	{
		v_line[i].step = 1;
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			v_line[i].r_value[j] = 0;
	}

	for (i = 1; i <= sparse->row_num; i++)		//Clean
	{
		c_line[i].step = 1;
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
			c_line[i].q_value[j] = 0;
	}


	for (i = 1; i <= sparse->col_num; i++)
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
		{
			k = v_line[i].link[j]->step;
			v_line[i].link[j]->q_value[k] = lp0[i];
			v_line[i].link[j]->step++;
		}

	for (k1 = 0; k1<max_iter; k1++)
	{
		iter++;
		for (i = 1; i <= sparse->row_num; i++) //Step Initiate
			c_line[i].step = 1;
		for (i = 1; i <= sparse->col_num; i++)
			v_line[i].step = 1;

		for (i = sparse->row_num; i >= 1; i--) //Calculate r
		{
			temp = c_line[i].q_value[1];
			for (j = 2; j <= sparse->row_weigh[i - 1]; j++)
				temp = LLR_add_1(temp, c_line[i].q_value[j]);

			for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
			{
				k = c_line[i].link[j]->step;
				c_line[i].link[j]->link[k] = &c_line[i];
				c_line[i].link[j]->r_value[k] = LLR_sub_1(temp, c_line[i].q_value[j]);
				c_line[i].link[j]->step++;
			}
		}

		for (i = 1; i <= sparse->col_num; i++)
			v_line[i].step = 1;

		for (i = 1; i <= sparse->col_num; i++)  //Calculate lp
		{
			temp = 0;
			for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
				temp = temp + v_line[i].r_value[j];
			lp[i] = lp0[i] + temp;
		}

		for (i = 1; i <= sparse->row_num; i++)
			c_line[i].step = 1;

		for (i = sparse->col_num; i >= 1; i--)   //Calculate q
		{
			for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			{
				k = v_line[i].link[j]->step;
				v_line[i].link[j]->link[k] = &v_line[i];
				temp1 = lp[i] - v_line[i].r_value[j];
				v_line[i].link[j]->q_value[k] = temp1;
				v_line[i].link[j]->step++;
			}
		}

		check_result2(sparse, resultcode, lp, decode_valid);

		if (decode_valid)              //Parity Check
		{
			break;
		}
	}



}





//MIN-SUM

double sisym(double x)
{
	double y;
	if (x == 0)  //当x==0时，函数返回0
		y = 0;
	else if (x > 0)  //当x>0时，函数返回1
		y = 1;
	else  //当x<0时，函数返回－1
		y = -1;

	return y;
}


//MIN-SUM 非均匀量化
void decode_quasi_ms(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, double q_step, int q_bits, double d_quasi)
{
	int i, j, k, k1;
	double temp, temp1;
	double min_1, min_2;		//Minimum and the second minimum value

	int  min_indx;		//The sign of check-sum, the index of the minimum number

	for (i = 1; i <= sparse->col_num; i++)			//Initiate
		resultcode[i] = 0;

	for (i = 0; i <= sparse->col_num; i++)
		lp[i] = 0.0;

	for (i = 0; i <= sparse->col_num; i++)
		lp0[i] = 0.0;

	for (i = 1; i <= sparse->col_num; i++)
		lp0[i] = transcode[i];

	/*lp  y_i,quantization*/

	for (i = 1; i <= sparse->col_num; i++)                 //R(0)=q
	{
		lp0[i] = (int)(transcode[i] / q_step);
		if (transcode[i] - lp0[i] * q_step>(lp0[i] + 1)*q_step - transcode[i]) {
			lp0[i]++;
		}
		if (lp0[i]>(1 << (q_bits - 1)) - 1) {
			lp0[i] = (1 << (q_bits - 1)) - 1;
		}
		if (lp0[i]<-(1 << (q_bits - 1)) + 1) {
			lp0[i] = -(1 << (q_bits - 1)) + 1;
		}
	}

	for (i = 1; i <= sparse->col_num; i++) {
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++) {
			v_line[i].link[j] = &c_line[sparse->col[i - 1][j - 1]];
		}
	}

	for (i = 1; i <= sparse->row_num; i++) {
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++) {
			c_line[i].link[j] = &v_line[sparse->row[i - 1][j - 1]];
		}
	}

	for (i = 1; i <= sparse->col_num; i++)		//Clean
	{
		v_line[i].step = 1;
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++) {
			v_line[i].r_value[j] = 0;
			v_line[i].m_value[j] = 0;
		}
	}

	for (i = 1; i <= sparse->row_num; i++)		//Clean
	{
		c_line[i].step = 1;
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++) {
			c_line[i].q_value[j] = 0;
			c_line[i].m_value[j] = 0;
		}
	}


	for (i = 1; i <= sparse->col_num; i++) {
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
		{
			k = v_line[i].link[j]->step;
			v_line[i].link[j]->q_value[k] = lp0[i];
			v_line[i].link[j]->m_value[k] = sgn(lp0[i]);
			v_line[i].link[j]->step++;
		}
	}
	for (k1 = 0; k1<max_iter; k1++)
	{
		iter++;
		for (i = 1; i <= sparse->row_num; i++) { //Step Initiate
			c_line[i].step = 1;
		}
		for (i = 1; i <= sparse->col_num; i++) {
			v_line[i].step = 1;

		}
		for (i = sparse->row_num; i >= 1; i--) //Calculate r
		{
			min_1 = fabs(c_line[i].q_value[1]);
			min_indx = 1;

			min_2 = 1e10;
			temp = c_line[i].m_value[1];                //Calculate a
			for (j = 2; j <= sparse->row_weigh[i - 1]; j++)
			{
				temp = temp*c_line[i].m_value[j];

			}

			for (j = 2; j <= sparse->row_weigh[i - 1]; j++)
			{
				if (min_1>fabs(c_line[i].q_value[j]))
				{
					min_2 = min_1;
					min_1 = fabs(c_line[i].q_value[j]);
					min_indx = j;

				}
				else
				{
					if (min_2>fabs(c_line[i].q_value[j]))
					{
						min_2 = fabs(c_line[i].q_value[j]);
					}

				}
			}

			min_1 *= delta;
			min_2 *= delta;

			for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
			{
				k = c_line[i].link[j]->step;
				c_line[i].link[j]->link[k] = &c_line[i];
				c_line[i].link[j]->m_value[k] = (int)(temp*c_line[i].m_value[j]);
				if (j != min_indx) {
					c_line[i].link[j]->r_value[k] = min_1*c_line[i].link[j]->m_value[k];
				}
				else {
					c_line[i].link[j]->r_value[k] = min_2*c_line[i].link[j]->m_value[k];

				}
				c_line[i].link[j]->step++;
			}

		}

		for (i = 1; i <= sparse->col_num; i++)
			v_line[i].step = 1;

		for (i = 1; i <= sparse->col_num; i++)  //Calculate lp
		{
			temp = 0;
			for (j = 1; j <= sparse->col_weigh[i - 1]; j++) {
				temp = temp + v_line[i].r_value[j];
			}
			lp[i] = lp0[i] + temp;

		}

		for (i = 1; i <= sparse->row_num; i++)
			c_line[i].step = 1;

		for (i = sparse->col_num; i >= 1; i--)   //Calculate q
		{
			for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			{
				k = v_line[i].link[j]->step;
				v_line[i].link[j]->link[k] = &v_line[i];
				temp1 = lp[i] - v_line[i].r_value[j];
				if (temp1>(1 << (q_bits - 1)) - 1) {
					temp1 = (1 << (q_bits - 1)) - 1;
				}
				if (temp1<-(1 << (q_bits - 1)) + 1) {
					temp1 = -(1 << (q_bits - 1)) + 1;
				}
				v_line[i].link[j]->q_value[k] = temp1;
				v_line[i].link[j]->m_value[k] = sgn(temp1);
				v_line[i].link[j]->step++;
			}
		}

		check_result2(sparse, resultcode, lp, decode_valid);

		if (decode_valid)              //Parity Check
		{
			break;
		}
	}
}



void decode_ms(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, double q_step, int q_bits)
{
	int i, j, k, k1;
	double temp, temp1;
	double min_1, min_2;		//Minimum and the second minimum value

	int  min_indx;		//The sign of check-sum, the index of the minimum number

	for (i = 1; i <= sparse->col_num; i++) {			//Initiate
		resultcode[i] = 0;
	}

	for (i = 0; i <= sparse->col_num; i++) {
		lp[i] = 0.0;
	}

	for (i = 1; i <= sparse->col_num; i++) {
		lp0[i] = transcode[i];
	}

	for (i = 1; i <= sparse->col_num; i++) {
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++) {
			v_line[i].link[j] = &c_line[sparse->col[i - 1][j - 1]];
		}
	}

	for (i = 1; i <= sparse->row_num; i++) {
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++) {
			c_line[i].link[j] = &v_line[sparse->row[i - 1][j - 1]];
		}
	}

	for (i = 1; i <= sparse->col_num; i++)		//Clean
	{
		v_line[i].step = 1;
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++) {
			v_line[i].r_value[j] = 0;
			v_line[i].m_value[j] = 0;
		}
	}

	for (i = 1; i <= sparse->row_num; i++)		//Clean
	{
		c_line[i].step = 1;
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++) {
			c_line[i].q_value[j] = 0;
			c_line[i].m_value[j] = 0;
		}
	}


	for (i = 1; i <= sparse->col_num; i++) {
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
		{
			k = v_line[i].link[j]->step;
			v_line[i].link[j]->q_value[k] = lp0[i];
			v_line[i].link[j]->m_value[k] = sgn(lp0[i]);
			v_line[i].link[j]->step++;
		}
	}

	for (k1 = 0; k1<max_iter; k1++)
	{
		iter++;
		for (i = 1; i <= sparse->row_num; i++) //Step Initiate
			c_line[i].step = 1;
		for (i = 1; i <= sparse->col_num; i++)
			v_line[i].step = 1;

		for (i = sparse->row_num; i >= 1; i--) //Calculate r
		{
			min_1 = fabs(c_line[i].q_value[1]);
			min_indx = 1;

			min_2 = 1e10;
			temp = c_line[i].m_value[1];                //Calculate a
			for (j = 2; j <= sparse->row_weigh[i - 1]; j++)
			{
				temp = temp*c_line[i].m_value[j];

			}

			for (j = 2; j <= sparse->row_weigh[i - 1]; j++)
			{
				if (min_1>fabs(c_line[i].q_value[j]))
				{
					min_2 = min_1;
					min_1 = fabs(c_line[i].q_value[j]);
					min_indx = j;

				}
				else
				{
					if (min_2>fabs(c_line[i].q_value[j]))
					{
						min_2 = fabs(c_line[i].q_value[j]);
					}

				}
			}

			for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
			{
				k = c_line[i].link[j]->step;
				c_line[i].link[j]->link[k] = &c_line[i];
				c_line[i].link[j]->m_value[k] = (int)(temp*c_line[i].m_value[j]);
				if (j != min_indx) {
					c_line[i].link[j]->r_value[k] = min_1;
				}
				else {
					c_line[i].link[j]->r_value[k] = min_2;

				}
				c_line[i].link[j]->r_value[k] *= delta*c_line[i].link[j]->m_value[k];
				c_line[i].link[j]->step++;
			}

		}

		for (i = 1; i <= sparse->col_num; i++)
			v_line[i].step = 1;

		for (i = 1; i <= sparse->col_num; i++)  //Calculate lp
		{
			temp = 0;
			for (j = 1; j <= sparse->col_weigh[i - 1]; j++) {
				temp = temp + v_line[i].r_value[j];
			}
			lp[i] = lp0[i] + temp;
		}

		for (i = 1; i <= sparse->row_num; i++)
			c_line[i].step = 1;

		for (i = sparse->col_num; i >= 1; i--)   //Calculate q
		{
			for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			{
				k = v_line[i].link[j]->step;
				v_line[i].link[j]->link[k] = &v_line[i];
				temp1 = lp[i] - v_line[i].r_value[j];
				v_line[i].link[j]->q_value[k] = temp1;
				v_line[i].link[j]->m_value[k] = sgn(temp1);
				v_line[i].link[j]->step++;
			}
		}

		check_result2(sparse, resultcode, lp, decode_valid);

		if (decode_valid)              //Parity Check
		{
			break;
		}
	}
}


//downsampling min-sum
void decode_dsms(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, double q_step, int q_bits, int dsum)
{
	int i, j, k, k1;
	double temp, temp1;
	double min_1, min_2;		//Minimum and the second minimum value

	int  min_indx;		//The sign of check-sum, the index of the minimum number

	for (i = 1; i <= sparse->col_num; i++) {			//Initiate
		resultcode[i] = 0;
	}

	for (i = 0; i <= sparse->col_num; i++) {
		lp[i] = 0.0;
	}


	for (i = 1; i <= sparse->col_num; i++) {
		lp0[i] = transcode[i];
	}

	for (i = 1; i <= sparse->col_num; i++) {
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++) {
			v_line[i].link[j] = &c_line[sparse->col[i - 1][j - 1]];
		}
	}

	for (i = 1; i <= sparse->row_num; i++) {
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++) {
			c_line[i].link[j] = &v_line[sparse->row[i - 1][j - 1]];
		}
	}

	for (i = 1; i <= sparse->col_num; i++)		//Clean
	{
		v_line[i].step = 1;
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++) {
			v_line[i].r_value[j] = 0;
			v_line[i].m_value[j] = 0;
		}
	}

	for (i = 1; i <= sparse->row_num; i++)		//Clean
	{
		c_line[i].step = 1;
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++) {
			c_line[i].q_value[j] = 0;
			c_line[i].m_value[j] = 0;
		}
	}


	for (i = 1; i <= sparse->col_num; i++) {
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
		{
			k = v_line[i].link[j]->step;
			v_line[i].link[j]->q_value[k] = lp0[i];
			v_line[i].link[j]->m_value[k] = sgn(lp0[i]);
			v_line[i].link[j]->step++;
		}
	}

	for (k1 = 0; k1<max_iter; k1++)
	{
		iter++;
		for (i = 1; i <= sparse->row_num; i++) //Step Initiate
			c_line[i].step = 1;
		for (i = 1; i <= sparse->col_num; i++)
			v_line[i].step = 1;

		for (i = sparse->row_num; i >= 1; i--) //Calculate r
		{
			min_1 = fabs(c_line[i].q_value[1]);
			min_indx = 1;

			min_2 = 1e10;
			temp = c_line[i].m_value[1];                //Calculate a
			for (j = 2; j <= sparse->row_weigh[i - 1]; j++)
			{
				temp = temp*c_line[i].m_value[j];

			}

			for (j = 2; j <= sparse->row_weigh[i - 1]; j++)
			{
				if (min_1>fabs(c_line[i].q_value[j]))
				{
					min_2 = min_1;
					min_1 = fabs(c_line[i].q_value[j]);
					min_indx = j;

				}
				else
				{
					if (min_2>fabs(c_line[i].q_value[j]))
					{
						min_2 = fabs(c_line[i].q_value[j]);
					}

				}
			}

			for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
			{
				k = c_line[i].link[j]->step;
				c_line[i].link[j]->link[k] = &c_line[i];
				c_line[i].link[j]->m_value[k] = (int)(temp*c_line[i].m_value[j]);
				if (j != min_indx) {
					c_line[i].link[j]->r_value[k] = min_1;
				}
				else {
					c_line[i].link[j]->r_value[k] = min_2;

				}
				c_line[i].link[j]->r_value[k] *= delta*c_line[i].link[j]->m_value[k];
				c_line[i].link[j]->step++;
			}

		}

		for (i = 1; i <= sparse->col_num; i++)
			v_line[i].step = 1;

		for (i = 1; i <= sparse->col_num; i++)  //Calculate lp
		{
			temp = 0;
			for (j = 1; j <= sparse->col_weigh[i - 1]; j++) {
				temp = temp + v_line[i].r_value[j];
			}
			lp[i] = lp0[i] + temp;
		}

		for (i = 1; i <= sparse->row_num; i++)
			c_line[i].step = 1;

		for (i = sparse->col_num; i >= 1; i--)   //Calculate q
		{
			for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			{
				k = v_line[i].link[j]->step;
				v_line[i].link[j]->link[k] = &v_line[i];
				temp1 = lp[i] - v_line[i].r_value[j];
				v_line[i].link[j]->q_value[k] = temp1;
				v_line[i].link[j]->m_value[k] = sgn(temp1);
				v_line[i].link[j]->step++;
			}
		}

		check_result2(sparse, resultcode, lp, decode_valid);

		if (decode_valid)              //Parity Check
		{
			break;
		}
	}
}


// quantized min-sum
void decode_qua_ms(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, double q_step, int q_bits)
{
	int i, j, k, k1;
	double temp, temp1;
	double min_1, min_2;		//Minimum and the second minimum value

	int  min_indx;		//The sign of check-sum, the index of the minimum number

	for (i = 1; i <= sparse->col_num; i++)			//Initiate
		resultcode[i] = 0;

	for (i = 0; i <= sparse->col_num; i++)
		lp[i] = 0.0;

	for (i = 0; i <= sparse->col_num; i++)
		lp0[i] = 0.0;

	for (i = 1; i <= sparse->col_num; i++)
		lp0[i] = transcode[i];

	/*lp  y_i,quantization*/

	for (i = 1; i <= sparse->col_num; i++)                 //R(0)=q
	{
		lp0[i] = (int)(transcode[i] / q_step);
		if (transcode[i] - lp0[i] * q_step>(lp0[i] + 1)*q_step - transcode[i]) {
			lp0[i]++;
		}
		if (lp0[i]>(1 << (q_bits - 1)) - 1) {
			lp0[i] = (1 << (q_bits - 1)) - 1;
		}
		if (lp0[i]<-(1 << (q_bits - 1)) + 1) {
			lp0[i] = -(1 << (q_bits - 1)) + 1;
		}
	}

	for (i = 1; i <= sparse->col_num; i++) {
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++) {
			v_line[i].link[j] = &c_line[sparse->col[i - 1][j - 1]];
		}
	}

	for (i = 1; i <= sparse->row_num; i++) {
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++) {
			c_line[i].link[j] = &v_line[sparse->row[i - 1][j - 1]];
		}
	}

	for (i = 1; i <= sparse->col_num; i++)		//Clean
	{
		v_line[i].step = 1;
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++) {
			v_line[i].r_value[j] = 0;
			v_line[i].m_value[j] = 0;
		}
	}

	for (i = 1; i <= sparse->row_num; i++)		//Clean
	{
		c_line[i].step = 1;
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++) {
			c_line[i].q_value[j] = 0;
			c_line[i].m_value[j] = 0;
		}
	}


	for (i = 1; i <= sparse->col_num; i++) {
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
		{
			k = v_line[i].link[j]->step;
			v_line[i].link[j]->q_value[k] = lp0[i];
			v_line[i].link[j]->m_value[k] = sgn(lp0[i]);
			v_line[i].link[j]->step++;
		}
	}
	for (k1 = 0; k1<max_iter; k1++)
	{
		iter++;
		for (i = 1; i <= sparse->row_num; i++) { //Step Initiate
			c_line[i].step = 1;
		}
		for (i = 1; i <= sparse->col_num; i++) {
			v_line[i].step = 1;

		}
		for (i = sparse->row_num; i >= 1; i--) //Calculate r
		{
			min_1 = fabs(c_line[i].q_value[1]);
			min_indx = 1;

			min_2 = 1e10;
			temp = c_line[i].m_value[1];                //Calculate a
			for (j = 2; j <= sparse->row_weigh[i - 1]; j++)
			{
				temp = temp*c_line[i].m_value[j];

			}

			for (j = 2; j <= sparse->row_weigh[i - 1]; j++)
			{
				if (min_1>fabs(c_line[i].q_value[j]))
				{
					min_2 = min_1;
					min_1 = fabs(c_line[i].q_value[j]);
					min_indx = j;

				}
				else
				{
					if (min_2>fabs(c_line[i].q_value[j]))
					{
						min_2 = fabs(c_line[i].q_value[j]);
					}

				}
			}

			min_1 *= delta;
			min_2 *= delta;

			for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
			{
				k = c_line[i].link[j]->step;
				c_line[i].link[j]->link[k] = &c_line[i];
				c_line[i].link[j]->m_value[k] = (int)(temp*c_line[i].m_value[j]);
				if (j != min_indx) {
					c_line[i].link[j]->r_value[k] = min_1*c_line[i].link[j]->m_value[k];
				}
				else {
					c_line[i].link[j]->r_value[k] = min_2*c_line[i].link[j]->m_value[k];

				}
				c_line[i].link[j]->step++;
			}

		}

		for (i = 1; i <= sparse->col_num; i++)
			v_line[i].step = 1;

		for (i = 1; i <= sparse->col_num; i++)  //Calculate lp
		{
			temp = 0;
			for (j = 1; j <= sparse->col_weigh[i - 1]; j++) {
				temp = temp + v_line[i].r_value[j];
			}
			lp[i] = lp0[i] + temp;

		}

		for (i = 1; i <= sparse->row_num; i++)
			c_line[i].step = 1;

		for (i = sparse->col_num; i >= 1; i--)   //Calculate q
		{
			for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			{
				k = v_line[i].link[j]->step;
				v_line[i].link[j]->link[k] = &v_line[i];
				temp1 = lp[i] - v_line[i].r_value[j];
				if (temp1>(1 << (q_bits - 1)) - 1) {
					temp1 = (1 << (q_bits - 1)) - 1;
				}
				if (temp1<-(1 << (q_bits - 1)) + 1) {
					temp1 = -(1 << (q_bits - 1)) + 1;
				}
				v_line[i].link[j]->q_value[k] = temp1;
				v_line[i].link[j]->m_value[k] = sgn(temp1);
				v_line[i].link[j]->step++;
			}
		}

		check_result2(sparse, resultcode, lp, decode_valid);

		if (decode_valid)              //Parity Check
		{
			break;
		}
	}
}


//MIN-SUM by Mass
void vMSDecoder(struct matrix *S_Sparse, char *szResultCode, double **dppCNValue, int **ippVNSgn, double **dppVNValue, double *dpLLRTotal, double *dpChannelLLR, double *lp0, int &iIter, bool &bDecodeValid, int iMaxIter, double delta, char *szCode, int q_bits, double q_step, double d_quasi)
{
	/*	double **dppCNValue;
	double **dppVNValue;
	double *dpLLRTotal;*/
	int iMinIndex;
	int i2ndMinIndex;
	int iExtrSgn;
	int i, j, k, iter;
	//	DWORD dwSysClk;

	//	dwSysClk = GetTickCount();

	int d_0 = (1 << (q_bits - 2)) - 1;
	FILE *dp;
	FILE *fp;

	int iBE;
	iBE = 0;
	check_result2(S_Sparse, szResultCode, dpChannelLLR, bDecodeValid);

	for (i = 1; i <= S_Sparse->col_num; i++) {
		/*lp0[i] = (int)(((int)(abs(dpChannelLLR[i])/(q_step/2))+1)/2);
		if(lp0[i]>pow(2.0,q_bits-1)-1){
		lp0[i] =pow(2.0,q_bits-1)-1;
		}else;*/

		/*if(lp0[i]>d_0){
		if(lp0[i]<d_quasi*d_0)
		lp0[i]=d_0;
		else if(lp0[i]>=pow(d_quasi,d_0+1)*d_0)
		lp0[i]=pow(d_quasi,d_0+1)*d_0;
		else{
		for(int r=1;r<=d_0;r++){
		if(lp0[i]>=pow(d_quasi,r)*d_0&&dpLLRTotal[i]<pow(d_quasi,r+1)*d_0){
		lp0[i]=pow(d_quasi,r)*d_0;
		break;
		}
		}
		}
		}else;*/

		lp0[i] = (dpChannelLLR[i]/*>0?1:-1*/);
	}
	/*dp=fopen("quantization of 4095 code.txt","a");
	fprintf(dp,"bits:    transcode:     quantization:   \n");
	for(i=1;i<=S_Sparse->col_num;i++)
	{
	fprintf(dp,"%5d     ",i);
	fprintf(dp,"%6.3f        ",dpChannelLLR[i]);
	fprintf(dp,"%5.0f ",lp0[i]);
	fprintf(dp,"\n");
	}
	fclose(dp);*/

	//VN value initial.
	for (i = 0; i < S_Sparse->col_num; i++)
	{
		for (j = 0; j < S_Sparse->col_weigh[i]; j++)
		{
			dppVNValue[i][S_Sparse->col[i][j] - 1] = lp0[i + 1];
			//			ippVNSgn[i][S_Sparse->col[i][j] - 1] = ((dpChannelLLR[i + 1] >= 0) ? 0 : 1);		
		}
	}


	for (iter = 0; iter < iMaxIter; iter++)
	{
		iIter++;
		//CN value update.
		for (i = 0; i < S_Sparse->row_num; i++)
		{
			//Find the minimum and second-minimum LLRs' index.
			if (fabs(dppVNValue[S_Sparse->row[i][0] - 1][i]) < fabs(dppVNValue[S_Sparse->row[i][1] - 1][i]))
			{
				iMinIndex = S_Sparse->row[i][0] - 1;
				i2ndMinIndex = S_Sparse->row[i][1] - 1;
			}
			else
			{
				iMinIndex = S_Sparse->row[i][1] - 1;
				i2ndMinIndex = S_Sparse->row[i][0] - 1;
			}

			for (j = 2; j < S_Sparse->row_weigh[i]; j++)
			{
				if (fabs(dppVNValue[S_Sparse->row[i][j] - 1][i]) < fabs(dppVNValue[iMinIndex][i]))
				{
					i2ndMinIndex = iMinIndex;
					iMinIndex = S_Sparse->row[i][j] - 1;
				}
				else if (fabs(dppVNValue[S_Sparse->row[i][j] - 1][i]) < fabs(dppVNValue[i2ndMinIndex][i]))
				{
					i2ndMinIndex = S_Sparse->row[i][j] - 1;
				}
				else;
			}//End of the for(j=2;...) of minimum index
			 //Pass the extrinsic imformation one by one.
			for (j = 0; j < S_Sparse->row_weigh[i]; j++)
			{
				//Pass the minimum value
				if (S_Sparse->row[i][j] - 1 != iMinIndex)
				{
					dppCNValue[i][S_Sparse->row[i][j] - 1] = fabs(dppVNValue[iMinIndex][i]);
				}
				else
				{
					dppCNValue[i][S_Sparse->row[i][j] - 1] = fabs(dppVNValue[i2ndMinIndex][i]);
				}
				/*				//Pass the sign.
				iExtrSgn = 0;
				for(k = 0; k < S_Sparse->row_weigh[i]; k++)
				{
				iExtrSgn += ippVNSgn[S_Sparse->row[i][k] - 1][i];
				}//End of for(k=0;...) of calculate the sign.
				iExtrSgn += ippVNSgn[S_Sparse->row[i][j] - 1][i];
				iExtrSgn = ((iExtrSgn % 2 == 0) ? 1 : -1);
				*/
				///*
				iExtrSgn = 0;
				for (k = 0; k < S_Sparse->row_weigh[i]; k++)
				{
					if (dppVNValue[S_Sparse->row[i][k] - 1][i] < 0)
					{
						iExtrSgn++;
					}
				}//End of for(k=0;...) of calculate the sign.
				iExtrSgn += (dppVNValue[S_Sparse->row[i][j] - 1][i] >= 0 ? 0 : 1);
				iExtrSgn = ((iExtrSgn % 2 == 0) ? 1 : -1);
				//*/
				dppCNValue[i][S_Sparse->row[i][j] - 1] *= delta * iExtrSgn;

			}//End of for(j=0;..;..) of extrinsic imformation passing.
		}//End of for(i=0;..;..) of CN update.

		 //Total LLR calculation.
		for (i = 1; i <= S_Sparse->col_num; i++)
		{
			dpLLRTotal[i] = lp0[i];

			for (j = 0; j < S_Sparse->col_weigh[i - 1]; j++)
			{
				dpLLRTotal[i] += dppCNValue[S_Sparse->col[i - 1][j] - 1][i - 1];
			}
		}
		//Parity Check
		check_result2(S_Sparse, szResultCode, dpLLRTotal, bDecodeValid);

		iBE = 0;
		for (i = 1; i <= S_Sparse->col_num; i++)
		{
			if (szResultCode[i] != szCode[i])
			{
				iBE++;
			}
		}

		//if(iter>120){
		//	fp=fopen("trapping set of 4095 code.txt","a");
		//	fprintf(fp,"iter:  %d\n",iter);

		//	fprintf(fp,"VN:   CN:\n");
		//	for(i = 1; i <= S_Sparse->col_num; i++){
		//		if(szResultCode[i] != szCode[i]){
		//			fprintf(fp,"%5d   ",i);
		//			for(j = 0; j < S_Sparse->col_weigh[i]; j++)
		//			{
		//				fprintf(fp,"%5d  ",S_Sparse->col[i-1][j]);
		//			} 
		//			fprintf(fp,"\n");
		//		}else;


		//	}
		//	// fprintf(fp,"\n");
		//	fprintf(fp,"code:   ");
		//	for(i = 1; i <= S_Sparse->col_num; i++){
		//		if(szResultCode[i] != szCode[i]){
		//			fprintf(fp,"%5d ",szCode[i]);
		//		}else;
		//	}
		//	fprintf(fp,"\n");

		//	fprintf(fp,"transcode: ");
		//	for(i = 1; i <= S_Sparse->col_num; i++){
		//		if(szResultCode[i] != szCode[i]){
		//			fprintf(fp,"%6.3f ",dpChannelLLR[i]);
		//		}else;
		//	}
		//	fprintf(fp,"\n");

		//	fprintf(fp,"rj,i:\n");
		//	for(i =1; i <=S_Sparse->col_num; i++)
		//	{
		//		if(szResultCode[i] != szCode[i]){
		//			for(j = 0; j < S_Sparse->col_weigh[i]; j++)
		//			{
		//				fprintf(fp,"    %6.3f ",dppCNValue[S_Sparse->col[i-1][j]-1][i-1]);
		//			} 
		//			fprintf(fp,"\n");
		//		}else;

		//	}

		//	fprintf(fp,"qi,j,c:\n");
		//	for(i =1; i <=S_Sparse->col_num; i++)
		//	{
		//		if(szResultCode[i] != szCode[i]){
		//			for(j = 0; j<S_Sparse->col_weigh[i]; j++)
		//			{
		//				fprintf(fp,"    %6.3f ",dppVNValue[i-1][S_Sparse->col[i-1][j] - 1]);
		//			} 
		//			fprintf(fp,"\n");
		//		}else;

		//	}
		//	fprintf(fp,"Qi:  ");
		//	for(i = 1; i <= S_Sparse->col_num; i++){
		//		if(szResultCode[i] != szCode[i]){
		//			fprintf(fp,"%6.3f ",dpLLRTotal[i]);
		//		}else;
		//	}
		//	fprintf(fp,"\n");
		//	fclose(fp);
		//}

		if (bDecodeValid)
		{
			break;
		}
		//VN update
		for (i = 0; i < S_Sparse->col_num; i++)
		{
			for (j = 0; j < S_Sparse->col_weigh[i]; j++)
			{
				dppVNValue[i][S_Sparse->col[i][j] - 1] = dpLLRTotal[i + 1] - dppCNValue[S_Sparse->col[i][j] - 1][i];
				//				ippVNSgn[i][S_Sparse->col[i][j] - 1] = ((dppVNValue[i][S_Sparse->col[i][j] - 1] >= 0) ? 0 : 1);		
			}
		}
		//printf("\r%d",iter);
	}//End of for(iIter=0;..;..).
	 //	cout << GetTickCount() - ulSysClk << "ms" << endl;


	 /*
	 free(dppVNValue[0]);
	 free(dppVNValue);
	 free(dppCNValue[0]);
	 free(dppCNValue);
	 free(dpLLRTotal);*/
}



void decode_ddbmp(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, double q_step, int q_bits)
{
	int i, j, k, k1;
	double temp;
	double temp1;

	for (i = 1; i <= sparse->col_num; i++)			//Initiate
		resultcode[i] = 0;

	/*Parity Check*/
	check_result2(sparse, resultcode, transcode, decode_valid);

	if (decode_valid)              //Parity Check
	{
		return;
	}

	/*Initiate DD BMP Algorithm*/
	for (i = 0; i <= sparse->col_num; i++)
		lp[i] = 0.0;

	for (i = 0; i <= sparse->col_num; i++) {
		lp0[i] = 0.0;
	}

	//hard decision
	for (i = 1; i <= sparse->col_num; i++) {
		lp0[i] = transcode[i];
	}


	/*Quantizing Yi*/
	for (i = 1; i <= sparse->col_num; i++)
	{
		lp0[i] = (int)((transcode[i]) / q_step);

		if (transcode[i] - lp0[i] * q_step > (lp0[i] + 1)*q_step - transcode[i]) {
			lp0[i]++;
		}

		if (lp0[i]>(1 << (q_bits - 1)) - 1) {
			lp0[i] = (1 << (q_bits - 1)) - 1;
		}

		if (lp0[i]<-(1 << (q_bits - 1)) + 1) {
			lp0[i] = -(1 << (q_bits - 1)) + 1;
		}

	}
	//	printf("\n");

	//Set Up Nodes
	for (i = 1; i <= sparse->col_num; i++)
	{
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
		{
			v_line[i].link[j] = &c_line[sparse->col[i - 1][j - 1]];
		}
	}

	for (i = 1; i <= sparse->row_num; i++)
	{
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
		{
			c_line[i].link[j] = &v_line[sparse->row[i - 1][j - 1]];
		}
	}

	for (i = 1; i <= sparse->col_num; i++)		//Clean
	{
		v_line[i].step = 1;
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
		{
			v_line[i].r_value[j] = 0;
			v_line[i].m_value[j] = 0;
		}
	}

	for (i = 1; i <= sparse->row_num; i++)			//Clean
	{
		c_line[i].step = 1;
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
		{
			c_line[i].q_value[j] = 0;
			c_line[i].m_value[j] = 0;
		}
	}


	for (i = 1; i <= sparse->col_num; i++)
	{
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
		{
			k = v_line[i].link[j]->step;
			v_line[i].link[j]->q_value[k] = lp0[i];
			v_line[i].link[j]->m_value[k] = sgn(lp0[i]);
			v_line[i].link[j]->step++;
		}
	}


	for (k1 = 0; k1<max_iter; k1++)
	{
		iter++;
		for (i = 1; i <= sparse->row_num; i++) //Step Initiate
			c_line[i].step = 1;
		for (i = 1; i <= sparse->col_num; i++)
			v_line[i].step = 1;

		for (i = sparse->row_num; i >= 1; i--) //Calculate r
		{
			temp = c_line[i].m_value[1];
			for (j = 2; j <= sparse->row_weigh[i - 1]; j++) {
				temp = temp*c_line[i].m_value[j];
			}

			for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
			{
				k = c_line[i].link[j]->step;
				c_line[i].link[j]->link[k] = &c_line[i];
				c_line[i].link[j]->m_value[k] = (int)(temp*c_line[i].m_value[j]);                   //type(1)
				c_line[i].link[j]->step++;
			}
		}



		for (i = 1; i <= sparse->col_num; i++)
			v_line[i].step = 1;

		for (i = 1; i <= sparse->col_num; i++)  //Calculate (2)a
		{
			temp = 0;
			for (j = 1; j <= sparse->col_weigh[i - 1]; j++) {
				temp += v_line[i].m_value[j];
			}
			lp[i] = temp;

		}


		for (i = 1; i <= sparse->row_num; i++) {
			c_line[i].step = 1;
		}


		for (i = 1; i <= sparse->col_num; i++)
		{
			temp1 = 0;

			for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			{
				k = v_line[i].link[j]->step;
				v_line[i].link[j]->link[k] = &v_line[i];
				temp1 = v_line[i].link[j]->q_value[k] + lp[i] - v_line[i].m_value[j];
				if (temp1>(1 << (q_bits - 1)) - 1)
					temp1 = (1 << (q_bits - 1)) - 1;

				if (temp1<-(1 << (q_bits - 1)) + 1)
					temp1 = -(1 << (q_bits - 1)) + 1;
				v_line[i].link[j]->q_value[k] = temp1;
				v_line[i].link[j]->m_value[k] = sgn(temp1);

				v_line[i].link[j]->step++;

			}
		}

		for (i = 1; i <= sparse->row_num; i++)
			c_line[i].step = 1;

		for (i = 1; i <= sparse->col_num; i++)  //Calculate
		{
			temp = 0;
			for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			{
				k = v_line[i].link[j]->step;
				temp = temp + v_line[i].link[j]->m_value[k];
				v_line[i].link[j]->step++;
			}
			lp[i] = sgn(lp0[i]) + temp;
		}

		check_result2(sparse, resultcode, lp, decode_valid);

		if (decode_valid)              //Parity Check
		{

			break;
		}
	}
}

void decode_ddwmld(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, int q_bits)
{
	int i, j, k1;
	int temp, temp1;

	for (i = 1; i <= sparse->col_num; i++)			//Initiate
		resultcode[i] = 0;

	/*Parity Check*/
	check_result2(sparse, resultcode, transcode, decode_valid);

	if (decode_valid)              //Parity Check
	{
		return;
	}

	/*Initiate DD Weighted MLD Algorithm*/
	for (i = 0; i <= sparse->col_num; i++)
		lp[i] = 0.0;

	for (i = 0; i <= sparse->col_num; i++)
		lp0[i] = 0.0;

	/*Quantizing y_i*/
	for (i = 1; i <= sparse->col_num; i++)
	{
		lp0[i] = (int)((transcode[i]) / delta);

		if (transcode[i] - lp0[i] * delta > (lp0[i] + 1)*delta - transcode[i])
			lp0[i]++;

		if (lp0[i]>(1 << (q_bits - 1)) - 1)
			lp0[i] = (1 << (q_bits - 1)) - 1;

		if (lp0[i]<-(1 << (q_bits - 1)) + 1)
			lp0[i] = -(1 << (q_bits - 1)) + 1;

	}

	/*Iterative Updating*/
	for (k1 = 0; k1<max_iter; k1++)
	{
		iter++;


		/*Hard information*/
		for (j = 1; j <= sparse->col_num; j++)
		{
			code_sgn[j] = (1 + sgn(lp0[j])) >> 1;
		}

		for (i = 1; i <= sparse->row_num; i++)  //Calculate Check Sums
		{
			temp = 0;
			for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
				temp ^= code_sgn[sparse->row[i - 1][j - 1]];
			chk_sum[i] = temp;		//lp size more than requirement.
		}

		/*Update by MLD method*/
		for (j = 1; j <= sparse->col_num; j++)
		{
			for (i = 1; i <= sparse->col_weigh[j - 1]; i++)
			{
				temp1 = ((chk_sum[sparse->col[j - 1][i - 1]] ^ code_sgn[j]) << 1) - 1;
				lp0[j] += temp1;
			}

			//Protection
			if (lp0[j]>(1 << (q_bits - 1)) - 1)
				lp0[j] = (1 << (q_bits - 1)) - 1;

			if (lp0[j]<-(1 << (q_bits - 1)) + 1)
				lp0[j] = -(1 << (q_bits - 1)) + 1;
		}



		check_result2(sparse, resultcode, lp0, decode_valid);

		if (decode_valid)              //Parity Check
		{
			break;
		}
	}
}

// Weighted Majority Logic Decoding Algorithm
void decode_weighted_ddwmld(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, int q_bits, int chk_d)
{
	int i, j, k1;
	int temp, temp1;
	int chk_step[50] = { 0 };

	for (i = 1; i <= sparse->col_num; i++)			//Initiate
		resultcode[i] = 0;


	//Reliability Mapping
	for (i = 1; i <= chk_d; i++)
	{
		chk_step[i] = i* (1 << (q_bits - 1)) / chk_d;
	}

	/*Parity Check*/
	check_result2(sparse, resultcode, transcode, decode_valid);

	if (decode_valid)              //Parity Check
	{
		return;
	}

	/*Initiate DD Weighted MLD Algorithm*/
	for (i = 0; i <= sparse->col_num; i++)
		lp[i] = 0.0;

	for (i = 0; i <= sparse->col_num; i++)
		lp0[i] = 0.0;

	/*Quantizing y_i*/
	for (i = 1; i <= sparse->col_num; i++)
	{
		lp0[i] = (int)((transcode[i]) / delta);

		if (transcode[i] - lp0[i] * delta > (lp0[i] + 1)*delta - transcode[i])
			lp0[i]++;

		if (lp0[i]>(1 << (q_bits - 1)) - 1)
			lp0[i] = (1 << (q_bits - 1)) - 1;

		if (lp0[i]<-(1 << (q_bits - 1)) + 1)
			lp0[i] = -(1 << (q_bits - 1)) + 1;

		//		printf("%1.1f ",lp0[i]);
	}

	/*Iterative Updating*/
	for (k1 = 0; k1<max_iter; k1++)
	{
		iter++;


		/*Hard information*/
		for (j = 1; j <= sparse->col_num; j++)
		{
			code_sgn[j] = (1 + sgn(lp0[j])) >> 1;
		}

		for (i = 1; i <= sparse->row_num; i++)  //Calculate Check Sums
		{
			temp = 0;
			for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
				temp ^= code_sgn[sparse->row[i - 1][j - 1]];
			chk_sum[i] = temp;		//lp size more than requirement.
		}


		/*Update Reliability measure of Check nodes*/

		//		printf("\n");
		for (i = 1; i <= sparse->row_num; i++)
		{
			lp[i] = fabs(lp0[sparse->row[i - 1][0]]);
			for (j = 2; j <= sparse->row_weigh[i - 1]; j++)
			{
				if (lp[i]>fabs(lp0[sparse->row[i - 1][j - 1]]))
				{
					lp[i] = fabs(lp0[sparse->row[i - 1][j - 1]]);
				}
			}
			//			printf("%1.1f ",lp[i]);
		}

		// 		printf("\n");
		/*Quantize to Range [0,d]*/
		for (i = 1; i <= sparse->row_num; i++)
		{
			for (j = 1; j<chk_d; j++)
			{
				if (lp[i]<chk_step[j])
				{
					lp[i] = j;
					break;
				}
			}
		}

		/*Update by MLD method*/
		for (j = 1; j <= sparse->col_num; j++)
		{
			for (i = 1; i <= sparse->col_weigh[j - 1]; i++)
			{
				temp1 = (int)((((chk_sum[sparse->col[j - 1][i - 1]] ^ code_sgn[j]) << 1) - 1) * lp[sparse->col[j - 1][i - 1]]);
				lp0[j] += temp1;
			}

			//Protection
			if (lp0[j]>(1 << (q_bits - 1)) - 1)
				lp0[j] = (1 << (q_bits - 1)) - 1;

			if (lp0[j]<-(1 << (q_bits - 1)) + 1)
				lp0[j] = -(1 << (q_bits - 1)) + 1;
		}



		check_result2(sparse, resultcode, lp0, decode_valid);

		if (decode_valid)              //Parity Check
		{
			break;
		}
	}
}
// soft relative iterative majority logic decoding
void decode_srbmld(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *weight_1, int *weight_2, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, float scale, double q_step, int q_bits)
{
	int i, j, k1;
	int temp;
	double temp1;
	double min_1, min_2;		//Minimum and the second minimum value
	int min_indx, min_indx_2;
	int d = 3;
	for (i = 1; i <= sparse->col_num; i++)			//Initiate
		resultcode[i] = 0;

	/*Parity Check*/
	check_result2(sparse, resultcode, transcode, decode_valid);

	if (decode_valid)              //Parity Check
	{
		return;
	}

	//printf("%f\n",scale);

	/*Initiate  Weighted MLD Algorithm*/
	for (i = 0; i <= sparse->col_num; i++)
		lp[i] = 0.0;

	for (i = 0; i <= sparse->col_num; i++)
		lp0[i] = 0.0;

	for (i = 1; i <= sparse->col_num; i++)
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			v_line[i].link[j] = &c_line[sparse->col[i - 1][j - 1]];

	for (i = 1; i <= sparse->row_num; i++)
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
			c_line[i].link[j] = &v_line[sparse->row[i - 1][j - 1]];

	for (i = 1; i <= sparse->col_num; i++)		//Clean
	{
		v_line[i].step = 1;
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			v_line[i].r_value[j] = 0;
	}

	for (i = 1; i <= sparse->row_num; i++)		//Clean
	{
		c_line[i].step = 1;
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
			c_line[i].q_value[j] = 0;
	}
	/*lp  y_i,quantization*/

	for (i = 1; i <= sparse->col_num; i++)                 //R(0)=q
	{
		//lp[i]=transcode[i];
		lp[i] = (int)(transcode[i] / q_step);
		if (transcode[i] - lp[i] * q_step>(lp[i] + 1)*q_step - transcode[i])
			lp[i]++;
		if (lp[i]>(1 << (q_bits - 1)) - 1)
			lp[i] = (1 << (q_bits - 1)) - 1;
		if (lp[i]<-(1 << (q_bits - 1)) + 1)
			lp[i] = -(1 << (q_bits - 1)) + 1;

	}


	for (i = 1; i <= sparse->row_num; i++)  //Calculate Check Sums
	{
		min_1 = fabs(lp[sparse->row[i - 1][0]]);
		min_indx = sparse->row[i - 1][0];
		min_indx_2 = 0;
		min_2 = 1e10;
		for (j = 2; j <= sparse->row_weigh[i - 1]; j++)
		{
			if (min_1>fabs(lp[sparse->row[i - 1][j - 1]]))
			{
				min_2 = min_1;
				min_indx_2 = min_indx;
				min_1 = fabs(lp[sparse->row[i - 1][j - 1]]);
				min_indx = sparse->row[i - 1][j - 1];

			}
			else
			{
				if (min_2>fabs(lp[sparse->row[i - 1][j - 1]]))
				{
					min_2 = fabs(lp[sparse->row[i - 1][j - 1]]);
					min_indx_2 = sparse->row[i - 1][j - 1];
				}
			}
		}
		weight_1[i] = (int)fabs(lp[min_indx]);
		weight_2[i] = (int)fabs(lp[min_indx_2]);

	}


	for (i = 1; i <= sparse->col_num; i++)		//q
	{
		temp = (int)fabs(lp[i]);
		/*printf("%f %d\n",lp[i],temp);*/
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
		{
			if (temp>weight_1[sparse->col[i - 1][j - 1]])
				v_line[i].r_value[j] = weight_1[sparse->col[i - 1][j - 1]];
			else
				v_line[i].r_value[j] = weight_2[sparse->col[i - 1][j - 1]];
		}
	}


	/*Iterative Updating*/
	for (k1 = 0; k1<max_iter; k1++)
	{
		//for(i=30;i<50;i++){
		//	printf("%6.1f ",lp[i]);
		//}
		//printf("\n");
		iter++;
		/*Hard information*/
		for (j = 1; j <= sparse->col_num; j++)
		{
			code_sgn[j] = (1 + sgn(lp[j])) >> 1;         //Z(k)
		}

		for (i = 1; i <= sparse->row_num; i++)  //Calculate Check Sums
		{
			temp = 0;
			for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
				temp ^= code_sgn[sparse->row[i - 1][j - 1]];
			chk_sum[i] = temp;		//lp size more than requirement.
		}



		/*Update by MLD method*/
		for (j = 1; j <= sparse->col_num; j++)
		{

			for (i = 1; i <= sparse->col_weigh[j - 1]; i++)
			{

				temp1 = (((chk_sum[sparse->col[j - 1][i - 1]] ^ code_sgn[j]) << 1) - 1);
				lp[j] += temp1;
			}

		}

		//	printf("\n");

		check_result2(sparse, resultcode, lp, decode_valid);

		if (decode_valid)              //Parity Check
		{
			break;
		}

	}
}


void decode_improved_srbmld(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *weight_1, int *weight_2, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, float scale, double q_step, int q_bits)
{
	int i, j, k1;
	int temp;
	double temp1;
	double min_1, min_2;		//Minimum and the second minimum value
	int min_indx, min_indx_2;
	int d = 3;
	for (i = 1; i <= sparse->col_num; i++)			//Initiate
		resultcode[i] = 0;

	/*Parity Check*/
	check_result2(sparse, resultcode, transcode, decode_valid);

	if (decode_valid)              //Parity Check
	{
		return;
	}

	//printf("%f\n",scale);

	/*Initiate  Weighted MLD Algorithm*/
	for (i = 0; i <= sparse->col_num; i++)
		lp[i] = 0.0;

	for (i = 0; i <= sparse->col_num; i++)
		lp0[i] = 0.0;

	for (i = 1; i <= sparse->col_num; i++)
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			v_line[i].link[j] = &c_line[sparse->col[i - 1][j - 1]];

	for (i = 1; i <= sparse->row_num; i++)
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
			c_line[i].link[j] = &v_line[sparse->row[i - 1][j - 1]];

	for (i = 1; i <= sparse->col_num; i++)		//Clean
	{
		v_line[i].step = 1;
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			v_line[i].r_value[j] = 0;
	}

	for (i = 1; i <= sparse->row_num; i++)		//Clean
	{
		c_line[i].step = 1;
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
			c_line[i].q_value[j] = 0;
	}
	/*lp  y_i,quantization*/

	for (i = 1; i <= sparse->col_num; i++)                 //R(0)=q
	{
		//lp[i]=transcode[i];
		lp[i] = (int)(transcode[i] / q_step);
		if (transcode[i] - lp[i] * q_step>(lp[i] + 1)*q_step - transcode[i])
			lp[i]++;
		if (lp[i]>(1 << (q_bits - 1)) - 1)
			lp[i] = (1 << (q_bits - 1)) - 1;
		if (lp[i]<-(1 << (q_bits - 1)) + 1)
			lp[i] = -(1 << (q_bits - 1)) + 1;

	}


	for (i = 1; i <= sparse->row_num; i++)  //Calculate Check Sums
	{
		min_1 = fabs(lp[sparse->row[i - 1][0]]);
		min_indx = sparse->row[i - 1][0];
		min_indx_2 = 0;
		min_2 = 1e10;
		for (j = 2; j <= sparse->row_weigh[i - 1]; j++)
		{
			if (min_1>fabs(lp[sparse->row[i - 1][j - 1]]))
			{
				min_2 = min_1;
				min_indx_2 = min_indx;
				min_1 = fabs(lp[sparse->row[i - 1][j - 1]]);
				min_indx = sparse->row[i - 1][j - 1];

			}
			else
			{
				if (min_2>fabs(lp[sparse->row[i - 1][j - 1]]))
				{
					min_2 = fabs(lp[sparse->row[i - 1][j - 1]]);
					min_indx_2 = sparse->row[i - 1][j - 1];
				}
			}
		}
		weight_1[i] = (int)fabs(lp[min_indx]);
		weight_2[i] = (int)fabs(lp[min_indx_2]);

	}


	for (i = 1; i <= sparse->col_num; i++)		//q
	{
		temp = (int)fabs(lp[i]);
		/*printf("%f %d\n",lp[i],temp);*/
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
		{
			if (temp>weight_1[sparse->col[i - 1][j - 1]])
				v_line[i].r_value[j] = weight_1[sparse->col[i - 1][j - 1]];
			else
				v_line[i].r_value[j] = weight_2[sparse->col[i - 1][j - 1]];
		}
	}


	/*Iterative Updating*/
	for (k1 = 0; k1<max_iter; k1++)
	{
		//for(i=30;i<50;i++){
		//	printf("%6.1f ",lp[i]);
		//}
		//printf("\n");
		iter++;
		/*Hard information*/
		for (j = 1; j <= sparse->col_num; j++)
		{
			code_sgn[j] = (1 + sgn(lp[j])) >> 1;         //Z(k)
		}

		for (i = 1; i <= sparse->row_num; i++)  //Calculate Check Sums
		{
			temp = 0;
			for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
				temp ^= code_sgn[sparse->row[i - 1][j - 1]];
			chk_sum[i] = temp;		//lp size more than requirement.
		}



		/*Update by MLD method*/
		for (j = 1; j <= sparse->col_num; j++)
		{

			for (i = 1; i <= sparse->col_weigh[j - 1]; i++)
			{

				temp1 = (((chk_sum[sparse->col[j - 1][i - 1]] ^ code_sgn[j]) << 1) - 1)*delta*v_line[j].r_value[i];
				lp[j] += temp1;
			}

		}

		//	printf("\n");

		check_result2(sparse, resultcode, lp, decode_valid);

		if (decode_valid)              //Parity Check
		{
			break;
		}

	}
}


void decode_mrbmld(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *weight_1, int *weight_2, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, double q_step, int q_bits)
{
	int i, j, k1;
	int temp;
	double temp1;
	int d = 3;
	for (i = 1; i <= sparse->col_num; i++)			//Initiate
		resultcode[i] = 0;

	/*Parity Check*/
	check_result2(sparse, resultcode, transcode, decode_valid);

	if (decode_valid)              //Parity Check
	{
		return;
	}

	//printf("%f\n",scale);

	/*Initiate  Weighted MLD Algorithm*/
	for (i = 0; i <= sparse->col_num; i++)
		lp[i] = 0.0;

	for (i = 0; i <= sparse->col_num; i++)
		lp0[i] = 0.0;

	for (i = 1; i <= sparse->col_num; i++)
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			v_line[i].link[j] = &c_line[sparse->col[i - 1][j - 1]];

	for (i = 1; i <= sparse->row_num; i++)
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
			c_line[i].link[j] = &v_line[sparse->row[i - 1][j - 1]];

	for (i = 1; i <= sparse->col_num; i++)		//Clean
	{
		v_line[i].step = 1;
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			v_line[i].r_value[j] = 0;
	}

	for (i = 1; i <= sparse->row_num; i++)		//Clean
	{
		c_line[i].step = 1;
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
			c_line[i].q_value[j] = 0;
	}
	/*lp  y_i,quantization*/

	for (i = 1; i <= sparse->col_num; i++)                 //R(0)=q
	{
		//lp[i]=transcode[i];
		lp0[i] = (int)(transcode[i] / q_step);
		if (transcode[i] - lp0[i] * q_step>(lp[i] + 1)*q_step - transcode[i])
			lp0[i]++;
		if (lp0[i]>(1 << (q_bits - 1)) - 1)
			lp0[i] = (1 << (q_bits - 1)) - 1;
		if (lp0[i]<-(1 << (q_bits - 1)) + 1)
			lp0[i] = -(1 << (q_bits - 1)) + 1;
		/* if(lp[i]>(1<<(q_bits-2))-1){
		if(lp[i]<d*((1<<(q_bits-2))-1))
		lp[i]=(1<<(q_bits-2))-1;
		else if(lp[i]>=pow(d,(1<<(q_bits-2)))*((1<<(q_bits-2))-1))
		lp[i]=pow(d,(1<<(q_bits-2)))*((1<<(q_bits-2))-1);
		else{
		for(int r=1;r++;r<=(1<<(q_bits-2))-1){
		if(lp0[i]>=pow(d,r)*((1<<(q_bits-2))-1)&&lp0[i]<pow(d,r+1)*((1<<(q_bits-2))-1)){
		lp0[i]=pow(d,r)*((1<<(q_bits-2))-1);
		break;
		}
		}
		}
		}

		if(lp[i]<-(1<<(q_bits-2))+1){
		if(lp[i]>-d*((1<<(q_bits-2))-1))
		lp[i]=-(1<<(q_bits-2))+1;
		else if(lp[i]<=-pow(d,(1<<(q_bits-2)))*((1<<(q_bits-2))-1))
		lp[i]=-pow(d,(1<<(q_bits-2)))*((1<<(q_bits-2))-1);
		else{
		for(int r=1;r++;r<=(1<<(q_bits-2))-1){
		if(lp[i]>-pow(d,r+1)*((1<<(q_bits-2))-1)&&lp0[i]<=-pow(d,r)*((1<<(q_bits-2))-1)){
		lp[i]=-pow(d,r)*((1<<(q_bits-2))-1);
		break;
		}
		}
		}
		}	*/
	}

	/*Iterative Updating*/
	for (k1 = 0; k1<max_iter; k1++)
	{
		iter++;
		/*Hard information*/
		for (j = 1; j <= sparse->col_num; j++)
		{
			code_sgn[j] = (1 + sgn(lp[j])) >> 1;         //Z(k)
		}

		for (i = 1; i <= sparse->row_num; i++)  //Calculate Check Sums
		{
			temp = 0;
			for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
				temp ^= code_sgn[sparse->row[i - 1][j - 1]];
			chk_sum[i] = temp;		//lp size more than requirement.
		}

		for (i = 0; i <= sparse->col_num; i++)
			lp[i] = 0.0;
		/*Update by MLD method*/
		for (j = 1; j <= sparse->col_num; j++)
		{

			for (i = 1; i <= sparse->col_weigh[j - 1]; i++)
			{

				temp1 = (((chk_sum[sparse->col[j - 1][i - 1]] ^ code_sgn[j]) << 1) - 1)*delta;
				lp[j] += temp1;
			}
			lp[j] = lp0[j] + lp[j];
			/*   if(lp[j]>(1<<(q_bits-1))-1)
			lp[j]=(1<<(q_bits-1))-1;
			if(lp[j]<-(1<<(q_bits-1))+1)
			lp[j]=-(1<<(q_bits-1))+1;*/
		}
		//	printf("\n");
		check_result2(sparse, resultcode, lp, decode_valid);

		if (decode_valid)              //Parity Check
		{
			break;
		}

	}
}

void decode_rbi_mlgd(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *weight_1, int *weight_2, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, double q_step, int q_bits)
{
	int i, j, k1;
	int temp;
	double temp1;
	double d = 1.3;
	double min_1, min_2;		//Minimum and the second minimum value
	int min_indx, min_indx_2;

	for (i = 1; i <= sparse->col_num; i++)			//Initiate
		resultcode[i] = 0;

	/*Parity Check*/
	check_result2(sparse, resultcode, transcode, decode_valid);

	if (decode_valid)              //Parity Check
	{
		return;
	}

	//printf("%f\n",scale);

	/*Initiate  Weighted MLD Algorithm*/
	for (i = 0; i <= sparse->col_num; i++)
		lp[i] = 0.0;

	for (i = 0; i <= sparse->col_num; i++)
		lp0[i] = 0.0;

	for (i = 1; i <= sparse->col_num; i++)
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			v_line[i].link[j] = &c_line[sparse->col[i - 1][j - 1]];

	for (i = 1; i <= sparse->row_num; i++)
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
			c_line[i].link[j] = &v_line[sparse->row[i - 1][j - 1]];

	for (i = 1; i <= sparse->col_num; i++)		//Clean
	{
		v_line[i].step = 1;
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			v_line[i].r_value[j] = 0;
	}

	for (i = 1; i <= sparse->row_num; i++)		//Clean
	{
		c_line[i].step = 1;
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
			c_line[i].q_value[j] = 0;
	}
	/*lp  y_i,quantization*/

	for (i = 1; i <= sparse->col_num; i++)                 //R(0)=q
	{
		//lp[i]=transcode[i];
		lp0[i] = (int)(transcode[i] / q_step);
		if (transcode[i] - lp0[i] * q_step>(lp[i] + 1)*q_step - transcode[i])
			lp0[i]++;
		if (lp0[i]>(1 << (q_bits - 1)) - 1)
			lp0[i] = (1 << (q_bits - 1)) - 1;
		if (lp0[i]<-(1 << (q_bits - 1)) + 1)
			lp0[i] = -(1 << (q_bits - 1)) + 1;
		/* if(lp0[i]>(1<<(q_bits-2))-1){
		if(lp0[i]<d*((1<<(q_bits-2))-1))
		lp0[i]=(1<<(q_bits-2))-1;
		else if(lp0[i]>=powf(d,(1<<(q_bits-2)))*(1<<(q_bits-2)))
		lp0[i]=(int)(powf(d,(1<<(q_bits-2)))*(1<<(q_bits-2)));
		else{
		for(int r=1;r++;r<=(1<<(q_bits-2))-1){
		if(lp0[i]>=powf(d,r)*((1<<(q_bits-2))-1)&&lp0[i]<powf(d,r+1)*((1<<(q_bits-2))-1)){
		lp0[i]=(int)(powf(d,r)*((1<<(q_bits-2))-1));
		break;
		}
		}
		}
		}

		if(lp0[i]<-(1<<(q_bits-2))+1){
		if(lp0[i]>-d*((1<<(q_bits-2))-1))
		lp0[i]=-(1<<(q_bits-2))+1;
		else if(lp0[i]<=-powf(d,(1<<(q_bits-2)))*(1<<(q_bits-2)))
		lp0[i]=-(int)(powf(d,(1<<(q_bits-2)))*(1<<(q_bits-2)));
		else{
		for(int r=1;r++;r<=(1<<(q_bits-2))-1){
		if(lp0[i]>-powf(d,r+1)*((1<<(q_bits-2))-1)&&lp0[i]<=-powf(d,r)*((1<<(q_bits-2))-1)){
		lp0[i]=-(int)(powf(d,r)*((1<<(q_bits-2))-1));
		break;
		}
		}

		}

		}*/



		lp[i] = lp0[i];
	}



	/*Iterative Updating*/
	for (k1 = 0; k1<max_iter; k1++)
	{
		iter++;
		/*Hard information*/
		for (j = 1; j <= sparse->col_num; j++)
		{
			code_sgn[j] = (1 + sgn(lp[j])) >> 1;         //Z(k)
		}

		for (i = 1; i <= sparse->row_num; i++)  //Calculate Check Sums
		{
			temp = 0;
			for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
				temp ^= code_sgn[sparse->row[i - 1][j - 1]];
			chk_sum[i] = temp;		//lp size more than requirement.
		}

		for (i = 1; i <= sparse->row_num; i++)  //Calculate Check Sums
		{
			min_1 = fabs(lp[sparse->row[i - 1][0]]);
			min_indx = sparse->row[i - 1][0];
			min_indx_2 = 0;
			min_2 = 1e10;
			for (j = 2; j <= sparse->row_weigh[i - 1]; j++)
			{
				if (min_1>fabs(lp[sparse->row[i - 1][j - 1]]))
				{
					min_2 = min_1;
					min_indx_2 = min_indx;
					min_1 = fabs(lp[sparse->row[i - 1][j - 1]]);
					min_indx = sparse->row[i - 1][j - 1];

				}
				else
				{
					if (min_2>fabs(lp[sparse->row[i - 1][j - 1]]))
					{
						min_2 = fabs(lp[sparse->row[i - 1][j - 1]]);
						min_indx_2 = sparse->row[i - 1][j - 1];
					}
				}
			}
			weight_1[i] = (int)fabs(lp[min_indx]);
			weight_2[i] = (int)fabs(lp[min_indx_2]);

		}

		for (i = 1; i <= sparse->col_num; i++)		//q
		{
			temp = (int)fabs(lp[i]);
			/*printf("%f %d\n",lp[i],temp);*/
			for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			{
				if (temp>weight_1[sparse->col[i - 1][j - 1]])
					v_line[i].r_value[j] = weight_1[sparse->col[i - 1][j - 1]];
				else
					v_line[i].r_value[j] = weight_2[sparse->col[i - 1][j - 1]];
			}
		}


		for (i = 0; i <= sparse->col_num; i++)
			lp[i] = 0.0;
		/*Update by MLD method*/
		for (j = 1; j <= sparse->col_num; j++)
		{

			for (i = 1; i <= sparse->col_weigh[j - 1]; i++)
			{

				temp1 = (((chk_sum[sparse->col[j - 1][i - 1]] ^ code_sgn[j]) << 1) - 1)*delta*v_line[j].r_value[i];
				lp[j] += temp1;
			}
			lp[j] = lp0[j] + lp[j];

		}
		//	printf("\n");
		check_result2(sparse, resultcode, lp, decode_valid);

		if (decode_valid)              //Parity Check
		{
			break;
		}

	}
}



void decode_quai_rbi_mlgd(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *weight_1, int *weight_2, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, double q_step, int q_bits)
{
	int i, j, k1;
	int temp;
	double temp1;
	double d = 40.0;
	double k = d*q_step + 0.5*q_step;
	double k2 = k;
	double kmax = 4.025;
	double q_step_2 = 0.01;
	double min_1, min_2;		//Minimum and the second minimum value
	int min_indx, min_indx_2;

	for (i = 1; i <= sparse->col_num; i++)			//Initiate
		resultcode[i] = 0;

	/*Parity Check*/
	check_result2(sparse, resultcode, transcode, decode_valid);

	if (decode_valid)              //Parity Check
	{
		return;
	}

	//printf("%f\n",scale);

	/*Initiate  Weighted MLD Algorithm*/
	for (i = 0; i <= sparse->col_num; i++)
		lp[i] = 0.0;

	for (i = 0; i <= sparse->col_num; i++)
		lp0[i] = 0.0;

	for (i = 1; i <= sparse->col_num; i++)
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			v_line[i].link[j] = &c_line[sparse->col[i - 1][j - 1]];

	for (i = 1; i <= sparse->row_num; i++)
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
			c_line[i].link[j] = &v_line[sparse->row[i - 1][j - 1]];

	for (i = 1; i <= sparse->col_num; i++)		//Clean
	{
		v_line[i].step = 1;
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			v_line[i].r_value[j] = 0;
	}

	for (i = 1; i <= sparse->row_num; i++)		//Clean
	{
		c_line[i].step = 1;
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
			c_line[i].q_value[j] = 0;
	}
	/*lp  y_i,quantization*/
	for (i = 1; i <= sparse->col_num; i++)                 //R(0)=q
	{
		//lp[i]=transcode[i];
		lp0[i] = (int)(transcode[i] / q_step);
		if (transcode[i] - lp0[i] * q_step>(lp[i] + 1)*q_step - transcode[i])
			lp0[i]++;
		if (lp0[i] >= d) {
			for (j = 1; j <= 19; j++) {
				if (k2 <= transcode[i] && transcode[i]<k2 + j*q_step_2) {
					lp0[i] = d + j;
					break;
				}
				k2 += j*q_step_2;
			}
		}

		if (lp0[i] <= -d) {
			for (j = 1; j <= 19; j++) {
				if (-(k2 + j*q_step_2)<transcode[i] && transcode[i] <= -k2) {
					lp0[i] = -d - j;
					break;
				}
				k2 += j*q_step_2;
			}
		}

		if (transcode[i] >= kmax)
			lp0[i] = 60.0;
		if (transcode[i] <= -kmax)
			lp0[i] = -60.0;

		//      if(lp0[i]>(1<<(q_bits-1))-1)
		//	lp0[i]=(1<<(q_bits-1))-1;
		//if(lp0[i]<-(1<<(q_bits-1))+1)
		//	lp0[i]=-(1<<(q_bits-1))+1;
		//if(lp0[i]>d){
		//	if(lp0[i]<d1)
		//		lp0[i]=61.0;
		//	else if(lp0[i]>d2)
		//		lp0[i]=63.0;
		//	else
		//		lp0[i]=62.0;
		//}	
		//if(lp0[i]<-d){
		//	if(lp0[i]<d1)
		//		lp0[i]=-61.0;
		//	else if(lp0[i]<d2)
		//		lp0[i]=-63.0;
		//	else
		//		lp0[i]=-62.0;
		//}
		//if(lp0[i]>d)
		//           lp0[i]=d;
		//      if(lp0[i]<-d)
		//	 lp0[i]=-d;

		//   if(lp0[i]>(1<<(q_bits-2))-1){
		//    if(lp0[i]<d*((1<<(q_bits-2))-1))
		//        lp0[i]=(1<<(q_bits-2))-1;
		//    else if(lp0[i]>=pow(d,(1<<(q_bits-2)))*(1<<(q_bits-2)))
		//         lp0[i]=pow(d,(1<<(q_bits-2)))*(1<<(q_bits-2));
		//    else{
		//          for(int r=1;r++;r<=(1<<(q_bits-2))-1){
		//               if(lp0[i]>=pow(d,r)*((1<<(q_bits-2))-1)&&lp0[i]<pow(d,r+1)*((1<<(q_bits-2))-1)){
		//                    lp0[i]=pow(d,r)*((1<<(q_bits-2))-1);

		//                     break;
		//               }
		//          }  
		//      }
		//}
		//if(lp0[i]<-(1<<(q_bits-2))+1){
		//     if(lp0[i]>-d*((1<<(q_bits-2))-1))
		//         lp0[i]=-(1<<(q_bits-2))+1;
		//     else if(lp0[i]<=-pow(d,(1<<(q_bits-2)))*(1<<(q_bits-2)))
		//         lp0[i]=-pow(d,(1<<(q_bits-2)))*(1<<(q_bits-2));
		//     else{
		//         for(int r=1;r++;r<=(1<<(q_bits-2))-1){
		//             if(lp0[i]>-pow(d,r+1)*((1<<(q_bits-2))-1)&&lp0[i]<=-pow(d,r)*((1<<(q_bits-2))-1)){
		//                  lp0[i]=-pow(d,r)*((1<<(q_bits-2))-1);
		//                  break;
		//              }
		//          }
		//          
		//	 }
		//      
		//}		
		lp[i] = lp0[i];
	}

	/*Iterative Updating*/
	for (k1 = 0; k1<max_iter; k1++)
	{
		iter++;
		/*Hard information*/
		for (j = 1; j <= sparse->col_num; j++)
		{
			code_sgn[j] = (1 + sgn(lp[j])) >> 1;         //Z(k)
		}

		for (i = 1; i <= sparse->row_num; i++)  //Calculate Check Sums
		{
			temp = 0;
			for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
				temp ^= code_sgn[sparse->row[i - 1][j - 1]];
			chk_sum[i] = temp;		//lp size more than requirement.
		}

		for (i = 1; i <= sparse->row_num; i++)  //Calculate Check Sums
		{
			min_1 = fabs(lp[sparse->row[i - 1][0]]);
			min_indx = sparse->row[i - 1][0];
			min_indx_2 = 0;
			min_2 = 1e10;
			for (j = 2; j <= sparse->row_weigh[i - 1]; j++)
			{
				if (min_1>fabs(lp[sparse->row[i - 1][j - 1]]))
				{
					min_2 = min_1;
					min_indx_2 = min_indx;
					min_1 = fabs(lp[sparse->row[i - 1][j - 1]]);
					min_indx = sparse->row[i - 1][j - 1];

				}
				else
				{
					if (min_2>fabs(lp[sparse->row[i - 1][j - 1]]))
					{
						min_2 = fabs(lp[sparse->row[i - 1][j - 1]]);
						min_indx_2 = sparse->row[i - 1][j - 1];
					}
				}
			}
			weight_1[i] = (int)fabs(lp[min_indx]);
			weight_2[i] = (int)fabs(lp[min_indx_2]);

		}

		for (i = 1; i <= sparse->col_num; i++)		//q
		{
			temp = (int)fabs(lp[i]);
			for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			{
				if (temp>weight_1[sparse->col[i - 1][j - 1]])
					v_line[i].r_value[j] = weight_1[sparse->col[i - 1][j - 1]];
				else
					v_line[i].r_value[j] = weight_2[sparse->col[i - 1][j - 1]];
			}
		}


		for (i = 0; i <= sparse->col_num; i++)
			lp[i] = 0.0;
		/*Update by MLD method*/
		double tmp = 0;
		for (j = 1; j <= sparse->col_num; j++)
		{

			for (i = 1; i <= sparse->col_weigh[j - 1]; i++)
			{

				temp1 = (((chk_sum[sparse->col[j - 1][i - 1]] ^ code_sgn[j]) << 1) - 1)*delta*v_line[j].r_value[i];
				lp[j] += temp1;
			}
			lp[j] = lp0[j] + lp[j];

		}

		check_result2(sparse, resultcode, lp, decode_valid);

		if (decode_valid)              //Parity Check
		{
			break;
		}
	}
}


/*Iterative One Step Majority Logic Decoding*/
void decode_iosmld(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, int q_bits)
{
	int i, j, k1;
	int temp, temp1;

	for (i = 1; i <= sparse->col_num; i++)			//Initiate
		resultcode[i] = 0;

	/*Parity Check*/
	check_result2(sparse, resultcode, transcode, decode_valid);

	if (decode_valid)              //Parity Check
	{
		return;
	}

	/*Initiate DD Weighted MLD Algorithm*/
	for (i = 0; i <= sparse->col_num; i++)
		lp[i] = 0.0;

	/*lp  y_i*/
	for (i = 1; i <= sparse->col_num; i++)
	{
		lp[i] = sgn(transcode[i])*sparse->max_colweigh;
	}

	/*Iterative Updating*/
	for (k1 = 0; k1<max_iter; k1++)
	{
		iter++;

		/*Hard information*/
		for (j = 1; j <= sparse->col_num; j++)
		{
			code_sgn[j] = (1 + sgn(lp[j])) >> 1;
		}

		for (i = 1; i <= sparse->row_num; i++)  //Calculate Check Sums
		{
			temp = 0;
			for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
				temp ^= code_sgn[sparse->row[i - 1][j - 1]];
			chk_sum[i] = temp;		//lp size more than requirement.
		}

		// 		for(i=0;i<=sparse->col_num;i++)
		// 			lp[i]=0.0;

		/*Update by MLD method*/
		for (j = 1; j <= sparse->col_num; j++)
		{
			for (i = 1; i <= sparse->col_weigh[j - 1]; i++)
			{
				temp1 = ((chk_sum[sparse->col[j - 1][i - 1]] ^ code_sgn[j]) << 1) - 1;
				lp[j] += temp1;
			}

			if (lp[j]>(sparse->max_colweigh))
			{
				lp[j] = (sparse->max_colweigh);
			}

			if (lp[j]<-(sparse->max_colweigh))
			{
				lp[j] = -(sparse->max_colweigh);
			}

			if (lp[j] == 0)
			{
				lp[j] = (code_sgn[j] << 1) - 1;
			}
		}

		check_result2(sparse, resultcode, lp, decode_valid);

		if (decode_valid)              //Parity Check
		{
			break;
		}

	}
}

/*One Step Majority Logic Decoding*/
void decode_osmld(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, int q_bits)
{
	int i, j, k1;
	int temp, temp1;
	int cnt = 0;

	for (i = 1; i <= sparse->col_num; i++)			//Initiate
		resultcode[i] = 0;

	/*Parity Check*/
	check_result2(sparse, resultcode, transcode, decode_valid);

	if (decode_valid)              //Parity Check
	{
		return;
	}

	/*Initiate DD Weighted MLD Algorithm*/
	for (i = 0; i <= sparse->col_num; i++)
		lp[i] = 0.0;

	/*lp  y_i*/
	for (i = 1; i <= sparse->col_num; i++)
	{
		lp[i] = sgn(transcode[i]);
	}

	/*Iterative Updating*/
	for (k1 = 0; k1<1; k1++)
	{
		iter++;
		cnt = 0;
		/*Hard information*/
		for (j = 1; j <= sparse->col_num; j++)
		{
			code_sgn[j] = (1 + sgn(lp[j])) >> 1;

			if (code_sgn[j] == 1)
			{
				cnt++;
			}
		}
		//		printf("\n");

		for (i = 1; i <= sparse->row_num; i++)  //Calculate Check Sums
		{
			temp = 0;
			for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
				temp ^= code_sgn[sparse->row[i - 1][j - 1]];
			chk_sum[i] = temp;		//lp size more than requirement.
		}

		for (i = 0; i <= sparse->col_num; i++)
			lp[i] = 0.0;

		/*Update by MLD method*/
		for (j = 1; j <= sparse->col_num; j++)
		{
			for (i = 1; i <= sparse->col_weigh[j - 1]; i++)
			{
				temp1 = (chk_sum[sparse->col[j - 1][i - 1]] << 1) - 1;
				lp[j] += temp1;
			}

			//		printf("%f , %d ;", lp[j], j);


			if (lp[j] <= 0)
			{
				lp[j] = (code_sgn[j] << 1) - 1;
			}
			else
			{
				lp[j] = -(code_sgn[j] << 1) + 1;
			}
		}
		//	printf("\n");

		check_result2(sparse, resultcode, lp, decode_valid);

		if (decode_valid)              //Parity Check
		{
			break;
		}
		else
		{
			if (cnt <= (sparse->col_weigh[1]) / 2)
				printf("*************\n");
		}

	}
}

/*BF*/
void decode_bf(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, int q_bits)
{
	int i, j, k1;
	int temp, temp1;
	double lrgst_sum;
	int lrgst_indx;

	for (i = 1; i <= sparse->col_num; i++)			//Initiate
		resultcode[i] = 0;

	/*Parity Check*/
	check_result2(sparse, resultcode, transcode, decode_valid);

	if (decode_valid)              //Parity Check
	{
		return;
	}

	/*Initiate DD Weighted MLD Algorithm*/
	for (i = 0; i <= sparse->col_num; i++)
		lp[i] = 0.0;

	/*lp  y_i*/
	for (i = 1; i <= sparse->col_num; i++)
	{
		lp[i] = sgn(transcode[i]);
	}


	/*Iterative Updating*/
	for (k1 = 0; k1<max_iter; k1++)
	{
		iter++;


		/*Hard information*/
		for (j = 1; j <= sparse->col_num; j++)
		{
			code_sgn[j] = (1 + sgn(lp[j])) >> 1;
			//		printf("%d %d ",code_sgn[j],j);
		}
		//		printf("\n");

		for (i = 1; i <= sparse->row_num; i++)  //Calculate Check Sums
		{
			temp = 0;
			for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
				temp ^= code_sgn[sparse->row[i - 1][j - 1]];
			chk_sum[i] = temp;		//lp size more than requirement.
		}

		for (i = 0; i <= sparse->col_num; i++)
			lp0[i] = 0.0;

		lrgst_indx = 0;
		lrgst_sum = -1;

		/*Update by MLD method*/
		for (j = 1; j <= sparse->col_num; j++)
		{
			for (i = 1; i <= sparse->col_weigh[j - 1]; i++)
			{
				temp1 = chk_sum[sparse->col[j - 1][i - 1]];
				lp0[j] += temp1;
			}

			if (lrgst_sum<lp0[j])
			{
				lrgst_indx = j;
				lrgst_sum = lp0[j];
			}
			//			printf("%d ",lrgst_sum);
		}
		//	printf("\n");

		/*Flip Multiple Bits*/
		for (j = 1; j <= sparse->col_num; j++)
		{
			if (lrgst_sum == lp0[j])
			{
				lp[j] = -lp[j];
			}
		}




		check_result2(sparse, resultcode, lp, decode_valid);

		if (decode_valid)              //Parity Check
		{
			break;
		}

	}
}


/*WBF*/
void decode_wbf(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, int q_bits)
{
	int i, j, k1;
	int temp;
	double temp1;
	double lrgst_sum;
	int	lrgst_indx;
	double min_y;			//min{|y_i|}

	double tmp_E;			//temp E

	for (i = 1; i <= sparse->col_num; i++)			//Initiate
		resultcode[i] = 0;

	/*Parity Check*/
	check_result2(sparse, resultcode, transcode, decode_valid);

	if (decode_valid)              //Parity Check
	{
		return;
	}

	/*Initiate WBF Algorithm*/
	for (i = 0; i <= sparse->col_num; i++)
		lp[i] = 0.0;

	/*lp  y_i*/
	for (i = 1; i <= sparse->col_num; i++)
	{
		lp[i] = sgn(transcode[i]);
	}


	/*Iterative Updating*/
	for (k1 = 0; k1<max_iter; k1++)
	{
		iter++;


		/*Hard information*/
		for (j = 1; j <= sparse->col_num; j++)
		{
			code_sgn[j] = (1 + sgn(lp[j])) >> 1;
			//		printf("%d %d ",code_sgn[j],j);
		}
		//		printf("\n");

		for (i = 1; i <= sparse->row_num; i++)  //Calculate Check Sums
		{
			temp = 0;
			for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
				temp ^= code_sgn[sparse->row[i - 1][j - 1]];
			chk_sum[i] = temp;		//lp size more than requirement.

									//Calculate the minimum y_i
			min_y = fabs(transcode[sparse->row[i - 1][0]]);

			for (j = 2; j <= sparse->row_weigh[i - 1]; j++)
			{
				if (min_y>fabs(transcode[sparse->row[i - 1][j - 1]]))
					min_y = fabs(transcode[sparse->row[i - 1][j - 1]]);
			}

			lp0[i] = min_y;
		}

		lrgst_indx = 0;
		lrgst_sum = -1000;
		/*Update by MLD method*/
		for (j = 1; j <= sparse->col_num; j++)
		{
			tmp_E = 0;
			for (i = 1; i <= sparse->col_weigh[j - 1]; i++)
			{
				temp1 = (chk_sum[sparse->col[j - 1][i - 1]] << 1) - 1;
				tmp_E += temp1*lp0[sparse->col[j - 1][i - 1]];
			}

			if (lrgst_sum<tmp_E)
			{
				lrgst_indx = j;
				lrgst_sum = tmp_E;
			}
			//			printf("%d ",lrgst_sum);
		}
		//	printf("\n");

		/*Flip One Bit*/
		lp[lrgst_indx] = -lp[lrgst_indx];



		check_result2(sparse, resultcode, lp, decode_valid);

		if (decode_valid)              //Parity Check
		{
			break;
		}

	}
}

/*MWBF*/
void decode_mwbf(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, int q_bits, double alpha)
{
	int i, j, k1;
	int temp;
	double temp1;
	double lrgst_sum;
	int	lrgst_indx;
	double min_y;			//min{|y_i|}

	double tmp_E;			//temp E

	for (i = 1; i <= sparse->col_num; i++)			//Initiate
		resultcode[i] = 0;

	/*Parity Check*/
	check_result2(sparse, resultcode, transcode, decode_valid);

	if (decode_valid)              //Parity Check
	{
		return;
	}

	/*Initiate WBF Algorithm*/
	for (i = 0; i <= sparse->col_num; i++)
		lp[i] = 0.0;

	/*lp  y_i*/
	for (i = 1; i <= sparse->col_num; i++)
	{
		lp[i] = sgn(transcode[i]);
	}


	/*Iterative Updating*/
	for (k1 = 0; k1<max_iter; k1++)
	{
		iter++;


		/*Hard information*/
		for (j = 1; j <= sparse->col_num; j++)
		{
			code_sgn[j] = (1 + sgn(lp[j])) >> 1;
			//		printf("%d %d ",code_sgn[j],j);
		}
		//		printf("\n");

		for (i = 1; i <= sparse->row_num; i++)  //Calculate Check Sums
		{
			temp = 0;
			for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
				temp ^= code_sgn[sparse->row[i - 1][j - 1]];
			chk_sum[i] = temp;		//lp size more than requirement.

									//Calculate the minimum y_i
			min_y = fabs(transcode[sparse->row[i - 1][0]]);

			for (j = 2; j <= sparse->row_weigh[i - 1]; j++)
			{
				if (min_y>fabs(transcode[sparse->row[i - 1][j - 1]]))
					min_y = fabs(transcode[sparse->row[i - 1][j - 1]]);
			}

			lp0[i] = min_y;
		}

		lrgst_indx = 0;
		lrgst_sum = -1000;
		/*Update by MLD method*/
		for (j = 1; j <= sparse->col_num; j++)
		{
			tmp_E = 0;
			for (i = 1; i <= sparse->col_weigh[j - 1]; i++)
			{
				temp1 = (chk_sum[sparse->col[j - 1][i - 1]] << 1) - 1;
				tmp_E += temp1*lp0[sparse->col[j - 1][i - 1]];
			}

			tmp_E = tmp_E - alpha*fabs(transcode[j]);

			if (lrgst_sum<tmp_E)
			{
				lrgst_indx = j;
				lrgst_sum = tmp_E;
			}
			//			printf("%d ",lrgst_sum);
		}
		//	printf("\n");

		/*Flip One Bit*/
		lp[lrgst_indx] = -lp[lrgst_indx];



		check_result2(sparse, resultcode, lp, decode_valid);

		if (decode_valid)              //Parity Check
		{
			break;
		}

	}
}


/*IMWBF*/
void decode_imwbf(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int *chk_sum, int *code_sgn, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, float delta, int q_bits, double alpha, double *lp1, int *min_indx)
{
	int i, j, k1;
	int temp;
	double temp1;
	double lrgst_sum;
	int	lrgst_indx;
	double min_y;			//min{|y_i|}
	double submin_y;

	double tmp_E;			//temp E

	for (i = 1; i <= sparse->col_num; i++)			//Initiate
		resultcode[i] = 0;

	/*Parity Check*/
	check_result2(sparse, resultcode, transcode, decode_valid);

	if (decode_valid)              //Parity Check
	{
		return;
	}

	/*Initiate WBF Algorithm*/
	for (i = 0; i <= sparse->col_num; i++)
		lp[i] = 0.0;

	/*lp  y_i*/
	for (i = 1; i <= sparse->col_num; i++)
	{
		lp[i] = sgn(transcode[i]);
	}


	/*Iterative Updating*/
	for (k1 = 0; k1<max_iter; k1++)
	{
		iter++;


		/*Hard information*/
		for (j = 1; j <= sparse->col_num; j++)
		{
			code_sgn[j] = (1 + sgn(lp[j])) >> 1;
			//printf("%d %d ",code_sgn[j],j);
		}
		//printf("\n");

		for (i = 1; i <= sparse->row_num; i++)  //Calculate Check Sums
		{
			temp = 0;
			for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
				temp ^= code_sgn[sparse->row[i - 1][j - 1]];
			chk_sum[i] = temp;		//lp size more than requirement.

									//Calculate the minimum y_i
			min_y = fabs(transcode[sparse->row[i - 1][0]]);
			submin_y = 1000000;

			min_indx[i] = sparse->row[i - 1][0];

			for (j = 2; j <= sparse->row_weigh[i - 1]; j++)
			{
				if (min_y>fabs(transcode[sparse->row[i - 1][j - 1]]))
				{
					submin_y = min_y;
					min_y = fabs(transcode[sparse->row[i - 1][j - 1]]);
					min_indx[i] = sparse->row[i - 1][j - 1];
				}
			}
			if (submin_y == 0)
			{
				for (j = 2; j <= sparse->row_weigh[i - 1]; j++)
				{
					if (submin_y>fabs(transcode[sparse->row[i - 1][j - 1]]))
					{
						submin_y = fabs(transcode[sparse->row[i - 1][j - 1]]);
					}
				}
			}

			lp0[i] = min_y;
			lp1[i] = submin_y;
		}

		lrgst_indx = 0;
		lrgst_sum = -1000000;
		/*Update by MLD method*/
		for (j = 1; j <= sparse->col_num; j++)
		{
			tmp_E = 0;
			for (i = 1; i <= sparse->col_weigh[j - 1]; i++)
			{
				temp1 = (chk_sum[sparse->col[j - 1][i - 1]] << 1) - 1;

				if (j != min_indx[sparse->col[j - 1][i - 1]])
					tmp_E += temp1*lp0[sparse->col[j - 1][i - 1]];
				else
					tmp_E += temp1*lp1[sparse->col[j - 1][i - 1]];
			}

			tmp_E = tmp_E - alpha*fabs(transcode[j]);

			if (lrgst_sum<tmp_E)
			{
				lrgst_indx = j;
				lrgst_sum = tmp_E;
			}
			//			printf("%d ",lrgst_sum);
		}
		//	printf("\n");

		/*Flip One Bit*/
		lp[lrgst_indx] = -lp[lrgst_indx];



		check_result2(sparse, resultcode, lp, decode_valid);

		if (decode_valid)              //Parity Check
		{
			break;
		}

	}

}



//Show the subgraph if a< a_thrsd,  b< b_thrsd,  thrsd (threshold)
void check_graph(struct matrix *sparse, vector<int>& ts_v, vector<int>& ts_c)
{
	int i, j;
	int a, b;

	a = ts_v.size();
	b = ts_c.size();


	printf("Subgraph(%d,%d) \n", a, b);
	for (j = 0; j<a; j++)
	{
		printf("%d :  ", ts_v[j]);
		for (i = 0; i<sparse->col_weigh[ts_v[j]]; i++)
		{
			printf("%d ", sparse->col[ts_v[j]][i] - 1);
		}
		printf("\n");
	}

	printf("Odd degree check-sums\n");
	for (i = 0; i<b; i++)
	{
		printf("%d: ", ts_c[i]);
		for (j = 0; j<sparse->row_weigh[ts_c[i]]; j++)
		{
			printf("%d ", sparse->row[ts_c[i]][j] - 1);
		}
		printf("\n");
	}
	printf("\n");

}


//return the location of nonzero entries and unsatisfied check sums
//void check_result(struct matrix *sparse, char *resultcode,double *lp,bool &decode_valid, vector<int>& ts_v, vector<int>& ts_c)
//{
//	int i,j,temp;
//	int a,b;
//
//	//     ts_v: nonzero entries;  ts_c: unsatisfied check-sums
//	ts_v.clear();
//	ts_c.clear();
//
//	decode_valid=true;
//
//	for (i=1; i<=sparse->col_num; i++)
//	{
//		if(lp[i]>=0.0)
//		{
//			resultcode[i]=1;
//			ts_v.push_back(i-1);
//		}
//		else resultcode[i]=0;
//	}
//
//	for (i=0; i<sparse->row_num; i++)
//	{
//		temp=0;
//		for (j=0; j<sparse->row_weigh[i]; j++)
//		{
//			temp=(temp+resultcode[sparse->row[i][j]]);
//			if (temp==2)
//				temp=0;
//		}
//
//		if (temp==0)  {}
//		else
//		{
//			decode_valid=false;
//			ts_c.push_back(i);
//		}
//	}
//
//	a=ts_v.size();
//	b=ts_c.size();
//
//	if (a!=0&&a<A_THRSD && b< B_THRSD)
//	{
//		check_graph(sparse, ts_v, ts_c);
//	}
//
//}



//Boost Decoding, hope it works. FROM CYC2QC draft!

void decode_bst(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter)
{
	int i, j, k, k1;
	double temp, temp1;


	vector<int> ts_v;		//count non-zero entries

	vector<int> ts_c;		//count un-sa check-sums





	for (i = 1; i <= sparse->col_num; i++)			//Initiate
		resultcode[i] = 0;

	for (i = 0; i <= sparse->col_num; i++)
		lp[i] = 0.0;

	for (i = 0; i <= sparse->col_num; i++)
		lp0[i] = 0.0;

	for (i = 1; i <= sparse->col_num; i++)
		lp0[i] = transcode[i] * 2 / noise;


	//	check_result(sparse, resultcode,lp0, decode_valid, ts_v, ts_c);


	for (i = 1; i <= sparse->col_num; i++)
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			v_line[i].link[j] = &c_line[sparse->col[i - 1][j - 1]];

	for (i = 1; i <= sparse->row_num; i++)
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
			c_line[i].link[j] = &v_line[sparse->row[i - 1][j - 1]];

	for (i = 1; i <= sparse->col_num; i++)		//Clean
	{
		v_line[i].step = 1;
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			v_line[i].r_value[j] = 0;
	}

	for (i = 1; i <= sparse->row_num; i++)		//Clean
	{
		c_line[i].step = 1;
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
			c_line[i].q_value[j] = 0;
	}


	for (i = 1; i <= sparse->col_num; i++)
	{
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
		{
			k = v_line[i].link[j]->step;
			v_line[i].link[j]->q_value[k] = lp0[i];
			v_line[i].link[j]->step++;
		}
	}

	for (k1 = 0; k1<max_iter; k1++)
	{
		iter++;
		for (i = 1; i <= sparse->row_num; i++) //Step Initiate
			c_line[i].step = 1;
		for (i = 1; i <= sparse->col_num; i++)
			v_line[i].step = 1;

		for (i = sparse->row_num; i >= 1; i--) //Calculate r
		{
			temp = c_line[i].q_value[1];
			for (j = 2; j <= sparse->row_weigh[i - 1]; j++)
				temp = LLR_add(temp, c_line[i].q_value[j]);

			for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
			{
				k = c_line[i].link[j]->step;
				c_line[i].link[j]->link[k] = &c_line[i];
				c_line[i].link[j]->r_value[k] = LLR_sub(temp, c_line[i].q_value[j]);
				c_line[i].link[j]->step++;
			}
		}

		for (i = 1; i <= sparse->col_num; i++)
			v_line[i].step = 1;

		for (i = 1; i <= sparse->col_num; i++)  //Calculate lp
		{
			temp = 0;
			for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
				temp = temp + v_line[i].r_value[j];
			lp[i] = lp0[i] + temp;
		}

		for (i = 1; i <= sparse->row_num; i++)
			c_line[i].step = 1;

		for (i = sparse->col_num; i >= 1; i--)   //Calculate q
		{
			for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			{
				k = v_line[i].link[j]->step;
				v_line[i].link[j]->link[k] = &v_line[i];
				temp1 = lp[i] - v_line[i].r_value[j];

				if (temp1>10)
					temp1 = 10;
				if (temp1<-10)
					temp1 = -10;


				v_line[i].link[j]->q_value[k] = temp1;
				v_line[i].link[j]->step++;
			}
		}

		check_result2(sparse, resultcode, lp, decode_valid);

		if (decode_valid)              //Parity Check
		{
			break;
		}
	}


	//	check_result(sparse, resultcode,lp, decode_valid, ts_v, ts_c);

}


void spa_qua(struct matrix *sparse, char *resultcode, struct v_node *v_line, struct c_node *c_line, double *lp, double *lp0, int &iter, double *transcode, double noise, bool &decode_valid, int max_iter, double q_step, int q_bits)
{
	int i, j, k, k1;
	double temp, temp1;

	for (i = 1; i <= sparse->col_num; i++)			//Initiate
		resultcode[i] = 0;

	for (i = 0; i <= sparse->col_num; i++)
		lp[i] = 0.0;

	for (i = 0; i <= sparse->col_num; i++)
		lp0[i] = 0.0;

	for (i = 1; i <= sparse->col_num; i++)
		lp0[i] = transcode[i];

	/*lp  y_i,quantization*/

	for (i = 1; i <= sparse->col_num; i++)                 //R(0)=q
	{
		lp0[i] = (int)(transcode[i] / q_step);
		if (transcode[i] - lp0[i] * q_step>(lp0[i] + 1)*q_step - transcode[i]) {
			lp0[i]++;
		}
		if (lp0[i]>(1 << (q_bits - 1)) - 1) {
			lp0[i] = (1 << (q_bits - 1)) - 1;
		}
		if (lp0[i]<-(1 << (q_bits - 1)) + 1) {
			lp0[i] = -(1 << (q_bits - 1)) + 1;
		}
	}

	for (i = 1; i <= sparse->col_num; i++) {
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++) {
			v_line[i].link[j] = &c_line[sparse->col[i - 1][j - 1]];
		}
	}

	for (i = 1; i <= sparse->row_num; i++) {
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++) {
			c_line[i].link[j] = &v_line[sparse->row[i - 1][j - 1]];
		}
	}


	for (i = 1; i <= sparse->col_num; i++)		//Clean
	{
		v_line[i].step = 1;
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			v_line[i].r_value[j] = 0;
	}

	for (i = 1; i <= sparse->row_num; i++)		//Clean
	{
		c_line[i].step = 1;
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
			c_line[i].q_value[j] = 0;
	}


	for (i = 1; i <= sparse->col_num; i++)
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
		{
			k = v_line[i].link[j]->step;
			v_line[i].link[j]->q_value[k] = lp0[i];
			v_line[i].link[j]->step++;
		}

	for (k1 = 0; k1<max_iter; k1++)
	{
		iter++;
		for (i = 1; i <= sparse->row_num; i++) //Step Initiate
			c_line[i].step = 1;
		for (i = 1; i <= sparse->col_num; i++)
			v_line[i].step = 1;

		for (i = sparse->row_num; i >= 1; i--) //Calculate r
		{
			temp = c_line[i].q_value[1];
			for (j = 2; j <= sparse->row_weigh[i - 1]; j++)
				temp = int(LLR_add(temp, c_line[i].q_value[j]));
			for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
			{
				k = c_line[i].link[j]->step;
				c_line[i].link[j]->link[k] = &c_line[i];
				c_line[i].link[j]->r_value[k] =int( LLR_sub(temp, c_line[i].q_value[j]));
				c_line[i].link[j]->step++;
			}
		}

		for (i = 1; i <= sparse->col_num; i++)
			v_line[i].step = 1;

		for (i = 1; i <= sparse->col_num; i++)  //Calculate lp
		{
			temp = 0;
			for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
				temp = temp + v_line[i].r_value[j];
			lp[i] = lp0[i] + temp;
		}

		for (i = 1; i <= sparse->row_num; i++)
			c_line[i].step = 1;

		for (i = sparse->col_num; i >= 1; i--)   //Calculate q
		{
			for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			{
				k = v_line[i].link[j]->step;
				v_line[i].link[j]->link[k] = &v_line[i];
				temp1 = lp[i] - v_line[i].r_value[j];
				if (temp1>(1 << (q_bits - 1)) - 1) {
					temp1 = (1 << (q_bits - 1)) - 1;
				}
				if (temp1<-(1 << (q_bits - 1)) + 1) {
					temp1 = -(1 << (q_bits - 1)) + 1;
				}
				v_line[i].link[j]->q_value[k] = temp1;
				v_line[i].link[j]->step++;
			}
		}

		check_result2(sparse, resultcode, lp, decode_valid);

		if (decode_valid)              //Parity Check
		{
			break;
		}
	}



}
double LLR_add_1(double in1, double in2)
{
	double	result;
	double	tp1 = 0, tp2 = 0;
	double	max1, max2;
	double  peak = 7.5;

	if (in1 > peak)
		in1 = peak;
	else if (in1 < -peak)
		in1 = -peak;
	if (in2 > peak)
		in2 = peak;
	else if (in2 < -peak)
		in2 = -peak;


	if (in1<in2)	max1 = in2;
	else		max1 = in1;

	if ((in1 + in2)<0)	max2 = 0;
	else			max2 = in1 + in2;

	//tp1 = exp(-fabs(in1-in2));	//if 0 modulate -1 and 1modulate 1
	tp1 = exp(-fabs(in1 + in2));		//if 0 modulate 1 and 1modulate -1

	tp1 = log(1 + tp1);
	//if 0 modulate -1 and 1modulate 1
	//tp2 = exp(-fabs(in1+in2));
	tp2 = exp(-fabs(in1 - in2));		//if 0 modulate 1 and 1modulate -1

	tp2 = log(1 + tp2);

	//result = max1 - max2 + tp1 - tp2;//if 0 modulate -1 and 1modulate 1
	result = -max1 + max2 + tp1 - tp2;//if 0 modulate 1 and 1 modulate -1

	return(result);
}

double LLR_sub_1(double in1, double in2)
{
	double	result;
	double	tp1 = 0, tp2 = 0;
	double	max1, max2;

	double  peak = 7.5;

	if (in1 > peak)
		in1 = peak;
	else if (in1 < -peak)
		in1 = -peak;
	if (in2 > peak)
		in2 = peak;
	else if (in2 < -peak)
		in2 = -peak;

	if (fabs(in1) == fabs(in2))
	{
		if (in1 >= 0) in1 = in1 + 0.000000000000000000000000000000000000000000000000000000000000000001;
		else in1 = in1 - 0.000000000000000000000000000000000000000000000000000000000000000001;
	}


	if (in1<in2)	max1 = in2;
	else		max1 = in1;

	if ((in1 + in2)<0)	max2 = 0;
	else			max2 = in1 + in2;

	//tp1 = exp(-fabs(in1-in2));	//if 0 modulate -1 and 1modulate 1
	tp1 = exp(-fabs(in1 + in2));		//if 0 modulate 1 and 1modulate -1
	if (tp1>0.9999999999999999)
		tp1 = 0.9999999999999999;
	tp1 = log(1 - tp1);

	//tp2 = exp(-fabs(in1+in2));	//if 0 modulate -1 and 1modulate 1
	tp2 = exp(-fabs(in1 - in2));		//if 0 modulate 1 and 1modulate -1
	if (tp2>0.9999999999999999)
		tp2 = 0.9999999999999999;
	tp2 = log(1 - tp2);

	//result = max1 - max2 + tp1 - tp2;//if 0 modulate -1 and 1modulate 1
	result = -max1 + max2 + tp1 - tp2;//if 0 modulate 1 and 1modulate -1

	return(result);
}
