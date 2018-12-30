#include "function.h"




void init_sparse(struct matrix *sparse, char* filename)
{
	int i, j;
	int s1, s2;
	FILE *fp;


	fopen_s(&fp, filename, "r");
	fscanf_s(fp, "%d", &sparse->col_num);

	fscanf_s(fp, "%d", &sparse->row_num);
	sparse->message_num = sparse->col_num - sparse->row_num;

	printf("\r%s: %d \n", "row_num",sparse->row_num);
	printf("\r%s: %d \n", "col_num", sparse->col_num);

	fscanf_s(fp, "%d", &sparse->max_colweigh);
	fscanf_s(fp, "%d", &sparse->max_rowweigh);

	//	sparse->col_weigh=new int[sparse->col_num];
	sparse->col_weigh = (int *)malloc(sizeof(int)*(sparse->col_num));

	for (i = 0; i<abs(sparse->col_num); i++)
		fscanf_s(fp, "%d", &sparse->col_weigh[i]);

	//	sparse->row_weigh=new int[sparse->row_num];
	sparse->row_weigh = (int *)malloc(sizeof(int)*(sparse->row_num));

	for (i = 0; i<abs(sparse->row_num); i++)
	{
		fscanf_s(fp, "%d", &sparse->row_weigh[i]);
	}


	//	sparse->col=new int *[sparse->col_num];
	sparse->col = (int **)malloc(sizeof(int *)*(sparse->col_num));

	for (i = 0; i<abs(sparse->col_num); i++)
	{
		//	sparse->col[i]=new int[sparse->col_weigh[i]];
		sparse->col[i] = (int *)malloc(sizeof(int)*(abs(sparse->col_weigh[i])));
		for (j = 0; j<abs(sparse->col_weigh[i]); j++)
			fscanf_s(fp, "%d", &sparse->col[i][j]);
		//printf("%d ",sparse->col[i][j]);}printf("\n");
	}

	//	sparse->row=new int *[sparse->row_num];
	sparse->row = (int **)malloc(sizeof(int *)*(sparse->row_num));

	for (i = 0; i<abs(sparse->row_num); i++)
	{
		//		sparse->row[i]=new int[sparse->row_weigh[i]];
		sparse->row[i] = (int *)malloc(sizeof(int)*(abs(sparse->row_weigh[i])));
		fscanf_s(fp, "%d", &s1);
		fscanf_s(fp, "%d", &s2);
		if (s1<s2)
		{
			sparse->row[i][0] = s1;
			sparse->row[i][1] = s2;
			for (j = 2; j<abs(sparse->row_weigh[i]); j++)
				fscanf_s(fp, "%d", &sparse->row[i][j]);
		}
		else
		{
			sparse->row[i][sparse->row_weigh[i] - 1] = s1;
			sparse->row[i][sparse->row_weigh[i] - 2] = s2;
			for (j = sparse->row_weigh[i] - 3; j >= 0; j--)
				fscanf_s(fp, "%d", &sparse->row[i][j]);
		}
	}
}



void init_gauss(struct matrix *sparse, char **gauss, int *exchange)
{
	int i, j;

	for (i = 0; i <= sparse->row_num; i++)
		for (j = 0; j <= sparse->col_num; j++)
			gauss[i][j] = 0;

	for (i = 1; i <= sparse->row_num; i++)
	{
		for (j = 0; j<sparse->row_weigh[i - 1]; j++)
			gauss[i][sparse->row[i - 1][j]] = 1;
	}

}

void init_decode(struct matrix *sparse, struct v_node *v_line, struct c_node *c_line)
{
	int i, j;

	//Construct Graph
	for (i = 1; i <= sparse->col_num; i++)
	{
		//	v_line[i].m_value=new int [sparse->col_weigh[i-1]+1];
		v_line[i].m_value = (int *)malloc(sizeof(int)*(sparse->col_weigh[i - 1] + 1));

		//	v_line[i].r_value=new double [sparse->col_weigh[i-1]+1];
		v_line[i].r_value = (double *)malloc(sizeof(double)*(sparse->col_weigh[i - 1] + 1));

		//	v_line[i].link=new c_node*[sparse->col_weigh[i-1]+1];
		v_line[i].link = (c_node **)malloc(sizeof(c_node *)*(sparse->col_weigh[i - 1] + 1));

		v_line[i].sybol = i;
		v_line[i].step = 1;
	}

	for (i = 1; i <= sparse->row_num; i++)
	{
		//	c_line[i].m_value=new int [sparse->row_weigh[i-1]+1];
		c_line[i].m_value = (int *)malloc(sizeof(int)*(sparse->row_weigh[i - 1] + 1));

		//	c_line[i].q_value=new double [sparse->row_weigh[i-1]+1];
		c_line[i].q_value = (double *)malloc(sizeof(double)*(sparse->row_weigh[i - 1] + 1));

		//	c_line[i].link=new v_node*[sparse->row_weigh[i-1]+1];
		c_line[i].link = (v_node **)malloc(sizeof(v_node *)*(sparse->row_weigh[i - 1] + 1));

		c_line[i].sybol = i;
		c_line[i].step = 1;
	}

	for (i = 1; i <= sparse->col_num; i++)
		for (j = 1; j <= sparse->col_weigh[i - 1]; j++)
			v_line[i].link[j] = &c_line[sparse->col[i - 1][j - 1]];

	for (i = 1; i <= sparse->row_num; i++)
		for (j = 1; j <= sparse->row_weigh[i - 1]; j++)
			c_line[i].link[j] = &v_line[sparse->row[i - 1][j - 1]];
}