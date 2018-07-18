#include "function.h"

void count_gauss(struct matrix *sparse, char **gauss, int *exchange, int& rank)
{
	int  n, m, i, j, k, temp, col, change_col, swap;
	bool unrank;

	n = sparse->col_num;
	m = sparse->row_num;

	for (col = n; col >= n - m + 1; col--)
	{
		temp = 0;
		unrank = false;  //Not Full Rank

		for (i = 1; i <= m - n + col; i++) //Search the 'one in first column' row from bottom to up
		{
			if (gauss[i][col] == 1)
			{
				temp = i;
				break;
			}
		}

		if (temp == 0)  //If fail, skip to another column
		{
			change_col = 0; //The Index of Changing column

			for (j = col - 1; j >= 1; j--) //Search the changing column
			{
				for (k = 1; k <= m - n + col; k++)
				{
					if (gauss[k][j] == 1)
					{
						change_col = j;
						break;
					}
				}
				if (change_col != 0)break;
			}

			if (change_col == 0)
				unrank = true;  //No valid column, end the cycle
			else              //Change Column
			{
				for (j = 1; j <= m; j++)
				{
					swap = gauss[j][change_col];
					gauss[j][change_col] = gauss[j][col];
					gauss[j][col] = swap;
				}
				swap = exchange[change_col];
				exchange[change_col] = exchange[col];
				exchange[col] = swap;
				col++;     //Next column
			}
		}
		else if (temp == m - n + col) {} //Find the last column
		else                //GE
		{
			for (i = temp + 1; i <= m - n + col - 1; i++) //Search the 'one in first column' row from bottom to up
			{
				if (gauss[i][col] == 1)
				{
					for (j = 1; j <= n; j++)
						gauss[i][j] = (gauss[i][j] + gauss[temp][j]) % 2;
				}
			}

			if (gauss[m - n + col][col] == 1)   //Process the last column
			{
				for (j = 1; j <= n; j++)
					gauss[temp][j] = (gauss[m - n + col][j] + gauss[temp][j]) % 2;
			}
			else
			{
				for (j = 1; j <= n; j++)
				{
					swap = gauss[temp][j];
					gauss[temp][j] = gauss[m - n + col][j];
					gauss[m - n + col][j] = swap;
				}
			}
		}

		if (unrank)
		{
			rank = n - col;
			printf("%s", "Rank: ");
			printf("%d", rank);
			break;
		}
		else rank = sparse->row_num;
	}
	sparse->message_num = sparse->message_num + sparse->row_num - rank;
	printf("\r%s%d", "Information Bits: ", sparse->message_num);
}


void check_gauss(struct matrix *sparse, char **gauss, int& rank)
{

	int i, j;
	bool one, zero;
	printf("\n");

	if (rank == sparse->row_num)
	{
		one = true;
		for (j = sparse->col_num - sparse->row_num + 1; j <= sparse->col_num; j++)
		{
			i = j - sparse->col_num + sparse->row_num;
			if (gauss[i][j] == 0) one = false;
		}

		if (one) printf("matrix 1 right; ");
		else printf("matrix 1 false; ");
		printf("\n");
		zero = true;

		for (j = sparse->col_num - sparse->row_num + 2; j <= sparse->col_num; j++)
		{
			for (i = 1; i<(j - sparse->col_num + sparse->row_num); i++)
				if (gauss[i][j] == 1) zero = false;
		}

		if (zero) printf("matrix 0 right; ");
		else printf("matrix 0 false; ");
		printf("\n");
	}
	else
	{
		one = true;
		for (j = sparse->col_num - rank + 1; j <= sparse->col_num; j++)
		{
			i = j - sparse->col_num + sparse->row_num;
			if (gauss[i][j] == 0) one = false;//cout<<i<<j;
		}
		if (one) printf("matrix 1 right; ");
		else printf("matrix 1 false; ");

		zero = true;
		for (j = 1; j <= sparse->col_num; j++)
		{
			for (i = 1; i <= (sparse->row_num - rank); i++)
			{
				if (gauss[i][j] == 1) zero = false;//cout<<i<<j;
			}
		}

		if (zero) printf("matrix 0 right; ");
		else printf("matrix 0 false; ");
	}
}
