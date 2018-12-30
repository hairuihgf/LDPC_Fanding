#include "function.h"

void encode(struct matrix *sparse, char **gauss, int *exchange, char *message, char *code)
{
	int i, j, col, temp;
	char *codeswap;

	codeswap = (char*)malloc(sizeof(char)*(sparse->col_num + 1));

	for (i = 1; i <= sparse->col_num; i++)
		codeswap[i] = 0;

	for (i = 1; i <= sparse->message_num; i++)
	{
		if (message[i] == 0) code[i] = 0;
		else code[i] = 1;
	}


	for (j = (sparse->message_num + 1); j <= sparse->col_num; j++)
	{
		temp = 0;
		i = j - sparse->col_num + sparse->row_num;
		for (col = 1; col <= j - 1; col++)
			temp = (temp + gauss[i][col] * code[col]) % 2;
		code[j] = temp;
	}

	for (i = 1; i <= sparse->col_num; i++)
		codeswap[exchange[i]] = code[i];

	for (i = 1; i <= sparse->col_num; i++)
		code[i] = codeswap[i];

	free(codeswap);
}


bool check_code(struct matrix *sparse, char*code)
{
	int i, j, temp;
	bool code_valid;

	code_valid = true;

	for (i = 0; i<sparse->row_num; i++)
	{
		temp = 0;
		for (j = 0; j<sparse->row_weigh[i]; j++)
		{
			temp = (temp + code[sparse->row[i][j]]) % 2;
		}

		if (temp == 0) {}
		else {
			code_valid = false;
		}
	}
	return code_valid;
}



/*
The purpose of this function is to convert an unsigned
binary number to reflected binary Gray code.
*/
unsigned short bin2gray(unsigned short num)
{
	return (num >> 1) ^ num;
}


/*
The purpose of this function is to convert a reflected binary
Gray code number to a binary number.
*/
unsigned short gray2bin(unsigned short num)
{
	unsigned short temp = num ^ (num >> 8);
	temp ^= (temp >> 4);
	temp ^= (temp >> 2);
	temp ^= (temp >> 1);
	return temp;
}

// convert a d-bit positive decimal number into a d-bit binary sequence
// mapping: bin[0]<- MSB, bin[d-1]<- LSB, ex: 4 ->{bin[0],bin[1],bin[2]}={100}
void dec2bin(int dec, int d, char* bin)
{
	for (int idx = d - 1; idx >= 0; idx--) {
		bin[idx] = dec & 1;
		dec >>= 1;
	}
}

// convert a d-bit binary sequence into a d-bit un-signed decimal number
// mapping:  bin[0] -> MSB, bin[d-1] -> LSB
int bin2dec(int d, char* bin)
{
	int dec = 0;
	for (int i = 0; i<d; i++)
		dec += (bin[i] << (d - 1 - i));
	return dec;
}