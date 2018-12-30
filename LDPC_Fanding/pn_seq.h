#pragma once
void pn_init(char *a, char *reg, int length, bool reset)
{
	int i, j;
	char sum;

	if (reset)
	{
		for (i = 0; i <= 23; i++)
		{
			reg[i] = 0;
		}

		reg[0] = 1;
	}

	for (i = 1; i <= length; i++)
	{
		a[i] = reg[23];
		sum = reg[1] ^ reg[2] ^ reg[3] ^ reg[4] ^ reg[5] ^ reg[6]
			^ reg[7] ^ reg[8] ^ reg[9] ^ reg[10] ^ reg[11] ^ reg[12] ^ reg[13] ^ reg[14] ^ reg[15] ^ reg[16]
			^ reg[17] ^ reg[18] ^ reg[19] ^ reg[20] ^ reg[21] ^ reg[22] ^ reg[23];
		for (j = 23; j >= 1; j--)
		{
			reg[j] = reg[j - 1];
		}
		reg[0] = sum;
	}
}

