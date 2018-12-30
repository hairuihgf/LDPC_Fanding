#define _CRT_SECURE_NO_DEPRECATE
#include <string.h>
#include <stdlib.h>
#include "function.h"

void ReadProfile(struct profile *profiles)
{
	FILE	*fin;
	char	ch;
	float	readtemp;

	fopen_s(& fin, "profile.txt", "r");

	/*Read Profile*/
	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%s", profiles->matrixfile);


	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &profiles->test_num);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%f", &readtemp);
	profiles->snr = (double)readtemp;

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%f", &readtemp);
	profiles->snr_step = (double)readtemp;

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d ", &profiles->seed1);
	fscanf(fin, "%d ", &profiles->seed2);
	fscanf(fin, "%d ", &profiles->seed3);
	


	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &profiles->leastframe);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &profiles->leasterror);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &profiles->codemethod);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &profiles->max_iter);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &profiles->err_era);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &profiles->burstlen);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &profiles->burstnum);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &profiles->de_alg);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%f", &readtemp);
	profiles->delta = readtemp;

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%f", &profiles->scale);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &profiles->codeword);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &profiles->q_bits);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%f", &readtemp);
	profiles->alpha = readtemp;


	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%f", &readtemp);
	profiles->q_step = readtemp;
	fscanf(fin, "%f", &readtemp);
	profiles->d_quasi = readtemp;
	fscanf(fin, "%d", &profiles->max_q_level);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &profiles->dsum);


	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &profiles->modulate_mode);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%f", &readtemp);
	profiles->DOPPLER = readtemp;

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%f", &readtemp);
	profiles->T_SAMPLE = readtemp;

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &profiles->RAYLEIGH_M);

	while ((ch = getc(fin)) != ':') {}
	fscanf(fin, "%d", &profiles->channel_style);

	fclose(fin);

	fin = NULL;
}