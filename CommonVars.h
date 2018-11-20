#ifndef _COMMON_VARS_H_
#define _COMMON_VARS_H_

#include <stdint.h>

#define PI 3.1415926
#define SIN_TABLE_LEN 2048
#define WINDOW_SIZE 512

extern float SIN_TABLE[];
extern float HAMMING_WINDOW[];

void initSinTable();
void initHammingWindow();
#endif