#ifndef _COMMON_VARS_H_
#define _COMMON_VARS_H_

#include <stdint.h>

#define PI 3.1415926
#define SIN_TABLE_LEN 1024

extern float SIN_TABLE[];

void initSinTable();
#endif