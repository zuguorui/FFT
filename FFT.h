#ifndef _FFT_H_
#define _FFT_H_

/*
This file aimed at processing FFT.
All of the input and output use a float array to present a complex array, and the format is: 
{real, imag, real, imag....}
All the size and len is half of array's length.

*/

#include <stdint.h>

extern float *W;

void fftInit(int32_t fftLen);
int32_t bitReverse(int32_t index, int32_t fftSize);
void fft(float *in, int32_t size, float *out, int8_t flag);
void fftClean();

void oddEvenSplite(float *in, float *outOdd, float *outEven, int32_t size);

#endif