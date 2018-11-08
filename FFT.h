#ifndef _FFT_H_
#define _FFT_H_

#include <stdint.h>

int32_t bitReverse(int32_t index, int32_t fftSize);
void fft(float *in, int32_t size, float *out);
void ifft(float *in, int32_t size, float *out);

#endif