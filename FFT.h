#ifndef _FFT_H_
#define _FFT_H_

/*
This file aimed at processing FFT.
All of the input and output use a float array to present a complex array, and the format is: 
{real, imag, real, imag....}
All the size and len is half of array's length.

*/

#include <stdint.h>
/*
Holding some temp vars during FFT to avoid cost of allocating and
freeing arrays.
*/
typedef struct{
    int32_t fftSize;
    float *array1;
    float *array2;
    float *W;
} FFT_instance;


void fftInit(int32_t fftLen, FFT_instance *instance);
int32_t bitReverse(int32_t index, int32_t fftSize);
void cfft(float *in, int32_t size, FFT_instance *instance, float *out, int8_t flag);

void fftClean(FFT_instance *instance);

void rfft(float *in, int32_t size, FFT_instance *instance, float *out);

void irfft(float *in, int32_t size, FFT_instance *instance, float *out);

void oddEvenSplite(float *in, float *outReal, float *outImag, int32_t size);

void computeW(uint32_t N, float* W, int8_t flag);

#endif