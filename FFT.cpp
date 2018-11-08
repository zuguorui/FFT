#include "FFT.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>


using namespace std;

//fftSize must be power of 2.
int32_t bitReverse(int32_t index, int32_t fftSize)
{
    int32_t result = 0;
    int mPtr = 0;
    while(fftSize >> mPtr != 1)
    {
        result = (result << 1) | ((index >> mPtr) & 0x01);
        mPtr++;
    }
    return result;
}

void fft(float *in, int32_t size, float *out)
{

}

void ifft(float *in, int32_t size, float *out)
{

}