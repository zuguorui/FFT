#include "FFT.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include "CommonVars.h"


using namespace std;

float *W = NULL;

void fftInit(int32_t fftLen)
{
    W = (float *)calloc(fftLen * 2, sizeof(float));
}

void fftClean()
{
    if(W != NULL)
    {
        free(W);
        W = NULL;
    }
}

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

/*
in: it presents a complex, {real, imag, real, imag....};
size: count of complex numbers, it means in.length = 2 * size, and the same with out. size must be power of 2.
out: it is the same with in.
flag: 1 means fft, -1 means ifft.
*/
void fft(float *in, int32_t size, float *out, int8_t flag)
{
    //bit reverse
    int index;
    float tempReal, tempImag;
    out[0] = in[0];
    out[1] = in[1];
    out[2 * (size - 1)] = in[2 * (size - 1)];
    out[2 * (size - 1) + 1] = in[2 * (size - 1) + 1];
    for (int i = 1; i < size - 1; i++)
    {
        index = bitReverse(i, size);
        out[2 * i] = in[2 * index];
        out[2 * i + 1] = in[2 * index + 1];
        
    }

    uint32_t butterFlySize = 2;
    uint32_t cosOffset = SIN_TABLE_LEN / 4;
    uint32_t sinPos = 0, cosPos = 0;

    uint32_t step = 0;
    float aReal, aImag, bReal, bImag, wReal, wImag;

    while(butterFlySize <= size)
    {
        cout << "butterFlySize = " << butterFlySize << endl;
        //compute W factor
        for(int i = 0; i < butterFlySize >> 1; i++)
        {
            sinPos = ((SIN_TABLE_LEN * i) / butterFlySize) % SIN_TABLE_LEN;
            cosPos = (sinPos + cosOffset) % SIN_TABLE_LEN;
            W[2 * i] = SIN_TABLE[cosPos];
            W[2 * i + 1] = -1 * flag * SIN_TABLE[sinPos]; //be careful at this, W of fft is cos() - jsin(). W of ifft is cos() + jsin()
            // cout << "W[" << 2 * i << "] = " << W[2 * i] << endl;
            // cout << "W[" << 2 * i + 1 << "] = " << W[2 * i + 1] << endl;
        }
        int32_t halfWSize = butterFlySize >> 1;
        // for(int i = butterFlySize >> 1; i < butterFlySize; i++)
        // {
        //     W[2 * i] = -W[2 * (i - halfWSize)];
        //     W[2 * i + 1] = -W[2 * (i - halfWSize) + 1];
        // }
        
        
        for(int i = 0; i < size / butterFlySize; i++)
        {
            step = 0;
            while(step < halfWSize)
            {
                index = i * butterFlySize + step;
                aReal = out[2 * index];
                aImag = out[2 * index + 1];
                bReal = out[2 * (index + halfWSize)];
                bImag = out[2 * (index + halfWSize) + 1];
                wReal = W[2 * step];
                wImag = W[2 * step + 1];
                out[2 * index] = aReal + bReal * wReal - bImag * wImag;
                out[2 * index + 1] = aImag + bImag * wReal + bReal * wImag;
                wReal = -wReal;
                wImag = -wImag;
                out[2 * (index + halfWSize)] = aReal + bReal * wReal - bImag * wImag;
                out[2 * (index + halfWSize) + 1] = aImag + bImag * wReal + bReal * wImag;
                step++;
                
            }
        }

        butterFlySize *= 2;
        // cout << "out = [";
        // for(int i = 0; i < 2 * size; i++)
        // {
        //     cout << out[i] << ", ";
        // }
        // cout << "]" << endl;

    }
    //ifft
    if(flag == -1)
    {
        for(int i = 0; i < 2 * size; i++)
        {
            out[i] = out[i] / size;
        }
    }

}
/*
in: signal will be splited
outReal: the odd part of in, 
outImag: the even part of in
size: these three arrays must have the same length, and the length = 2 * size

This is a FFT future, a signal which only has real part through FFT, the real part of FFT result is even symmetry and the imag part is
odd symmetry. A signal which only has imag part through FFT,  the real part of FFT result is odd symmetry and the imag part is even 
symmetry.
We can splite a signal to odd-symmetry part and even-symmetry part by follow way:
Xo[n] = (X[n] - X[N-n]) / 2
Xe[n] = (X[n] + X[N-n]) / 2
Attention that FFT output is indexed [0, N-1], and it is periodic. It means when n = 0, X[N] = X[0].

The outputs of this function are outReal and outImag.
outReal: the FFT result of real part of the input signal. This equals you pass the real part as a real signal to FFT directly.
outImag: the FFT result of imag part of the input signal, but it dose not equals pass the imag part as real signal. If you want
treat this output as that, you must multiply a complex unit j to every element.
*/
void oddEvenSplite(float *in, float *outReal, float *outImag, int32_t size)
{
    //for real signal, real part is even symmetry, and imag part is odd symmetry
    outReal[0] = in[0];
    outReal[1] = 0;
    //for imag signal, real part is odd symmetry, and imag part is even symmetry
    outImag[0] = 0;
    outImag[1] = in[1];
    float r1, i1, r2, i2;
    for(int i = 1; i < size; i++)
    {
        r1 = in[2 * i];
        i1 = in[2 * i + 1];

        r2 = in[2 * (size - i)];
        i2 = in[2 * (size - i) + 1];

        outReal[2 * i] = (r1 + r2) / 2;
        outReal[2 * i + 1] = (i1 - i2) / 2;

        outImag[2 * i] = (r1 - r2) / 2;
        outImag[2 * i + 1] = (i1 + i2) / 2;

    }
}
