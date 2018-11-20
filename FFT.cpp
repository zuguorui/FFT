#include "FFT.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include "CommonVars.h"


using namespace std;


void fftInit(int32_t fftLen, FFT_instance *instance)
{
    instance->fftSize = fftLen;
    instance->array1 = (float*)calloc(fftLen * 2, sizeof(float));
    instance->array2 = (float*)calloc(fftLen * 2, sizeof(float));
    instance->W = (float*)calloc(fftLen, sizeof(float));
}

void fftClean(FFT_instance *instance)
{
    if(instance == NULL)
    {
        return;
    }
    instance->fftSize = 0;
    if(instance->array1 != NULL)
    {
        free(instance->array1);
        instance->array1 = NULL;
    }
    if (instance->array2 != NULL)
    {
        free(instance->array2);
        instance->array2 = NULL;
    }
    if(instance->W != NULL)
    {
        free(instance->W);
        instance->W = NULL;
    }
}

//fftSize must be power of 2.
int32_t bitReverse(int32_t index, int32_t fftSize)
{
    int32_t result = 0;
    int mPtr = 0;
    while((fftSize >> mPtr) != 1)
    {
        result = (result << 1) | ((index >> mPtr) & 0x01);
        mPtr++;
    }
    return result;
}

/*
Compute W factories of FFT or IFFT
N: butterfly-compute size
flag: 1: compute W for FFT. -1: compute for IFFT
*/
void computeW(uint32_t N, float *W, int8_t flag)
{
    uint32_t cosOffset = SIN_TABLE_LEN / 4;
    uint32_t sinPos = 0, cosPos = 0;
    int halfSize = N >> 1;
    //compute W factor, cause W[k + N/2] = -W[k], we only need to compute half of W.
    for (int i = 0; i < halfSize; i++)
    {
        sinPos = ((SIN_TABLE_LEN * i) / N) % SIN_TABLE_LEN;
        cosPos = (sinPos + cosOffset) % SIN_TABLE_LEN;
        W[2 * i] = SIN_TABLE[cosPos];
        W[2 * i + 1] = -1 * flag * SIN_TABLE[sinPos]; //be careful at this, W of fft is cos() - jsin(). W of ifft is cos() + jsin()
        // cout << "W[" << 2 * i << "] = " << W[2 * i] << endl;
        // cout << "W[" << 2 * i + 1 << "] = " << W[2 * i + 1] << endl;
    }
}

/*
This is standard real FFT, and will be convert to a complex FFT.
in: a real sequence in length = size
size: FFT size, must be power of 2
out: a float array presents a complex sequence in format: {real, imag, real, imag ....}. length = 2 * size

*/
void rfft(float *in, int32_t size, FFT_instance *instance, float *out)
{
    /*
    According to the format of input and output of cfft, we don't need to do anything on input to convert a N point real FFT to a 
    N/2 point complex FFT. Just make size = size / 2.
    */
    
    int32_t cSize = size >> 1;
    float *W = instance->W;
    float *realOut = instance->array1;
    float *imagOut = instance->array2;
    //float *buffer = (float *)calloc(size, sizeof(float));

    cfft(in, cSize, instance, out, 1);

    oddEvenSplite(out, realOut, imagOut, cSize);

    // cout << "realOut: " << endl;
    // for(int i = 0; i < cSize; i++)
    // {
    //     cout << realOut[2 * i] << "+" << realOut[2 * i + 1] << "i, ";
    // }
    // cout << endl;
    // cout << "imagOut: " << endl;
    // for(int i = 0; i < cSize; i++)
    // {
    //     cout << imagOut[2 * i + 1] << "+" << -imagOut[2 * i] << "i, ";
    // }
    // cout << endl;
    computeW(size, W, 1);
    float aReal, aImag, bReal, bImag, wReal, wImag;

    // at this step, we need to merge two domien outputs of real and imag part.
    // Attention to the operate on imagOut. According to the outImag of oddEvenSplite(), we must convert it to normal real FFT result.
    
    for(int i = 0; i < cSize; i++)
    {
        aReal = realOut[2 * i];
        aImag = realOut[2 * i + 1];

        bReal = imagOut[2 * i + 1];
        bImag = -imagOut[2 * i];

        wReal = W[2 * i];
        wImag = W[2 * i + 1];

        out[2 * i] = aReal + bReal * wReal - bImag * wImag;
        out[2 * i + 1] = aImag + bImag * wReal + bReal * wImag;

        wReal = -wReal;
        wImag = -wImag;

        out[2 * (i + cSize)] = aReal + bReal * wReal - bImag * wImag;
        out[2 * (i + cSize) + 1] = aImag + bImag * wReal + bReal * wImag;
    }
    cout << "out:" << endl;
    for (int i = 0; i < size; i++)
    {
        cout << out[2 * i] << "+" << out[2 * i + 1] << "i, ";
    }
    cout << endl;

}

/*
This is standard real IFFT, and will be convert to complex IFFT.
in: a float array presents a complex sequence in format: {real, imag, real, imag....}, length = 2*size
size: FFT size, must be power of 2
out: a real sequence in length = size.

In fact, this equals cfft. Because the frequence data is always complex. We just pick the real part of output of cfft.
*/
void irfft(float *in, int32_t size, FFT_instance *instance, float *out)
{
    float *tempOut = instance->array1;
    // cout << "in:" << endl;
    // for(int i = 0; i < size; i++)
    // {
    //     cout << in[2 * i] << "+" << in[2 * i + 1] << "i, ";
    // }
    // cout << endl;
    cfft(in, size, instance, tempOut, -1);
    for(int i = 0; i < size; i++)
    {
        out[i] = tempOut[2 * i];
    }
    // cout << "tempOut:" << endl;
    // for (int i = 0; i < size; i++)
    // {
    //     cout << tempOut[2 * i] << "+" << tempOut[2 * i + 1] << "i, ";
    // }
    // cout << endl;
}

/*
This is standard complex FFT
in: it presents a complex, {real, imag, real, imag....};
size: count of complex numbers, it means in.length = 2 * size, and the same with out. size must be power of 2.
out: it is the same with in.
flag: 1 means fft, -1 means ifft.
*/
void cfft(float *in, int32_t size, FFT_instance *instance, float *out, int8_t flag)
{
    //bit reverse
    int index;
    float *W = instance->W;
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
        //cout << "butterFlySize = " << butterFlySize << endl;
        computeW(butterFlySize, W, flag);
        int32_t halfWSize = butterFlySize >> 1;
        
        
        
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
outImag: the FFT result of imag part of the input signal, but it dose not equals pass the imag part as real signal to FFT. If you want
treat this output as that, you must multiply a -1j to every element. That means you need to exchange elements in outImag:
temp = outImag[2 * n];
outImag[2 * n] = outImag[2 * n + 1];
outImag[2 * n + 1] = -temp;
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
