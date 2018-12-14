#include <iostream>
#include <stdlib.h>
#include "FFT.h"
#include "CommonVars.h"
#include <fstream>
#include <random>
#include <memory.h>
#include <typeinfo.h>

#define MIN(x, y) (x > y ? y : x)
#define MAX(x, y) (x > y ? x : y)

using namespace std;



void getSinTable()
{

    ofstream fout;
    fout.open("./sin_table.txt", ios_base::out);
    initSinTable();
    cout << SIN_TABLE << endl;
    for (int i = 0; i < SIN_TABLE_LEN; i++)
    {
        //cout << SIN_TABLE[i] << ",";
        printf("%ff, ", SIN_TABLE[i]);
        fout << SIN_TABLE[i] << "f,";
    }
    fout.close();
}

void getHammingWindow()
{
    ofstream fout;
    fout.open("./hamming_window.txt", ios_base::out);
    initHammingWindow();
    
    for (int i = 0; i < WINDOW_SIZE; i++)
    {
        cout << HAMMING_WINDOW[i] << ",";

        fout << HAMMING_WINDOW[i] << ",";
    }
    fout.close();
}

void getHanningWindow()
{
    ofstream fout;
    fout.open("./hanning_window.txt", ios_base::out);
    initHanningWindow();

    for (int i = 0; i < WINDOW_SIZE; i++)
    {
        cout << HANNING_WINDOW[i] << ",";

        fout << HANNING_WINDOW[i] << ",";
    }
    fout.close();
}

void cfftTest()
{
    ofstream fout;
    fout.open("./cfft_result.txt", ios_base::out);

    int fftSize = 256;
    float data[2 * fftSize] = {0};
    float out[2 * fftSize] = {0};
    float ifftOut[2 * fftSize] = {0};
    FFT_instance instance;
    fftInit(fftSize, &instance);
    default_random_engine e;
    uniform_real_distribution<float> u(-10, 10);
    fout << "data: " << endl;
    int ptr = 0;
    while(ptr < 2 * fftSize)
    {
        float real = u(e);
        float imag = u(e);
        fout << real << "+" << imag << "i, ";
        data[ptr] = real;
        data[ptr + 1] = imag;
        ptr += 2;
    }
    fout << endl;
    cfft(data, fftSize, &instance, out, 1);
    fout << "out: " << endl;
    ptr = 0;
    while (ptr < 2 * fftSize)
    {
        float real = out[ptr];
        float imag = out[ptr + 1];
        fout << real << "+" << imag << "i, ";
        ptr += 2;
    }
    fout << endl;
    cfft(out, fftSize, &instance, ifftOut, -1);
    fout << "ifft: " << endl;
    ptr = 0;
    while (ptr < 2 * fftSize)
    {
        float real = ifftOut[ptr];
        float imag = ifftOut[ptr + 1];
        fout << real << "+" << imag << "i, ";
        ptr += 2;
    }
    fout << endl;

    float outImag[2 * fftSize] = {0};
    float outReal[2 * fftSize] = {0};
    oddEvenSplite(out, outReal, outImag, fftSize);
    fout << "outReal: " << endl;
    ptr = 0;
    while (ptr < 2 * fftSize)
    {
        float real = outReal[ptr];
        float imag = outReal[ptr + 1];
        fout << real << "+" << imag << "i, ";
        ptr += 2;
    }
    fout << endl;

    fout << "outImag: " << endl;
    ptr = 0;
    while (ptr < 2 * fftSize)
    {
        float real = outImag[ptr];
        float imag = outImag[ptr + 1];
        fout << real << "+" << imag << "i, ";
        ptr += 2;
    }
    fout << endl;

    fout.close();
    fftClean(&instance);
}

void rfftTest()
{
    ofstream fout;
    fout.open("./rfft_result.txt", ios_base::out);

    int fftSize = 8;
    float rData[fftSize] = {0};
    
    float rfftOut[2 * fftSize] = {0};
    
    float irfftOut[fftSize] = {0};
    FFT_instance *instance = (FFT_instance*)calloc(1, sizeof(FFT_instance));
    fftInit(fftSize, instance);
    default_random_engine e;
    uniform_real_distribution<float> u(-10, 10);
    
    for(int i = 0; i < fftSize; i++)
    {
        float real = u(e);
        rData[i] = real;
    }

    rfft(rData, fftSize, instance, rfftOut);
    irfft(rfftOut, fftSize, instance, irfftOut);

    fout << "rData:" << endl;
    for(int i = 0; i < fftSize; i++)
    {
        fout << rData[i] << ", ";
    }
    fout << endl;
    fout << "realPart:" << endl;
    for(int i = 0; i < fftSize; i+= 2)
    {
        fout << rData[i] << ", ";
    }
    fout << endl;
    fout << "imagPart:" << endl;
    for(int i = 1; i < fftSize; i+= 2)
    {
        fout << rData[i] << ", ";
    }
    fout << endl;
    fout << "rfftOut:" << endl;
    for(int i = 0; i < fftSize; i++)
    {
        fout << rfftOut[2 * i] << "+" << rfftOut[2 * i + 1] << "i, ";
    }
    fout << endl;

    

    fout << "irfftOut:" << endl;
    for(int i = 0; i < fftSize; i++)
    {
        fout << irfftOut[i] << ", ";
    }
    fout << endl;

    fout.close();
    fftClean(instance);
    free(instance);
}


void mmecpyTest()
{
    float a[8] = {1, 2, 3, 4, 5, 6, 7, 8};
    float b[8] = {10, 11, 12, 13, 14, 15, 16, 17};
    memcpy(a, a + 4, 4 * sizeof(float));
    memmove(a + 4, b, 4 * sizeof(float));
    
    for(int i = 0; i < 8; i++)
    {
        cout << a[i] << " ";
    }
}

void minMaxTest()
{
    int a = 2;
    int b = 9;
    cout << MIN(a, b) << endl;
    cout << MAX(a, b);
}

void roundFloatTest()
{
    float a = 5.99999;
    int b = (int)a;
    cout << b;
}



int main()
{
    //minMaxTest();
    //getHammingWindow();
    //mmecpyTest();
    getSinTable();
    //cfftTest();
    //getHanningWindow();
    //roundFloatTest();
    cout << endl << "over" << endl;
    getchar();
}