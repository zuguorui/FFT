#include <iostream>
#include <stdlib.h>
#include "FFT.h"
#include "CommonVars.h"
#include <fstream>
#include <random>

using namespace std;



void getSinTable()
{

    ofstream out;
    out.open("./sin_table.txt", ios_base::out);
    initSinTable();
    cout << SIN_TABLE << endl;
    for (int i = 0; i < SIN_TABLE_LEN; i++)
    {
        //cout << SIN_TABLE[i] << ",";
        printf("%f, ", SIN_TABLE[i]);
        out << SIN_TABLE[i] << ",";
    }
    out.close();
}

void fftTest()
{
    ofstream fout;
    fout.open("./fft_result.txt", ios_base::out);

    int fftSize = 256;
    float data[2 * fftSize] = {0};
    float out[2 * fftSize] = {0};
    float ifftOut[2 * fftSize] = {0};
    fftInit(fftSize);
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
    cfft(data, fftSize, out, 1);
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
    cfft(out, fftSize, ifftOut, -1);
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

    float outEven[2 * fftSize] = {0};
    float outOdd[2 * fftSize] = {0};
    oddEvenSplite(out, outOdd, outEven, fftSize);
    fout << "outOdd: " << endl;
    ptr = 0;
    while (ptr < 2 * fftSize)
    {
        float real = outOdd[ptr];
        float imag = outOdd[ptr + 1];
        fout << real << "+" << imag << "i, ";
        ptr += 2;
    }
    fout << endl;

    fout << "outEven: " << endl;
    ptr = 0;
    while (ptr < 2 * fftSize)
    {
        float real = outEven[ptr];
        float imag = outEven[ptr + 1];
        fout << real << "+" << imag << "i, ";
        ptr += 2;
    }
    fout << endl;

    fout.close();

}


int main()
{
    fftTest();
    cout << endl << "over" << endl;
    getchar();
}