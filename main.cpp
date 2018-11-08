#include <iostream>
#include <stdlib.h>
#include "FFT.h"
#include "CommonVars.h"
#include <fstream>

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



int main()
{
    getSinTable();
    cout << "over" << endl;
    getchar();
}