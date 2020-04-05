#pragma once
#include <complex>
#include <iostream>
#include <valarray>
#include "data.h"

const double PI = 3.141592653589793238460;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

struct MyArray {
    double X; // X axes
    double Y; // Y axes
};

int Time2Freq(TData  * Data, TData *FreqData);
void fft(CArray& x);

//template <class ValType>
std::complex<double>* FillComp(const double* Data, std::complex<double>* CData, int Len, int N);

