#pragma once
#include <complex>
#include <iostream>
#include <valarray>
#include <vector>
#include "data.h"


using namespace std;
const double PI = 3.141592653589793238460;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

struct MyArray {
    double X; // X axes
    double Y; // Y axes
};

int Time2Freq(TData  * Data, TData *FreqData);
void MyFFT(TData* Data, TData* FreqData, int Start, int End);
void fft(CArray& x);

template <class ValType>
complex<ValType>* Num2Complex(const ValType * Data, complex<ValType>* CData, int Len)
{
    for (int I = 0; I < Len; I++) 
    {
        CData[I] = Data[I];
    }
    return CData;
}


/*---------------------------------------------------------------------------------------------------
FUNCTION: Time2Freq
INPUT:
DESCRIPTION:
RETURN:
----------------------------------------------------------------------------------------------------*/
int Time2Freq(TData * Data, TData *FreqData) {

    int FSize = Data->Len();
    int Start = 0;
    int End = pow(2, (int)log2(FSize)) - 1; // Get the correct chunk size for the fft function

        while (FSize != pow(2, (int)log2(FSize))) {
        
        cout << "Start at: " << Start << " Up to: " << End << " = " << End - Start + 1 << endl;

        int CSize = End - Start + 1;
        double* SizedData = new double[CSize];
        complex<double>* CompData = new Complex[CSize];
        valarray<Complex> Cdata(Num2Complex(Data->Y(SizedData, Start, End), CompData, CSize), CSize);
        delete[] SizedData;
        delete[] CompData;

        fft(Cdata);
        cout << endl << "Output from FFT" << endl;
        for (int I = 0; I < CSize / 2; I++) {
            if (Cdata[I].real() > 0.01) {
                cout << I << " : " << Cdata[I] << endl;
            }
        }

        float Factor = (Data->Len() - 1) / (Data->X(Data->Len() - 1) - Data->X(0)) / CSize; // Factor = [Sample Freqeuncy]/[fft data size]
        for (int I = 0; I <= CSize / 2; I++) {
            FreqData->SetX(I, Factor * I);
            FreqData->SetY(I, Cdata[I].real());
            //      cout << I << " " ;
        }

        cout << FreqData->X(0) << " " << FreqData->X(1) << " " << FreqData->X(2) << endl;
        
        Start = 1 + End; // Move the Start position to the next value in Data
        FSize = FSize - pow(2, (int)log2(FSize)); // Get the Size of the remaining data;
        End = Start + pow(2, (int)log2(FSize))-1; // Get the correct End position in Data according to FSize

        if (FSize < 15) {
            cout << "  Done with last set of 2^N. " << endl;
            break;
        }
        cout << "Next Data Chunk size: " << pow(2, (int)log2(FSize)) << endl;
    }
    cout << "FSize " << FSize << "    power " << pow(2, (int)log2(FSize)) << endl;
    if (FSize > 15) {
        cout << endl << "Start at " << Start << " Up to: " << Start + FSize << endl;
    }
    return 0;
}
/*---------------------------------------------------------------------------------------------------
FUNCTION: MyFFT
INPUT:
DESCRIPTION:
RETURN:
----------------------------------------------------------------------------------------------------*/
void MyFFT(TData* Data, TData *FreqData, int Start, int End) {

    int FSize = End - Start + 1;
    double* SizedData = new double[FSize];
    complex<double>* CompData = new Complex[FSize];
    valarray<Complex> Cdata(Num2Complex(Data->Y(SizedData, Start, End), CompData, FSize), FSize);
    delete[] SizedData;
    delete[] CompData;

    /*for (int I = 0; I < FSize; I++) {
        cout << I << " : " << Cdata[I] << endl;
    }*/
 //   cout << "Data in MyFFT Cdata[0]:" << Cdata[0] << " Cdata[end]: " << Cdata[FSize-1] << endl;
    fft(Cdata);
    cout << endl << "Output from FFT" << endl;
    for (int I = 0; I < FSize /2; I++) {
        if (Cdata[I].real() > 0.01) {
            cout << I << " : " << Cdata[I] << endl;
        }
    }

    float Factor = (Data->Len() - 1) / (Data->X(Data->Len() - 1) - Data->X(0)) / FSize; // Factor = [Sample Freqeuncy]/[fft data size]
    for (int I = 0; I <= FSize / 2; I++) {
        FreqData->SetX(I, Factor * I);
        FreqData->SetY(I, Cdata[I].real());
 //      cout << I << " " ;
    }
    
}

/*---------------------------------------------------------------------------------------------------
FUNCTION: fft
INPUT:
DESCRIPTION: !!! WARNING  !!!!
             This only works for datastes with NumberOfPoints = 2^n.
             Therefore massage the data to be available in chunks of 2^n.
RETURN:
----------------------------------------------------------------------------------------------------*/
// !!! WARNING  !!!! 
// This only works for datastes with NumberOfPoints = 2^n.   Therefore massage the data to be available in chunks of 2^n.
void fft(CArray& x)
{
    // conjugate the complex numbers
    x = x.apply(std::conj);

    // DFT
    unsigned int N = x.size(), k = N, n;
    double thetaT = 3.14159265358979323846264338328L / N;
    Complex phiT = Complex(cos(thetaT), -sin(thetaT)), T;
    while (k > 1)
    {
        n = k;
        k >>= 1;
        phiT = phiT * phiT;
        T = 1.0L;
        for (unsigned int l = 0; l < k; l++)
        {
            for (unsigned int a = l; a < N; a += n)
            {
                unsigned int b = a + k;
                Complex t = x[a] - x[b];
                x[a] += x[b];
                x[b] = t * T;
            }
            T *= phiT;
        }
    }
    // Decimate
    unsigned int m = (unsigned int)log2(N);
    for (unsigned int a = 0; a < N; a++)
    {
        unsigned int b = a;
        // Reverse bits
        b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
        b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
        b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
        b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
        b = ((b >> 16) | (b << 16)) >> (32 - m);
        if (b > a)
        {
            Complex t = x[a];
            x[a] = x[b];
            x[b] = t;
        }
    }

    // conjugate the complex numbers again
    x = x.apply(std::conj);

    // scale the numbers
    x /= x.size();
}

