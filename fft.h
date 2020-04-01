#pragma once
#include <complex>
#include <iostream>
#include <valarray>
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
void fft(CArray& x);

//template <class ValType>
complex<double>* FillComp(const double* Data, complex<double>* CData, int Len, int N);


/*---------------------------------------------------------------------------------------------------
FUNCTION: Time2Freq
INPUT:
DESCRIPTION:
RETURN:
----------------------------------------------------------------------------------------------------*/
int Time2Freq(TData * Data, TData *FreqData) {

    // First determine the biggest 2^N value that fits inside this dataset
    int FSize = pow(2, floor(log2(Data->Len())));
    int NyqSize = FSize / 2; // since FSize will always be even this will always be an integer
    int RemainingSize = Data->Len();
    int N = log2(FSize); // Get the value of N for 2^N of the Intermediate Dataset
    
    // Add one for the last value e.g for 100Hz it would display 99.9Hz if you don't add this.
    FreqData->ReSize(NyqSize + 1);

    complex<double>* CompData = new complex<double>[FSize]; 
    // Factor = [Sample Frequency]/[fft data size]
    float Factor = ((Data->Len() - 1) / (Data->X(Data->Len() - 1) - Data->X(0))) / FSize;
    //  float Freq = (Data->Len() - 1) / (Data->X(Data->Len() - 1) - Data->X(0));

    // Populate FreqData's Freq collumn
    for (int I = 0; I < NyqSize + 1; I++) {
        FreqData->SetX(I, Factor * I);
    }

    int Start = 0;
    int End = FSize - 1; // 
    int CSize = End - Start + 1; // Chunk Size

    while (CSize > 8) 
    {
        double* SizedData = new double[CSize];
        FillComp(Data->Y(SizedData, Start, End), CompData, CSize, N);
        delete[] SizedData; // cleaning up

        CArray Cdata(CompData, FSize);
        fft(Cdata);

        int ScaleN = pow(2, N - (int)log2(CSize));

        for (int I = 0; I <= NyqSize; I++) {
            
            FreqData->SetY(I, FreqData->Y(I) + Cdata[I].real() / ScaleN);
        }
        
        Start = 1 + End; // Move the Start position to the next value in Data
        RemainingSize = RemainingSize - CSize;
        CSize = pow(2, floor(log2(RemainingSize)));
        End = Start + CSize - 1; // Get the correct End position in Data according to FSize
     }
    delete[] CompData; // cleaning up

    return 0;
}

/*---------------------------------------------------------------------------------------------------
FUNCTION: FillComp
INPUT:
DESCRIPTION:
RETURN:
----------------------------------------------------------------------------------------------------*/
complex<double>* FillComp(const double* Data, complex<double>* CData, int Len, int N)
{
    int Repeat = pow(2, N - (int)log2(Len));

    int Pos = 0;
    for (int I = 0; I < Repeat; I++)
    {
        for (int J = 0; J < Len; J++)
        {
            CData[Pos++] = Data[J];
        }
    }
    return CData;
}

/*---------------------------------------------------------------------------------------------------
FUNCTION: fft
INPUT:
DESCRIPTION: !!! WARNING  !!!!
             This only works for datastes with NumberOfPoints = 2^n.
             Therefore massage the data to be available in chunks of 2^n.
RETURN:
----------------------------------------------------------------------------------------------------*/
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

