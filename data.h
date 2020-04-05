#pragma once
#include <fstream>
#include <iostream>

class TData
{
private:
    double* Xval; // X axes
    double* Yval; // Y axes
    int Chunk = 10000;
    int DataSize = Chunk;
    
    template <class TType>
    void SetNull(TType * Xyz, int Size) //void TData::SetNull(<TType>* Xyz, int Size)
    {
        for (int I = 0; I < Size; I++) // !!!!!! FIX THIS !!!!!!!!!!
        {
            Xyz[I] = 0;
        }
    }

public:
    //    CONSTRUCTORS
    TData(void) { Xval = NULL; Yval = NULL; }; // CONSTRUCTOR: Default
    TData(int Size); // CONSTRUCTOR: Create empty object with correct Array sizes for Xval and Yval
    TData(std::string FileName, std::string* ErrMsg); // CONSTRUCTOR: Construct and read in Data form file
    TData(std::string FileName); // CONSTRUCTOR: Construct and read in Data form file
    TData(TData& copy); // CONSTRUCTOR: Copy Constructor
    ~TData(void); // DESTRUCTOR: Destructor

    // FUNCTIONS
    void ReadData(std::string FileName, std::string* ErrMsg);


    // FUNCTIONS - ARRAY


    // FUNCTIONS - TRIVIAL
    double X(int Pos); // Returns the value of 'Xval' at position 'Pos'
    double Y(int Pos); // Returns the value of 'Yval' at position 'Pos'
    double* X(double* Output, int Start, int End);
    double* Y(double* Output, int Start, int End);
    const double *X(void); // Returns the value of 'Xval' at position 'Pos'
    const double *Y(void); // Returns the value of 'Yval' at position 'Pos'
    int Len(void); // Returns DataSize - The length of the array
    void SetX(int Pos,  double Val);
    void SetY(int Pos, double Val);
    void ReSize(int Val);
    
};
