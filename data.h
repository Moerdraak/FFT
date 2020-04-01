#pragma once
#include <fstream>

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

//****************************************************************************
// CONSTRUCTOR: Construct and read in Data form file
//****************************************************************************
TData::TData(int Size) // CONSTRUCTOR: Default 
{
    Xval = new double[Size]; 
    Yval = new double[Size];
    DataSize = Size;

    SetNull(Xval, Size);
    SetNull(Yval, Size);
}

//****************************************************************************
// CONSTRUCTOR: Construct and read in Data form file
//****************************************************************************
TData::TData(std::string FileName) 
{
    std::string Ignore;
    ReadData(FileName, &Ignore);
}

//****************************************************************************
// CONSTRUCTOR: Construct and read in Data form file
//****************************************************************************
TData::TData(std::string FileName, std::string* ErrMsg)
{
    ReadData(FileName, ErrMsg);
}

//****************************************************************************
// CONSTRUCTOR: Construct and read in Data form file
//****************************************************************************
TData::TData(TData& copy) // CONSTRUCTOR: Copy Constructor
{
    DataSize = copy.DataSize;
    Xval = new double[copy.DataSize];
    Yval = new double[copy.DataSize];

    std::copy(copy.Xval, copy.Xval + DataSize, Xval);
    std::copy(copy.Yval, copy.Yval + DataSize, Yval);

}
//****************************************************************************
// CONSTRUCTOR: Construct and read in Data form file
//****************************************************************************
TData::~TData(void) // DESTRUCTOR: Destructor
{
    delete[] Xval;
    delete[] Yval;
}

/*---------------------------------------------------------------------------------------------------
  FUNCTION: ReadData
  INPUT:
  DESCRIPTION:
  RETURN:
  ----------------------------------------------------------------------------------------------------*/
void TData::ReSize(int Val)
{
    if (Xval == NULL)
    {
        Xval = new double[Val];
        Yval = new double[Val];
        SetNull(Xval, Val);
        SetNull(Yval, Val);
    }
    else if (DataSize == Val)
    {
        // Do Nothing
    }
    else if (DataSize < Val)
    {
        double* Tmp = new double[Val];
        SetNull(Tmp, Val);

        std::copy(Xval, Xval + DataSize, Tmp);
        delete[] Xval;
        Xval = Tmp;

        Tmp = new double[Val];
        SetNull(Tmp, Val);
        std::copy(Yval, Yval + DataSize, Tmp);
        delete[] Yval;
        Yval = Tmp;
    }
    else
    {
        double* Tmp = new double[Val];

        std::copy(Xval, Xval + Val, Tmp);
        delete[] Xval;
        Xval = Tmp;

        Tmp = new double[Val];
        std::copy(Yval, Yval + Val, Tmp);
        delete[] Yval;
        Yval = Tmp;
    }
    DataSize = Val;
}

/*---------------------------------------------------------------------------------------------------
  FUNCTION: ReadData
  INPUT:
  DESCRIPTION:
  RETURN:
  ----------------------------------------------------------------------------------------------------*/
void TData::ReadData(std::string FileName, std::string* ErrMsg)
{
    std::ifstream DataFile;
    DataFile.open(FileName);

    if (!DataFile)
    {
        *ErrMsg = "No such File";
    }

    Xval = new double[DataSize];
    Yval = new double[DataSize];

    int I = 0;
    while (DataFile >> Xval[I] >> Yval[I])
    {
        I++;
        static int Cnt = 1;

        // If the data is more than our first guess, Increase the size.       
        if (I >= DataSize)
        {
            DataSize = Chunk * ++Cnt;

            // Increase and copy Xval
            double* Tmp = new double[DataSize]; // Create new memory space
            std::copy(Xval, Xval + DataSize, Tmp); // Copy Xval into new space
            delete[] Xval; // Get rid of the old allocated memory - otherwise this is a memory leak
            Xval = Tmp; // Point Xval to new space

            // Increase and copy Yval                
            Tmp = new double[DataSize];
            std::copy(Yval, Yval + DataSize, Tmp);
            delete[] Yval; // Get rid of the old allocated memory - otherwise this is a memory leak
            Yval = Tmp;
        }
    }
    DataSize = I; // Set DataSize to the correct size of the data that wa read in

    // Now properly resize the data for Xval
    double* Tmp = new double[DataSize];
    std::copy(Xval, Xval + DataSize, Tmp);
    delete[] Xval; // Get rid of the old allocated memory - otherwise this is a memory leak
    Xval = Tmp;

    // Now properly resize the data for Yval
    Tmp = new double[DataSize];
    std::copy(Yval, Yval + DataSize, Tmp);
    delete[] Yval; // Get rid of the old allocated memory - otherwise this is a memory leak
    Yval = Tmp;

    DataFile.close();
}
/*---------------------------------------------------------------------------------------------------
  FUNCTION: TData::X(Position)
  INPUT: Position
  DESCRIPTION: Returns the value of 'Xval' at position 'Pos'
  RETURN: X[Pos]
  ----------------------------------------------------------------------------------------------------*/
  // Function to give X value
double TData::X(int Pos) {
    return Xval[Pos];
}

/*---------------------------------------------------------------------------------------------------
FUNCTION: Range of X
INPUT:
DESCRIPTION:
RETURN:
----------------------------------------------------------------------------------------------------*/
double * TData::X(double *Output, int Start, int End)
{
    int Cnt = 0;
    for (int I = Start; I < End; I++)
    {
        Output[Cnt] = Xval[I];
        Cnt++;
    }
    return Output;
}

/*---------------------------------------------------------------------------------------------------
FUNCTION: Range of Y
INPUT:
DESCRIPTION:
RETURN:
----------------------------------------------------------------------------------------------------*/
double* TData::Y(double* Output, int Start, int End)
{
    int Cnt = 0;
    for (int I = Start; I <= End; I++)
    {
        Output[Cnt] = Yval[I];
        Cnt++;
    }
    return Output;
}
/*---------------------------------------------------------------------------------------------------
FUNCTION: TData::Y(Position)
INPUT: Position
DESCRIPTION: Returns the value of 'Yval' at position 'Pos'
RETURN: Y[Pos]
----------------------------------------------------------------------------------------------------*/
double TData::Y(int Pos) {
    return Yval[Pos];
}

/*---------------------------------------------------------------------------------------------------
  FUNCTION: TData::X(void)
  INPUT: Position
  DESCRIPTION: Returns a copy of X
  RETURN: X
  ----------------------------------------------------------------------------------------------------*/
  // Function to give X value
const double  * TData::X(void) {

    return (const double*)Xval; 
}

/*---------------------------------------------------------------------------------------------------
  FUNCTION: TData::Y(void)
  INPUT: Position
  DESCRIPTION: Returns a copy of Y
  RETURN: Y
  ----------------------------------------------------------------------------------------------------*/
  // Function to give Y value
const double * TData::Y(void) {

    return (const double*)Yval;
}

/*---------------------------------------------------------------------------------------------------
FUNCTION: TData::Len(void)
INPUT: Position
DESCRIPTION: Give the length of the X array which is obviously the same as the Y array
RETURN: DataSize
----------------------------------------------------------------------------------------------------*/
int TData::Len(void) {
    return DataSize;
}

/*---------------------------------------------------------------------------------------------------
FUNCTION: SetX
INPUT: 
DESCRIPTION: 
RETURN: 
----------------------------------------------------------------------------------------------------*/
void TData::SetX(int Pos, double Val)
{
    Xval[Pos] = Val;
}

/*---------------------------------------------------------------------------------------------------
FUNCTION: SetY
INPUT:
DESCRIPTION:
RETURN:
----------------------------------------------------------------------------------------------------*/
void TData::SetY(int Pos, double Val)
{
    Yval[Pos] = Val;
}
