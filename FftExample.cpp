#include "fft.h"

// START OF MAIN FUNCTION
int main()
{

    string MyErr = "All Good"; // Default Value

    // Read in data from a txt file
    TData BaseData("Data.txt", &MyErr); // Or just use TData BaseData("Data.txt");
    
    // Create pointer of TData to store answer
    TData FreqData; 

    // Call the function that converts time based data to Frequency domain
    Time2Freq(&BaseData, &FreqData); 
    
    // Screen dump of the significant points
    for (int I = 0; I < FreqData.Len(); I++) {
        if (FreqData.Y(I) > 0.5) {
        cout << "Freq: " << FreqData.X(I) << " Ampl: " << FreqData.Y(I) << endl;
        }
    }
    return 0;
}

