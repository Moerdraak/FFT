#include "fft.h"

// START OF MAIN FUNCTION
int main()
{

    string MyErr = "All Good"; // Default Value
    TData BaseData("Data.txt", &MyErr); // Or just use TData BaseData("Data.txt");

    // Create a FreqData array for holding the data
    // It must be half the size of the Time data (Nyquist)
    TData FreqData;

    // Call the function that performs the FFT and Massages the data into 
    // an array ready to plot Amplitude versus Frequency
    Time2Freq(&BaseData, &FreqData);

    cout << endl <<  "NOW BACK TO THE REAL WORLD" << endl << endl;
    for (int I = 0; I < 8192 / 2; I++) {
        if (FreqData.Y(I) > 0.5) {
        cout << "Freq: " << FreqData.X(I) << " Ampl: " << FreqData.Y(I) << endl;
        }
    }
    return 0;
}

