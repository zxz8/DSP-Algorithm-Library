# DSP-Algorithm-Library
A C++ library of not so common signal processing algorithms. The library is programmed in platform independent ANSI C++ and can be used directly from source code. I use the [DISLIN](https://www.dislin.de) library for the visualization of the 1d and 2d data, generated in the test routines for each single algorithm. DISLIN is a scientific plotting software developed at the Max Planck Institute for Solar System Research for displaying data as curves, surfaces, 3-D color plots, bar graphs, pie charts, contours and maps. In order to run all the test routines for the algorithms, you have first make sure that DISLIN is installed properly for your operating system. I have compiled all algorithms on the Raspberry Pi and prepared all test routines including DISLIN for this platform. But you can simply take the algorithms from the source code you are interested in, without the use of DISLIN. If DISLIN is available for the Apple ARM64 M architecture, I will recompile all test routines for MacOS ARM64 platform. The test routines are separated from the algorithms, and if you are only interested in a single algorithm, you can simply copy your algorithm source code files. Each folder contains the encapsulated algorithm class files and a usage example in the test.cpp file. A simple make file can be used to build the algorithm test environment on any platform. Please note that all the algorithms here are in the clauer::math namespace resp. clauer::io. All functions and classes are written with the same coding conventions for the class/function/variable names, and the formatting rules defined in the clauer coding styles for c++. Most of the header files have a rich documentation with an explanation of the algorithm, so please have always a look into the header files for the corresponding algorithm.  **Please note that I cannot guarantee the correctness of the algorithms, and some minor algorithms are not yet fully developed.**

## Data Plotter
The class collects all the useful function for the plotting of data structures into image files, on the screen and onto the printer. Therefor we use the dislin programming library from the Max Planck Institute for Solar System Research in Lindau/Karlsruhe and would linke to say thanks to the **author Helmut Michels**. This class provides the caller with functions for the plotting of one and more dimensional functions and plots onto images, the screen and the printer.
 
## Wave File Handler
This template class supports the caller with a basic interface to extract and write samples from MS wave format files. The easiest way to extract samples is the **autoReadWaveFile** function. Note that all reading vectors are generated inside the reading function and must be deleted outside the function explicitly. Note that the error handling is accumulative, which means that the error flag can bi given to more than one function, and it will only set to true if an error occurs.
 
## Cepstrum
Like most of the mathematical classes this class also has only static members so no explicite instanciation is necessary. The cepstrum class implement like the name says, the  cepstral analysis methods. This function calcualtes the cepstrum for the given time domain signal. The real cepstrum will be calculated corresponding folowing formular: Cepstrum = ABS ( IFFT ( LOG ( ABS ( FFT ( x) ) ) ). Note that the resulting cepstrum will be calculated via an zero-padded fast-fourier-transformation, which has in the back-transformation a symetry and therefore the resulting cepstrum will be have at the end an zero zone.
  
## Cohen Class Distributions

### Smoothed Pseudo Wigner Ville Distribution
This function implements the algorithm for the calculation of a "Discrete Smoothed Pseudo Wigner-Ville-Distribution" for the given double vector input signal. The function calculates the PWVT corresponding to the equation:

```
            /        /                                  -j2pi f tau`
SPWV(t,f) = | h(tau) | g(tau-t) x(mu+tau/2)x*(mu-tau/2)e          dmu dtau`
            /        /`
```

To prevent aliasing, the real input signal must be band limited with the nyquist frequency SR/2. The discrete wigner ville distribution has a band limitation of SR/4, so normally for complex signals the input time domain vector must be up-sampled to the double length. In our case we have a real input signal, so we use the complex transformation with the analytic signal as input which prevents aliasing at the frequency SR/4 and the mirroring of the image. For this case, we can switch between the transformation with the analytic signal and the real input signal. The upsampling method only required for complex input signals, which is not part of our implementation.
 
### Pseudo Margenau Hill Distribution
See the description of the Smoothed Pseudo Wigner Ville Distribution which is also valid for this cohen class distribution apart from the formula.
 
### Choi Williams Distribution
See the description of the Smoothed Pseudo Wigner Ville Distribution which is also valid for this cohen class distribution apart from the formula.

### Zaho Atlas Marks Distribution
See the description of the Smoothed Pseudo Wigner Ville Distribution which is also valid for this cohen class distribution apart from the formula.

## Convolution

## Correlation

## Damping Constant

## Digital Filter

## Envelope

## Fourier Transformation

## Harmonic Analyze 

## Impulse Response Extraction 

## Jitter Analysis

## Linear Predictive Coding

## Median Frequency

## Normalization

## Polygonal Chain

## Sample rate Converter

## Short Time Fourier Transformation

## Smooth

## Tools

## Wavelet Packet Decomposition

## Wavelet Transformation

## Zero Crossing Rate

## Math Utilities
