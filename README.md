# DSP-Algorithm-Library
A C++ library of not so common signal processing algorithms. The library is programmed in platform independent ANSI C++ and can be used directly from source code. I use the [DISLIN](https://www.dislin.de) library for the visualization of the 1d and 2d data, generated in the test routines for each single algorithm. DISLIN is a scientific plotting software developed at the Max Planck Institute for Solar System Research for displaying data as curves, surfaces, 3-D color plots, bar graphs, pie charts, contours and maps. In order to run all the test routines for the algorithms, you have first make sure that DISLIN is installed properly for your operating system. I have compiled all algorithms on the Raspberry Pi and prepared all test routines including DISLIN for this platform. But you can simply take the algorithms from the source code you are interested in, without the use of DISLIN. If DISLIN is available for the Apple ARM64 M architecture, I will recompile all test routines for MacOS ARM64 platform. The test routines are separated from the algorithms, and if you are only interested in a single algorithm, you can simply copy your algorithm source code files. Each folder contains the encapsulated algorithm class files and a usage example in the test.cpp file. A simple make file can be used to build the algorithm test environment on any platform. Please note that all the algorithms here are in the clauer::math namespace resp. clauer::io. All functions and classes are written with the same coding conventions for the class/function/variable names, and the formatting rules defined in the clauer coding styles for c++. Most of the header files have a rich documentation with an explanation of the algorithm, so please have always a look into the header files for the corresponding algorithm. I developed this allgorithms in 2005 for an europan union founded research project. **Please note that I cannot guarantee the correctness of the algorithms, and some minor algorithms are not yet fully developed.**

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
This class collects functions for the calculation of the convolution. One static function is for  the "normal" convolution calculation, corresponding to the well known convolution equation in the time domain. Another function calculates the so-called fast convolution in the frequency domain with help of the FFT. Like all other mathematic functions and classes, both functions are implemented as static functions, so no instantiation of the class is necessary. 

## Correlation
This class file defines the discrete cross-correlation and autocorrelation functions. Like most of the math classes, this functions are defined as static class functions, so no instantiation is needed. This very short algorithm is defined directly into the header file of the class definition as template implementations. As defined in the clauer coding styles, this class template header file has the suffix HPP instead of H for class files. NOTE: There are much more-faster implementations available, e.G. the FFT Autocorrelation which is also implemented here. This algorithm is not tested properly, so I cannot totally guaranty the correctness of the algorithm.

## Damping Constant
This class collects functions to find a suitable exponential expression for an exponential function which matches best the envelope of a given impulse response from a resonance analysis. The class implements only two functions, the main function which steers everything, and the exponential regression function which try to find a suitable exponential function. This function makes heavily use of the envelope function class-collection which was defined outside of this class file. The main function is implemented to find the damping constant autonomously and try to match the best values. The higher the damping constant, the faster signal converges to the zero line => Higher Damping. The function try to match the best physical units and the result should return a value with unit 1/s which can be guaranteed by the sampling rate. Another core function of this class-collection is the exponential regression, which here implemented with a fallback to a linear regression where the solution can found in every statistical book. The exponential regression returns the amplitude and damping factor of the function ```f(t)=A*exp(-bt)```, where b is the damping factor and the amplitude A can be ignored.

## Digital Filter
This class collects a set of useful digital filter functions, which are available for the four base types of digital filters (CriticalDamping, Bessel, Butterwirth and Tschebyscheff) implemnted with a low pass, a high pass, a band-pass and a bandstop filter. The filter coefficients defined into this class are predefined for the IIR filter. Higher filter orders are realized with multiple biquad filters in concatenation. A digital biquad filter is a second-order recursive linear filter, containing two poles and two zeros. "Biquad" is an abbreviation of "biquadratic", which refers to the fact that in the Z domain, its transfer function is the ratio of two quadratic functions. 
 
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
