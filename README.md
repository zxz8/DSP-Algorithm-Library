# DSP-Algorithm-Library
A C++ library of not so common signal processing algorithms. The library is programmed in platform independent ANSI C++ and can be used directly from source code. I use the [DISLIN](https://www.dislin.de) library for the visualization of the 1d and 2d data, generated in the test routines for each single algorithm. DISLIN is a scientific plotting software developed at the Max Planck Institute for Solar System Research for displaying data as curves, surfaces, 3-D color plots, bar graphs, pie charts, contours and maps. In order to run all the test routines for the algorithms, you have first make sure that DISLIN is installed properly for your operating system. I have compiled all algorithms on the Raspberry Pi and prepared all test routines including DISLIN for this platform. But you can simply take the algorithms from the source code you are interested in, without the use of DISLIN. If DISLIN is available for the Apple ARM64 M architecture, I will recompile all test routines for MacOS ARM64 platform. The test routines are separated from the algorithms, and if you are only interested in a single algorithm, you can simply copy your algorithm source code files. Each folder contains the encapsulated algorithm class files and a usage example in the test.cpp file. A simple make file can be used to build the algorithm test environment on any platform. Please note that all the algorithms here are in the clauer::math namespace resp. clauer::io. All functions and classes are written with the same coding conventions for the class/function/variable names, and the formatting rules defined in the clauer coding styles for c++. **Please note that I cannot guarantee the correctness of the algorithms, and some minor algorithms are not yet fully developed.**

## Data Plotter

## Wave File Handler

## Cepstrum

## Cohen Class Distributions

### Smoothed Pseudo Wigner Ville Distribution

### Pseudo Margenau Hill Distribution

### Choi Williams Distribution

### Zaho Atlas Marks Distribution

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
