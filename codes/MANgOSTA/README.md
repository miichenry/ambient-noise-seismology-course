# MANgOSTA

TEST A Python package for Multiscale Ambient NOiSe TomogrAphy.

## Description

MANgOSTA is a software package for performing ambient noise tomography from cross-correlations calculated by frequency bands using wavelets. It's developed by Iván Cabrera-Pérez (_INVOLCAN, Consejo Insular de la Energía de Gran Canaria_), Luca D'Auria (_INVOLCAN, ITER_) and Jean Soubestre (_INVOLCAN_).

## Structure

The package is divided into modules:

1. **Finder :** Find all available data from seismic stations pairs.
2. **Censor :** Remove coherent sources (_e.g._ earthquakes) using network covariance matrix analysis.
3. **DWTO :** Decompose cleaned signals on particular frequency bands using wavelets.
4. **Correlator :* Performs a normalization through one-bit or spectral whitening and performs the cross-correlation.
5. **Stacker :* Performs the stack of the cross-correlation result.

Parameters to be modified:

1. **Wavelet :** Select Wavelet Families (Symlet, Daubechies,...).
2. **fmin :** Select the minimum frequency.
3. **fmax :** Select the maximum frequency.
4. **wpwr :** Number of points used for the window.
5. **datadir :** Directory of the folder where is the data. By default is data.
6. **dataext :** Type of extension that the program has to search.
7. **Stafile :** Station data in extension .csv. (Name: NAME; Longitude: LON; Latitud: LAT; Altitude: ALT)
8. **Finder :** Find all available data from seismic stations pairs.
9. **Censor :** Remove coherent sources (_e.g._ earthquakes) using network covariance matrix analysis. (1.minimum frequency, 2.maximum frequency, 3.coherence limit).
10. **Whitening :** Is a Spectral Whitening softened. Select the type of spectral whitening(between 0 and 1).
11. **Whitening uniform :** This is the traditional Spectral Withening. Change the whitening to 1.
12. **One-bit :** One-bit normalization. If one-bit is 0 it's the traditional one-bit and if it's 1 the one-bit is not done.
13. **ccdir :** Folder where the cross-correlation will be stored.
14. **figdir :** Folder where the figure will be stored.
15. **correlator :** Perform the cross-correlation by making the Discrete Wavelet Transform. You can select if you want to make it: Transversal-Transversal (TT), Traversal-Radial (TR), etc.
16. **correlator_trad :** Perform the traditional cross-correlation without doing the Discrete Wavelet Transform. Of course, you can select the component TT, RT, etc. 

            


TODO
* procesamiento DWT
* SQLite

