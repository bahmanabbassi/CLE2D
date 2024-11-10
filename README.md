# **CLE2D - Curvilinear Lineament Extraction 2D: A MATLAB GUI for 2D Curvilinear Lineament Extraction**

## **What is this repository for?**

CLE2D is a MATLAB GUI tool for 2D curvilinear lineament extraction. It harnesses advanced unsupervised source separation techniques, including Principal Component Analysis (PCA), Continuous Wavelet Transforms (CWT), and Bayesian Optimization, to provide a solution for curvilinear lineament  extraction in geophysical and geoscientific data analysis.

Please cite the software as:

Abbassi, B., Cheng, L.Z., 2024. Curvilinear Lineament Extraction: Bayesian Optimization of Principal Component Wavelet Analysis and Hysteresis Thresholding, Computers & Geosciences, 2024,105768, ISSN 0098-3004, https://doi.org/10.1016/j.cageo.2024.105768

## **How do I get set up?**

This program is designed to run on any Windows-based personal computer with at least 8 GB of RAM. Increasing the RAM size allows larger images to be processed at once. Since large matrices are operating in this program, the read/write speed of the storage is also essential. A solid-state drive (SSD) with a non-volatile memory express interface (NVMe) is recommended.
CLE2D is provided in MATLAB M-File format. M-Files require MATLAB 2024a version 24.1 (and later) to run. To use the program, copy the downloaded repository folder (M-files and data sets) to your desired directory. Locate the M-Files in the current folder of MATLAB and then type CLE2D in the MATLAB Command Window. This interface offers comprehensive tools for curvilinear lineament extraction with a user-friendly design, ensuring that users can easily manipulate and visualize geoscientific data.

Required MATLAB Toolboxes:

  - Statistics and Machine Learning Toolbox
  - Wavelet Toolbox
  - Optimization Toolbox

## **Usage**

The CLE2D interface comprises several key functionalities, summarized below:

1. Coordinates
  - Max Lat/Min Lat: Set the maximum and minimum latitude.
  - Max Lon/Min Lon: Set the maximum and minimum longitude.
  - Method: Choose the method for coordinate input (Rectangular coordinates).
  - 
Coordinates can be read from a prepared text file with the following format:
  - First line: 		Min Longitude (Min Lon)
  - Second line: 		Max Longitude (Max Lon)
  - Third line: 		Min Latitude (Min Lat)
  - Fourth line: 		Max Latitude (Max Lat)

For example, in the case above:

- -77.8
- -77.2
- 53.15
- 53.45

2. Spacing

  - Spacing: Define the spacing value in arcseconds for both longitude and latitude. In Quebec, each arcsecond of latitude is approximately 33 meters, and each arcsecond of longitude is about 17 meters. The program automatically adjusts these ratios for any location on Earth to ensure optimal image scaling.
  - Data Filter: Apply a filter to the input data to smooth it further.

3. Input Point Data

  - 2D Interpolation: Select the interpolation method for the data points, which can be regularly spaced or irregularly sampled. The symbol '#' indicates the number of input data points.
  - Xn / Yn: Retrieve the number of pixels in the X and Y coordinates after interpolation.

4. Digitized Lineaments

  - Targets Type: Choose the type of targets (digitized faults) either as point data in .csv format (fault densities) or as images of the digitized lines.
  - Spacing: Define the spacing for digitized lineaments, which should be equal to or less than the main spacing for input data.
  - Cut-off 1: Automatically set the cut-off value to confine a buffer zone for the lineaments.
  - Cut-off 2: Set a second cut-off value to narrow down the automatically generated buffer zone, which is useful if the target fault density map is too thick for upcoming optimizations.
  - Target Filter: Apply a filter to the target data.

5.  Spectral Feature Extraction

This window is divided into two main sections. First, a 2D CWT (Continuous Wavelet Transform) breaks down the desired inputs into raw spectral features. Users can adjust parameters such as the number of scales (na), scale dilation, and the angles at which the mother wavelet will operate. Each mother wavelet structure then automatically computes the scales vector, corresponding frequencies, and angles vector. Isotropic and anisotropic mother wavelets are available depending on the study's objectives. After the CWT process, the program automatically determines the number of CWT features. Second, S-PCA (Spectral Principal Component Analysis) is provided for spectral source separation. Users must decide if dimensionality reduction is necessary.
  - CWT Inputs: Choose to merge the uploaded point datasets for CWT and S-PCA.
  - Mother Wavelet: Select the type of mother wavelet.
  - OrderX (n) / OrderY (m): Define the order of X and Y for the Difference of Gaussian (DOG) mother wavelet.
  - Change order by scales: Option to change the orders automatically by scales, setting the values inversely proportional to the scales to detect higher frequency features.
  - Number of Scales (na): Set the number of wavelet scales, with Bayesian optimization available to fine-tune this parameter.
  - Scale Dilation: Define the scale dilation factor, allowing access to longer scales while keeping the number of scales constant. For example, with na = 4 and Scale Dilation = 1, scales = {1, 2, 3, 4}; with Scale Dilation = 2, scales = {1, 3, 5, 7}.
  - WSFR: Set the Wavelet Smoothness Filter Ratio, which is necessary to avoid interpolation artifacts in the wavelet coefficient features.
  - CWT number of Angles: Set the number of directions for which the mother wavelet can surf in 2D space.
  - Scales (a): Define the scale vector automatically or manually.
  - CWT Angles: Define the angle vector automatically or manually.
  - β: Set the angle of symmetry for CWT, indicating the wavelet's rotational invariance or symmetry.
  - CWT: Run the continuous wavelet transform.
  - S-PCA: Run the spectral principal component analysis after CWT.
  - Number of CWT Fs: Set the number of CWT features.
  - DR to: Specify the dimension reduction target, with Bayesian optimization available to fine-tune this parameter.

6. Select Features

  - Point Data: Use point data for lineament extraction.
  - CWT: Use CWT results for lineament extraction.
  - S-PCA: Use results of spectral principal component analysis for lineament extraction.

7. Lineament Extraction

  - Line Res: Define the resolution of the lineament extraction output. Larger resolutions result in crisper and smoother results but at a higher computational cost.
  - SF # of Angles: Define the number of angles for calculating Aspect in the hysteresis thresholding procedure.
  - Bayesian Opt MaxIter: Set the maximum iterations for Bayesian optimization.
  - w: Define the width of the step filter, typically 10 percent of the largest pixel numbers, with Bayesian optimization available to fine-tune this parameter.
  - VSFW: Variability of the Step Filtering Widths, set by default to 0.25, to determine the variability of w in correlation with the complexity of the extracted spectral features. Bayesian optimization helps fine-tune this parameter.
  - Lineaments: Extract lineaments.
  - AutoLine: Automatically fine-tune the hyperparameters (na, WSFR, DR, w, and VSFW).

8. Plot

  - Plot Results: Option to plot the results.
  - Filter: Apply a filter to the plotted results.
  - Close All: Close all plots.
  - Clear All: Clear all settings.

## **Input/output formats**

CLE2D supports input data in point datasets formatted as .csv files. The data must follow an XYZ-style format where:

  - X column: Represents the longitudes.
  - Y column: Represents the latitudes.
  - Z column: Represents the values of the geoscientific image, such as reflectance, magnetic field intensities, etc.
  
To use these inputs, users must first define the minimum and maximum values of the geographic coordinate system in decimal degrees for both latitude and longitude. Additionally, the spacing in arc seconds for both latitude and longitude needs to be specified. This setup ensures that the program accurately interprets the spatial extent and resolution of the input data.
The outputs generated by CLE2D are primarily in MATLAB .fig format. This format includes:

  - Extracted Spectral Features: The spectral features extracted from the input data using the 2D Continuous Wavelet Transform (CWT) and Spectral Principal Component Analysis (S-PCA).
  - Retrieved Curvilinear Lineaments: The detected curvilinear lineaments were identified through the lineament extraction process.

These .fig files visually represent the processed data, enabling users to analyze and interpret the results directly within the MATLAB environment.

## **Who do I talk to?**

Bahman Abbassi, 
Université du Québec en Abitibi-Témiscamingue, 
bahman.abbassi@uqat.ca


## **License**

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
