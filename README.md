# **CLE2D – Curvilinear Lineament Extraction 2D: A MATLAB GUI for 2D Curvilinear Lineament Extraction**

## **What is this repository for?**

CLE2D is a MATLAB GUI designed for the extraction of curvilinear lineaments in 2D.  
It relies on advanced unsupervised source separation techniques, including Principal Component Analysis (PCA), 2D Continuous Wavelet Transforms (CWT), and Bayesian optimization, to provide a complete workflow for analysing and extracting lineaments from geophysical and geoscientific data.

Please cite the software as:

Abbassi, B. (2024). CLE2D: Curvilinear Lineament Extraction: Bayesian Optimization of Principal Component Wavelet Analysis and Hysteresis Thresholding. GitHub. https://github.com/bahmanabbassi/CLE2D

## **System requirements and setup**

This program is designed to run on any Windows-based personal computer with at least 8 GB of RAM.  
Increasing the RAM allows larger images to be processed at once. Since the program operates on large matrices, disk read/write speed is also important. A fast NVMe SSD is recommended.

CLE2D is provided as MATLAB files (`.m` and `.mlapp`). It requires MATLAB 2024a (version 24.1) or later.  
To use the program:

1. Copy the repository folder (MATLAB files and datasets) to your preferred directory.
2. In MATLAB, set the CLE2D folder as the **Current Folder**.
3. Launch the interface by typing `CLE2D` in the MATLAB Command Window.

The CLE2D interface offers a complete set of tools for curvilinear lineament extraction, with a layout oriented toward geoscientific data analysis.

Required MATLAB Toolboxes:

  - Statistics and Machine Learning Toolbox  
  - Wavelet Toolbox  
  - Optimization Toolbox

### Dependencies and acknowledgments

This project uses components of **Yet Another Wavelet Toolbox (YAWTB)**,  
copyright (C) 2001–2002 by the YAWTB team.
The original license headers have been preserved in the relevant `.m` files,  
in accordance with the YAWTB licensing conditions.

CLE2D also uses:

- A MATLAB implementation of **PCA / ICA** (including `fastICA.m`) developed by **Brian Moore**:  
  Brian Moore (2026). *PCA and ICA Package* (`https://www.mathworks.com/matlabcentral/fileexchange/38300-pca-and-ica-package`),  
  MATLAB Central File Exchange. Retrieved March 16, 2026.  

- A **geological fault detection** code (hysteresis thresholding / step filtering)  
  developed by **Costas Panagiotakis**:  
  Costas Panagiotakis (2026). *Detection of Geological Faults*  
  (`https://www.mathworks.com/matlabcentral/fileexchange/64693-detection-of-geological-faults`),  
  MATLAB Central File Exchange. Retrieved March 16, 2026.  

- The toolbox **YAWTB: Yet Another Wavelet ToolBox**  
  (`https://github.com/jacquesdurden/yawtb`).  

In addition to citing CLE2D itself, we recommend citing these packages when the corresponding
modules (PCA/ICA, fault detection, or YAWTB-based wavelets) are explicitly used.

## **Using the interface**

The CLE2D interface comprises several main windows and functionalities, summarised below.

1. Coordinates  
  - Max Lat / Min Lat: set the maximum and minimum latitude.  
  - Max Lon / Min Lon: set the maximum and minimum longitude.  
  - Method: choose the coordinate input method (rectangular coordinates).
  
Coordinates can be read from a prepared text file with the following format:  
  - 1st line: minimum longitude (Min Lon)  
  - 2nd line: maximum longitude (Max Lon)  
  - 3rd line: minimum latitude (Min Lat)  
  - 4th line: maximum latitude (Max Lat)

For example:

- -77.8  
- -77.2  
- 53.15  
- 53.45  

2. Spacing  

  - Spacing: define the spacing value in arcseconds for longitude and latitude. In Québec, one arcsecond of latitude corresponds to approximately 33 metres, and one arcsecond of longitude to about 17 metres. The program automatically adjusts these ratios for any location on Earth to ensure appropriate image scaling.  
  - Filter: apply an additional filter to the input data to further smooth it.

3. Input point data  

  - 2D Interpolation: select the interpolation method for the input points, whether regularly spaced or irregularly sampled. The symbol `#` indicates the number of input data points.  
  - Xn / Yn: retrieve the number of pixels in X and Y after interpolation.

4. Digitised lineaments  

  - Targets Type: choose the type of targets (digitised faults) either as point data in `.csv` format (fault densities) or as images of the digitised lineaments.  
  - Spacing: define the spacing for the digitised lineaments, which should be less than or equal to the main spacing used for the input data.  
  - Cut-off 1: automatically set a first cut-off value to define a buffer zone around the lineaments.  
  - Cut-off 2: adjust a second cut-off value to refine this buffer zone, useful when the fault-density map is too wide for subsequent optimisations.  
  - Filter: apply a filter to the target data.

5. Spectral feature extraction  

This window is divided into two main parts.  
First, a 2D CWT (2D Continuous Wavelet Transform) decomposes the inputs into raw spectral features. The user can adjust parameters such as the number of scales (`na`), scale dilation, and the angles at which the mother wavelet operates.  
Each wavelet structure automatically computes the scale vector, corresponding frequencies, and angle vector. Isotropic and anisotropic mother wavelets are available depending on the study objectives.  
After the CWT, the program automatically determines the number of CWT features.

Second, spectral PCA (S‑PCA) is used for spectral source separation. The user decides whether dimensionality reduction is required.

  - CWT Inputs: option to merge loaded point datasets for use in both CWT and S‑PCA.  
  - OrderX (n) / OrderY (m): order in X and Y for the “Difference of Gaussian” (DOG) mother wavelet.  
  - Change orders by scales: option to automatically vary the orders as a function of scale (values inversely proportional to the scales) to emphasise high‑frequency features.  
  - Number of scales (na): number of wavelet scales, with optional Bayesian optimisation.  
  - Scale dilation: dilation factor on the scales, allowing access to longer scales while keeping the same number of scales (e.g. `na = 4`, dilation = 2 → scales {1, 3, 5, 7}).  
  - WSFR (Wavelet Smoothness Filter Ratio): smoothing ratio for wavelet coefficients to avoid interpolation artefacts.  
  - CWT number of angles: number of directions along which the mother wavelet propagates in 2D.  
  - Scales (a): scale vector (automatic or user‑defined).  
  - CWT Angles: angle vector (automatic or user‑defined).  
  - β: symmetry angle for the CWT, indicating rotational invariance or symmetry of the wavelet.  
  - CWT: run the 2D continuous wavelet transform.  
  - S‑PCA: run spectral principal component analysis after the CWT.  
  - Number of CWT Fs: number of CWT features.  
  - DR to: dimensionality‑reduction target (optionally tuned via Bayesian optimisation).

6. Feature selection  

  - Point Data: use point data as features for lineament extraction.  
  - CWT: use the CWT features.  
  - S‑PCA: use the spectral PCA features.

7. Lineament extraction  

  - Line resolution: resolution of the lineament‑extraction output. Higher resolutions produce smoother and sharper results at increased computational cost.  
  - SF # of angles: number of angles used to compute Aspect in the hysteresis‑thresholding procedure.  
  - Bayesian Opt MaxIter: maximum number of iterations for Bayesian optimisation.  
  - w: width of the step filter, typically 10% of the largest image dimension; can be tuned via Bayesian optimisation.  
  - VSFW (Variability of the Step Filtering Widths): default 0.25, controlling how `w` varies with the complexity of spectral features (also tuned by Bayesian optimisation).  
  - Lineaments: run the lineament‑extraction procedure.  
  - AutoLine: automatically tune the hyperparameters (`na`, `WSFR`, `DR`, `w`, and `VSFW`).

8. Plot  

  - Plot Results: display the results.  
  - Filter: apply a filter to the plotted results.  
  - Close All: close all figures.  
  - Clear All: reset the settings.

## **Input / output formats**

CLE2D supports input data as `.csv` point datasets in XYZ format, where:

  - X column: longitudes.  
  - Y column: latitudes.  
  - Z column: values of the geoscientific image (e.g. reflectance, magnetic field intensity, etc.).
  
To use these inputs, the user must first define the minimum and maximum values of the geographic coordinate system (in decimal degrees) for both latitude and longitude, as well as the spacing (in arcseconds) for each direction. This ensures that the program correctly interprets the spatial extent and resolution of the data.

Outputs generated by CLE2D are mainly in MATLAB `.fig` format, including:

  - Extracted spectral features: obtained from the 2D CWT and spectral PCA (S‑PCA).  
  - Retrieved curvilinear lineaments: lineaments detected by the extraction procedure.

These `.fig` files allow direct visualisation of the results in MATLAB for analysis and interpretation.

## **Contact**

Bahman Abbassi  
Université du Québec en Abitibi‑Témiscamingue  
bahman.abbassi@uqat.ca

## **License**

This program is free software: you can redistribute it and/or modify it  
under the terms of the GNU General Public License as published by  
the Free Software Foundation, either version 3 of the License, or  
(at your option) any later version.

This program is distributed in the hope that it will be useful,  
but WITHOUT ANY WARRANTY; without even the implied warranty of  
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the  
GNU General Public License for more details.

You should have received a copy of the GNU General Public License  
along with this program. If not, see `<https://www.gnu.org/licenses/>`.

Lead developer: Bahman Abbassi  
Principal investigator: Li‑Zhen Cheng  
Affiliation: Institut de Recherche en Mines et en Environnement (IRME), Université du Québec en Abitibi‑Témiscamingue (UQAT)
