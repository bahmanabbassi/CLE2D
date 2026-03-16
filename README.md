# CLE2D – Curvilinear Lineament Extraction 2D
**A MATLAB GUI for 2D Curvilinear Lineament Extraction**

### What is this repository for?
CLE2D is a MATLAB GUI designed for the extraction of curvilinear lineaments in 2D. It relies on advanced unsupervised source separation techniques, including Principal Component Analysis (PCA), 2D Continuous Wavelet Transforms (CWT), and Bayesian optimization, to provide a complete workflow for analysing and extracting lineaments from geophysical and geoscientific data.

**Please cite the software as:**
> Abbassi, B. (2024). CLE2D: Curvilinear Lineament Extraction: Bayesian Optimization of Principal Component Wavelet Analysis and Hysteresis Thresholding. GitHub. https://github.com/bahmanabbassi/CLE2D

---

### System Requirements and Setup
* **Hardware:** This program is designed to run on any Windows-based personal computer with at least 8 GB of RAM. Increasing the RAM allows larger images to be processed at once. Since the program operates on large matrices, disk read/write speed is also important. A fast NVMe SSD is recommended.
* **Software:** CLE2D is provided as MATLAB files (`.m` and `.mlapp`). It requires **MATLAB 2024a (version 24.1) or later** (fully optimized for R2025b).

**Required MATLAB Toolboxes:**
* Statistics and Machine Learning Toolbox
* Wavelet Toolbox
* Optimization Toolbox

**Installation:**
1. Copy the repository folder (MATLAB files and datasets) to your preferred directory.
2. In MATLAB, set the CLE2D folder as the **Current Folder**.
3. Launch the interface by typing `CLE2D` in the MATLAB Command Window.

---

### Dependencies and Acknowledgments
This project uses components of the **Yet Another Wavelet Toolbox (YAWTB)**, copyright (C) 2001–2002 by the YAWTB team. The original license headers have been preserved in the relevant `.m` files, in accordance with the YAWTB licensing conditions. [YAWTB GitHub](https://github.com/jacquesdurden/yawtb).

**CLE2D also uses:**
* A MATLAB implementation of **PCA / ICA** (including `fastICA.m`) developed by Brian Moore:
  *Brian Moore (2026). PCA and ICA Package, MATLAB Central File Exchange. Retrieved March 16, 2026.* [Link](https://www.mathworks.com/matlabcentral/fileexchange/38300-pca-and-ica-package)
* A **geological fault detection** code (hysteresis thresholding / step filtering) developed by Costas Panagiotakis:
  *Costas Panagiotakis (2026). Detection of Geological Faults, MATLAB Central File Exchange. Retrieved March 16, 2026.* [Link](https://www.mathworks.com/matlabcentral/fileexchange/64693-detection-of-geological-faults)

*Note: In addition to citing CLE2D itself, we recommend citing these packages when the corresponding modules are explicitly used.*

---

### Using the Interface
The CLE2D interface offers a complete set of tools with a layout oriented toward geoscientific data analysis.

#### Coordinates & Spacing
* **Max/Min Lat & Lon:** Set the geographic boundaries.
* **Method:** Choose the coordinate input method (rectangular coordinates). Coordinates can be read from a prepared text file (Min Lon, Max Lon, Min Lat, Max Lat).
* **Spacing:** Define the spacing value in arcseconds for longitude and latitude. (e.g., in Québec, 1 arcsecond lat ≈ 33m, 1 arcsecond lon ≈ 17m). The program automatically adjusts these ratios globally.
* **Filter:** Apply an additional smoothing filter to the input data.

#### Input Point Data & Digitised Lineaments
* **2D Interpolation:** Select the interpolation method for input points (regular or irregular). 
* **Xn / Yn:** Retrieve the number of pixels in X and Y after interpolation.
* **Targets Type:** Choose targets (digitised faults) as point data `.csv` or as images.
* **Cut-off 1 & 2:** Automatically set/adjust cut-off values to define buffer zones around lineaments.

#### Spectral Feature Extraction
1. **2D CWT:** Decomposes inputs into raw spectral features. Adjust scales (`na`), dilation, and angles. Isotropic and anisotropic mother wavelets are available.
2. **S‑PCA:** Used for spectral source separation. The user decides whether dimensionality reduction is required.
   * **WSFR (Wavelet Smoothness Filter Ratio):** Smoothing ratio for wavelet coefficients to avoid interpolation artefacts.
   * **β:** Symmetry angle for the CWT.

#### Lineament Extraction
* **Line Resolution:** Resolution of the extraction output. Higher resolutions produce smoother/sharper results at an increased computational cost.
* **SF # of angles:** Number of angles used to compute Aspect in the hysteresis-thresholding procedure.
* **Bayesian Opt MaxIter:** Maximum iterations for tuning.
* **w & VSFW:** Controls the step filter width and its variability with spectral complexity.
* **AutoLine:** Automatically tune hyperparameters (`na`, `WSFR`, `DR`, `w`, and `VSFW`).

---

### Input / Output Formats
CLE2D supports input data as `.csv` point datasets in XYZ format:
* **X column:** Longitudes
* **Y column:** Latitudes
* **Z column:** Geoscientific image values (e.g., reflectance, magnetic field intensity).

Outputs generated by CLE2D are mainly in MATLAB `.fig` format (extracted spectral features and retrieved curvilinear lineaments), allowing direct visualisation in MATLAB.

---

### Contact & License
**Lead Developer:** Bahman Abbassi (bahman.abbassi@uqat.ca)
**Principal Investigator:** Li-Zhen Cheng
**Affiliation:** Institut de Recherche en Mines et en Environnement (IRME), Université du Québec en Abitibi-Témiscamingue (UQAT)

This program is free software: you can redistribute it and/or modify it under the terms of the **GNU General Public License (v3 or later)** as published by the Free Software Foundation. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY. See the [GNU General Public License](https://www.gnu.org/licenses/) for more details.
