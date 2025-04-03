# 2DMatSuite

**2DMatSuite** is a MATLAB-based toolkit designed for advanced material characterization and electrode design for 2D materials. This repository contains scripts that enable you to perform peak fitting for material characterization, point-by-point ellipsometry data inversion, and automated electrode layout generation (for both two-terminal and Hall bar device configurations) on exfoliated/CVD-grown 2D materials.

## Contents

- **Peakfit.m**  
  Performs Gaussian/Lorentzian peak fitting on spectroscopy data (e.g., Raman or PL spectra).  
  **Requirements:** MATLAB Curve Fitting Toolbox (for `lsqcurvefit` or other optimization functions).

- **pointElli.m**  
  Conducts a point-by-point inversion of ellipsometry data to extract the refractive index (n) and extinction coefficient (k) of thin films. The substrate is assumed to be 280nm SOI. There will be a later update on more sophisticated ellipsometry data modeling.
  **Requirements:** MATLAB Optimization Toolbox (using `fmincon`).

- **Two_terminal_electrodes_for_flakes.m**  
  Generates electrode designs for two-terminal devices. This script extracts coordinates from images of exfoliated 2D flakes (using pre-printed alignment markers) and creates coordinates that can be imported in custom GDS layouts for electrode fabrication.
  **Requirements:** MATLAB Image Processing Toolbox.

- **Hall_bar_electrodes_for_flakes.m**  
  Automates the design of Hall bar devices (6-electrode configuration) by using a marker-based coordinate system. It calculates the electrode geometry based on user input and predefined alignment markers, outputting coordinates ready for GDS conversion.
  **Requirements:** MATLAB Image Processing Toolbox.

## Getting Started

1. **Installation**  
   - Clone or download the repository.
   - Ensure MATLAB is installed along with the necessary toolboxes:
     - **Curve Fitting Toolbox** (for `peakfit.m`)
     - **Optimization Toolbox** (for `pointElli.m`)
     - **Image Processing Toolbox** (if using interactive image functions for electrode design)

2. **Usage**  
   - Open the MATLAB scripts from the repository.
   - Follow the on-screen instructions in each script to load your experimental data and images.
   - The output (e.g., electrode coordinates, fitted parameters) will be saved in text files or plotted for further analysis.

3. **Data Requirements**  
   - For **peakfit.m**, ensure you have a text file containing your spectroscopy data with proper column formatting.
   - For **pointElli.m**, prepare your ellipsometry data file (wavelength, Ψ, Δ) and the necessary refractive index data for the substrate.
   - For the electrode design scripts (**twoterminal.m** and **hallbar.m**), use high-resolution images of your SOI chip with pre-printed markers to allow accurate coordinate transformation.



Happy coding and good luck with your 2D materials research!
