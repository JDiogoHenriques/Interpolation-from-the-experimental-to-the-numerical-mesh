# Interpolation-from-the-experimental-to-the-numerical-mesh
MATLAB code for the interpolation from experimental DIC data points to the numerical data points developed in MATLAB

## Features
- Interpolation between different meshes.
- Multiple interpolation methods.
- Possibility of auto-alignment of both meshes.
- Includes the interpolation from the DIC results of 2 software: MatchID and Aramis.
- Allows the graphical visualisation of the results.

## Installation

1. Clone this repository.
2. Download the following dataset for a working example:
  2.1 MatchID: https://uapt33090-my.sharepoint.com/:f:/g/personal/joaodiogofh_ua_pt/EhGPX_zm7aFAjWWc4Od5E38BPFYYku_Z9zp1oUI8w4RzQQ?e=oTd7Wn
   ##
  2.2 Aramis: https://uapt33090-my.sharepoint.com/:f:/g/personal/joaodiogofh_ua_pt/EjuOupgm15tFntUHusApXBoBILyqLz8HQweHLb2e0zO_ZQ?e=b0f4rc
4. Put the folder inside the A00 and overwrite the contents of the preexisting folder.
5. Run the code "Exp2Num_MatchID.m" or "Exp2Num_Aramis.m", depending on the DIC processor used.

## Usage
1. Make sure the experimental data is in the correct data format as shown in the provided example.
2. Update specimen name and specimen-specific data in the header at the beginning of the script.
3. Run the "Exp2Num_*.m" code, it does the interpolation of the experimental data to the numerical mesh and plots the results. Additionally it performs the averaging of the boundary conditions for the input on a FE software as well as writes an "*.inp" with the displacement and strain interpolated results.
