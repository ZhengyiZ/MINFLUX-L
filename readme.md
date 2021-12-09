<h1 align="center">
  <br>MINFLUX-L<br>
</h1>

<h4 align="center">Simultaneous estimation of single-molecular position and orientation with ultra-high precision</h4>

## Introduction
As an variant of [MINFLUX](https://doi.org/10.1126/science.aak9913), this method enables the simultaneous estimation of the position and orientation of a single molecule with ultra-high precision.

The datasets and part of the code are provided here. For some privacy reasons, the code for the Debye diffraction integral is not available for now, and an approximate expression (numCal.m) is used instead.

## Instruction
Currently, four types of Gaussian or Doughnut EBPs ([3/4](https://doi.org/10.1126/science.aak9913)/[6/7](https://doi.org/10.1038/s41467-021-21652-z) exposures) and six types of line-shaped dark spot (LDS) EBPs are supported.

Four examples are provided, showing in detail how to call the functions:
- example_of_CRB.m shows how to calculate various CRBs
- example_of_Field_Simu.m shows how to simulate the RMSE in the XY plane
- example_of_N_Simu.m shows how to simulate the dependence of the RMSE on the total number of photons *N*
- example_of_Tracking_Simu.m shows how to simulate tracking.

To view the datasets of paper, please run Fig_x_xxxxxx_Display.m.

## Environment
This code has been tested on:
- MATLAB R2020b on Windows 10
- MATLAB R2019b on Windows 7
