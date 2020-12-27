# Cell-shape

MATLAB and Python source codes for calculating volume, surface area and other geometric characteristics of various geometric shapes. 

File "../data/CellSamples.xlsx" contains data samples with linear dimensions of plankton cells with different geometry. The Shape ID in this file correspond to the shape ID in the shape atlas http://phytobioimaging.unisalento.it/Products/AtlasOfShapes.aspx?ID_Tipo=0
and also in CalculationsOfCellVolume&Surface.pdf

To derive formulas for the surface area and volume we use a Maple script **CellVolumeand Surface.mw**.  This script was also saved as **CalculationsOfCellVolume&Surface.pdf** in case you do not have Maple. 

**If you use this file please cite**

Alexey Ryabov, Onur Kerimoglu, Elena Litchman, Irina Olenina, Leonilde Roselli, Alberto Basset, Elena Stanca, Bernd Blasius. 
"Shape matters: cell geometry determines phytoplankton diversity"
**bioRxiv** 2020.02.06.937219; doi: https://doi.org/10.1101/2020.02.06.937219

## Description

## Main Scripts:
### matlab/cellgeom.m
is a MATLAB script for calculating mean parameters. It takes data from "..\data\CellSamples.xlsx" and saves results into  "..\data\CellSamples_VA.xlsx". There is a flag, which allows to average linear dimensions for each species+location before calculating volume, area, etc. 
### python/main.py & calc_cell_geom.py
are the equivalent Python (3.6) scripts, but they work at the moment for the datasets from Olenina et al., and not from Roselli et al.

## Parameters calculated:
* cell volume, *V*
* cell surface area, *A*
* equivalent radius of a sphere with the same volume, 
   <a href="https://www.codecogs.com/eqnedit.php?latex=R_v&space;=&space;\left\(\frac{3}{4}&space;\frac{V}{\pi}\right\)^{1/3}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?R_v&space;=&space;\left\(\frac{3}{4}&space;\frac{V}{\pi}\right\)^{1/3}" title="R_v = \left\(\frac{3}{4} \frac{V}{\pi}\right\)^{1/3}" /></a>
* equivalent radius of a sphere with the same area 
<a href="https://www.codecogs.com/eqnedit.php?latex=R_s&space;=&space;\sqrt{\frac{A}{4\pi}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?R_s&space;=&space;\sqrt{\frac{A}{4\pi}}" title="R_s = \sqrt{\frac{A}{4\pi}}" /></a>
* surface extension (inverse sphericity),
<a href="https://www.codecogs.com/eqnedit.php?latex=\varepsilon&space;=&space;\frac{A}{4\pi&space;R^2_v}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\varepsilon&space;=&space;\frac{A}{4\pi&space;R^2_v}" title="\varepsilon = \frac{A}{4\pi R^2_v}" /></a>

* minimal, middle and maximal linear cell dimensions  (L_min, L_mid, L_max)
* aspect ratio (r> for prolate cells and r<1 for oblate cells)
* cell elongation class
  * prolate: r> 3/2, maximal cell dimensions exceeds the minimal dimension by more than 50% 
  * compact: 2/3<r<3/2  (difference between linear dimensions is less than 50%)
  * oblate : r<2/3, maximal cell dimensions exceeds the minimal dimension by more than 50%
