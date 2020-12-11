# Cell-shape
if you use this file please cite 

Alexey Ryabov, Onur Kerimoglu, Elena Litchman, Irina Olenina, Leonilde Roselli, Alberto Basset, Elena Stanca, Bernd Blasius. 
"Shape matters: cell geometry determines phytoplankton diversity"
bioRxiv 2020.02.06.937219; doi: https://doi.org/10.1101/2020.02.06.937219

MATLAB source code for calculating volume, surface area and other geometric characteristis of various geometric shapes

File ../data/CellSamples.xlsx contains samples of linear dimentsions of plankton cells with different geometry 
the Shape ID in this file correcpond to the shape ID in the shape atals http://phytobioimaging.unisalento.it/Products/AtlasOfShapes.aspx?ID_Tipo=0
and also in CalculationsOfCellVolume&Surface.pdf

To get surface area and volume we use the formulas 
CellVolumeand Surface.mw  a Maple Script for deriving formulas for cell surface area and volumes
CalculationsOfCellVolume&Surface.pdf is a pdf file of this script, in case zou do no have Maple. 

## cellgeom.m  
  the main purpose to **calcualte volume and surface area** based on "..\data\CellSamples.xlsx"
  the script can  also averages linear dimensions for each species+location
  it saves data into  "..\data\CellSamples_VA.xlsx"

### The paramters of cell geometry calculated:
* cell volume, *V*
* cell surface area, *A*
* equvivalent radius of a sphere with the same volume, 
   <a href="https://www.codecogs.com/eqnedit.php?latex=R_v&space;=&space;\left\(\frac{3}{4}&space;\frac{V}{\pi}\right\)^{1/3}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?R_v&space;=&space;\left\(\frac{3}{4}&space;\frac{V}{\pi}\right\)^{1/3}" title="R_v = \left\(\frac{3}{4} \frac{V}{\pi}\right\)^{1/3}" /></a>
* equvivalent radius of a sphere with the same area 
<a href="https://www.codecogs.com/eqnedit.php?latex=R_s&space;=&space;\sqrt{\frac{A}{4\pi}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?R_s&space;=&space;\sqrt{\frac{A}{4\pi}}" title="R_s = \sqrt{\frac{A}{4\pi}}" /></a>

* surface extension (inverse sphericity),
<a href="https://www.codecogs.com/eqnedit.php?latex=\varepsilon&space;=&space;\frac{A}{4\pi&space;R^2_v}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\varepsilon&space;=&space;\frac{A}{4\pi&space;R^2_v}" title="\varepsilon = \frac{A}{4\pi R^2_v}" /></a>

* minimal, middle and maximal linear cell dimensions  (L_min, L_mid, L_max)
* aspect ratio (r> for prolate cells and r<1 for oblate cells)
* classify the cell as prolate, oblate and compact
  * prolate: r> 3/2, maximal cell dimensions exceeds the minimal dimension by more than 50% 
  * compact: 2/3<r<3/2  (difference between linear dimensions is less than 50%)
  * oblate : r<2/3, maximal cell dimensions exceeds the minimal dimension by more than 50%

