# Cell-shape
Data and source code (MATLAB) for calculating volume and surface area of various geometric shapes

File ../data/CellSamples.xlsx contains samples of linear dimentsion of plankton cells from of different geometry 
the id of phytoplankton shapes can be found in the shape atals http://phytobioimaging.unisalento.it/Products/AtlasOfShapes.aspx?ID_Tipo=0

To get surface area and volume we use the formulas 
CellVolumeand Surface.mw  a Maple Script for deriving formulas for cell surface area and volumes
CalculationsOfCellVolume&Surface.pdf is a pdf file of this script, in case zou do no have Maple. 

cellgeom.m  
  the main purpose to **calcualte volume and surface area** based on "..\data\CellSamples.xlsx"
  the script can  also averages linear dimensions for each species+location
  it saves data into  "..\data\CellSamples_VA.xlsx"

The paramters of cell geometry calculated:
cell volume
cell surface area
surface extension (inverse sphericity)
equvivalent spherical radius, volume based (Eqv_Rad_v = (3/4 * V/pi))^(1/3))
equvivalent spherical radius, surface area based (Eqv_Rad_s = sqrt(A/(4 * pi))) ;
minimal, middle and maximal linear cell dimensions  (L_min, L_mid, L_max)
aspect ratio (r> for prolate cells and r<1 for oblate cells)
classify the cell as prolate, oblate and compact
prolate: r> 3/2, maximal cell dimensions exceeds the minimal dimension by more than 50% 
compact: 2/3<r<3/2  (difference between linear dimensions is less than 50%)
oblate : r<2/3, maximal cell dimensions exceeds the minimal dimension by more than 50%

