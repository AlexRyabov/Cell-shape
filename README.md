# Cell-shape
Data and source code (MATLAB) for calculating volume and surface area of various geometric shapes

File ./Roselli/ALL linear dimensions.xlsx contains linear dimentsion of planton cells from LifeWatch ERIC 
the id of phytoplankton shapes are in the atals http://phytobioimaging.unisalento.it/Products/AtlasOfShapes.aspx?ID_Tipo=0

To get surface area and volume we use formulas 
CellVolumeand Surface.mw  a Maple Script for deriving formulas for cell surface area and volumes
CalculationsOfCellVolume&Surface.pdf is a pdf file of this script, in case zou do no have Maple. 

MATLAB scripts 
process_RoselliData.m  
  the main purpose to **calcualte volume and surface area** based on "..\Roselli\ALL linear dimensions.xlsx"
  it also averages linear dimensions for each species+location+orientation(Podosira stelligera (girdle view) and Podosira stelligera (girdle view) are averaged separately)
  it also corrects genus names using names from https://www.algaebase.org/
  it saves data into  "..\Roselli\ALL linear dimensions_VA.xlsx"
  
CleanSpeciesNames.m
tTable = CleanSpeciesNames(tTable, AllowedLostRows)
changes the genus name according to data from  https://www.algaebase.org/


