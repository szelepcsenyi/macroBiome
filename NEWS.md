# macroBiome 0.2.0

Released 2023-02-12

## Added

-   function `cliKoppenPoints()` designate the Köppen-Geiger classification (KGC) type and calculate the associated bioclimiatic indicies using the monthly time series of temperature and precipitation
-   function `cliKoppenGrid()` to apply the above algorithm to the appropriate raster datasets
-   in the functions `cliBioCliIdxPoints()` and `cliBioCliIdxGrid()`, the range of selectable bioclimatic indices has been expanded by those used in the KGC system

## Changed

-   the classification schemes based on the Holdridge Life Zone (HLZ) system are now vectorized

## Fixed

-   function `cliHoldridgeGrid()` ignored the class "BaSl" (Bare soil and no vegetation) marked with a value of 39
