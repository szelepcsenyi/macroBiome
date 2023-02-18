\name{macroBiome-package}

\alias{macroBiome-package}

\docType{package}

\title{
A Tool for Mapping the Distribution of the Biomes and Bioclimate
}

\description{
Procedures for simulating biomes by equilibrium vegetation models, with a special focus on paleoenvironmental applications.

Three widely used equilibrium biome models are currently implemented in the package: the Holdridge Life Zone (HLZ) system (Holdridge 1947, 1967), the Köppen-Geiger classification (KGC) system (Köppen 1936) and the BIOME model (Prentice et al. 1992). Three climatic forest-steppe models are also implemented.

In the BIOME model, the water balance module used to calculate daily radiation, evapotranspiration and soil moisture was replaced by the SPLASH v.1.0 model (Davis et al. 2017). The methodology used in this model was modified in two aspects: the 'bucket size' was made freely adjustable, and variations of the Earth's orbital parameters were taken into account by implementing the procedure proposed by Berger and Loutre (1991).

The application of process-based models, besides temperature and precipitation data, requires also a meteorological variable directly related to radiation. Sunshine duration or cloud cover can be such a variable. Unfortunately, paleoclimatic datasets developed to support paleoecological studies using correlative species distribution models do not include such variables. For example, see the CHELSA-TraCE21k v1.0 dataset (\url{https://chelsa-climate.org/chelsa-trace21k/}).

However, estimating relative sunshine duration data from commonly available meteorological variables may be a solution to overcome the lack of data. To achieve this goal, here, the procedure described by Yin (1999) is recommended with some minor modifications. Evidences of validity of the procedures implemented in the package were shown by Szelepcsényi et al. (2022).

The use of 'Grid' functions allows to generate distribution maps of biomes (and bioclimatic indices used in the given procedure) in a fast and controlled way given the appropriate raster datasets.
}

\references{

\cite{Berger A, Loutre MF (1991) Insolation values for the climate of the last 10 million years. Quat Sci Rev
    10(4):297-317. \doi{10.1016/0277-3791(91)90033-Q}}

\cite{Davis TW, Prentice IC, Stocker BD, Thomas RT, Whitley RJ, Wang H, Evans BJ, Gallego-Sala AV, Sykes MT,
    Cramer W (2017) Simple process-led algorithms for simulating habitats (SPLASH v.1.0): robust indices of
    radiation, evapotranspiration and plant-available moisture. Geosci Model Dev 10(2):689–708.
    \doi{10.5194/gmd-10-689-2017}}

\cite{Holdridge LR (1947) Determination of World Plant Formations From Simple Climatic Data. Science
    105(2727):367–368. \doi{ 10.1126/science.105.2727.367}}

\cite{Holdridge LR (1967) Life zone ecology. Tropical Science Center, San Jose, Costa Rica}

\cite{Köppen W (1936) Das geographische System der Klimate. In: Köppen W, Geiger R (eds) Handbuch der Klimatologie.        Verlag von Gebrüder Borntraeger, Berlin, Germany, pp 1–44}

\cite{Prentice IC, Cramer W, Harrison SP, Leemans R, Monserud RA, Solomon AM (1992) A Global Biome Model Based on
    Plant Physiology and Dominance, Soil Properties and Climate. J Biogeogr 19(2):117–134. \doi{10.2307/2845499}}

\cite{Szelepcsényi Z, Breuer H, Fodor N (2022) Estimating relative sunshine duration from commonly available
    meteorological variables for simulating biome distribution in the Carpathian Region. Hung Geogr Bull 71(1):3–19.
    \doi{10.15201/hungeobull.71.1.1}}

\cite{Yin X (1999) Bright Sunshine Duration in Relation to Precipitation, Air Temperature and Geographic Location.
    Theor Appl Climatol 64(1–2):61–68. \doi{10.1007/s007040050111}}

}
