#' Mandatory Data for the Example 'Points'
#'
#' @name inp_exPoints
#' @aliases inp_exPoints
#' @docType data
#'
#' @description Data needed to demonstrate the working of functions with a point mode. Measured monthly time series
#'     of three basic climate variables retrieved from the CarpatClim database (Spinoni et al. 2015), at a grid cell
#'     near Szeged, Hungary (46.3° N, 20.2° E), for the period 1981-2010; supplemented with the associated
#'     geographical data. Monthly mean relative sunshine duration has been obtained as a ratio of monthly total of
#'     sunshine duration and maximum potential number of sunshine hours under clear-sky conditions. Daylength
#'     (accumulated hours of daylight) was calculated via Eq 1.6.11 in Duffie and Beckman (1991), according to the
#'     SPLASH radiation algorithm (see \code{\link{cliAvgDlySolIrrPoints}}).
#'
#' @format A list with three data frames and three numeric scalars.
#'
#'     Data frames contain one-year monthly time series for the period 1981-2010, for the following climate variables:
#'     \itemize{
#'       \item{\code{temp}: mean air temperature (in °C)}
#'       \item{\code{prec}: precipitation sum (in mm)}
#'       \item{\code{bsdf}: mean relative sunshine duration (dimensionless)}
#'     }
#'
#'     The following geographical parameters can be extracted from the numeric scalars:
#'     \itemize{
#'       \item{\code{lat}: latitude coordinate (in decimal degrees)}
#'       \item{\code{lon}: longitude coordinate (in decimal degrees)}
#'       \item{\code{elv}: elevation (in meters above sea level)}
#'     }
#'
#' @references
#'
#' \emph{Duffie JA, Beckman WA (1991) Solar Engineering of Thermal Processes. Second Edition. Wiley-Interscience,
#'     New York, NY}
#'
#' \emph{Spinoni J, Szalai S, Szentimrey T, Lakatos M, Bihari Z, Nagy A, Németh Á, Kovács T, Mihic D,
#'     Dacic M, Petrovic P, Kržič A, Hiebl J, Auer I, Milkovic J, Štepánek P, Zahradnícek P, Kilar P, Limanowka D,
#'     Pyrc R, Cheval S, Birsan M-V, Dumitrescu A, Deak G, Matei M, Antolovic I, Nejedlík P, Štastný P, Kajaba P,
#'     Bochnícek O, Galo D,  Mikulová K, Nabyvanets Y, Skrynyk O, Krakovska S, Gnatiuk N, Tolasz R, Antofie T, Vogt J
#'     (2015) Climate  of  the  Carpathian  Region  in  the  period  1961–2010:  climatologies  and  trends  of  10
#'     variables. Int J Climatol 35(7):1322-1341. \doi{10.1002/joc.4059}}
#'
#' @examples
#' data(inp_exPoints)
#' str(inp_exPoints)
#'
NULL

#' Mandatory Data for the Example 'Single-Year Grid'
#' @name inp_exSglyGrid
#' @aliases inp_exSglyGrid
#' @docType data
#'
#' @description Data needed to demonstrate the working of functions with a grid mode. Measured monthly time series of
#'     the year 2010 for three basic climate variables retrieved from the the CarpatClim database
#'     (Spinoni et al. 2015), for Csongrád-Csanád County; supplemented with a digital elevation model. Monthly mean
#'     relative sunshine duration has been obtained as a ratio of monthly total of sunshine duration and maximum
#'     potential number of sunshine hours under clear-sky conditions. Daylength (accumulated hours of daylight) was
#'     calculated via Eq 1.6.11 in Duffie and Beckman (1991), according to the SPLASH radiation model
#'     (see \code{\link{cliAvgDlySolIrrPoints}}). For all grids, the WGS84 (EPSG:4326) coordinate system was used
#'     with a horizontal resolution of 0.1°.
#'
#' @format A list with three RasterBricks and one RasterLayer.
#'
#'     Each of the RasterBricks contains a one-year monthly time series for the year 2010,
#'     for the following climate variables:
#'     \itemize{
#'       \item{\code{temp}: mean air temperature (in °C)}
#'       \item{\code{prec}: precipitation sum (in mm)}
#'       \item{\code{bsdf}: mean relative sunshine duration (dimensionless)}
#'     }
#'
#'     Elevation data can be extracted from a single RasterLayer: \code{elv}.
#'
#' @references
#'
#' \emph{Duffie JA, Beckman WA (1991) Solar Engineering of Thermal Processes. Second Edition. Wiley-Interscience,
#'     New York, NY}
#'
#' \emph{Spinoni J, Szalai S, Szentimrey T, Lakatos M, Bihari Z, Nagy A, Németh Á, Kovács T, Mihic D,
#'     Dacic M, Petrovic P, Kržič A, Hiebl J, Auer I, Milkovic J, Štepánek P, Zahradnícek P, Kilar P, Limanowka D,
#'     Pyrc R, Cheval S, Birsan M-V, Dumitrescu A, Deak G, Matei M, Antolovic I, Nejedlík P, Štastný P, Kajaba P,
#'     Bochnícek O, Galo D,  Mikulová K, Nabyvanets Y, Skrynyk O, Krakovska S, Gnatiuk N, Tolasz R, Antofie T, Vogt J
#'     (2015) Climate  of  the  Carpathian  Region  in  the  period  1961–2010:  climatologies  and  trends  of  10
#'     variables. Int J Climatol 35(7):1322-1341. \doi{10.1002/joc.4059}}
#'
#' @examples
#' data(inp_exSglyGrid)
#' str(inp_exSglyGrid)
#'
NULL

#' Mandatory Data for the Example 'Climate Normal Grid'
#' @name inp_exClnrGrid
#' @aliases inp_exClnrGrid
#' @docType data
#'
#' @description Data needed to demonstrate the working of functions with a grid mode. Average monthly time series of
#'     three measured basic climate variables for the normal period 1981-2010 retrieved from the the CarpatClim
#'     database (Spinoni et al. 2015), for Csongrád-Csanád County; supplemented with a digital elevation model.
#'     Monthly mean relative sunshine duration has been obtained as a ratio of monthly total of sunshine duration and
#'     maximum potential number of sunshine hours under clear-sky conditions. Daylength (accumulated hours of
#'     daylight) was calculated via Eq 1.6.11 in Duffie and Beckman (1991), according to the SPLASH radiation model
#'     (see \code{\link{cliAvgDlySolIrrPoints}}). For all grids, the WGS84 (EPSG:4326) coordinate system was used
#'     with a horizontal resolution of 0.1°.
#'
#' @format A list with three RasterBricks and one RasterLayer.
#'
#'     Each of the RasterBricks contains a one-year average monthly time series for the normal period 1981-2010,
#'     for the following climate variables:
#'     \itemize{
#'       \item{\code{temp}: mean air temperature (in °C)}
#'       \item{\code{prec}: precipitation sum (in mm)}
#'       \item{\code{bsdf}: mean relative sunshine duration (dimensionless)}
#'     }
#'
#'     Elevation data can be extracted from a single RasterLayer: \code{elv}.
#'
#' @references
#'
#' \emph{Duffie JA, Beckman WA (1991) Solar Engineering of Thermal Processes. Second Edition. Wiley-Interscience,
#'     New York, NY}
#'
#' \emph{Spinoni J, Szalai S, Szentimrey T, Lakatos M, Bihari Z, Nagy A, Németh Á, Kovács T, Mihic D,
#'     Dacic M, Petrovic P, Kržič A, Hiebl J, Auer I, Milkovic J, Štepánek P, Zahradnícek P, Kilar P, Limanowka D,
#'     Pyrc R, Cheval S, Birsan M-V, Dumitrescu A, Deak G, Matei M, Antolovic I, Nejedlík P, Štastný P, Kajaba P,
#'     Bochnícek O, Galo D,  Mikulová K, Nabyvanets Y, Skrynyk O, Krakovska S, Gnatiuk N, Tolasz R, Antofie T, Vogt J
#'     (2015) Climate  of  the  Carpathian  Region  in  the  period  1961–2010:  climatologies  and  trends  of  10
#'     variables. Int J Climatol 35(7):1322-1341. \doi{10.1002/joc.4059}}
#'
#' @examples
#' data(inp_exClnrGrid)
#' str(inp_exClnrGrid)
#'
NULL

#' Supplemental Data Frame for Decoding Outputs of Vegetation Classifiers
#' @name vegClsNumCodes
#' @aliases vegClsNumCodes
#' @docType data
#'
#' @description The key to the classes used by climate-based vegetation classifiers implemented here. Currently, three
#'     bioclimatic vegetation classification approaches are implemented:
#'     \itemize{
#'       \item{\code{HLZ}: a version with no altitudinal belts of the Holdridge life zone (HLZ) system (Holdridge
#'       1947, 1967), in accordance with works of Szelepcsényi et al. (2014, 2018)}
#'       \item{\code{KGC}: the Köppen-Geiger classification (KGC) system (Köppen 1936) with some modifications
#'       suggested by Peel et al. (2007)}
#'       \item{\code{BIOME}: the initial version of the BIOME model developed by Prentice et al. (1992)}
#'     }
#'
#' @format A data frame that allows for interpreting the return objects provided by climate-based vegetation
#'     classifiers implemented here. Two columns belong to each vegetation classification approach. Columns whose
#'     names begin with the string \code{'Name.'} contain the full names of the vegetation classes. While columns
#'     whose names begin with the string \code{'Code.'} summarize the abbreviations used by functions with a point
#'     mode. Row numbers of the data frame have a special role because they are the same as the numbers returned by
#'     the functions with a grid mode.
#'
#' @references
#'
#' \emph{Holdridge LR (1947) Determination of World Plant Formations From Simple Climatic Data. Science
#'     105(2727):367–368. \doi{10.1126/science.105.2727.367}}
#'
#' \emph{Holdridge LR (1967) Life zone ecology. Tropical Science Center, San Jose, Costa Rica}
#'
#' \cite{Köppen W (1936) Das geographische System der Klimate. In: Köppen W, Geiger R (eds) Handbuch der
#'     Klimatologie. Verlag von Gebrüder Borntraeger, Berlin, Germany, pp 1–44}
#'
#' \cite{Peel MC, Finlayson BL, McMahon TA (2007) Updated world map of the Köppen-Geiger climate classification.
#'     Hydrol Earth Syst Sci 11(5):1633–1644. \doi{10.5194/hess-11-1633-2007}}
#'
#' \emph{Prentice IC, Cramer W, Harrison SP, Leemans R, Monserud RA, Solomon AM (1992) A global biome model based on
#'     plant physiology and dominance, soil properties and climate. J Biogeogr 19(2):117–134. \doi{10.2307/2845499}}
#'
#' \emph{Szelepcsényi Z, Breuer H, Sümegi P (2014) The climate of Carpathian Region in the 20th century based on the
#'     original and modified Holdridge life zone system. Cent Eur J Geosci 6(3):293–307.
#'     \doi{10.2478/s13533-012-0189-5}}
#'
#' \emph{Szelepcsényi Z, Breuer H, Kis A, Pongrácz R, Sümegi P (2018) Assessment of projected climate change in the
#'     Carpathian Region using the Holdridge life zone system. Theor Appl Climatol 131(1–2):593–610.
#'     \doi{10.1007/s00704-016-1987-3}}
#'
#' @examples
#' data(vegClsNumCodes)
#' str(vegClsNumCodes)
#'
NULL
