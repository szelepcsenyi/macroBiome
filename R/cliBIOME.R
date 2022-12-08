#' Vegetation Classifier Using the BIOME Model
#'
#' @description Calculates the values of bioclimatic indices used in the BIOME model developed by Prentice et al.
#'     (1992), and designates the biome type using these values, for a given geographical location (latitude and
#'     elevation) and year/epoch, by using the monthly time series of temperature, precipitation and relative
#'     sunshine duration.
#'
#' @param temp 'numeric' R object with one-year time series of monthly mean air temperature (in °C)
#' @param prec 'numeric' R object with one-year time series of monthly precipitation sum (in mm)
#' @param bsdf 'numeric' R object with one-year time series of monthly mean relative sunshine duration (dimensionless)
#' @param lat 'numeric' vector with the latitude coordinates (in decimal degrees)
#' @param elv 'numeric' vector with the elevation values (in meters above sea level)
#' @param year 'numeric' vector with values of the year (using astronomical year numbering)
#' @param MSMC 'numeric' vector with values of the maximum soil moisture capacity (aka 'bucket size') (in mm)
#' @param aprchTEMP 'character' vector of length 1 that indicates the scheme used to generate daily values of the
#'     daily mean air temperature for a specific year. Valid values are as follows: \cr
#'     (a) \code{'hip'} -
#'     this scheme applies the mean-preserving 'harmonic' interpolation method of Epstein (1991) to the values of
#'     monthly mean air temperature in order to generate daily values; \cr
#'     (b) \code{'tsi'} -
#'     this scheme uses an iterative interpolation technique (Lüdeke et al. 1994) to time series of the monthly mean
#'     air temperature, in order to generate a synthetic time series of the selected meteorological variable at a
#'     temporal resolution that is higher than the daily scale; finally, this synthetic time series is upscaled to a
#'     daily resolution; \cr
#'     (c) \code{'const'} -
#'     this scheme is assumed that values of the daily mean air temperature are constant within each month.
#' @param aprchPREC 'character' vector of length 1 that indicates the scheme to generate daily values of the
#'     daily precipitation sum. Valid values are as follows: \cr
#'     (a) \code{'tsi'} -
#'     this scheme uses an iterative interpolation technique (Lüdeke et al. 1994) to time series of the monthly mean
#'     precipitation intensity, in order to generate a synthetic time series of the selected meteorological variable
#'     at a temporal resolution that is higher than the daily scale; finally, this synthetic time series is upscaled
#'     to a daily resolution; \cr
#'     (b) \code{'hip'} -
#'     this scheme applies the mean-preserving 'harmonic' interpolation method of Epstein (1991) to the values of
#'     monthly mean precipitation intensity in order to generate daily values; \cr
#'     (c) \code{'const'} -
#'     this scheme is assumed that values of the daily precipitation sum are constant within each month (the monthly
#'     precipitation sum is divided equally across each day of the month).
#' @param aprchBSDF 'character' vector of length 1 that indicates the scheme used to generate daily values of the
#'     daily fractional sunshine duration for a specific year. Valid values are as follows: \cr
#'     (a) \code{'hip'} -
#'     this scheme applies the mean-preserving 'harmonic' interpolation method of Epstein (1991) to the values of
#'     monthly mean relative sunshine duration in order to generate daily values; \cr
#'     (b) \code{'const'} -
#'     this scheme is assumed that values of the daily relative sunshine duration are constant within each month.
#' @param dvTEMP 'numeric' vector of length 12 with monthly values of the damping variable for the air temperature
#'     data.
#' @param dvPREC 'numeric' vector of length 12 with monthly values of the damping variable for the precipitation data.
#' @param verbose 'logical' scalar that indicates whether or not values of the bioclimatic indices used should be
#'     added to the output.
#'
#' @details To classify vegetation, the BIOME model developed by Prentice et al. (1992) uses the values of the
#'     following 5 bioclimatic indices:
#'
#'     \itemize{
#'       \item{\code{tc}: Mean Temperature of the Coldest Month (in °C)}
#'       \item{\code{tw}: Mean Temperature of the Warmest Month (in °C)}
#'       \item{\code{gdd5}: Growing Degree-Days on 5°C base (in °C day)}
#'       \item{\code{gdd0}: Growing Degree-Days on 0°C base (in °C day)}
#'       \item{\code{ptc}: Priestley–Taylor Coefficient (dimensionless)}
#'     }
#'
#'     For details about calculating bioclimatic indices, see the function
#'     \code{\link[macroBiome]{cliBioCliIdxPoints}}. The Priestley–Taylor coefficient (\code{'ptc'}, dimensionless)
#'     is exceptional because its computation requires a simulation of evapotranspiration at daily time step via the
#'     implementation of the SPLASH algorithm (Davis et al. 2017) (see
#'     \code{\link[macroBiome]{dlyEngWtrFluxPoints}}). The application of the SPLASH algorithm requires, among other
#'     things, one-year time series of the climate variables at daily scale, which are generated from average monthly
#'     values using the function \code{\link[macroBiome]{dlyWeaGenPoints}}. \cr
#'     The designation of the biome type is implemented as a two-step procedure. First, the presence of each plant
#'     functional type (PFT) is estimated under the given climatic conditions. Subsequently, the biome type is
#'     designated by combining PFTs occurring at the maximal dominance level with each other (see Table 5 in Prentice
#'     et al. (1992)). Each PFT is described by constraints of bioclimatic variables associated with their climatic
#'     tolerances and requirements (see Table 1 in Prentice et al. (1992)). In the initial version of the BIOME
#'     model, a total of 17 biome types are distinguished (see \code{\link[macroBiome]{vegClsNumCodes}}).
#'
#' @return Depending on the setting, a data frame with one or more columns where the biome types are stored in the
#'     last (character) column, while the additional columns contain the values of bioclimatic indices used. The
#'     abbreviations of biome types can be found in the data frame \code{\link[macroBiome]{vegClsNumCodes}}. If
#'     \code{verbose = FALSE}, the return object is a one-column data frame with the biome types.
#'
#' @note As with any function with a point mode, a set of basic input data is defined here. In this case, they are as
#'     follows: \code{'temp'} (one-year time series of monthly mean air temperature), \code{'prec'} (one-year time
#'     series of monthly precipitation sum), and \code{'bsdf'} (one-year time series of monthly mean relative
#'     sunshine duration). The objects \code{'temp'}, \code{'prec'} and \code{'bsdf'} must be either vectors of
#'     length 12 or 12-column matrices. The first dimensions of these matrices have to be the same length. The
#'     function automatically converts vectors into single-row matrices during the error handling, and then uses these
#'     matrices. The first dimensions of these matrices determines the number of rows in the result matrix. In the
#'     case of arguments that do not affect the course of the calculation procedure or the structure of the return
#'     object, scalar values (i.e., 'numeric' vector of length 1) may also be allowed. In this case, they are as
#'     follows: \code{'lat'} (latitude coordinates in decimal degrees), \code{'elv'} (elevation in meters above sea
#'     level), \code{'year'} (year using astronomical year numbering), and \code{'MSMC'} ('bucket size' in mm). These
#'     scalars are converted to vectors by the function during the error handling, and these vectors are applied in
#'     the further calculations. If these data are stored in vectors of length at least 2, their length must be the
#'     same size of first dimension of the matrices containing the basic data.
#'
#' @references
#'
#' \cite{Davis TW, Prentice IC, Stocker BD, Thomas RT, Whitley RJ, Wang H, Evans BJ, Gallego-Sala AV, Sykes MT,
#'     Cramer W (2017) Simple process-led algorithms for simulating habitats (SPLASH v.1.0): robust indices of
#'     radiation, evapotranspiration and plant-available moisture. Geosci Model Dev 10(2):689–708.
#'     \doi{10.5194/gmd-10-689-2017}}
#'
#' \cite{Epstein ES (1991) On Obtaining Daily Climatological Values from Monthly Means. J Clim 4(3):365–368.
#'     \doi{10.1175/1520-0442(1991)004<0365:OODCVF>2.0.CO;2}}
#'
#' \cite{Lüdeke MKB, Badeck FW, Otto RD, Häger C, Dönges S, Kindermann J, Würth G, Lang T, Jäkel U, Klaudius A,
#'     Ramge P, Habermehl S, Kohlmaier GH (1994) The Frankfurt Biosphere Model: A global process-oriented model of
#'     seasonal and long-term CO2 exchange between terrestrial ecosystems and the atmosphere. I. Model description
#'     and illustrative results for cold deciduous and boreal forests. Clim Res 4(2):143-166. \doi{10.3354/cr004143}}
#'
#' \cite{Prentice IC, Cramer W, Harrison SP, Leemans R, Monserud RA, Solomon AM (1992) A Global Biome Model Based on
#'     Plant Physiology and Dominance, Soil Properties and Climate. J Biogeogr 19(2):117–134. \doi{10.2307/2845499}}
#'
#' @examples
#' \donttest{
#' # Loading mandatory data for the Example 'Points'
#' data(inp_exPoints)
#'
#' # Designate the biome type (using the related biolcimatic indices), with default settings,
#' # at a grid cell near Szeged, Hungary (46.3N, 20.2E) (for the normal period 1981-2010)
#' with(inp_exPoints, {
#' year <- trunc(mean(seq(1981, 2010)))
#' BIOME <- cliBIOMEPoints(colMeans(temp), colMeans(prec), colMeans(bsdf), lat, elv,
#'     year = year, verbose = TRUE)
#' numCode <- which(sapply(vegClsNumCodes$Code.BIOME, identical, BIOME[, "vegCls"]))
#' cbind(BIOME[,-c(6)], vegClsNumCodes[numCode, c("Name.BIOME", "Code.BIOME")])
#' })
#' }
#'
#' @importFrom stats setNames
#' @importFrom strex match_arg
#'
#' @export
#'
cliBIOMEPoints <- function(temp, prec, bsdf, lat, elv, year = 2000, MSMC = 150., aprchTEMP = c("hip", "tsi", "const"),
                           aprchPREC = c("tsi", "hip", "const"), aprchBSDF = c("hip", "const"), dvTEMP = rep(0.7, 12),
                           dvPREC = rep(0.7, 12), verbose = FALSE) {

  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  aprchTEMP <- strex::match_arg(aprchTEMP)
  aprchPREC <- strex::match_arg(aprchPREC)
  aprchBSDF <- strex::match_arg(aprchBSDF)

  errorChecking(year = year, MSMC = MSMC, dvTEMP = dvTEMP, dvPREC = dvPREC)

  # Vectorization of scalar variables
  cv.scl <- c("lat", "elv", "year", "MSMC")
  if (any(sapply(cv.scl, function(x) { (length(get(x)) == 1 & is.numeric(get(x))) | identical(get(x), NA) }))) {
    lgth <- errorHandling(temp = temp, prec = prec, bsdf = bsdf)$lgth
    list2env(sapply(cv.scl, function(x) {
      if ((length(get(x)) == 1 & is.numeric(get(x))) | identical(get(x), NA)) {
        assign(x, rep(get(x), lgth)) } else { assign(x, get(x)) } },
      simplify = FALSE), envir = environment())
  }
  err_han <- errorHandling(temp = temp, prec = prec, bsdf = bsdf, lat = lat, elv = elv, year = year, MSMC = MSMC)
  list2env(Filter(Negate(is.null), err_han), envir = environment())

  cv.arg <- c("temp", "prec", "bsdf", "lat", "elv")
  for (i in 1 : length(cv.arg)) {
    if (is.null(get(cv.arg[i]))) { stop("Invalid argument: '", cv.arg[i], "' is missing, with no default.") }
  }

  cv.bci <- c("tc", "gdd5", "gdd0", "tw", "ptc")


  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 01. Calculate values of each bioclimatic index required to classify vegetation
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  bioCliIdx <- cliBioCliIdxPoints(temp, prec, bsdf = bsdf, lat = lat, elv = elv, year = year, MSMC = MSMC,
                                  aprchTEMP = aprchTEMP, aprchPREC = aprchPREC, aprchBSDF = aprchBSDF,
                                  dvTEMP = dvTEMP, dvPREC = dvPREC, bciOpVar = cv.bci, argCkd = T)
  list2env(setNames(split(bioCliIdx, col(bioCliIdx)), colnames(bioCliIdx)), envir = environment())

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 02. Set the result object containing vegetation class
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  vegCls <- rep(NA, length = lgth)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 03. Determine the vegetation class by using the BIOME model
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  presBds <- array(NA, dim = c(nrow(bioPFTDefinitions), length(cv.bci), lgth),
                   dimnames = list(rownames(bioPFTDefinitions), cv.bci, NULL))
  for (i_bci in 1 : length(cv.bci)) {
    bci <- get(cv.bci[i_bci])
    lwrLim <- bioPFTDefinitions[, paste0(cv.bci[i_bci], "_", "lwr")]
    uprLim <- bioPFTDefinitions[, paste0(cv.bci[i_bci], "_", "upr")]
    for (i_vcl in 1 : nrow(bioPFTDefinitions)) {
      presBds[i_vcl, i_bci, ] <- ifelse(is.na(bci), NA, ifelse(bci >= lwrLim[i_vcl] & bci < uprLim[i_vcl],
                                                               as.integer(1), as.integer(0)))
    }
  }
  presVal <- sapply(1 : lgth, function(i) { apply(presBds[, , i], 1, prod) })

  for (i_dom in 1 : max(bioPFTDefinitions$Dominance)) {
    slctd <- which(colSums(presVal[bioPFTDefinitions$Dominance == i_dom, , drop = FALSE]) != 0)
    presVal[bioPFTDefinitions$Dominance != i_dom, slctd] <- 0
  }

  for (i_vcl in 1 : nrow(bioBiomeDefinitions)) {
    psblPFTComp <- as.numeric(subset(bioBiomeDefinitions, select = -c(Vegetation.class, Numeric.code))[i_vcl, ])
    slctd <- which(apply(presVal, 2, function(x, want) identical(as.numeric(x), want), psblPFTComp))
    vegCls[slctd] <- rownames(bioBiomeDefinitions)[i_vcl]
  }

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  if (verbose) {
    rslt <- data.frame(do.call(cbind, mget(cv.bci)), vegCls = vegCls)
  } else {
    rslt <- data.frame(vegCls = vegCls)
  }
  return(rslt)

}


#' Vegetation Classifier Using the BIOME Model
#'
#' @description Calculates the values of bioclimatic indices used in the BIOME model developed by Prentice et al.
#'     (1992), and designates the biome type using these values, for a given region and year/epoch, by using the
#'     monthly time series of temperature, precipitation and relative sunshine duration, and the elevation data.
#'
#' @param rs.temp multi-layer Raster* object with one-year time series of monthly mean air temperature (in °C)
#' @param rs.prec multi-layer Raster* object with one-year time series of monthly precipitation sum (in mm)
#' @param rs.bsdf multi-layer Raster* object with one-year time series of monthly mean relative sunshine duration
#'     (dimensionless)
#' @param rl.elv single-layer Raster* object with the elevation values (in meters above sea level)
#' @param sc.year 'numeric' scalar with the value of the year (using astronomical year numbering)
#' @param rl.MSMC 'numeric' scalar or single-layer Raster* object with the value/values of the maximum soil moisture
#'     capacity (aka 'bucket size') (in mm)
#' @param aprchTEMP 'character' vector of length 1 that indicates the scheme used to generate daily values of the
#'     daily mean air temperature for a specific year. Valid values are as follows: \cr
#'     (a) \code{'hip'} -
#'     this scheme applies the mean-preserving 'harmonic' interpolation method of Epstein (1991) to the values of
#'     monthly mean air temperature in order to generate daily values; \cr
#'     (b) \code{'tsi'} -
#'     this scheme uses an iterative interpolation technique (Lüdeke et al. 1994) to time series of the monthly mean
#'     air temperature, in order to generate a synthetic time series of the selected meteorological variable at a
#'     temporal resolution that is higher than the daily scale; finally, this synthetic time series is upscaled to a
#'     daily resolution; \cr
#'     (c) \code{'const'} -
#'     this scheme is assumed that values of the daily mean air temperature are constant within each month.
#' @param aprchPREC 'character' vector of length 1 that indicates the scheme to generate daily values of the
#'     daily precipitation sum. Valid values are as follows: \cr
#'     (a) \code{'tsi'} -
#'     this scheme uses an iterative interpolation technique (Lüdeke et al. 1994) to time series of the monthly mean
#'     precipitation intensity, in order to generate a synthetic time series of the selected meteorological variable
#'     at a temporal resolution that is higher than the daily scale; finally, this synthetic time series is upscaled
#'     to a daily resolution; \cr
#'     (b) \code{'hip'} -
#'     this scheme applies the mean-preserving 'harmonic' interpolation method of Epstein (1991) to the values of
#'     monthly mean precipitation intensity in order to generate daily values; \cr
#'     (c) \code{'const'} -
#'     this scheme is assumed that values of the daily precipitation sum are constant within each month (the monthly
#'     precipitation sum is divided equally across each day of the month).
#' @param aprchBSDF 'character' vector of length 1 that indicates the scheme used to generate daily values of the
#'     daily fractional sunshine duration for a specific year. Valid values are as follows: \cr
#'     (a) \code{'hip'} -
#'     this scheme applies the mean-preserving 'harmonic' interpolation method of Epstein (1991) to the values of
#'     monthly mean relative sunshine duration in order to generate daily values; \cr
#'     (b) \code{'const'} -
#'     this scheme is assumed that values of the daily relative sunshine duration are constant within each month.
#' @param dvTEMP 'numeric' vector of length 12 with monthly values of the damping variable for the air temperature
#'     data.
#' @param dvPREC 'numeric' vector of length 12 with monthly values of the damping variable for the precipitation data.
#' @param verbose 'logical' scalar that indicates whether or not values of the bioclimatic indices used should be
#'     added to the output.
#' @param filename output filename
#' @param ... additional arguments passed on to \code{\link[raster]{writeRaster}}
#'
#' @details See \code{\link[macroBiome]{cliBIOMEPoints}}.
#'
#' @return Depending on the setting, a RasterStack with one or more layers where the numeric integers encoding the
#'     biome type are stored at the last layer, while the additional layers contain the values of bioclimatic indices
#'     used. The meaning of integers is given in the data frame \code{\link[macroBiome]{vegClsNumCodes}}. If
#'     \code{verbose = FALSE}, the return object is a single-layer RasterStack with numeric integers encoding the
#'     biome type.
#'
#' @note The objects \code{'rs.temp'}, \code{'rs.prec'} and \code{'rs.bsdf'} must be 12-layer Raster* objects, while
#'     the object \code{'rl.elv'} has to be a single-layer Raster* object. The object \code{'rl.MSMC'} must be either
#'     a single positive number (a universal bucket size) or a single-layer Raster* object (a regionally-specified
#'     bucket size). These Raster* objects must have the same bounding box, projection, and resolution. The object
#'     \code{'sc.year'} has to be a single integer number.
#'
#' @references
#'
#' \cite{Epstein ES (1991) On Obtaining Daily Climatological Values from Monthly Means. J Clim 4(3):365–368.
#'     \doi{10.1175/1520-0442(1991)004<0365:OODCVF>2.0.CO;2}}
#'
#' \cite{Lüdeke MKB, Badeck FW, Otto RD, Häger C, Dönges S, Kindermann J, Würth G, Lang T, Jäkel U, Klaudius A,
#'     Ramge P, Habermehl S, Kohlmaier GH (1994) The Frankfurt Biosphere Model: A global process-oriented model of
#'     seasonal and long-term CO2 exchange between terrestrial ecosystems and the atmosphere. I. Model description
#'     and illustrative results for cold deciduous and boreal forests. Clim Res 4(2):143-166. \doi{10.3354/cr004143}}
#'
#' \cite{Prentice IC, Cramer W, Harrison SP, Leemans R, Monserud RA, Solomon AM (1992) A Global Biome Model Based on
#'     Plant Physiology and Dominance, Soil Properties and Climate. J Biogeogr 19(2):117–134. \doi{10.2307/2845499}}
#'
#' @examples
#' \donttest{
#' # Loading mandatory data for the Example 'Climate Normal Grid'
#' data(inp_exClnrGrid)
#' inp_exClnrGrid <- lapply(inp_exClnrGrid, crop, extent(20.15, 20.25, 46.25, 46.35))
#'
#' # Designate the biome type (using the related bioclimatic indices), with default settings,
#' # at a grid cell near Szeged, Hungary (46.3N, 20.2E) (for the normal period 1981-2010)
#' with(inp_exClnrGrid, {
#' year <- trunc(mean(seq(1981, 2010)))
#' rs.BIOME <- cliBIOMEGrid(temp, prec, bsdf, elv, sc.year = year, verbose = TRUE)
#' rs.BIOME
#' })
#' }
#'
#' @importFrom strex match_arg
#' @import raster
#'
#' @export
#'
cliBIOMEGrid <- function(rs.temp, rs.prec, rs.bsdf, rl.elv, sc.year = 2000, rl.MSMC = 150.,
                         aprchTEMP = c("hip", "tsi", "const"), aprchPREC = c("tsi", "hip", "const"),
                         aprchBSDF = c("hip", "const"), dvTEMP = rep(0.7, 12), dvPREC = rep(0.7, 12),
                         verbose = FALSE, filename = "", ...) {

  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  aprchTEMP <- strex::match_arg(aprchTEMP)
  aprchPREC <- strex::match_arg(aprchPREC)
  aprchBSDF <- strex::match_arg(aprchBSDF)

  errorChecking(dvTEMP = dvTEMP, dvPREC = dvPREC)

  cv.arg <- c("rs.temp", "rs.prec", "rs.bsdf", "rl.elv")
  for (i in 1 : length(cv.arg)) {
    if (is.null(get(cv.arg[i]))) { stop("Invalid argument: '", cv.arg[i], "' is missing, with no default.") }
  }

  if (length(sc.year) == 1L & is.numeric(sc.year)) {
    if (sc.year %% 1 != 0) {
      stop("Invalid argument: 'sc.year' has to be a single integer number.")
    }
  } else {
    stop("Invalid argument: 'sc.year' has to be a single integer number.")
  }

  cv.rstr_cls <- c("RasterLayer", "RasterBrick", "RasterStack")
  if (length(which(!is.na(match(class(rl.MSMC), cv.rstr_cls)))) < 1L) {
    if (length(rl.MSMC) == 1L & is.numeric(rl.MSMC)) {
      rl <- raster(rs.temp, layer = 1)
      values(rl) <- rl.MSMC
      rl.MSMC <- mask(rl, raster(rs.temp, layer = 1))
    } else {
      stop("Invalid argument: 'rl.MSMC' has to be a single number or ",
           "a RasterLayer, RasterBrick, or a RasterStack with one layer.")
    }
  } else {
    if (nlayers(rl.MSMC) != 1) {
      stop("Invalid argument: 'rl.MSMC' has to be a single number or ",
           "a RasterLayer, RasterBrick, or a RasterStack with one layer.")
    }
  }

  errorCheckingGrid(rs.temp = rs.temp, rs.prec = rs.prec, rs.bsdf = rs.bsdf, rl.elv = rl.elv, rl.MSMC = rl.MSMC)


  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

  rl.lat <- getGeogrCoord(raster(rs.temp, layer = 1), "lat")

  n_lyr <- ifelse(verbose, 6, 1)

  rs.rslt <- brick(rs.temp, nl = n_lyr)

  small <- canProcessInMemory(rs.rslt, 3)
  filename <- trim(filename)
  if (!small & filename == '') {
    filename <- rasterTmpFile()
  }
  if (filename != '') {
    rs.rslt <- writeStart(rs.rslt, filename, overwrite = TRUE)
    todisk <- TRUE
  } else {
    arr <- array(dim = c(ncol(rs.rslt), nrow(rs.rslt), nlayers(rs.rslt)))
    todisk <- FALSE
  }
  bs <- blockSize(rs.temp)
  pb <- pbCreate(bs$n, ...)

  if (todisk) {
    for (i in 1 : bs$n) {
      cv.mly_var <- c("rs.temp", "rs.prec", "rs.bsdf")
      cv.loc_dta <- c("rl.lat", "rl.elv", "rl.MSMC")
      cv.arg <- c(cv.mly_var, cv.loc_dta)
      for (i_arg in 1 : length(cv.arg)) {
        assign(substring(cv.arg[i_arg], 4), getValues(get(cv.arg[i_arg]), row = bs$row[i], nrows = bs$nrows[i]))
      }
      df.rslt <- cliBIOMEPoints(temp, prec, bsdf, lat, elv, sc.year, MSMC, aprchTEMP, aprchPREC, aprchBSDF,
                                dvTEMP, dvPREC, verbose)
      numCode <- bioBiomeDefinitions$Numeric.code[match(df.rslt[["vegCls"]], rownames(bioBiomeDefinitions))]
      df.rslt[["vegCls"]] <- numCode

      rs.rslt <- writeValues(rs.rslt, as.matrix(df.rslt), bs$row[i])
      pbStep(pb, i)
    }
    rs.rslt <- writeStop(rs.rslt)
  } else {
    for (i in 1 : bs$n) {
      cv.mly_var <- c("rs.temp", "rs.prec", "rs.bsdf")
      cv.loc_dta <- c("rl.lat", "rl.elv", "rl.MSMC")
      cv.arg <- c(cv.mly_var, cv.loc_dta)
      for (i_arg in 1 : length(cv.arg)) {
        assign(substring(cv.arg[i_arg], 4), getValues(get(cv.arg[i_arg]), row = bs$row[i], nrows = bs$nrows[i]))
      }
      df.rslt <- cliBIOMEPoints(temp, prec, bsdf, lat, elv, sc.year, MSMC, aprchTEMP, aprchPREC, aprchBSDF,
                                dvTEMP, dvPREC, verbose)
      numCode <- bioBiomeDefinitions$Numeric.code[match(df.rslt[["vegCls"]], rownames(bioBiomeDefinitions))]
      df.rslt[["vegCls"]] <- numCode

      cols <- bs$row[i] : (bs$row[i] + bs$nrows[i] - 1)
      arr[, cols, ] <- array(as.matrix(df.rslt), dim = c(bs$nrows[i], ncol(rs.rslt), nlayers(rs.rslt)))
      pbStep(pb, i)
    }
    for (lyr in 1 : nlayers(rs.rslt)) {
      rs.rslt <- setValues(rs.rslt, as.vector(arr[, , lyr]), layer = lyr)
    }
  }

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  names(rs.rslt) <- colnames(df.rslt)
  return(stack(rs.rslt))

}
