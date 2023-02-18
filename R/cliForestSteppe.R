#' Forest-Steppe Models
#'
#' @description Calculates the values of bioclimatic indices used in forest-steppe models with different theoretical
#'     backgrounds, and estimates the presence/absence of 'forest-steppe' ecotone, for a given geographical location
#'     (latitude and elevation) and year/epoch, by using the monthly time series of climate variables.
#'
#' @param temp 'numeric' R object with one-year time series of monthly mean air temperature (in °C)
#' @param prec 'numeric' R object with one-year time series of monthly precipitation sum (in mm)
#' @param bsdf 'numeric' R object with one-year time series of monthly mean relative sunshine duration (dimensionless)
#' @param lat 'numeric' vector with the latitude coordinates (in decimal degrees)
#' @param elv 'numeric' vector with the elevation values (in meters above sea level)
#' @param year 'numeric' vector with values of the year (using astronomical year numbering)
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
#' @param aprchBSDF 'character' vector of length 1 that indicates the scheme used to generate daily values of the
#'     daily fractional sunshine duration for a specific year. Valid values are as follows: \cr
#'     (a) \code{'hip'} -
#'     this scheme applies the mean-preserving 'harmonic' interpolation method of Epstein (1991) to the values of
#'     monthly mean relative sunshine duration in order to generate daily values; \cr
#'     (b) \code{'const'} -
#'     this scheme is assumed that values of the daily relative sunshine duration are constant within each month.
#' @param dvTEMP 'numeric' vector of length 12 with monthly values of the damping variable for the air temperature
#'     data.
#' @param verbose 'logical' scalar that indicates whether or not values of the bioclimatic indices used should be
#'     added to the output.
#'
#' @details Here, three forest-steppe models with different theoretical backgrounds are implemented:
#'
#'     \itemize{
#'       \item{\code{fsp_hlz}: A modified variant of the widely used Holdridge life zone (HLZ) system (see for the
#'       basic concept Holdridge 1947, 1967; for a proposed variant Szelepcsényi et al. 2014).}
#'       \item{\code{fsp_fai}: A clarified version of the forestry climate classification (see for the basic concept
#'       Führer et al. 2011; for a proposed variant Mátyás et al. 2018)}
#'       \item{\code{fsp_svm}: An initial version of the Siberian Vegetation Model (see Monserud et al. 1993)}
#'     }
#'
#'     The HLZ system classifies the vegetation type based on the distance from the ideal (theoretical) point in the
#'     3-dimensional space of the following bioclimatic indices:
#'
#'     \itemize{
#'       \item{\code{abt}: Mean Annual Biotemperature (Eq 1 in Szelepcsényi et al. (2014); in °C)}
#'       \item{\code{tap}: Total Annual Precipitation (in mm)}
#'       \item{\code{per}: Potential Evapotranspiration Ratio (Eq 4 in Szelepcsényi et al. (2014); dimensionless)}
#'     }
#'
#'     The plotting of thresholds of the above-mentioned bioclimatic indices in the HLZ chart leads to emerge a set
#'     of hexagons and triangles. The hexagons indicate the so-called core HLZ types, while the so-called
#'     transitional HLZ types are circumscribed by equilateral triangles in the HLZ chart (see Szelepcsényi et al.
#'     2014). However, in contrast to this study, here, the transitional types are defined as separate zones
#'     designated by the centres of the triangles. As a result, hexagons appear around the triangles in the HLZ
#'     chart, and in parallel, the size of the hexagons denoting the core types also decreases. Thus, the size of the
#'     core and transitional types are the same in this approach. During the classification, all forest-steppe types
#'     designated by Szelepcsényi et al. (2014) (and redefined by us) are aggregated into one class. \cr
#'     The forestry climate classification developed by Führer et al. (2011) was reworked by Mátyás et al. (2018). In
#'     the context of assessing the effects of future climate change, the 'forest-steppe' climate class was
#'     introduced in the model. In the work of Mátyás et al. (2018), this type is characterized by the Forestry
#'     Aridity Index (\code{fai}, dimensionless) values between 7.25 and 8. This definition is used here. \cr
#'     The Siberian Vegetation Model (Monserud et al. 1993) defines numerous types of forest-steppe on the basis of
#'     values of the Growing Degree-Days above 5°C (\code{gdd5}, in °C day), the Budyko's Dryness Index
#'     (\code{bdi}, dimensionless), and the Condrad's Continentality Index (\code{cci}, in per cent). Here, all
#'     such ecotone types are aggregated into one class, in order to estimate the presence/absence of the
#'     ‘forest-steppe’ ecotone.
#'
#' @return Depending on the setting, a data frame with three or more columns where the presence/absence data are
#'     stored in the last three columns labelled \code{'fsp_hlz'}, \code{'fsp_fai'} and \code{'fsp_svm'}, while the
#'     additional columns contain the values of bioclimatic indices used. If \code{verbose = FALSE}, the return
#'     object is a two- or three-column data frame with the presence/absence data, depending on the available data.
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
#'     level), and \code{'year'} (year using astronomical year numbering). These scalars are converted to vectors by
#'     the function during the error handling, and these vectors are applied in the further calculations. If these
#'     data are stored in vectors of length at least 2, their length must be the same size of first dimension of the
#'     matrices containing the basic data.
#'
#' @references
#'
#' \cite{Epstein ES (1991) On Obtaining Daily Climatological Values from Monthly Means. J Clim 4(3):365–368.
#'     \doi{10.1175/1520-0442(1991)004<0365:OODCVF>2.0.CO;2}}
#'
#' \cite{Führer E, Horváth L, Jagodics A, Machon A, Szabados I (2011) Application of a new aridity index in Hungarian
#'     forestry practice. Időjárás 115(3):205–216}
#'
#' \cite{Holdridge LR (1947) Determination of World Plant Formations From Simple Climatic Data. Science
#'     105(2727):367–368. \doi{10.1126/science.105.2727.367}}
#'
#' \cite{Holdridge LR (1967) Life zone ecology. Tropical Science Center, San Jose, Costa Rica}
#'
#' \cite{Lüdeke MKB, Badeck FW, Otto RD, Häger C, Dönges S, Kindermann J, Würth G, Lang T, Jäkel U, Klaudius A,
#'     Ramge P, Habermehl S, Kohlmaier GH (1994) The Frankfurt Biosphere Model: A global process-oriented model of
#'     seasonal and long-term CO2 exchange between terrestrial ecosystems and the atmosphere. I. Model description
#'     and illustrative results for cold deciduous and boreal forests. Clim Res 4(2):143-166. \doi{10.3354/cr004143}}
#'
#' \cite{Mátyás Cs, Berki I, Bidló A, Csóka Gy, Czimber K, Führer E, Gálos B, Gribovszki Z, Illés G, Hirka A, Somogyi
#'     Z (2018) Sustainability of Forest Cover under Climate Change on the Temperate-Continental Xeric Limits.
#'     Forests 9(8):489. \doi{10.3390/f9080489}}
#'
#' \cite{Monserud RA, Denissenko OV, Tchebakova NM (1993) Comparison of Siberian paleovegetation to current and
#'     future vegetation under climate change. Clim Res 3(3):143–159. \doi{10.3354/cr003143}}
#'
#' \cite{Szelepcsényi Z, Breuer H, Sümegi P (2014) The climate of Carpathian Region in the 20th century based on the
#'     original and modified Holdridge life zone system. Cent Eur J Geosci 6(3):293–307.
#'     \doi{10.2478/s13533-012-0189-5}}
#'
#' @examples
#' # Loading mandatory data for the Example 'Points'
#' data(inp_exPoints)
#'
#' # Predict the 'forest-steppe' ecotone (using the related bioclimatic indices),
#' # with default settings, at a grid cell near Szeged, Hungary (46.3N, 20.2E)
#' # (for the normal period 1981-2010)
#' with(inp_exPoints, {
#' year <- trunc(mean(seq(1981, 2010)))
#' fsp <- cliForestSteppePoints(colMeans(temp), colMeans(prec), colMeans(bsdf), lat, elv,
#'     year = year, verbose = TRUE)
#' fsp
#' })
#'
#' @importFrom stats complete.cases setNames
#' @importFrom strex match_arg
#'
#' @export
#'
cliForestSteppePoints <- function(temp, prec, bsdf = NULL, lat = NULL, elv = NULL, year = 2000,
                                  aprchTEMP = c("hip", "tsi", "const"), aprchBSDF = c("hip", "const"),
                                  dvTEMP = rep(0.7, 12), verbose = FALSE) {

  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  aprchTEMP <- strex::match_arg(aprchTEMP)
  aprchBSDF <- strex::match_arg(aprchBSDF)

  errorChecking(year = year, dvTEMP = dvTEMP)

  # Vectorization of scalar variables
  cv.scl <- c("lat", "elv", "year")
  if (any(sapply(cv.scl, function(x) { (length(get(x)) == 1 & is.numeric(get(x))) }))) {
    lgth <- errorHandling(temp = temp, prec = prec, bsdf = bsdf)$lgth
    list2env(sapply(cv.scl, function(x) {
      if (length(get(x)) == 1 & is.numeric(get(x))) { assign(x, rep(get(x), lgth)) } else { assign(x, get(x)) } },
      simplify = FALSE), envir = environment())
  }
  err_han <- errorHandling(temp = temp, prec = prec, bsdf = bsdf, lat = lat, elv = elv, year = year)
  list2env(Filter(Negate(is.null), err_han), envir = environment())

  cv.arg <- c("temp", "prec")
  for (i in 1 : length(cv.arg)) {
    if (is.null(get(cv.arg[i]))) { stop("Invalid argument: '", cv.arg[i], "' is missing, with no default.") }
  }

  # Create a list of bioclimatic indices that will be used in vegetation classification
  # The availability of input variables, avblIpVar (cls: mtx)
  reqdIpVar <- setdiff(colnames(ipVarRequirements), c("MSMC", "PREC"))
  avblIpVar <- matrix(NA, ncol = length(reqdIpVar), dimnames = list(NULL, reqdIpVar))
  for (i in 1 : 3) {
    if (!is.null(get(c(cv.arg, "bsdf")[i]))) {
      assign(toupper(c(cv.arg, "bsdf")[i]), NA)
    } else {
      assign(toupper(c(cv.arg, "bsdf")[i]), NULL)
    }
  }
  for (i in 1 : length(avblIpVar)) { avblIpVar[i] <- !is.null(get(reqdIpVar[i])) }

  # The number of missing input variables, numMisIpVar (cls: mtx)
  ipVarReq <- subset(ipVarRequirements, select = -c(MSMC, PREC))
  numMisIpVar <- rowSums(ipVarReq) - avblIpVar %*% t(ipVarReq)

  intMtx <- bciRequirements

  # The eligibility of the vegetation classification schemes, elgVegCls (cls: mtx)
  elgVegCls <- matrix(c(T, T, T, F, F), ncol = ncol(intMtx), dimnames = list(NULL, colnames(intMtx)))

  # The availability of the vegetation classification schemes, avblVegCls (cls: mtx)
  avblVegCls <- colSums(intMtx) - (numMisIpVar[1, ] == 0) %*% as.matrix(intMtx) == 0

  # A list of selected vegetation classification schemes, cv.svc (cls: cv)
  cv.svc <- as.character(apply(elgVegCls * avblVegCls, 1, function(i) paste(names(i[i == 1]))))

  # A list of bioclimatic indices to be calculated, cv.bci (cls: mtx)
  cv.bci <- rownames(intMtx)[rowSums(intMtx[cv.svc]) != 0]


  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 01. Calculate values of each bioclimatic index required to classify vegetation
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  bioCliIdx <- cliBioCliIdxPoints(temp, prec, bsdf = bsdf, lat = lat, elv = elv, year = year,
                                  aprchTEMP = aprchTEMP, aprchBSDF = aprchBSDF, dvTEMP = dvTEMP, bciOpVar = cv.bci)
  list2env(unclass(as.data.frame(bioCliIdx)), envir = environment())

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 02. Set the result object containing presence/absence data of 'forest-steppe'
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  presAbseDt <- matrix(nrow = lgth, ncol = length(cv.svc), dimnames = list(NULL, paste0("fsp_", cv.svc)))

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 03. Use each classification scheme to the values of available bioclimatic indices
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # A. Transitional HLZ type
  if (any("hlz" %in% cv.svc)) {
    cv.bci <- rownames(intMtx)[rowSums(intMtx["hlz"]) != 0]
    distBds <- matrix(nrow = lgth, ncol = nrow(fspChlzDefSubset), dimnames = list(NULL, rownames(fspChlzDefSubset)))
    for (i_vcl in 1 : nrow(fspChlzDefSubset)) {
      optVabt <- log2(fspChlzDefSubset$abt[i_vcl]) + 0.5
      tmpVabt <- ifelse(abt == 0., NA, (log2(abt) - optVabt) ** 2.)
      optVtap <- log2(fspChlzDefSubset$tap[i_vcl]) + 0.5
      tmpVtap <- ifelse(tap == 0., NA, (log2(tap) - optVtap) ** 2.)
      optVper <- log2(fspChlzDefSubset$per[i_vcl]) + 0.5
      tmpVper <- ifelse(per == 0., NA, (log2(per) - optVper) ** 2.)
      distBds[, i_vcl] <- sqrt(tmpVabt + tmpVtap + tmpVper)
    }
    psblVegCls <- apply(distBds, 1, function(x) ifelse(any(is.na(x)), NA, names(which.min(x))))
    fsp_hlz <- ifelse(is.na(psblVegCls), NA, ifelse(grepl("Fs", psblVegCls), as.integer(1), as.integer(0)))
  }

  # B. Forestry climate class
  if (any("fai" %in% cv.svc)) {
    fsp_fai <- ifelse(is.na(fai), NA, ifelse(fai >= 7.25 & fai < 8.5, as.integer(1), as.integer(0)))
  }

  # C. Biome type in the Siberian Vegetation Model
  if (any("svm" %in% cv.svc)) {
    cv.bci <- rownames(intMtx)[rowSums(intMtx["svm"]) != 0]
    presBds <- array(NA, dim = c(nrow(svmDefinitions), length(cv.bci), lgth),
                     dimnames = list(rownames(svmDefinitions), cv.bci, NULL))
    for (i_bci in 1 : length(cv.bci)) {
      bci <- get(cv.bci[i_bci])
      lwrLim <- svmDefinitions[, paste0(cv.bci[i_bci], "_", "lwr")]
      uprLim <- svmDefinitions[, paste0(cv.bci[i_bci], "_", "upr")]
      for (i_vcl in 1 : nrow(svmDefinitions)) {
        presBds[i_vcl, i_bci, ] <- ifelse(is.na(bci), NA, ifelse(bci >= lwrLim[i_vcl] & bci < uprLim[i_vcl],
                                                                 as.integer(1), as.integer(0)))
      }
    }
    vld <- complete.cases(list(gdd5, bdi, cci))
    presVal <- fsp_svm <- rep(NA, length = lgth)
    presVal[vld] <- sapply(seq(1, lgth)[vld], function(i) { as.numeric(which(apply(presBds[, , i], 1, prod) == 1)) })
    fsp_svm[vld] <- ifelse(is.na(presVal[vld]), NA,
                           ifelse(svmDefinitions$Biome.type[presVal[vld]] == "Forest-Steppe",
                                  as.integer(1), as.integer(0)))
  }

  presAbseDt <- do.call(cbind, as.list(setNames(mget(paste0("fsp_", cv.svc)), paste0("fsp_", cv.svc))))

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  if (verbose) {
    rslt <- cbind(bioCliIdx, presAbseDt)
  } else {
    rslt <- cbind(presAbseDt)
  }

  return(rslt)

}


#' Forest-Steppe Models
#'
#' @description Calculates the values of bioclimatic indices used in forest-steppe models with different theoretical
#'     backgrounds, and estimates the presence/absence of 'forest-steppe' ecotone, for a given region and year/epoch,
#'     by using the monthly time series of climate variables, and the elevation data.
#'
#' @param rs.temp multi-layer Raster* object with one-year time series of monthly mean air temperature (in °C)
#' @param rs.prec multi-layer Raster* object with one-year time series of monthly precipitation sum (in mm)
#' @param rs.bsdf multi-layer Raster* object with one-year time series of monthly mean relative sunshine duration
#'     (dimensionless)
#' @param rl.elv single-layer Raster* object with the elevation values (in meters above sea level)
#' @param sc.year 'numeric' scalar with the value of the year (using astronomical year numbering)
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
#' @param aprchBSDF 'character' vector of length 1 that indicates the scheme used to generate daily values of the
#'     daily fractional sunshine duration for a specific year. Valid values are as follows: \cr
#'     (a) \code{'hip'} -
#'     this scheme applies the mean-preserving 'harmonic' interpolation method of Epstein (1991) to the values of
#'     monthly mean relative sunshine duration in order to generate daily values; \cr
#'     (b) \code{'const'} -
#'     this scheme is assumed that values of the daily relative sunshine duration are constant within each month.
#' @param dvTEMP 'numeric' vector of length 12 with monthly values of the damping variable for the air temperature
#'     data.
#' @param verbose 'logical' scalar that indicates whether or not values of the bioclimatic indices used should be
#'     added to the output.
#' @param filename output filename
#' @param ... additional arguments passed on to \code{\link[raster]{writeRaster}}
#'
#' @details See \code{\link[macroBiome]{cliForestSteppePoints}}.
#'
#' @return Depending on the settings, a RasterStack with two or more layers where the presence/absence data are
#'     stored in layers labelled \code{'fsp_hlz'}, \code{'fsp_fai'} and \code{'fsp_svm'}, while the additional layers
#'     contain the values of bioclimatic indices used. If \code{verbose = FALSE}, the return object is a two- or
#'     three-layer RasterStack with presence/absence data, depending on the available data.
#'
#' @note The objects \code{'rs.temp'}, \code{'rs.prec'} and \code{'rs.bsdf'} must be 12-layer Raster* objects, while
#'     the object \code{'rl.elv'} has to be a single-layer Raster* object. These Raster* objects must have the same
#'     bounding box, projection, and resolution. The object \code{'sc.year'} has to be a single integer number.
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
#' @examples
#' # Loading mandatory data for the Example 'Climate Normal Grid'
#' data(inp_exClnrGrid)
#'
#' # Predict the 'forest-steppe' ecotone (using the related bioclimatic indices),
#' # with default settings, for Csongrad-Csanad County (for the normal period 1981-2010)
#' with(inp_exClnrGrid, {
#' year <- trunc(mean(seq(1981, 2010)))
#' rs.fsp1 <- cliForestSteppeGrid(temp, prec, verbose = TRUE)
#' rs.fsp1
#' })
#'
#' \donttest{
#' # Predict the 'forest-steppe' ecotone (using the related bioclimatic indices),
#' # with default settings, for Csongrad-Csanad County (for the normal period 1981-2010)
#' with(inp_exClnrGrid, {
#' year <- trunc(mean(seq(1981, 2010)))
#' rs.fsp2 <- cliForestSteppeGrid(temp, prec, bsdf, elv, sc.year = year, verbose = TRUE)
#' rs.fsp2
#' })
#' }
#'
#' @importFrom strex match_arg
#' @import raster
#'
#' @export
#'
cliForestSteppeGrid <- function(rs.temp, rs.prec, rs.bsdf = NULL, rl.elv = NULL, sc.year = 2000,
                                aprchTEMP = c("hip", "tsi", "const"), aprchBSDF = c("hip", "const"),
                                dvTEMP = rep(0.7, 12), verbose = FALSE, filename = "", ...) {

  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  aprchTEMP <- strex::match_arg(aprchTEMP)
  aprchBSDF <- strex::match_arg(aprchBSDF)

  errorChecking(dvTEMP = dvTEMP)

  cv.arg <- c("rs.temp", "rs.prec")
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

  errorCheckingGrid(rs.temp = rs.temp, rs.prec = rs.prec, rs.bsdf = rs.bsdf, rl.elv = rl.elv)


  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

  rl.lat <- getGeogrCoord(raster(rs.temp, layer = 1), "lat")

  if (verbose) {
    n_lyr <- ifelse(is.null(rs.bsdf), 6, 10)
  } else {
    n_lyr <- ifelse(is.null(rs.bsdf), 2, 3)
  }

  rs.rslt <- stack(brick(rs.temp, nl = n_lyr))

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
      cv.loc_dta <- c("rl.lat", "rl.elv")
      cv.arg <- c(cv.mly_var, cv.loc_dta)
      for (i_arg in 1 : length(cv.arg)) {
        if (!is.null(get(cv.arg[i_arg]))) {
          assign(substring(cv.arg[i_arg], 4), getValues(get(cv.arg[i_arg]), row = bs$row[i], nrows = bs$nrows[i]))
        } else {
          assign(substring(cv.arg[i_arg], 4), NULL)
        }
      }
      mx.rslt <- cliForestSteppePoints(temp, prec, bsdf, lat, elv, sc.year, aprchTEMP, aprchBSDF, dvTEMP, verbose)

      rs.rslt <- writeValues(rs.rslt, mx.rslt, bs$row[i])
      pbStep(pb, i)
    }
    rs.rslt <- writeStop(rs.rslt)
  } else {
    for (i in 1 : bs$n) {
      cv.mly_var <- c("rs.temp", "rs.prec", "rs.bsdf")
      cv.loc_dta <- c("rl.lat", "rl.elv")
      cv.arg <- c(cv.mly_var, cv.loc_dta)
      for (i_arg in 1 : length(cv.arg)) {
        if (!is.null(get(cv.arg[i_arg]))) {
          assign(substring(cv.arg[i_arg], 4), getValues(get(cv.arg[i_arg]), row = bs$row[i], nrows = bs$nrows[i]))
        } else {
          assign(substring(cv.arg[i_arg], 4), NULL)
        }
      }
      mx.rslt <- cliForestSteppePoints(temp, prec, bsdf, lat, elv, sc.year, aprchTEMP, aprchBSDF, dvTEMP, verbose)

      cols <- bs$row[i] : (bs$row[i] + bs$nrows[i] - 1)
      arr[, cols, ] <- array(mx.rslt, dim = c(bs$nrows[i], ncol(rs.rslt), nlayers(rs.rslt)))
      pbStep(pb, i)
    }
    for (lyr in 1 : nlayers(rs.rslt)) {
      rs.rslt <- setValues(rs.rslt, as.vector(arr[, , lyr]), layer = lyr)
    }
  }


  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  names(rs.rslt) <- colnames(mx.rslt)
  return(rs.rslt)

}
