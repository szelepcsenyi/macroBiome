#' Calculator for Bioclimatic Indices
#'
#' @description Calculates the values of selected bioclimatic indices, for a given geographical location (latitude
#'     and elevation) and year/epoch, by using the monthly time series of temperature, precipitation and relative
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
#' @param bciOpVar 'character' vector of at least one length that indicates which of the bioclimatic indices is/are
#'     to be computed. Valid values are as follows: \cr
#'     (a) \code{'abt'} - Mean Annual Biotemperature (in °C); \cr
#'     (b) \code{'tap'} - Total Annual Precipitation (in mm); \cr
#'     (c) \code{'per'} - Potential Evapotranspiration Ratio (dimensionless); \cr
#'     (d) \code{'fai'} - Forestry Aridity Index (dimensionless); \cr
#'     (e) \code{'gdd0'} - Growing Degree-Days on 0°C base (in °C day); \cr
#'     (f) \code{'gdd5'} - Growing Degree-Days on 5°C base (in °C day);  \cr
#'     (g) \code{'bdi'} - Budyko's Dryness Index (dimensionless); \cr
#'     (h) \code{'cci'} - Condrad's Continentality Index (in per cent);  \cr
#'     (i) \code{'mat'} - Mean Annual Temperature (in °C); \cr
#'     (j) \code{'tc'} - Mean Temperature of the Coldest Month (in °C); \cr
#'     (k) \code{'tw'} - Mean Temperature of the Warmest Month (in °C);  \cr
#'     (l) \code{'tm10'} - Number of Months with Mean Temperature above 10°C (dimensionless);  \cr
#'     (m) \code{'pdry'} - Precipitation Sum of the Driest Month (in mm); \cr
#'     (n) \code{'psdry'} - Precipitation Sum of the Driest Month in the Summer Half-Year (in mm); \cr
#'     (o) \code{'pwdry'} - Precipitation Sum of the Driest Month in the Winter Half-Year (in mm); \cr
#'     (p) \code{'pswet'} - Precipitation Sum of the Wettest Month in the Summer Half-Year (in mm); \cr
#'     (q) \code{'pwwet'} - Precipitation Sum of the Wettest Month in the Winter Half-Year (in mm); \cr
#'     (r) \code{'ps'} - Precipitation Sum of the Summer Half-Year (in mm); \cr
#'     (s) \code{'pw'} - Precipitation Sum of the Winter Half-Year (in mm); \cr
#'     (t) \code{'ptc'} - Priestley–Taylor Coefficient (dimensionless).
#' @param argCkd 'logical' scalar that indicates whether or not the checking and correction of arguments can be
#'    omitted.
#'
#' @details Taking into account all implemented bioclimatic indices, the following five require only temperature data:
#'
#'     \itemize{
#'       \item{\code{abt}: Mean Annual Biotemperature (Eq 1 in Szelepcsényi et al. (2014); in °C)}
#'       \item{\code{mat}: Mean Annual Temperature (in °C)}
#'       \item{\code{tc}: Mean Temperature of the Coldest Month (in °C)}
#'       \item{\code{tw}: Mean Temperature of the Warmest Month (in °C)}
#'       \item{\code{tm10}: Number of Months with Mean Temperature above 10°C (dimensionless)}
#'       \item{\code{gdd5}: Growing Degree-Days on 5°C base (in °C day)}
#'       \item{\code{gdd0}: Growing Degree-Days on 0°C base (in °C day)}
#'     }
#'
#'     Monthly data are sufficient to calculate values of the mean temperatures of the coldest and warmest months,
#'     the mean annual temperature/biotemperature and the number of months with temperature > 10°C, while daily
#'     values are needed to compute values of the growing degree-days. If only a set of these bioclimatic indices has
#'     to be calculated, the setting \code{prec = NULL} must be used. \cr
#'     An important bioclimatic index for both the Holdridge life zone system and Köppen-Geiger climate
#'     classification system is the total annual precipitation, for the calculation of which requires only monthly
#'     precipitation data. If only this bioclimatic index has to be computed, the setting \code{temp = NULL} must be
#'     used. The same setting has to be used for calculation of the precipitation sum of the driest month. \cr
#'     In addition to monthly temperature data, latitude coordinates are also required to calculate the Condrad's
#'     Continentality Index (\code{cci}: Eq 4 in Conrad (1946); in per cent). \cr
#'     For calculating values of the Potential Evapotranspiration Ratio used in the Holdridge life zone system
#'     (\code{per}: Eq 4 in Szelepcsényi et al. (2014); dimensionless) and the Forestry Aridity Index introduced by
#'     the forestry climate classification (\code{fai}: Eq 1 in Führer et al. (2011); dimensionless), both
#'     temperature and precipitation data at a monthly timescale are also required. Same data are needed to calculate
#'     most precipitation statistics used by the Köppen-Geiger climate classification system:
#'
#'     \itemize{
#'       \item{\code{psdry}: Precipitation Sum of the Driest Month in the Summer Half-Year (in mm)}
#'       \item{\code{pwdry}: Precipitation Sum of the Driest Month in the Winter Half-Year (in mm)}
#'       \item{\code{pswet}: Precipitation Sum of the Wettest Month in the Summer Half-Year (in mm)}
#'       \item{\code{pwwet}: Precipitation Sum of the Wettest Month in the Winter Half-Year (in mm)}
#'       \item{\code{ps}: Precipitation Sum of the Summer Half-Year (in mm)}
#'       \item{\code{pw}: Precipitation Sum of the Winter Half-Year (in mm)}
#'     }
#'
#'     For these bioclimatic indices, summer (winter) half-year is defined as the warmer (cooler) six month period of
#'     AMJJAS (from April to September) and ONDJFM (from October to March). \cr
#'     The computation of the Budyko's Dryness Index (\code{bdi}, dimensionless) and the Priestley–Taylor Coefficient
#'     (\code{ptc}, dimensionless) requires a simulation of evapotranspiration at daily time step via the
#'     implementation of the SPLASH algorithm (Davis et al. 2017) (see
#'     \code{\link[macroBiome]{dlyEngWtrFluxPoints}}). In addition to one-year time series of daily temperature and
#'     precipitation data, the application of the SPLASH algorithm requires values of the relative sunshine duration
#'     at a daily timescale, latitude coordinate, altitude, year/epoch, and the so-called 'bucket size'. The Dryness
#'     Index is a ratio of annual potential evapotranspiration to precipitation (see Monserud et al. 1993). The value
#'     of the Priestley–Taylor coefficient is calculated as the ratio of actual evapotranspiration to equilibrium
#'     evapotranspiration, which represents the fraction of plant-available surface moisture (see Prentice et al.
#'     1992, Davis et al. 2017). \cr
#'     The function applies only monthly time series to compute values of each bioclimatic index, considering the
#'     field of application of the package. However, as we can see, in some cases there is a need for daily time
#'     series that are here generated by using the function \code{\link[macroBiome]{dlyWeaGenPoints}}.
#'
#' @return A matrix with one or more columns where each column contain the values of a given bioclimatic index.
#'
#' @note As with any function with a point mode, a set of basic input data is defined here. In this case, they are as
#'     follows: \code{'temp'} (one-year time series of monthly mean air temperature), \code{'prec'} (one-year time
#'     series of monthly precipitation sum), and \code{'bsdf'} (one-year time series of monthly mean relative sunshine
#'     duration). The objects \code{'temp'}, \code{'prec'} and \code{'bsdf'} must be either vectors of length 12 or
#'     12-column matrices. The first dimensions of these matrices have to be the same length. The function
#'     automatically converts vectors into single-row matrices during the error handling, and then uses these
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
#' \cite{Conrad V (1964) Usual formulas of continentality and their limits of validity. EOS, Trans Am Geophys Union
#'     27(5):663-664. \doi{10.1029/TR027i005p00663}}
#'
#' \cite{Davis TW, Prentice IC, Stocker BD, Thomas RT, Whitley RJ, Wang H, Evans BJ, Gallego-Sala AV, Sykes MT,
#'     Cramer W (2017) Simple process-led algorithms for simulating habitats (SPLASH v.1.0): robust indices of
#'     radiation, evapotranspiration and plant-available moisture. Geosci Model Dev 10(2):689–708.
#'     \doi{10.5194/gmd-10-689-2017}}
#'
#' \cite{Epstein ES (1991) On Obtaining Daily Climatological Values from Monthly Means. J Clim 4(3):365–368.
#'     \doi{10.1175/1520-0442(1991)004<0365:OODCVF>2.0.CO;2}}
#'
#' \cite{Führer E, Horváth L, Jagodics A, Machon A, Szabados I (2011) Application of a new aridity index in Hungarian
#'     forestry practice. Időjárás 115(3):205–216}
#'
#' \cite{Lüdeke MKB, Badeck FW, Otto RD, Häger C, Dönges S, Kindermann J, Würth G, Lang T, Jäkel U, Klaudius A,
#'     Ramge P, Habermehl S, Kohlmaier GH (1994) The Frankfurt Biosphere Model: A global process-oriented model of
#'     seasonal and long-term CO2 exchange between terrestrial ecosystems and the atmosphere. I. Model description
#'     and illustrative results for cold deciduous and boreal forests. Clim Res 4(2):143-166. \doi{10.3354/cr004143}}
#'
#' \cite{Monserud RA, Denissenko OV, Tchebakova NM (1993) Comparison of Siberian paleovegetation to current and
#'     future vegetation under climate change. Clim Res 3(3):143–159. \doi{10.3354/cr003143}}
#'
#' \cite{Prentice IC, Cramer W, Harrison SP, Leemans R, Monserud RA, Solomon AM (1992) A Global Biome Model Based on
#'     Plant Physiology and Dominance, Soil Properties and Climate. J Biogeogr 19(2):117–134. \doi{10.2307/2845499}}
#'
#' \cite{Szelepcsényi Z, Breuer H, Sümegi P (2014) The climate of Carpathian Region in the 20th century based on the
#'     original and modified Holdridge life zone system. Cent Eur J Geosci 6(3):293–307.
#'     \doi{10.2478/s13533-012-0189-5}}
#'
#' @examples
#' # Loading mandatory data for the Example 'Points'
#' data(inp_exPoints)
#'
#' # Calculate values of all default bioclimatic indices with default settings,
#' # at a grid cell near Szeged, Hungary (46.3N, 20.2E) (for the normal period 1981-2010)
#' with(inp_exPoints, {
#' bci1 <- cliBioCliIdxPoints(colMeans(temp), colMeans(prec))
#' bci1
#' })
#'
#' \donttest{
#' # Calculate values of all selected bioclimatic indices with default settings,
#' # at a grid cell near Szeged, Hungary (46.3N, 20.2E) (for the normal period 1981-2010)
#' with(inp_exPoints, {
#' year <- trunc(mean(seq(1981, 2010)))
#' bciOpVar <- c("gdd5", "bdi", "cci", "tc", "gdd0", "tw", "ptc")
#' bci2 <- cliBioCliIdxPoints(colMeans(temp), colMeans(prec), colMeans(bsdf), lat, elv,
#'    year = year, bciOpVar = bciOpVar)
#' bci2
#' })
#' }
#'
#' @importFrom stats complete.cases setNames
#' @importFrom strex match_arg
#'
#' @export
#'
cliBioCliIdxPoints <- function(temp, prec, bsdf = NULL, lat = NULL, elv = NULL, year = 2000, MSMC = 150.,
                               aprchTEMP = c("hip", "tsi", "const"), aprchPREC = c("tsi", "hip", "const"),
                               aprchBSDF = c("hip", "const"), dvTEMP = rep(0.7, 12), dvPREC = rep(0.7, 12),
                               bciOpVar = c("abt", "tap", "per", "fai"), argCkd = FALSE) {

  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  bciOpVar <- strex::match_arg(bciOpVar, rownames(ipVarRequirements), several_ok = TRUE)
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
  err_han <- errorHandling(temp = temp, prec = prec, bsdf = bsdf, lat = lat, elv = elv, year = year, MSMC = MSMC,
                           onlyLgthChecking = argCkd)
  list2env(Filter(Negate(is.null), err_han), envir = environment())

  # Check the availability of inputs to calculate bioclimatic indices
  slctdBciOpVar <- matrix(F, ncol = nrow(ipVarRequirements), dimnames = list(NULL, rownames(ipVarRequirements)))
  for (i in 1 : ncol(slctdBciOpVar)) {
    for (j in 1 : length(bciOpVar)) {
      if (rownames(ipVarRequirements)[i] == bciOpVar[j]) slctdBciOpVar[i] <- TRUE
    }
  }
  reqIpVar <- slctdBciOpVar %*% as.matrix(ipVarRequirements)
  cv.arg <- c("temp", "prec", "bsdf", "lat", "elv", "year", "MSMC")
  reqBasIpVar <- subset(reqIpVar, select = cv.arg)
  for (i in 1 : ncol(reqBasIpVar)) {
    if (reqBasIpVar[1, i] != 0 & is.null(get(names(reqBasIpVar[1, i])))) {
      stop("Invalid argument: '", names(reqBasIpVar[, i]), "' is missing, with no default. ",
           "That is required for calculating the following bioclimatic index/indices: ",
           paste(rownames(ipVarRequirements)[ipVarRequirements[, i] == 1], collapse = ", "), ".")
    }
  }


  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  # n_moy: the number of months in the year (12)
  n_moy <- 12

  # amjjas: the serial numbers of months in a six-month period (from April to September)
  amjjas <- seq(4, 9)

  # ondjfm: the serial numbers of months in a six-month period (from October to March)
  ondjfm <- setdiff(seq(1, 12), amjjas)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 01. Calculate values of the bioclimatic indices for which monthly time series of climate variables are required
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Mean Annual Biotemperature, deg C (abt) | ref: Eq 1 in Szelepcsényi et al. (2014)
  if (any(c("abt", "per") %in% bciOpVar)) {
    abt <-  rowSums(ifelse(temp >= 0. & temp <= 30., temp, 0.)) / n_moy
  }

  # Total Annual Precipitation, mm (tap) | ref: Eq 3 in Szelepcsényi et al. (2014)
  if (any(c("tap", "per") %in% bciOpVar)) {
    tap <- rowSums(prec)
  }

  # Potential Evapotranspiration Ratio, dimensionless (per) | ref: Eq 4 in Szelepcsényi et al. (2014)
  if ("per" %in% bciOpVar) {
    numer <- (58.93 * abt)
    denom <- tap
    per <- ifelse(!is.finite(numer / denom), NA, (numer / denom))
  }

  # Forestry Aridity Index, dimensionless (fai) | ref: Eq 1 in Führer et al. (2011)
  if ("fai" %in% bciOpVar) {
    numer <- rowMeans(subset(temp, select = seq(7, 8)))
    denom <- rowSums(subset(prec, select = seq(5, 7))) + rowSums(subset(prec, select = seq(7, 8)))
    fai <- ifelse(!is.finite(numer / denom), NA, (numer / denom) * 100.)
  }

  # Mean Temperature of the Coldest Month, deg C (tc)
  if (any(c("tc", "cci") %in% bciOpVar)) {
    tc <- do.call(pmin, as.data.frame(temp))
  }

  # Mean Temperature of the Warmest Month, deg C (tw)
  if (any(c("tw", "cci") %in% bciOpVar)) {
    tw <- do.call(pmax, as.data.frame(temp))
  }

  # Condrad's Continentality Index, per cent (cci) | ref: Eq 4 in Conrad (1946)
  if ("cci" %in% bciOpVar) {
    numer <- 1.7 * (tw  - tc)
    denom <- sin(deg2rad(abs(lat + 10.)))
    cci <- ifelse(!is.finite(numer / denom), NA, (numer / denom) - 14.)
  }

  # Mean Annual Temperature, deg C (mat)
  if ("mat" %in% bciOpVar) {
    mat <- rowSums(temp) / n_moy
  }

  # Number of Months with Mean Temperature above 10 deg C, dimensionless (tm10)
  if ("tm10" %in% bciOpVar) {
    tm10 <- rowSums(temp > 10.)
  }

  # Precipitation Sum of the Driest Month, mm (pdry)
  if ("pdry" %in% bciOpVar) {
    pdry <- do.call(pmin, as.data.frame(prec))
  }

  # Monthly Precipitation Sums in the Summer Half-Year, mm (smrPrec)
  if (any(c("psdry", "pswet", "ps") %in% bciOpVar)) {
    smr <- rowMeans(temp[, amjjas, drop = FALSE]) >= rowMeans(temp[, ondjfm, drop = FALSE])
    mskOps <- matrix(c(rep(NA, 6), amjjas, ondjfm), nrow = 3, byrow = TRUE)
    smrMsk <- mskOps[ifelse(is.na(smr), 1, ifelse(smr, 2, 3)), , drop = FALSE]
    smrPrec <- t(sapply(1 : nrow(prec), function(i) { prec[i, smrMsk[i, ]] }))
  }

  # Precipitation Sum of the Driest Month in the Summer Half-Year, mm (psdry)
  if ("psdry" %in% bciOpVar) {
    psdry <- do.call(pmin, as.data.frame(smrPrec))
  }

  # Precipitation Sum of the Wettest Month in the Summer Half-Year, mm (pswet)
  if ("pswet" %in% bciOpVar) {
    pswet <- do.call(pmax, as.data.frame(smrPrec))
  }

  # Precipitation Sum of the Summer Half-Year, mm (ps)
  if ("ps" %in% bciOpVar) {
    ps <- rowSums(smrPrec)
  }

  # Monthly Precipitation Sums in the Winter Half-Year, mm (winPrec)
  if (any(c("pwdry", "pwwet", "pw") %in% bciOpVar)) {
    win <- rowMeans(temp[, amjjas, drop = FALSE]) < rowMeans(temp[, ondjfm, drop = FALSE])
    mskOps <- matrix(c(rep(NA, 6), amjjas, ondjfm), nrow = 3, byrow = TRUE)
    winMsk <- mskOps[ifelse(is.na(smr), 1, ifelse(win, 2, 3)), , drop = FALSE]
    winPrec <- t(sapply(1 : nrow(prec), function(i) { prec[i, winMsk[i, ]] }))
  }

  # Precipitation Sum of the Driest Month in the Winter Half-Year, mm (pwdry)
  if ("pwdry" %in% bciOpVar) {
    pwdry <- do.call(pmin, as.data.frame(winPrec))
  }

  # Precipitation Sum of the Wettest Month in the Winter Half-Year, mm (pwwet)
  if ("pwwet" %in% bciOpVar) {
    pwwet <- do.call(pmax, as.data.frame(winPrec))
  }

  # Precipitation Sum of the Winter Half-Year, mm (pw)
  if ("pw" %in% bciOpVar) {
    pw <- rowSums(winPrec)
  }

  dlyBsdIpvReq <- ipVarRequirements[rowSums(ipVarRequirements[, c("TEMP", "PREC", "BSDF")]) != 0,
                                    c("TEMP", "PREC", "BSDF")]
  dlyBsdBioCli <- rownames(dlyBsdIpvReq)
  if (any(dlyBsdBioCli %in% bciOpVar)) {

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 02. Generate the daily values of climate variables
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (min(reqIpVar[1, c("TEMP", "BSDF")]) != 0) {
      dlyCli <- dlyWeaGenPoints(temp, prec, bsdf, year, aprchTEMP, aprchPREC, aprchBSDF, dvTEMP, dvPREC, argCkd = T)
      list2env(dlyCli, envir = environment())
    } else {
      if (reqIpVar[, "TEMP"] != 0) {
        if (aprchTEMP == "tsi") {
          TEMP <- dwgTmpScLinIntrpl(temp, "temp", year, dvTEMP, argCkd = T)
        } else if (aprchTEMP == "const") {
          TEMP <- dwgTmpScLinIntrpl(temp, "temp", year, rep(1., n_moy), argCkd = T)
        } else {
          TEMP <- dwgHarmonicIntrpl(temp, "temp", year, argCkd = T)
        }
      }
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 03. Set the result vectors of daily-based bioclimatic indices
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (i_var in 1 : length(dlyBsdBioCli)) {
      if (dlyBsdBioCli[i_var] %in% bciOpVar) { assign(dlyBsdBioCli[i_var], rep(NA, lgth)) }
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 04. Select the calendar(s) applied
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # n_moy: the number of months in the year (12)
    # n_doy: the number of days in the year (365 or 366, depending on the leap year)
    # IDNM: the number of days in the specified month

    n_moy <- 12
    IDNM <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    CDNM <- c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)

    MvIDNM <- MvCDNM <- matrix(nrow = 2, ncol = n_moy, dimnames = list(c("nrml", "leap"), month.abb))
    MvIDNM[c("nrml", "leap"), ] <- matrix(rep(IDNM, 2), nrow = 2, byrow = TRUE)
    MvIDNM[c("leap"), 2] <- 29
    MvCDNM[c("nrml", "leap"), ] <- matrix(rep(CDNM, 2), nrow = 2, byrow = TRUE)
    for (i_moy in 3 : n_moy) {
      MvCDNM[("leap"), i_moy] <- MvCDNM[("leap"), i_moy] + 1
    }

    n_doy <- rep(NA, length(year))
    n_doy[!is.na(year)] <- ifelse(is.leap.year(year[!is.na(year)]), 366, 365)

    N_doy <- c(365, 366)
    cv.dly <- vector(length = 0)
    if (any(c("gdd5", "gdd0") %in% bciOpVar)) { cv.dly <- c(cv.dly, "heatSum") }
    if (any(c("bdi", "ptc") %in% bciOpVar)) { cv.dly <- c(cv.dly, "dlyFlux") }
    cal <- setNames(lapply(vector("list", length(cv.dly)), function(x) vector()), cv.dly)
    for (i_dly in 1 : length(cv.dly)) {
      cal[[i_dly]] <- setNames(lapply(vector("list", 2), function(x) vector()), c("nrml", "leap"))

      if (cv.dly[i_dly] == "heatSum") {
        vld <- Reduce("&", list(do.call(complete.cases, list(year, temp))))
      } else {
        vld <- Reduce("&", list(do.call(complete.cases, list(year, temp, prec, bsdf))))
      }

      cal[[i_dly]]$nrml <- as.numeric(
        which(Reduce("&", list(vld, mapply(function(x) { x == N_doy[1] & !is.na(x) }, n_doy)))))
      cal[[i_dly]]$leap <- as.numeric(
        which(Reduce("&", list(vld, mapply(function(x) { x == N_doy[2] & !is.na(x) }, n_doy)))))
    }

    applCal <- match(sort(unique(n_doy[!is.na(n_doy)])), N_doy)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 05. Run algorithms to calculate values of each daily-based bioclimatic index for each calendar type
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Growing Degree-Days above 5 deg C, deg C day (gdd5)
    # Growing Degree-Days above 0 deg C, deg C day (gdd0)
    if (any(c("gdd5", "gdd0") %in% bciOpVar)) {

      if (any(sapply(cal$heatSum, function(x) { !identical(x, numeric(0)) }))) {

        for (i_cal in applCal) {
          slctd <- cal$heatSum[[c("nrml", "leap")[i_cal]]]
          tstp <- seq(1, N_doy[i_cal])

          slctdTEMP <- TEMP[slctd, tstp, drop = FALSE]
          if ("gdd5" %in% bciOpVar) {
            gdd5[slctd] <- rowSums(ifelse(slctdTEMP >= 5., slctdTEMP - 5., 0.), na.rm = FALSE)
          }
          if ("gdd0" %in% bciOpVar) {
            gdd0[slctd] <- rowSums(ifelse(slctdTEMP >= 0., slctdTEMP, 0.), na.rm = FALSE)
          }
        }
      }

    }

    # Dryness Index, dimensionless (bdi) | ref: 'text' in Monserud et al. (1993)
    # Priestley-Taylor Coefficient, dimensionless (ptc) | ref: 'text' in Prentice et al. (1992)
    if (any(c("bdi", "ptc") %in% bciOpVar)) {

      if (any(sapply(cal$dlyFlux, function(x) { !identical(x, numeric(0)) }))) {

        for (i_cal in applCal) {

          slctd <- cal$dlyFlux[[c("nrml", "leap")[i_cal]]]
          tstp <- seq(1, N_doy[i_cal])

          slctdTEMP <- TEMP[slctd, tstp, drop = FALSE]
          slctdPREC <- PREC[slctd, tstp, drop = FALSE]
          slctdBSDF <- BSDF[slctd, tstp, drop = FALSE]

          if ("ptc" %in% bciOpVar) {

            mlyBciVar <- dlyEngWtrFluxPoints(slctdTEMP, slctdPREC, slctdBSDF, lat[slctd], elv[slctd], year[slctd],
                                             MSMC[slctd], daily = FALSE)

            numer <- rowSums(mlyBciVar$AET_m)
            denom <- rowSums(mlyBciVar$EET_m)
            ptc[slctd] <- ifelse(!is.finite(numer / denom), NA, (numer / denom))

            if ("bdi" %in% bciOpVar) {
              numer <- rowSums(mlyBciVar$PET_m)
              denom <- rowSums(slctdPREC)
              bdi[slctd] <- ifelse(!is.finite(numer / denom), NA, (numer / denom))
            }

          } else {

            PET_mo_mm.mo1 <- matrix(0., nrow = length(slctd), ncol = n_moy, dimnames = list(NULL, month.abb))

            for (i_moy in 1 : n_moy) { # monthly

              for (i_dom in 1 : MvIDNM[c("nrml", "leap")[i_cal], i_moy]) { # daily

                i_doy <- MvCDNM[c("nrml", "leap")[i_cal], i_moy] + i_dom

                # Calculate daily radiation, condensation, and evaporation fluxes for a given day
                dlyQty <- calcDlyEvapot(year[slctd], rep(i_doy, length(slctd)), lat[slctd], elv[slctd],
                                        slctdBSDF[, i_doy], slctdTEMP[, i_doy], argCkd = TRUE)

                # Update the monthly totals of evapotranspiration
                PET_mo_mm.mo1[, i_moy] <- PET_mo_mm.mo1[, i_moy] + dlyQty$PET_mm

              } # daily

            } # monthly

            numer <- rowSums(PET_mo_mm.mo1)
            denom <- rowSums(slctdPREC)
            bdi[slctd] <- ifelse(!is.finite(numer / denom), NA, (numer / denom))

          }

        }

      }

    }

  }

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  rslt <- do.call(cbind, as.list(setNames(mget(bciOpVar), bciOpVar)))
  return(rslt)

}


#' Calculator for Bioclimatic Indices
#'
#' @description Calculates the values of selected bioclimatic indices, for a given region and year/epoch, by using the
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
#' @param bciOpVar 'character' vector of at least one length that indicates which of the bioclimatic indices is/are
#'     to be computed. Valid values are as follows: \cr
#'     (a) \code{'abt'} - Mean Annual Biotemperature (in °C); \cr
#'     (b) \code{'tap'} - Total Annual Precipitation (in mm); \cr
#'     (c) \code{'per'} - Potential Evapotranspiration Ratio (dimensionless); \cr
#'     (d) \code{'fai'} - Forestry Aridity Index (dimensionless); \cr
#'     (e) \code{'gdd0'} - Growing Degree-Days on 0°C base (in °C day); \cr
#'     (f) \code{'gdd5'} - Growing Degree-Days on 5°C base (in °C day);  \cr
#'     (g) \code{'bdi'} - Budyko's Dryness Index (dimensionless); \cr
#'     (h) \code{'cci'} - Condrad's Continentality Index (in per cent);  \cr
#'     (i) \code{'mat'} - Mean Annual Temperature (in °C); \cr
#'     (j) \code{'tc'} - Mean Temperature of the Coldest Month (in °C); \cr
#'     (k) \code{'tw'} - Mean Temperature of the Warmest Month (in °C);  \cr
#'     (l) \code{'tm10'} - Number of Months with Mean Temperature above 10°C (dimensionless);  \cr
#'     (m) \code{'pdry'} - Precipitation Sum of the Driest Month (in mm); \cr
#'     (n) \code{'psdry'} - Precipitation Sum of the Driest Month in the Summer Half-Year (in mm); \cr
#'     (o) \code{'pwdry'} - Precipitation Sum of the Driest Month in the Winter Half-Year (in mm); \cr
#'     (p) \code{'pswet'} - Precipitation Sum of the Wettest Month in the Summer Half-Year (in mm); \cr
#'     (q) \code{'pwwet'} - Precipitation Sum of the Wettest Month in the Winter Half-Year (in mm); \cr
#'     (r) \code{'ps'} - Precipitation Sum of the Summer Half-Year (in mm); \cr
#'     (s) \code{'pw'} - Precipitation Sum of the Winter Half-Year (in mm); \cr
#'     (t) \code{'ptc'} - Priestley–Taylor Coefficient (dimensionless).
#' @param filename output filename
#' @param ... additional arguments passed on to \code{\link[raster]{writeRaster}}
#'
#' @details See \code{\link[macroBiome]{cliBioCliIdxPoints}}.
#'
#' @return A RasterStack with one or more layers where each layer contain the values of a given bioclimatic index.
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
#' @examples
#' # Loading mandatory data for the Example 'Climate Normal Grid'
#' data(inp_exClnrGrid)
#'
#' # Calculate values of all default bioclimatic indices with default settings
#' # for Csongrad-Csanad County (for the normal period 1981-2010)
#' with(inp_exClnrGrid, {
#' rs.bci1 <- cliBioCliIdxGrid(temp, prec)
#' rs.bci1
#' })
#'
#' \donttest{
#' # Calculate values of all selected bioclimatic indices with default settings
#' # for Csongrad-Csanad County (for the normal period 1981-2010)
#' with(inp_exClnrGrid, {
#' year <- trunc(mean(seq(1981, 2010)))
#' bciOpVar <- c("gdd5", "bdi", "cci", "tc", "gdd0", "tw", "ptc")
#' rs.bci2 <- cliBioCliIdxGrid(temp, prec, bsdf, elv, sc.year = year, bciOpVar = bciOpVar)
#' rs.bci2
#' })
#' }
#'
#' @importFrom strex match_arg
#' @import raster
#'
#' @export
#'
cliBioCliIdxGrid <- function(rs.temp, rs.prec, rs.bsdf = NULL, rl.elv = NULL, sc.year = 2000, rl.MSMC = 150.,
                             aprchTEMP = c("hip", "tsi", "const"), aprchPREC = c("tsi", "hip", "const"),
                             aprchBSDF = c("hip", "const"), dvTEMP = rep(0.7, 12), dvPREC = rep(0.7, 12),
                             bciOpVar = c("abt", "tap", "per", "fai"), filename = "", ...) {

  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  bciOpVar <- strex::match_arg(bciOpVar, rownames(ipVarRequirements), several_ok = TRUE)
  aprchTEMP <- strex::match_arg(aprchTEMP)
  aprchPREC <- strex::match_arg(aprchPREC)
  aprchBSDF <- strex::match_arg(aprchBSDF)

  errorChecking(dvTEMP = dvTEMP, dvPREC = dvPREC)

  rs.aux <- if (is.null(rs.temp)) { rs.prec } else { rs.temp }

  # Check the availability of inputs to calculate bioclimatic indices
  slctdBciOpVar <- matrix(F, ncol = nrow(ipVarRequirements), dimnames = list(NULL, rownames(ipVarRequirements)))
  for (i in 1 : ncol(slctdBciOpVar)) {
    for (j in 1 : length(bciOpVar)) {
      if (rownames(ipVarRequirements)[i] == bciOpVar[j]) slctdBciOpVar[i] <- TRUE
    }
  }
  reqIpVar <- slctdBciOpVar %*% as.matrix(ipVarRequirements)
  cv.arg <- c("rs.temp", "rs.prec", "rs.bsdf", "rl.elv", "sc.year", "rl.MSMC")
  reqBasIpVar <- subset(reqIpVar, select = substring(cv.arg, 4))
  for (i in 1 : ncol(reqBasIpVar)) {
    if (reqBasIpVar[1, i] != 0 & is.null(get(cv.arg[which(colnames(reqBasIpVar) == substring(cv.arg[i], 4))]))) {
      stop("Invalid argument: '", cv.arg[which(colnames(reqBasIpVar) == substring(cv.arg[i], 4))],
           "' is missing, with no default. ",
           "That is required for calculating the following bioclimatic index/indices: ",
           paste(rownames(ipVarRequirements)[ipVarRequirements[, i] == 1], collapse = ", "), ".")
    }
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
    if (length(rl.MSMC) == 1L & is.numeric(rl.MSMC) & rl.MSMC > 0.) {
      rl <- raster(rs.aux, layer = 1)
      values(rl) <- rl.MSMC
      rl.MSMC <- mask(rl, raster(rs.aux, layer = 1))
    } else {
      stop("Invalid argument: 'rl.MSMC' has to be a single positive number or ",
           "a RasterLayer, RasterBrick, or a RasterStack with one layer.")
    }
  } else {
    if (nlayers(rl.MSMC) != 1) {
      stop("Invalid argument: 'rl.MSMC' has to be a single positive number or ",
           "a RasterLayer, RasterBrick, or a RasterStack with one layer.")
    }
  }

  errorCheckingGrid(rs.temp = rs.temp, rs.prec = rs.prec, rs.bsdf = rs.bsdf, rl.elv = rl.elv, rl.MSMC = rl.MSMC)


  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

  rl.lat <- getGeogrCoord(raster(rs.aux, layer = 1), "lat")

  n_lyr <- length(bciOpVar)

  rs.rslt <- brick(rs.aux, nl = n_lyr)

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
  bs <- blockSize(rs.aux)
  pb <- pbCreate(bs$n, ...)

  if (todisk) {
    for (i in 1 : bs$n) {
      cv.mly_var <- c("rs.temp", "rs.prec", "rs.bsdf")
      cv.loc_dta <- c("rl.lat", "rl.elv", "rl.MSMC")
      cv.arg <- c(cv.mly_var, cv.loc_dta)
      for (i_arg in 1 : length(cv.arg)) {
        if (!is.null(get(cv.arg[i_arg]))) {
          assign(substring(cv.arg[i_arg], 4), getValues(get(cv.arg[i_arg]), row = bs$row[i], nrows = bs$nrows[i]))
        } else {
          assign(substring(cv.arg[i_arg], 4), NULL)
        }
      }
      mx.rslt <- cliBioCliIdxPoints(temp, prec, bsdf, lat, elv, sc.year, MSMC, aprchTEMP, aprchPREC, aprchBSDF,
                                    dvTEMP, dvPREC, bciOpVar)

      rs.rslt <- writeValues(rs.rslt, mx.rslt, bs$row[i])
      pbStep(pb, i)
    }
    rs.rslt <- writeStop(rs.rslt)
  } else {
    for (i in 1 : bs$n) {
      cv.mly_var <- c("rs.temp", "rs.prec", "rs.bsdf")
      cv.loc_dta <- c("rl.lat", "rl.elv", "rl.MSMC")
      cv.arg <- c(cv.mly_var, cv.loc_dta)
      for (i_arg in 1 : length(cv.arg)) {
        if (!is.null(get(cv.arg[i_arg]))) {
          assign(substring(cv.arg[i_arg], 4), getValues(get(cv.arg[i_arg]), row = bs$row[i], nrows = bs$nrows[i]))
        } else {
          assign(substring(cv.arg[i_arg], 4), NULL)
        }
      }
      mx.rslt <- cliBioCliIdxPoints(temp, prec, bsdf, lat, elv, sc.year, MSMC, aprchTEMP, aprchPREC, aprchBSDF,
                                    dvTEMP, dvPREC, bciOpVar)

      cols <- bs$row[i] : (bs$row[i] + bs$nrows[i] - 1)
      arr[, cols, ] <- array(mx.rslt, dim = c(bs$nrows[i], ncol(rs.rslt), nlayers(rs.rslt)))
      pbStep(pb, i)
    }
    for (lyr in 1 : nlayers(rs.rslt)) {
      rs.rslt <- setValues(rs.rslt, as.vector(arr[, , lyr]), layer = lyr)
    }
  }

  pbClose(pb)

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  names(rs.rslt) <- bciOpVar
  return(stack(rs.rslt))

}
