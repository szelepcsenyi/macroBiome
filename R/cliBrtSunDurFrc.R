#' Estimator for Fraction of Bright Sunshine Duration
#'
#' @description Estimates monthly averages for daily fraction of bright sunshine duration, for a given geographical
#'     location (latitude, longitude, and elevation) and year, by using the monthly time series of temperature and
#'     precipitation.
#'
#' @param temp 'numeric' R object with one-year time series of monthly mean air temperature (in °C)
#' @param prec 'numeric' R object with one-year time series of monthly precipitation sum (in mm)
#' @param lat 'numeric' vector with the latitude coordinates (in decimal degrees)
#' @param lon 'numeric' vector with the longitude coordinates (in decimal degrees)
#' @param elv 'numeric' vector with the elevation values (in meters above sea level)
#' @param year 'numeric' vector with values of the year (using astronomical year numbering)
#' @param aprchSIM 'character' vector of length 1 that indicates the formula used to estimate the value of solar
#'     irradiance/irradiation for a specific day. Valid values are as follows: \cr
#'     (a) \code{'Solar123'} -
#'     in this approach, first, the mean hourly solar irradiance under cloudless-sky conditions is calculated as
#'     proposed by Yin (1997b), with a minor modification, using the daytime means of optical air mass and
#'     cosine zenith; the former is computed as recommended by Yin (1997b), while the latter is estimated by using
#'     Eq 5 of Yin (1997a); however, in contrast to the original approach, where the solar constant was fixed at
#'     \eqn{4.9212 MJ m^{-2} hr^{-1}}{4.9212 MJ m^{-2} hr^{-1}}, according to Yin (1999), its value is corrected by
#'     calendar day for the variable ellipticity of the Earth's orbit, by using the scheme of Brock (1981); in the
#'     calculations, the values of solar declination and daylength are derived by using the approach of Brock (1981);
#'     \cr
#'     (b) \code{'SPLASH'} -
#'     in this approach, first, under varying orbital parameters, the daily solar radiation at the top of the
#'     atmosphere is calculated (\eqn{H_{0}}{H_{0}}, Eq 7 in Davis et al. (2017)), and then this value is
#'     multiplied by the atmospheric transmittivity to obtain the value of daily surface radiation; in this case as
#'     well, cloudless conditions are assumed, i.e., the transmission coefficient is taken into account with an
#'     universal value of 0.75, however, its value is modified as a function of elevation, by using the scheme of
#'     Allen (1996); the daylength is calculated via Eq 1.6.11 in Duffie and Beckman (1991), using the sunset hour
#'     angle (\eqn{h_{s}}{h_{s}}, Eq 8. in Davis et al. (2017)); finally, the mean hourly surface radiation is
#'     derived as the quotient of the daily surface radiation and the daylength.
#'
#' @details To estimate the monthly averages of relative sunlight duration, the approach presented by Yin (1999) is
#'     implemented here. Many variables in this estimation scheme can be easily and unambiguously determined, but the
#'     approach uses two important quantities, the calculation method of which can be chosen here depending on the
#'     purpose of the investigations. One of them is the estimated value of the mean hourly solar irradiance under
#'     cloudless-sky conditions. This quantity can be estimated in this implementation of the approach with the
#'     original method (\code{aprchSIM = 'Solar123'}) or with the solar radiation model used in the SPLASH algorithm,
#'     considering the variability of orbital parameters of the Earth over time (\code{aprchSIM = 'SPLASH'}). The
#'     latter is recommended for paleo-climatological and paleo-environmental studies. These solar radiation models
#'     is also applied to calculate the daylength, whose monthly averages are used to estimate monthly averages of
#'     daily potential evapotranspiration (Eqs. A10 and A11 in Yin (1998)). \cr
#'     The procedure proposed by Yin (1999) requires the calculation of several regional factors (see Eq 3.3 in Yin
#'     (1999)). Each regional factor is activated as a function of latitude and longitude. However, it is important
#'     to note that in this implementation, these factors are activated with the current configuration of continents
#'     and islands. Continents and regions are classified using the high-resolution world map of the
#'     \code{\link[rworldmap]{rworldmap-package}}. In checking whether or not a given geographic location can be
#'     defined as an island, the high-resolution world map of the \code{\link[rnaturalearth]{rnaturalearth}}
#'     is applied.
#'
#' @return A 12-column matrix with monthly averages of relative sunshine duration.
#'
#' @note As with any function with a point mode, a set of basic input data is defined here. In this case, they are as
#'    follows: \code{'temp'} (one-year time series of monthly mean air temperature), and \code{'prec'} (one-year time
#'    series of monthly precipitation sum). The objects \code{'temp'} and \code{'prec'} must be either 12-length
#'    vectors or 12-column matrices. The first dimensions of these matrices have to be the same length. The function
#'    automatically converts vectors into single-row matrices during the error handling, and then uses these
#'    matrices. The first dimensions of these matrices determines the number of rows in the result matrix. In the
#'    case of arguments that do not affect the course of the calculation procedure or the structure of the return
#'    object, scalar values (i.e., 'numeric' vector of length 1) may also be allowed. In this case, they are as
#'    follows: \code{'lat'} (latitude coordinates in decimal degrees), \code{'lon'} (longitude coordinates in decimal
#'    degrees), \code{'elv'} (elevation in meters above sea level), and \code{'year'} (year using astronomical year
#'    numbering). These scalars are converted to vectors by the function during the error handling, and these vectors
#'    are applied in the further calculations. If these data are stored in vectors of length at least 2, their length
#'    must be the same size of first dimension of the matrices containing the basic data.
#'
#' @references
#'
#' \cite{Allen RG (1996) Assessing integrity of weather data for reference evapotranspiration estimation. J Irrig
#'     Drain Eng 122(2):97–106. \doi{10.1061/(ASCE)0733-9437(1996)122:2(97)}}
#'
#' \cite{Berger A, Loutre MF (1991) Insolation values for the climate of the last 10 million years. Quat Sci Rev
#'     10(4):297-317. \doi{10.1016/0277-3791(91)90033-Q}}
#'
#' \cite{Brock TD (1981) Calculating solar radiation for ecological studies. Ecol Model 14(1–2):1-19.
#'     \doi{10.1016/0304-3800(81)90011-9}}
#'
#' \cite{Davis TW, Prentice IC, Stocker BD, Thomas RT, Whitley RJ, Wang H, Evans BJ, Gallego-Sala AV, Sykes MT,
#'     Cramer W (2017) Simple process-led algorithms for simulating habitats (SPLASH v.1.0): robust indices of
#'     radiation, evapotranspiration and plant-available moisture. Geosci Model Dev 10(2):689–708.
#'     \doi{10.5194/gmd-10-689-2017}}
#'
#' \cite{Duffie JA, Beckman WA (1991) Solar Engineering of Thermal Processes. Second Edition. Wiley-Interscience,
#'     New York, NY}
#'
#' \cite{Yin X (1997a) Calculating daytime mean relative air mass. Agric For Meteorol 87(2-3):85-90.
#'     \doi{10.1016/S0168-1923(97)00029-4}}
#'
#' \cite{Yin X (1997b) Optical Air Mass: Daily Integration and its Applications. Meteorol Atmos Phys 63(3-4):227-233.
#'     \doi{10.1007/BF01027387}}
#'
#' \cite{Yin X (1998) The Albedo of Vegetated Land Surfaces: Systems Analysis and Mathematical Modeling. Theor Appl
#'     Climatol 60(1–4):121–140. \doi{10.1007/s007040050038}}
#'
#' \cite{Yin X (1999) Bright Sunshine Duration in Relation to Precipitation, Air Temperature and Geographic Location.
#'     Theor Appl Climatol 64(1–2):61–68. \doi{10.1007/s007040050111}}
#'
#' @examples
#' library (graphics)
#'
#' # Loading mandatory data for the Example 'Points'
#' data(inp_exPoints)
#'
#' # Measured and estimated one-year time series of the monthly mean relative sunshine duration,
#' # at a grid cell near Szeged, Hungary (46.3N, 20.2E), in the year 2010
#' with(inp_exPoints, {
#' bsdf01 <- matrix(nrow = 0, ncol = 12, dimnames = list(NULL, month.abb))
#' bsdf01 <- rbind(bsdf01, "Measured" = bsdf["2010", ])
#' bsdf01 <- rbind(bsdf01, "Solar123" = cliBrtSunDurFrcPoints(temp["2010", ], prec["2010", ],
#'     lat, lon, elv, year = 2010))
#' bsdf01 <- rbind(bsdf01, "SPLASH" = cliBrtSunDurFrcPoints(temp["2010", ], prec["2010", ],
#'     lat, lon, elv, year = 2010, aprchSIM = "SPLASH"))
#' cols <- c("black", "green", "blue")
#' matplot(t(bsdf01), type = "l", lwd = 2, col = cols, xaxt = "n", xlab = "Month",
#'     ylab = "Average relative sunshine duration (unitless)")
#' axis(1, at = seq(1, ncol(bsdf01)), labels = colnames(bsdf01))
#' legend(1, 0.7, legend = rownames(bsdf01), col = cols, lty = 1 : 2, lwd = 2, xpd = TRUE)
#' })
#'
#' \donttest{
#' # Relative root mean square error between measured and estimated values for the 'bsdf',
#' # at a grid cell near Szeged, Hungary (46.3N, 20.2E), in the period 1981-2010
#' with(inp_exPoints, {
#' years <- seq(1981, 2010)
#' bsdf02 <- cliBrtSunDurFrcPoints(temp, prec, lat, lon, elv, year = years)
#' rrmse <- function(pre, obs) { (sqrt(mean((pre - obs) ^ 2.)) / mean(obs)) * 100. }
#' rrmse_bsdf <- sapply(1 : 12, function(i) { rrmse(bsdf02[, i], bsdf[, i])  })
#' cols <- c("black", "green")
#' plot(rrmse_bsdf, type = "l", lwd = 2, col = cols, xaxt = "n", xlab = "Month",
#'     ylab = "Relative root mean square error (%)")
#' axis(1, at = 1 : 12, labels = month.abb)
#' })
#' }
#'
#' @importFrom rworldmap getMap
#' @importFrom sp CRS over proj4string SpatialPoints
#' @importFrom stats complete.cases setNames
#' @importFrom strex match_arg
#' @import rworldxtra
#'
#' @export
#'
cliBrtSunDurFrcPoints <- function(temp, prec, lat, lon, elv, year = 2000, aprchSIM = c("Solar123", "SPLASH")) {

  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  aprchSIM <- strex::match_arg(aprchSIM)

  errorChecking(year = year)

  # Vectorization of scalar variables
  cv.scl <- c("lat", "lon", "elv", "year")
  if (any(sapply(cv.scl, function(x) { (length(get(x)) == 1 & is.numeric(get(x))) | identical(get(x), NA) }))) {
    lgth <- errorHandling(temp = temp, prec = prec)$lgth
    list2env(sapply(cv.scl, function(x) {
      if ((length(get(x)) == 1 & is.numeric(get(x))) | identical(get(x), NA)) {
        assign(x, rep(get(x), lgth)) } else { assign(x, get(x)) } },
      simplify = FALSE), envir = environment())
  }
  err_han <- errorHandling(temp = temp, prec = prec, lat = lat, lon = lon, elv = elv, year = year)
  list2env(Filter(Negate(is.null), err_han), envir = environment())

  cv.arg <- c("temp", "prec", "lat", "lon", "elv")
  for (i in 1 : length(cv.arg)) {
    if (is.null(get(cv.arg[i]))) { stop("Invalid argument: '", cv.arg[i], "' is missing, with no default.") }
  }


  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 01. Select the calendar(s) applied
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # n_moy: the number of months in the year (12)
  # n_doy: the number of days in the year (365 or 366, depending on the leap year)
  # IDNM: the number of days in the specified month

  n_moy <- 12
  IDNM <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

  MvIDNM <- matrix(nrow = 2, ncol = n_moy, dimnames = list(c("nrml", "leap"), month.abb))
  MvIDNM[c("nrml", "leap"), ] <- matrix(rep(IDNM, 2), nrow = 2, byrow = TRUE)
  MvIDNM[c("leap"), 2] <- 29

  n_doy <- rep(NA, length(year))
  n_doy[!is.na(year)] <- ifelse(is.leap.year(year[!is.na(year)]), 366, 365)

  N_doy <- c(365, 366)
  cal <- setNames(lapply(vector("list", 2), function(x) vector()), c("nrml", "leap"))
  vld <- Reduce("&", list(do.call(complete.cases, list(temp, prec, lat, lon, elv, year))))
  cal$nrml <- as.numeric(which(Reduce("&", list(vld, sapply(n_doy, function(x) { x == N_doy[1] & !is.na(x) })))))
  cal$leap <- as.numeric(which(Reduce("&", list(vld, sapply(n_doy, function(x) { x == N_doy[2] & !is.na(x) })))))

  applCal <- match(sort(unique(n_doy[!is.na(n_doy)])), N_doy)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 02. Set the result object
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  bsdf <- matrix(nrow = lgth, ncol = n_moy, dimnames = list(NULL, month.abb))
  if (all(sapply(cal, function(x) { identical(x, numeric(0)) }))) { return(bsdf) }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 03. Estimate the monthly means of relative sunshine duration (bsdf), unitless
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Calculate monthly averages of daily solar irradiance for cloudless-sky conditions (R_E, MJ m-2 hr-1)
  # 'Solar123' - ref: 'text' in Yin (1999); Eq 3.2 in Yin (1997b)
  # 'SPLASH' - ref: Eqs 7, 10 and 11 in Davis et al. (2017); Eq 2 in Allen (1996)
  # Compute monthly means for daylength (DL, hr)
  # 'Solar123' - ref: 'text' in Yin (1997a); Eq 4 in Brock (1981); Eqs 5 and 6 in Forsythe et al. (1995)
  # 'SPLASH' - ref: Eq 8 in Davis et al. (2017); Eq 1.6.11 in Duffie and Beckman (1991)
  l.adsi <- cliAvgDlySolIrrPoints(lat, elv = elv, year = year, aprchSIM = aprchSIM, aprchTR = "hourly",
                                  mlyOpVar = c("R_E", "DL"))
  MvR_Emoa <- l.adsi$R_E_moa_MJ.m2.hr1
  MvDLmoa <- l.adsi$DL_moa_hr


  # Calculate the monthly average of daily potential evapotranspiration (EPmoa, mm day-1)
  # Saturated water vapor density, V_Ds (in g m-3)
  # ref: Eq A11 in Yin (1998); Burman and Pochop (1994)
  V_Ds <- 2170. / (temp + 273.) * (0.6108 * exp((17.27 * temp) / (temp + 237.3)))

  # Monthly mean for daily potential evapotranspiration, EPmoa (in mm day-1)
  # ref: Eq A10 in Yin (1998); Eq 2 in Hamon (1964)
  MvEPmoa <- 0.1651 * V_Ds * (MvDLmoa / 12.)


  # Compute the regional factors used by Eq 3.1 in Yin (1999) (f_i, unitless)
  # a multiplicative function reflective of regional idiosyncrasies
  # ref: Eq 3.3 in Yin (1999)
  f_i <- array(dim = c(lgth, 6, n_moy), dimnames = list(NULL, c("IA", "NA", "SA", "EA", "MA", "AF"), month.abb))

  # Regional data of the point to consider regional idiosyncrasies
  countriesSP <- getMap(resolution = 'high')
  vldPoints <- data.frame("lon" = lon[vld], "lat" = lat[vld])
  cntryCls <- over(SpatialPoints(vldPoints, proj4string = CRS(proj4string(countriesSP))), countriesSP)
  islCls <- over(SpatialPoints(vldPoints, proj4string = CRS(proj4string(islandsSP))), islandsSP)

  # island or Australia
  f_i[vld, "IA", ] <- ifelse(!is.na(islCls$featurecla) | cntryCls$continent == "Australia",
                          (15.6 - rowMeans(temp)) / 46.3, NA)

  # North America
  for (i_cal in applCal) {
    slctd <- cal[[c("nrml", "leap")[i_cal]]]
    MvDM <- MvIDNM[c("nrml", "leap")[i_cal], ]
    mbr2 <- 0.231 * apply(MvEPmoa[slctd, , drop = F], 1, function(x, y) { stats::weighted.mean(x, y) }, y = MvDM)
    mbr3 <- abs(rowMeans(temp[slctd, , drop = F])) / 72.5
    mbr4 <- do.call(pmin, as.data.frame(replace(temp[slctd, , drop = F], temp[slctd, , drop = F] < 0, 0))) / 30.
    f_i[slctd, "NA", ] <- ifelse(cntryCls$continent[slctd] == "North America", 0.331 - mbr2 - mbr3 + mbr4, NA)
  }

  # South America
  f_i[vld, "SA", ] <- ifelse(cntryCls$continent == "South America and the Caribbean", 0., NA)

  # Eurasia
  for (i_cal in applCal) {
    slctd <- cal[[c("nrml", "leap")[i_cal]]]
    MvDM <- MvIDNM[c("nrml", "leap")[i_cal], ]
    wrmst <- ifelse(lat >= 0., 7, 1)
    cldst <- ifelse(lat >= 0., 1, 7)
    numer <- (prec[slctd, wrmst, drop = F] / MvDM[wrmst]) - (prec[slctd, cldst, drop = F] / MvDM[cldst])
    denom <- rowSums(prec[slctd, , drop = F]) / N_doy[i_cal]
    f_i[slctd, "EA", ] <- ifelse(cntryCls$continent[slctd] == "Eurasia",
                                 (1. / 10.2) + (numer / (36.8 * denom)) * (2.3 - (numer / denom)), NA)
  }

  # monsoonal Asia
  patMA <- paste(c("Asia, South", "Asia, East", "Asia Pacific, High Income"), collapse="|")
  mbr2 <- 0.314 * MvR_Emoa
  mbr3 <- (31.8 - rowMeans(temp)) * (rowMeans(temp) / 140.) ** 2.
  slctd <- grepl(patMA, cntryCls$GBD)
  f_i[slctd, "MA", ] <- (-0.643 + mbr2 + mbr3)[slctd]

  # Africa
  df.temp <- as.data.frame(temp)[vld, ]
  f_i[vld, "AF", ] <- ifelse(cntryCls$continent == "Africa",
                             (7.26 - (do.call(pmax, df.temp) - do.call(pmin, df.temp))) / 32.3,
                          NA)

  # The regional factors used by Eq 3.1 in Yin (1999)
  # a multiplicative function reflective of regional idiosyncrasies (sumf_i, unitless)
  # ref: Eq 3.3 in Yin (1999)
  sumf_i <- matrix(apply(f_i, 3, function(x) { rowSums(x, na.rm = T) }), nrow = lgth, byrow = F)


  # Calculate the global factor used by Eq 3.1 in Yin (1999) (f_o, unitless)
  # Irradiation member in the equation representing the global trend (unitless)
  # ref: Eq 3.2 in Yin (1999) (see f_RE)
  f_RE <- (1. + 0.756 * MvR_Emoa * (3. - MvR_Emoa)) ** (-1.)

  # Precipitation member in the equation representing the global trend (unitless)
  # ref: Eq 3.2 in Yin (1999) (see f_P)
  f_P <- matrix(nrow = lgth, ncol = n_moy, dimnames = list(NULL, month.abb))
  for (i_cal in applCal) {
    slctd <- cal[[c("nrml", "leap")[i_cal]]]
    MvDM <- MvIDNM[c("nrml", "leap")[i_cal], ]
    numer <- 1. + 0.785 * t(apply(prec[slctd, , drop = F], 1, function(x, y) {x / y}, y = MvDM))
    denom <- 1. + 0.222 * t(apply(prec[slctd, , drop = F], 1, function(x, y) {x / y}, y = MvDM))
    f_P[slctd, ] <- t(sapply(1 : length(slctd), function(i) {numer[i, ] / denom[i, ]}))
  }

  # Potential evapotranspiration member in the equation representing the global trend (unitless)
  # ref: Eq 3.2 in Yin (1999) (see f_TEP)
  f_TEP <- 1. + (((7.66 * ifelse(temp <= 0., 1., 0.) - 4.98) * temp) + (MvEPmoa ** 2.)) / 184.

  # Latitude member in the equation representing the global trend (unitless)
  # ref: Eq 3.2 in Yin (1999) (see f_L)
  f_L <- 1. + 0.152 * sin((lat + 15.) * (2. * pi/ 77.))

  # The global factor used by Eq 3.1 in Yin (1999) (unitless)
  # a multiplicative function reflective of global generalities
  # ref: 3.2 in Yin (1999)
  f_o <- t(sapply(1 : lgth, function(i) {f_RE[i, ] * f_P[i, ] * f_TEP[i, ] * f_L[i]}))


  # Station atmospheric pressure, sap (unitless)
  # ref: Eq 2.11 in Yin (1997b); Linacre (1992)
  sap <- exp(-elv / 8000.)

  # Fraction of bright sunshine duration, n_N (unitless)
  # bsdf (or n_N): a ratio of the monthly mean of bright sunshine duration to the monthly mean for daylength
  # ref: 3.1 in Yin (1999)
  bsdf <- t(sapply(1 : lgth, function(i) { exp((-1.65 * sap * f_o)[i, ] * (1 + sumf_i[i, ])) }))

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  return(bsdf)

}


#' Estimator for Fraction of Bright Sunshine Duration
#'
#' @description Estimates monthly averages for daily fraction of bright sunshine duration, for a given region and
#'     year, by using the monthly time series of temperature and precipitation, and the elevation data.
#'
#' @param rs.temp multi-layer Raster* object with one-year time series of monthly mean air temperature (in °C)
#' @param rs.prec multi-layer Raster* object with one-year time series of monthly precipitation sum (in mm)
#' @param rl.elv single-layer Raster* object with the elevation values (in meters above sea level)
#' @param sc.year 'numeric' scalar with the value of the year (using astronomical year numbering)
#' @param aprchSIM 'character' vector of length 1 that indicates the formula used to estimate the value of solar
#'     irradiance/irradiation for a specific day. Valid values are as follows: \cr
#'     (a) \code{'Solar123'} -
#'     in this approach, first, the mean hourly solar irradiance under cloudless-sky conditions is calculated as
#'     proposed by Yin (1997b), with a minor modification, using the daytime means of optical air mass and
#'     cosine zenith; the former is computed as recommended by Yin (1997b), while the latter is estimated by using
#'     Eq 5 of Yin (1997a); however, in contrast to the original approach, where the solar constant was fixed at
#'     \eqn{4.9212 MJ m^{-2} hr^{-1}}{4.9212 MJ m^{-2} hr^{-1}}, according to Yin (1999), its value is corrected by
#'     calendar day for the variable ellipticity of the Earth's orbit, by using the scheme of Brock (1981); in the
#'     calculations, the values of solar declination and daylength are derived by using the approach of Brock (1981);
#'     \cr
#'     (b) \code{'SPLASH'} -
#'     in this approach, first, under varying orbital parameters, the daily solar radiation at the top of the
#'     atmosphere is calculated (\eqn{H_{0}}{H_{0}}, Eq 7 in Davis et al. (2017)), and then this value is
#'     multiplied by the atmospheric transmittivity to obtain the value of daily surface radiation; in this case as
#'     well, cloudless conditions are assumed, i.e., the transmission coefficient is taken into account with an
#'     universal value of 0.75, however, its value is modified as a function of elevation, by using the scheme of
#'     Allen (1996); the daylength is calculated via Eq 1.6.11 in Duffie and Beckman (1991), using the sunset hour
#'     angle (\eqn{h_{s}}{h_{s}}, Eq 8. in Davis et al. (2017)); finally, the mean hourly surface radiation is
#'     derived as the quotient of the daily surface radiation and the daylength.
#' @param filename output filename
#' @param ... additional arguments passed on to \code{\link[raster]{writeRaster}}
#'
#' @details See \code{\link[macroBiome]{cliBrtSunDurFrcPoints}}.
#'
#' @return A 12-layer RasterStack with one-year time series of monthly mean relative sunshine duration.
#'
#' @note The objects \code{'rs.temp'} and \code{'rs.prec'} must be 12-layer Raster* objects, while the object
#'     \code{'rl.elv'} has to be a single-layer Raster* object. These Raster* objects must have the same bounding
#'     box, projection, and resolution. The object \code{'sc.year'} has to be a single integer number.
#'
#' @references
#'
#' \cite{Allen RG (1996) Assessing integrity of weather data for reference evapotranspiration estimation. J Irrig
#'     Drain Eng 122(2):97–106. \doi{10.1061/(ASCE)0733-9437(1996)122:2(97)}}
#'
#' \cite{Berger A, Loutre MF (1991) Insolation values for the climate of the last 10 million years. Quat Sci Rev
#'     10(4):297-317. \doi{10.1016/0277-3791(91)90033-Q}}
#'
#' \cite{Brock TD (1981) Calculating solar radiation for ecological studies. Ecol Model 14(1–2):1-19.
#'     \doi{10.1016/0304-3800(81)90011-9}}
#'
#' \cite{Davis TW, Prentice IC, Stocker BD, Thomas RT, Whitley RJ, Wang H, Evans BJ, Gallego-Sala AV, Sykes MT,
#'     Cramer W (2017) Simple process-led algorithms for simulating habitats (SPLASH v.1.0): robust indices of
#'     radiation, evapotranspiration and plant-available moisture. Geosci Model Dev 10(2):689–708.
#'     \doi{10.5194/gmd-10-689-2017}}
#'
#' \cite{Duffie JA, Beckman WA (1991) Solar Engineering of Thermal Processes. Second Edition. Wiley-Interscience,
#'     New York, NY}
#'
#' \cite{Yin X (1997a) Calculating daytime mean relative air mass. Agric For Meteorol 87(2-3):85-90.
#'     \doi{10.1016/S0168-1923(97)00029-4}}
#'
#' \cite{Yin X (1997b) Optical Air Mass: Daily Integration and its Applications. Meteorol Atmos Phys 63(3-4):227-233.
#'     \doi{10.1007/BF01027387}}
#'
#' \cite{Yin X (1999) Bright Sunshine Duration in Relation to Precipitation, Air Temperature and Geographic Location.
#'     Theor Appl Climatol 64(1–2):61–68. \doi{10.1007/s007040050111}}
#'
#' @examples
#' # Loading mandatory data for the Example 'Single-Year Grid'
#' data(inp_exSglyGrid)
#' inp_exSglyGrid <- lapply(inp_exSglyGrid, crop, extent(20.15, 20.25, 46.25, 46.35))
#'
#' # Estimate values of the monthly mean relative sunshine duration
#' # at a grid cell near Szeged, Hungary (46.3N, 20.2E), in the year 2010
#' with(inp_exSglyGrid, {
#' rs.bsdf <- cliBrtSunDurFrcGrid(temp, prec, elv, sc.year = 2010)
#' rs.bsdf
#' })
#'
#' @importFrom strex match_arg
#' @import raster
#'
#' @export
#'
cliBrtSunDurFrcGrid <- function(rs.temp, rs.prec, rl.elv, sc.year = 2000, aprchSIM = c("Solar123", "SPLASH"),
                                filename = "", ...) {

  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  aprchSIM <- strex::match_arg(aprchSIM)

  if (length(sc.year) == 1L & is.numeric(sc.year)) {
    if (sc.year %% 1 != 0) {
      stop("Invalid argument: 'sc.year' has to be a single integer number.")
    }
  } else {
    stop("Invalid argument: 'sc.year' has to be a single integer number.")
  }

  errorCheckingGrid(rs.temp = rs.temp, rs.prec = rs.prec, rl.elv = rl.elv)

  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  rl.lat <- getGeogrCoord(raster(rs.temp, layer = 1), "lat")
  rl.lon <- getGeogrCoord(raster(rs.temp, layer = 1), "lon")

  rs.rslt <- stack(brick(rs.temp))

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
      cv.mly_var <- c("rs.temp", "rs.prec")
      cv.loc_dta <- c("rl.lat", "rl.lon", "rl.elv")
      cv.arg <- c(cv.mly_var, cv.loc_dta)
      for (i_arg in 1 : length(cv.arg)) {
        assign(substring(cv.arg[i_arg], 4), getValues(get(cv.arg[i_arg]), row = bs$row[i], nrows = bs$nrows[i]))
      }
      mx.rslt <- cliBrtSunDurFrcPoints(temp, prec, lat, lon, elv, year = sc.year, aprchSIM = aprchSIM)

      rs.rslt <- writeValues(rs.rslt, mx.rslt, bs$row[i])
      pbStep(pb, i)
    }
    rs.rslt <- writeStop(rs.rslt)
  } else {
    for (i in 1 : bs$n) {
      cv.mly_var <- c("rs.temp", "rs.prec")
      cv.loc_dta <- c("rl.lat", "rl.lon", "rl.elv")
      cv.arg <- c(cv.mly_var, cv.loc_dta)
      for (i_arg in 1 : length(cv.arg)) {
        assign(substring(cv.arg[i_arg], 4), getValues(get(cv.arg[i_arg]), row = bs$row[i], nrows = bs$nrows[i]))
      }
      mx.rslt <- cliBrtSunDurFrcPoints(temp, prec, lat, lon, elv, year = sc.year, aprchSIM = aprchSIM)

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
  names(rs.rslt) <- colnames(mx.rslt)
  return(rs.rslt)

}


getGeogrCoord <- function(rstr, varName = c("lon", "lat")) {
  # Set the variable name for the content of the result
  varName <- strex::match_arg(varName)
  # Create an array for the results
  out <- raster::raster(rstr)
  # Get the index of the blocks, every block has n rows, bigger the minblocks, smaller the chunk of rows
  bs <- raster::blockSize(rstr, minblocks = 200)
  # Check whether or not it is necessary to transform the coordinate system
  xform <- !(raster::compareCRS(rstr, sp::CRS("+init=epsg:4326")))
  filename <- tempfile()
  out <- raster::writeStart(out, filename, overwrite = TRUE)
  for (i in 1 : bs$n) {
    xncells <- raster::cellFromRow(rstr, (bs$row[i] : (bs$row[i] + bs$nrows[i] - 1)))
    xmat <- raster::getValues(rstr, bs$row[i], bs$nrows[i])
    if (xform) {
      coords <- data <- as.data.frame(raster::xyFromCell(rstr, xncells))
      xydata <- sp::spTransform(sp::SpatialPointsDataFrame(coords, data, proj4string = raster::crs(rstr)),
                                sp::CRS("+init=epsg:4326"))@coords
    } else {
      xydata <- raster::xyFromCell(rstr, xncells)
    }

    # Write the chunk of results, bs$row[i] is putting the results in the correct rows
    out <- raster::writeValues(out, xydata[, which(varName == c("lon", "lat"))], bs$row[i])
  }
  out <- raster::writeStop(out)
  out <- raster::mask(out, rstr)
  return(out)
}

# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script contains a vectorized function for converting cloud cover data to values of the sunshine duration:
#   cldn2bsdf(cldn)
#
#### DEFINE FUNCTIONS ################################################################################################

# ********************************************************************************************************************
# Name:     cldn2bsdf
# Inputs:   - double, monthly mean cloud cover fraction, in percentage (cldn)
# Returns:  - double, monthly mean relative sunshine duration, unitless (bsdf)
# Features: This function converts cloudiness values to relative sunshine duration data.
# Ref:      - Doorenbos J, Pruitt WO (1977) Guidelines for predicting crop water requirements. FAO Irrigation and
#             Drainage Paper No. 24, Technical Report. Food and Agriculture Organization of the United Nations,
#             Italy, Rome
# ********************************************************************************************************************
cldn2bsdf <- function(cldn) {
  bsdf <- rep(NA, length = length(cldn))
  vld <- !is.na(cldn) & cldn <= 100. & cldn >= 0.
  param <- cbind(upr = c(100., 80., 60., 50., 30.), lwr = c(80., 60., 50., 30., 10.),
                 a = c(-0.015, -0.01, -0.005, -0.01, -0.005), b = c(1.5, 1.1, 0.8, 1.05, 0.9))
  for (i_eq in 1 : nrow(param)) {
    slctd <- (cldn[vld] <= param[i_eq, "upr"] & cldn[vld] > param[i_eq, "lwr"])
    bsdf[vld][slctd] <- param[i_eq, "a"] * cldn[vld][slctd] + param[i_eq, "b"]
  }
  slctd <- (cldn[vld] <= 10.)
  bsdf[vld][slctd] <- -0.01 * cldn[vld][slctd] + 0.95
  return(bsdf)
}
