#' Daily Weather Generator
#'
#' @description Generates quasi-daily time series from the monthly mean values of temperature, precipitation and
#'     sunshine data.
#'
#' @param temp 'numeric' R object with one-year time series of monthly mean air temperature (in °C)
#' @param prec 'numeric' R object with one-year time series of monthly precipitation sum (in mm)
#' @param bsdf 'numeric' R object with one-year time series of monthly mean relative sunshine duration (dimensionless)
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
#' @param argCkd 'logical' scalar that indicates whether or not the checking and correction of arguments can be
#'    omitted.
#'
#' @details For relative sunshine duration and air temperature, it is recommended the 'harmonic' interpolation
#'     technique (\code{'hip'}) described by Epstein (1991), with a correction of physically impossible values. This
#'     technique can also be used to values of the mean precipitation intensity, however, in the case of
#'     precipitation, the temporal scaling (\code{'tsi'}) using an iterative interpolation technique described by
#'     Lüdeke et al. (1994) is recommended, with a damping variable of 0.7 for each month. The damping variable can
#'     be set separately for each month and must be between 0 and 1. This iterative scheme can also be applied to
#'     monthly mean data of air temperature. Furthermore, for all three climate variables, it is possible to ignore
#'     intra-month variability. For this, the setting \code{'const'} must be used.
#'
#' @return A list with three 365- or 366-column matrices that contain quasi-daily time series for the three basic
#'     climate variables:
#'
#'     \itemize{
#'       \item{\code{TEMP}: daily mean air temperature (in °C)}
#'       \item{\code{PREC}: daily precipitation sum (in mm)}
#'       \item{\code{BSDF}: daily fractional sunshine duration (dimensionless)}
#'     }
#'
#' @note As with any point function, a set of basic input data is defined here. In this case, they are as follows:
#'    \code{'temp'} (one-year time series of monthly mean air temperature), \code{'prec'} (one-year time series
#'    of monthly precipitation sum), and \code{'bsdf'} (one-year time series of monthly mean relative sunshine
#'    duration.) The objects \code{'temp'}, \code{'prec'} and \code{'bsdf'} must be either vectors of length 12 or
#'    12-column matrices. The first dimensions of these matrices have to be the same length. The function
#'    automatically converts vectors into single-row matrices during the error handling, and then uses these
#'    matrices. The first dimensions of these matrices determines the number of rows in the result matrix. In the
#'    case of arguments that do not affect the course of the calculation procedure or the structure of the return
#'    object, scalar values (i.e., 'numeric' vector of length 1) may also be allowed. In this case, it is as
#'    follow: \code{'year'} (year using astronomical year numbering). This scalar is converted to a vector by the
#'    function during the error handling, and this vector is applied in the further calculations. If these data are
#'    stored in vectors of length at least 2, their length must be the same size of first dimension of the matrices
#'    containing the basic data.
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
#' library(graphics)
#'
#' # Loading mandatory data for the Example 'Points'
#' data(inp_exPoints)
#'
#' with(inp_exPoints, {
#' # Generate quasi-daily time series for basic climate variables with default settings,
#' # at a grid cell near Szeged, Hungary (46.3N, 20.2E) (for the normal period 1981-2010)
#' year <- trunc(mean(seq(1981, 2010)))
#' wea01 <- dlyWeaGenPoints(colMeans(temp), colMeans(prec), colMeans(bsdf), year = year)
#'
#' # Modify the daily weather data generation techniques
#' # To temperature data, apply the iterative interpolation technique with basic settings
#' # To precipitation data, change the value of the damping variable, over the whole year
#' # To sunshine data, assume that its values are constant within each month
#' wea02 <- dlyWeaGenPoints(colMeans(temp), colMeans(prec), colMeans(bsdf), aprchTEMP = "tsi",
#'     aprchBSDF = "const", dvPREC = rep(0.6, 12), year = year)
#'
#' # Check the differences
#' vars <- c("TEMP", "PREC", "BSDF")
#' lbls <- list(expression(italic(T[a])~(~degree*C)), expression(italic(P[n])~(mm)),
#'     expression(italic(S[f])~(unitless)))
#' ys <- c(20, 2.5, 0.6)
#' ats <- list(seq(-4, 24, 4), seq(0, 3, 0.5), seq(0., 0.8, 0.2))
#' cols <- c("black", "green")
#' opar <- par(no.readonly = TRUE)
#' par(mfrow = c(3, 1))
#' for (i in 1 : length(vars)) {
#'   par(mar = c(2, 5, 1, 1))
#'   matplot(t(rbind(wea01[[vars[i]]], wea02[[vars[i]]])), type = "l", lwd = 2, col = cols,
#'       xaxt = "n", xlab = NA, ylab = NA, axes = FALSE)
#'   axis(side = 1, las = 1, tck = -0.03, labels = NA, at = seq(-60, 720, 30))
#'   axis(side = 2, las = 1, tck = -0.03, labels = NA, at = ats[[i]])
#'   axis(side = 2, las = 1, lwd = 0, line = -0.4, cex.axis = 1.6, at = ats[[i]])
#'   if (i == length(vars)) {
#'     axis(side = 1, las = 1, lwd = 0, line = -0.4, at = seq(-60, 720, 30), cex.axis = 1.6)
#'   }
#'   mtext(side = 2, lbls[[i]], line = 3, cex = 1.1)
#'   legend(1, ys[i], legend = c("default", "modified"), col = cols, lty = 1 : 2, lwd = 2, xpd = TRUE)
#' }
#' par(opar)
#' })
#'
#' @importFrom stats complete.cases setNames
#' @importFrom strex match_arg
#'
#' @export
#'
dlyWeaGenPoints <- function(temp, prec, bsdf, year = 2000, aprchTEMP = c("hip", "tsi", "const"),
                            aprchPREC = c("tsi", "hip", "const"), aprchBSDF = c("hip", "const"),
                            dvTEMP = rep(0.7, 12), dvPREC = rep(0.7, 12), argCkd = FALSE) {

  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  aprchTEMP <- strex::match_arg(aprchTEMP)
  aprchPREC <- strex::match_arg(aprchPREC)
  aprchBSDF <- strex::match_arg(aprchBSDF)

  errorChecking(year = year, dvTEMP = dvTEMP, dvPREC = dvPREC)

  # Vectorization of the scalar variable
  if ((length(year) == 1 & is.numeric(year))) {
    lgth <- errorHandling(bsdf = bsdf, temp = temp, prec = prec)$lgth
    year <- rep(year, lgth)
  }
  err_han <- errorHandling(bsdf = bsdf, temp = temp, prec = prec, year = year, onlyLgthChecking = argCkd)
  list2env(Filter(Negate(is.null), err_han), envir = environment())

  cv.arg <- c("temp", "prec", "bsdf")
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

  MvIDNM <- MvCDNM <- matrix(nrow = 2, ncol = n_moy, dimnames = list(c("nrml", "leap"), month.abb))
  MvIDNM[c("nrml", "leap"), ] <- matrix(rep(IDNM, 2), nrow = 2, byrow = TRUE)
  MvIDNM[c("leap"), 2] <- 29

  n_doy <- rep(NA, length(year))
  n_doy[!is.na(year)] <- ifelse(is.leap.year(year[!is.na(year)]), 366, 365)

  N_doy <- c(365, 366)
  cal <- setNames(lapply(vector("list", length(cv.arg)), function(x) vector()), cv.arg)
  for (i_arg in 1 : length(cv.arg)) {
    var <- get(cv.arg[i_arg])
    cal[[i_arg]] <- setNames(lapply(vector("list", 2), function(x) vector()), c("nrml", "leap"))
    vld <- Reduce("&", list(do.call(complete.cases, list(year, var))))
    cal[[i_arg]]$nrml <- as.numeric(
      which(Reduce("&", list(vld, mapply(function(x) { x == N_doy[1] & !is.na(x) }, n_doy)))))
    cal[[i_arg]]$leap <- as.numeric(
      which(Reduce("&", list(vld, mapply(function(x) { x == N_doy[2] & !is.na(x) }, n_doy)))))
  }

  applCal <- match(sort(unique(n_doy[!is.na(n_doy)])), N_doy)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 02. Set the result object
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rslt <- list()
  for (i_arg in 1 : length(cv.arg)) {
    if (all(sapply(cal[[i_arg]], function(x) { identical(x, numeric(0)) }))) {
      rslt[[toupper(cv.arg[i_arg])]] <-
        matrix(nrow = lgth, ncol = ifelse(!identical(integer(0), applCal), max(N_doy[applCal]), 365))
    }
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 03. Run the algorithm for each calendar type
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (i_arg in 1 : length(cv.arg)) {

    if (!all(sapply(cal[[i_arg]], function(x) { identical(x, numeric(0)) }))) {

      varLbl <- varName <- cv.arg[i_arg]
      var <- get(varLbl)
      aprchVAR <- get(paste0("aprch", toupper(varLbl)))
      if (varLbl != "bsdf") { dvVAR <- get(paste0("dv", toupper(varLbl))) }
      VAR <- matrix(nrow = lgth, ncol = ifelse(!identical(integer(0), applCal), max(N_doy[applCal]), 365))
      if (varName == "prec") { varName <- "pint" }

      for (i_cal in applCal) {

        slctd <- cal[[i_arg]][[c("nrml", "leap")[i_cal]]]

        if (varName == "pint") {
          slctdVar <- var[slctd, , drop = FALSE] %*% diag(1. / MvIDNM[c("nrml", "leap")[i_cal], ])
        } else {
          slctdVar <- var[slctd, , drop = FALSE]
        }

        tstp <- seq(1, N_doy[i_cal])
        if (aprchVAR == "tsi") {
          VAR[slctd, tstp] <- dwgTmpScLinIntrpl(slctdVar, varName, year[slctd], dvVAR, argCkd = TRUE)
        } else if (aprchVAR == "hip") {
          VAR[slctd, tstp] <- dwgHarmonicIntrpl(slctdVar, varName, year[slctd], argCkd = TRUE)
        } else {
          VAR[slctd, tstp] <- dwgTmpScLinIntrpl(slctdVar, varName, year[slctd], rep(1., n_moy), argCkd = TRUE)
        }

      }

      rslt[[toupper(varLbl)]] <- VAR
      rm(VAR)

    }

  }

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  rslt <- rslt[toupper(cv.arg)]
  return(rslt)

}


#### DEFINE FUNCTIONS ################################################################################################

# ********************************************************************************************************************
# Name:     dwgTmpScLinIntrpl
# Inputs:   - double, monthly time series of the given climate variable (var)
#           - character, name of the given climate variable (varName)
#           - double, year (using astronomical year numbering) (year)
#           - double, monthly time series of the damping variable for the given climate variable (dvVAR)
#           - logical, indicates whether or not the checking and correction of arguments can be omitted (argCkd)
# Returns:  - double, pseudo-daily time series of the given climate variable
# Features: This function interpolates monthly means of a climate variable to daily scale using an iterative
#           technique proposed by Lüdeke et al. (1994). This is an algorithm whose idea is to divide each doubled
#           value of the monthly mean (M) into two parts (p * M and (1 - p) * M) in the first step, as a function of
#           a ratio of adjacent monthly values. The value of p must be between 0 and 1. The value of p can be
#           modified via a damping variable (DV). If the value of DV is equal to 0, then the influence of adjacent
#           values is maximal during the decomposition, whereas for DV = 1, constant values are obtained. The typical
#           value for DV is 0.7. In accordance with the original procedure, six iteration steps are used. Then, the
#           resulting 768 values are redistributed over 365 or 366 days, depending on the astronomical year.
#           Different DV values can be set for each month (see dvVAR).
# Ref:      - Lüdeke MKB, Badeck FW, Otto RD, Häger C, Dönges S, Kindermann J, Würth G, Lang T, Jäkel U, Klaudius A,
#             Ramge P, Habermehl S, Kohlmaier GH (1994) The Frankfurt Biosphere Model: A global process-oriented
#             model of seasonal and long-term CO2 exchange between terrestrial ecosystems and the atmosphere. I.
#             Model description and illustrative results for cold deciduous and boreal forests. Clim Res 4(2):143-166.
#             DOI: 10.3354/cr004143
# ********************************************************************************************************************
dwgTmpScLinIntrpl <- function(var, varName = c("temp", "pint", "bsdf"), year = 2000, dvVAR = rep(0.7, 12),
                              argCkd = FALSE) {

  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  varName <- strex::match_arg(varName)

  if (!argCkd) { errorChecking(year = year, dvVAR = dvVAR) }

  arg <- as.list(formals(errorHandling))
  arg[[varName]] <- var
  arg$onlyLgthChecking <- argCkd
  if ((length(year) == 1 & is.numeric(year))) {
    arg$year <- rep(year, do.call(errorHandling, arg)$lgth)
  } else {
    arg$year <- year
  }
  list2env(Filter(Negate(is.null), do.call(errorHandling, arg)), envir = environment())
  var <- get(varName)

  if (is.null(var)) { stop("Invalid argument: 'var' is missing, with no default.") }


  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 01. Select the calendar(s) applied
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # n_moy: the number of months in the year (12)
  # n_doy: the number of days in the year (365 or 366, depending on the leap year)
  # i_moy: the month number of the year, ranging from 1 (January) to 12 (December)
  # i_doy: the day number of the year, ranging from 1 (1 January) to 365 (31 December)
  # IDNM: the number of days in the specified month
  # CDNM: the number of days in the months prior to the specified month

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
  cal <- stats::setNames(lapply(vector("list", 2), function(x) vector()), c("nrml", "leap"))
  vld <- do.call(stats::complete.cases, list(year, var))
  cal$nrml <- as.numeric(which(Reduce("&", list(vld, sapply(n_doy, function(x) { x == N_doy[1] & !is.na(x) })))))
  cal$leap <- as.numeric(which(Reduce("&", list(vld, sapply(n_doy, function(x) { x == N_doy[2] & !is.na(x) })))))

  applCal <- match(sort(unique(n_doy[!is.na(n_doy)])), N_doy)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 02. Set the result object
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  VAR <- matrix(nrow = lgth, ncol = ifelse(!identical(integer(0), applCal), max(N_doy[applCal]), 365))
  if (all(sapply(cal, function(x) { identical(x, numeric(0)) }))) { return(VAR) }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 03. Generate a high-resolution synthetic time series using an iterative interpolation technique
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  incl <- sort(as.numeric(unlist(cal)))
  synHgh <- var
  dvsHgh <- dvVAR
  n <- ncol(synHgh)
  n_lp <- 6
  for (i_lp in 1 : n_lp) {
    m <- n
    n <- n_moy * 2 ^ i_lp
    synLwr <- matrix(nrow = lgth, ncol = (m + 2))
    synLwr[incl, ] <- cbind(synHgh[incl, m, drop = FALSE], synHgh[incl, , drop = FALSE],
                            synHgh[incl, 1, drop = FALSE]) * 2.
    dvsLwr <- dvsHgh
    p <- matrix(nrow = lgth, ncol = m)
    for (i_ts in 1 : m) {
      numer <- synLwr[incl, (i_ts + 1) - 1, drop = FALSE]
      denom <- (synLwr[incl, (i_ts + 1) - 1, drop = FALSE] + synLwr[incl, (i_ts + 1) + 1, drop = FALSE])
      ratio <- ifelse(lapply(denom, function(x) {x == 0.}), 0.,
                      ifelse(lapply(numer / denom, function(x) {x > 1.}), 1.,
                             numer / denom))
      p[incl, i_ts] <- (1. - dvsLwr[i_ts]) * ratio + dvsLwr[i_ts] / 2.
    }
    synHgh <- matrix(nrow = lgth, ncol = n)
    for (i_ts in 1 : m) {
      synHgh[incl, (i_ts * 2) - 1] <- p[incl, i_ts, drop = FALSE] * synLwr[incl, i_ts + 1, drop = FALSE]
      synHgh[incl, (i_ts * 2)] <- (1. - p[incl, i_ts, drop = FALSE]) * synLwr[incl, i_ts + 1, drop = FALSE]
      dvsHgh[(i_ts * 2) - 1] <- dvsHgh[(i_ts * 2)] <- dvsLwr[i_ts]
    }
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 04. Upscale the synthetic time series to daily resolution (for each calendar type)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (i_cal in applCal) {

    slctd <- cal[[c("nrml", "leap")[i_cal]]]

    for (i_moy in 1 : n_moy) {

      synSbst <- synHgh[slctd, (((i_moy - 1) * 2 ^ n_lp) + 1) : (i_moy * 2 ^ n_lp), drop = FALSE]

      n_dom <- MvIDNM[c("nrml", "leap")[i_cal], i_moy]
      n_syn <- (2 ^ n_lp)
      mtx <- matrix(0, nrow = n_dom, ncol = n_syn)

      i_hor <- i_ver <- 1
      divd <- (n_syn / n_dom)
      while (i_hor <= n_syn & i_ver <= n_dom) {
        if (divd >= 1.) {
          extd <- 1. - sum(mtx[, i_hor])
          remd <- divd - extd
        } else {
          extd <- divd
          remd <- (n_syn / n_dom)
        }
        mtx[i_ver, i_hor] <- extd

        if (sum(mtx[, i_hor]) < 1.) {
          i_ver <- i_ver + 1
        } else {
          i_hor <- i_hor + 1
        }
        divd <- remd
      }

      dlySbst <- t(apply(synSbst, 1, function(x) {mtx %*% x})) * (n_dom / n_syn)

      for (i_dom in 1 : MvIDNM[c("nrml", "leap")[i_cal], i_moy]) {
        i_doy <- MvCDNM[c("nrml", "leap")[i_cal], i_moy] + i_dom
        VAR[slctd, i_doy] <- dlySbst[, i_dom]
      }

    }

  }

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  # Error checking of the result
  metVar <- data.frame(cep = c("temperature", "precipitation", "sunshine"), mly = c("temp", "pint", "bsdf"),
                       dly = c("TEMP", "PREC", "BSDF"))
  cepN <- metVar$cep[metVar$mly == varName]
  dlyN <- metVar$dly[metVar$mly == varName]

  pImp <- which(apply(VAR, 1, function(r) {
    any(r[!is.na(r)] < varFeatures[metVar$dly[metVar$mly == varName], "lwst_val"] |
          r[!is.na(r)] > varFeatures[metVar$dly[metVar$mly == varName], "hgst_val"]) }))

  if (!identical(integer(0), pImp)) {
    warning("For ", cepN , " data, the approach chosen to generate quasi-daily data leds to physically impossible ",
            "results for certain monthly time series. For this reason, to time series of the result object '", dlyN,
            "' that affected by this problem, the NA value is assigned.")
    VAR[pImp, ] <- NA
  }

  return(VAR)

}


# ********************************************************************************************************************
# Name:     dwgHarmonicIntrpl
# Inputs:   - double, monthly time series of the given climate variable (var)
#           - character, name of the given climate variable (varName)
#           - double, year (using astronomical year numbering) (year)
#           - logical, indicates whether or not the checking and correction of arguments can be omitted (argCkd)
# Returns:  - double, pseudo-daily time series of the given climate variable
# Features: This function interpolates the monthly data to pseudo-daily values by using the mean-preserving 'harmonic'
#           interpolation method of Epstein (1991). Currently, valid values are as follows: 'temp' - monthly mean air
#           temperature (in deg C); 'pint' - monthly mean precipitation intensity (in mm dy-1); 'bsdf' - monthly mean
#           relative sunshine duration (unitless).
#           This scheme can occasionally produce values that are physically impossible. It can happen if the method
#           is applied to variables like precipitation, which can have only non-negative values. Therefore, in the
#           last step, a correction procedure is applied for each climate variable that has any physical constraints.
#           Each physically impossible value is set to the maximal/minimal possible value, and then the resulting
#           difference are aggregated for each month. In the case of each month, finally, the monthly value of the
#           differences are distributed among the days of the given month. This procedure is performed iteratively
#           several times.
# Ref:      - Epstein ES (1991) On Obtaining Daily Climatological Values from Monthly Means. J Clim 4(3):365–368.
#             DOI: 10.1175/1520-0442(1991)004<0365:OODCVF>2.0.CO;2
# ********************************************************************************************************************
dwgHarmonicIntrpl <- function(var, varName = c("temp", "pint", "bsdf"), year = 2000, argCkd = FALSE) {

  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  varName <- strex::match_arg(varName)

  if (!argCkd) { errorChecking(year = year) }

  arg <- as.list(formals(errorHandling))
  arg[[varName]] <- var
  arg$onlyLgthChecking <- argCkd
  if ((length(year) == 1 & is.numeric(year))) {
    arg$year <- rep(year, do.call(errorHandling, arg)$lgth)
  } else {
    arg$year <- year
  }
  list2env(Filter(Negate(is.null), do.call(errorHandling, arg)), envir = environment())
  var <- get(varName)

  if (is.null(var)) { stop("Invalid argument: 'var' is missing, with no default.") }


  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 01. Select the calendar(s) applied
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # n_moy: the number of months in the year (12)
  # n_doy: the number of days in the year (365 or 366, depending on the leap year)
  # i_moy: the month number of the year, ranging from 1 (January) to 12 (December)
  # i_doy: the day number of the year, ranging from 1 (1 January) to 365 (31 December)
  # IDNM: the number of days in the specified month
  # CDNM: the number of days in the months prior to the specified month

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
  cal <- stats::setNames(lapply(vector("list", 2), function(x) vector()), c("nrml", "leap"))
  vld <- do.call(stats::complete.cases, list(year, var))
  cal$nrml <- as.numeric(which(Reduce("&", list(vld, sapply(n_doy, function(x) { x == N_doy[1] & !is.na(x) })))))
  cal$leap <- as.numeric(which(Reduce("&", list(vld, sapply(n_doy, function(x) { x == N_doy[2] & !is.na(x) })))))

  applCal <- match(sort(unique(n_doy[!is.na(n_doy)])), N_doy)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 02. Set the result object
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  VAR <- matrix(nrow = lgth, ncol = ifelse(!identical(integer(0), applCal), max(N_doy[applCal]), 365))
  if (all(sapply(cal, function(x) { identical(x, numeric(0)) }))) { return(VAR) }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 03. Calculate the calendar-independent coefficients
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  n_coef <- 6
  aCoef <- bCoef <- matrix(nrow = lgth, ncol = n_coef)
  # ref: Eq 6.1 in Epstein (1991)
  aCoef0 <- rowMeans(var)

  # ref: Eqs 6.2 and 6.3 in Epstein (1991)
  for (i_coef in 1 : (n_coef - 1)) {
    aSum <- bSum <- rep(0., lgth)
    cmpt1 <- pi * i_coef / n_moy
    for (i_moy in 1 : n_moy) {
      cmpt2 <- (2. * pi * i_coef * i_moy) / n_moy
      aSum <- aSum + (var[, i_moy] * cos(cmpt2)) / n_coef
      bSum <- bSum + (var[, i_moy] * sin(cmpt2)) / n_coef
    }
    aCoef[, i_coef] <- (cmpt1 / sin(cmpt1)) * aSum
    bCoef[, i_coef] <- (cmpt1 / sin(cmpt1)) * bSum
  }

  # ref: Eq 6.4 in Epstein (1991)
  aSum <- rep(0., lgth)
  for (i_moy in 1 : n_moy) {
    cmpt3 <- pi * i_moy
    aSum <- aSum + (var[, i_moy] * cos(cmpt3)) / n_moy
  }
  cmpt4 <- pi / 2.
  aCoef[, n_coef] <- (cmpt4 / sin(cmpt4)) * aSum

  # ref: Eq 6.5 in Epstein (1991)
  bCoef[, n_coef] <- 0.

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 03. Run the algorithm for each calendar type
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (i_cal in applCal) {

    slctd <- cal[[c("nrml", "leap")[i_cal]]]

    # ref: Eq 2 in Epstein (1991)
    i_doy <- 0
    for (i_moy in 1 : n_moy) {
      for (i_dom in 1 : MvIDNM[c("nrml", "leap")[i_cal], i_moy]) {
        i_doy <- MvCDNM[c("nrml", "leap")[i_cal], i_moy] + i_dom
        t <- (i_moy - 0.5) + (i_dom - 0.5) / MvIDNM[c("nrml", "leap")[i_cal], i_moy]
        VAR[slctd, i_doy] <- aCoef0[slctd]
        for (i_coef in 1 : n_coef) {
          cmpt5 <- ((2. * pi * i_coef * t) / n_moy)
          VAR[slctd, i_doy] <- VAR[slctd, i_doy] +
            aCoef[slctd, i_coef] * cos(cmpt5) + bCoef[slctd, i_coef] * sin(cmpt5)
        }
      }
    }

    if (varName == "pint" | varName == "bsdf") {
      # Set all negative values equal to zero
      VAR[slctd, ][VAR[slctd, ] < 0.] <- 0.

      # Set all values greater than one to one ('bsdf')
      if (varName == "bsdf") { VAR[slctd, ][VAR[slctd, ] > 1.] <- 1. }

      i_lp <- 0
      while(i_lp <= 30) {

        # Set each daily value equal to zero for which the monthly value is zero
        # Set each daily value equal to one for which the monthly value is one ('bsdf')
        i_doy <- 0
        for (i_moy in 1 : n_moy) {
          for (i_dom in 1 : MvIDNM[c("nrml", "leap")[i_cal], i_moy]) {
            i_doy <- MvCDNM[c("nrml", "leap")[i_cal], i_moy] + i_dom
            VAR[slctd[which(var[slctd, i_moy] == 0.)], i_doy] <- 0.
            if (varName == "bsdf") { VAR[slctd[which(var[slctd, i_moy] == 1.)], i_doy] <- 1. }
          }
        }

        # Calculate the correction factor for each month
        i_doy <- 0
        xdm <- vldNum <- matrix(0., nrow = lgth, ncol = n_moy)
        totDiff <- rep(NA, lgth)
        for (i_moy in 1 : n_moy) {
          for (i_dom in 1 : MvIDNM[c("nrml", "leap")[i_cal], i_moy]) {
            i_doy <- MvCDNM[c("nrml", "leap")[i_cal], i_moy] + i_dom
            xdm[slctd, i_moy] <- xdm[slctd, i_moy] + VAR[slctd, i_doy] / MvIDNM[c("nrml", "leap")[i_cal], i_moy]
            if (varName == "pint") {
              sls <- slctd[which(xdm[slctd, i_moy] < 0.)]
              vldNum[sls, i_moy] <- vldNum[sls, i_moy] + 1
            }
            if (varName == "bsdf") {
              sls <- slctd[which(xdm[slctd, i_moy] < 0. | xdm[slctd, i_moy] > 1.)]
              vldNum[sls, i_moy] <- vldNum[sls, i_moy] + 1
            }
          }
          totDiff[slctd] <- rowSums(cbind(totDiff[slctd], abs(var[slctd, i_moy] - xdm[slctd, i_moy])), na.rm = TRUE)
        }

        slctd <- slctd[which(totDiff[slctd] > 0.0001)]
        if (identical(integer(0), slctd)) { break }

        # Correct the value for each day
        i_doy <- 0
        for (i_moy in 1 : n_moy) {
          for (i_sls in slctd) {
            sls <- slctd[which(vldNum[i_sls, i_moy] != 0)]
            if (!identical(integer(0), sls)) {
              for (i_dom in 1 : MvIDNM[c("nrml", "leap")[i_cal], i_moy]) {
                i_doy <- MvCDNM[c("nrml", "leap")[i_cal], i_moy] + i_dom
                if (varName == "pint") {
                  sg <- sls[which(VAR[sls, i_doy] > 0.)]
                  VAR[sg, i_doy] <- VAR[sg, i_doy] + (var[sg, i_moy] - xdm[sg, i_moy])
                }
                if (varName == "bsdf") {
                  sg <- sls[which(VAR[sls, i_doy] > 0. & VAR[sls, i_doy] < 1.)]
                  VAR[sg, i_doy] <- VAR[sg, i_doy] + (var[sg, i_moy] - xdm[sg, i_moy])
                }
              }
            }
          }
        }

        i_lp <- i_lp + 1
      }

      slctd <- cal[[c("nrml", "leap")[i_cal]]]
      # Set each daily value equal to zero for which the monthly value is zero
      # Set each daily value equal to one for which the monthly value is one ('bsdf')
      i_doy <- 0
      for (i_moy in 1 : n_moy) {
        for (i_dom in 1 : MvIDNM[c("nrml", "leap")[i_cal], i_moy]) {
          i_doy <- MvCDNM[c("nrml", "leap")[i_cal], i_moy] + i_dom
          VAR[slctd[which(var[slctd, i_moy] == 0.)], i_doy] <- 0.
          if (varName == "bsdf") { VAR[slctd[which(var[slctd, i_moy] == 1.)], i_doy] <- 1. }
        }
      }

      # Set all negative values equal to zero
      VAR[slctd, ][VAR[slctd, ] < 0.] <- 0.

      # Set all values greater than one to one ('bsdf')
      if (varName == "bsdf") { VAR[slctd, ][VAR[slctd, ] > 1.] <- 1. }

    }

  }

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  # Error checking of the result
  # Error checking of the result
  metVar <- data.frame(cep = c("temperature", "precipitation", "sunshine"), mly = c("temp", "pint", "bsdf"),
                       dly = c("TEMP", "PREC", "BSDF"))
  cepN <- metVar$cep[metVar$mly == varName]
  dlyN <- metVar$dly[metVar$mly == varName]

  pImp <- which(apply(VAR, 1, function(r) {
    any(r[!is.na(r)] < varFeatures[metVar$dly[metVar$mly == varName], "lwst_val"] |
          r[!is.na(r)] > varFeatures[metVar$dly[metVar$mly == varName], "hgst_val"]) }))

  if (!identical(integer(0), pImp)) {
    warning("For ", cepN , " data, the approach chosen to generate quasi-daily data leds to physically impossible ",
            "results for certain monthly time series. For this reason, to time series of the result object '", dlyN,
            "' that affected by this problem, the NA value is assigned.")
    VAR[pImp, ] <- NA
  }

  return(VAR)

}
