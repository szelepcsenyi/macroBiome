#' Estimator for Daily Solar Irradiance/Irradiation
#'
#' @description Calculates monthly averages of daily solar irradiance/irradiation under cloudless-sky conditions, for
#'     a given latitude, elevation and year. In addition, it also optionally computes monthly means for daylength
#'     (accumulated hours of daylight), for a given latitude, elevation and year. In addition to monthly averages,
#'     daily values can also be directly access.
#'
#' @param lat 'numeric' vector with the latitude coordinates (in decimal degrees)
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
#' @param aprchTR 'character' vector of length 1 that indicates that the daily solar irradiance/irradiation is
#'     estimated with a daily or hourly time resolution. Valid values are as follows: \cr
#'     (a) \code{'daily'} - it is calculated as a daily amount (in \eqn{MJ m^{-2} dy^{-1}}{MJ m^{-2} dy^{-1}}); \cr
#'     (b) \code{'hourly'} - it is estimated as an hourly mean (in \eqn{MJ m^{-2} hr^{-1}}{MJ m^{-2} hr^{-1}}).
#' @param daily 'logical' scalar that indicates whether or not daily values should also be computed.
#' @param mlyOpVar 'character' vector of at least one length that indicates the variable(s) for which
#'     monthly time series are to be calculated. Valid values are as follows: \cr
#'     (a) \code{'R_E'} - monthly averages of daily solar irradiance/irradiation under cloudless-sky conditions
#'     (depending on the time resolution in \eqn{MJ m^{-2} dy^{-1}}{MJ m^{-2} dy^{-1}} or
#'     \eqn{MJ m^{-2} hr^{-1}}{MJ m^{-2} hr^{-1}}); \cr
#'     (b) \code{'DL'} - monthly means for daylength (accumulated hours of daylight) (in hours).
#'
#' @details To estimate the monthly averages of relative sunshine duration, Yin (1999) uses estimated values of the
#'     mean hourly solar irradiance under cloudless-sky conditions (\code{'R_E_MJ.m2.hr1'}), which is calculated via
#'     the scheme of Yin (1997b), with some minor modifications. Furthermore, in the approach proposed by Yin (1999),
#'     the monthly means of daylength (\code{'DL_moa_hr'}) is used to estimate monthly averages of daily potential
#'     evapotranspiration. The solar irradiance model presented by Yin (1997b) and slightly modified by Yin (1999) is
#'     labelled \code{'Solar123'} here. \cr
#'     The approach \code{'SPLASH'} uses a different radiation model, which is based on the procedure described by
#'     Davis et al. (2017). The approach used here differs from the base scheme only in that it takes into account
#'     changes in orbital parameters of the Earth over time. Temporal variability of orbital parameters is considered
#'     through the implementation of the procedure as proposed by Berger and Loutre (1991). In the SPLASH algorithm,
#'     the daily top-of-the-atmosphere solar radiation is computed as twice the integral of the extraterrestrial
#'     radiation flux realized from local solar noon to sunset (see Davis et al. 2017). Using this physical amount,
#'     the daily solar radiation at the surface of the Earth under clear-sky conditions (\code{'R_E_MJ.m2.dy1'})
#'     can be easily estimated. The hourly mean surface radiation is finally obtained as the quotient of this
#'     integrated value and the daylength. \cr
#'     In the both approach, the daylength is computed by using the sunset hour angle, however, due to differences in
#'     the calculation of the solar declination, the results obtained by the two procedures differ. The approach
#'     \code{'Solar123'} gives time-independent results, while using the approach \code{'SPLASH'}, varying results
#'     over time are obtained through the variable orbital parameters. \cr
#'     The two radiation models therefore differ significantly both in terms of conceptual frameworks and
#'     assumptions. The approach \code{'Solar123'} converts hourly mean values into daily amounts, while in the
#'     approach \code{'SPLASH'}, the integrated daily values are transformed to hourly averages. However, in terms of
#'     atmospheric transmittivity, both models consider the universal value of 0.75 as a basic value. To conclude,
#'     depending on the scope of application, both methods can give valid results, but for paleo-climatological and
#'     paleo-environmental studies, the approach \code{'SPLASH'} is clearly recommended.
#'
#' @return If daily values are also requested (\code{daily = TRUE}), the function returns a list of lists with daily
#'     and monthly data. If \code{daily = FALSE}, the return object is a list with the monthly means. The character
#'     vector \code{'mlyOpVar'} determines for which variables the monthly averages are calculated. Daily values also
#'     become available only for this/these variable(s). Each matrix in the list of daily data consists of 365 or 366
#'     columns, while monthly data are available, of course, as 12-column matrices. The former are accessible in the
#'     list \code{'dly'}, while the latter can be found in the list labelled as \code{'mly'}. Here, matrices with
#'     monthly data contain averages of the daily values. See the argument \code{'mlyOpVar'} for the content of
#'     matrices.
#'
#' @note As with any function with a point mode, a set of basic input data is defined here. In this case, it is as
#'     follow: \code{'lat'} (latitude coordinates in decimal degrees). The length of this vector determines the
#'     number of rows in matrices of the return list. In the case of arguments that do not affect the course of the
#'     calculation procedure or the structure of the return object, scalar values (i.e., 'numeric' vector of
#'     length 1) may also be allowed. In this case, they are as follows: \code{'elv'} (elevation in meters above sea
#'     level), and \code{'year'} (year using astronomical year numbering). These scalars are converted to vectors by
#'     the function during the error handling, and these vectors are applied in the further calculations. If these
#'     data are stored in vectors of length at least 2, their length must be the same size of first dimension of the
#'     matrix containing the basic data.
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
#' library(graphics)
#'
#' # Monthly average of the mean hourly solar irradiance under cloudless-sky conditions at sea level,
#' # in June 2000 along a meridian, according to the two different radiation models
#' lats <- seq(-90, 90, 10)
#' models <- c("Solar123", "SPLASH")
#' R_E_moa_MJ.m2.hr1 <- matrix(nrow = 2, ncol = length(lats), dimnames = list(models, lats))
#' R_E_moa_MJ.m2.hr1[1, ] <- cliAvgDlySolIrrPoints(lats, elv = 0,
#'     aprchTR = "hourly")$R_E_moa_MJ.m2.hr1[, "Jun"]
#' R_E_moa_MJ.m2.hr1[2, ] <- cliAvgDlySolIrrPoints(lats, elv = 0, aprchSIM = "SPLASH",
#'     aprchTR = "hourly")$R_E_moa_MJ.m2.hr1[, "Jun"]
#' cols <- c("black", "green")
#' matplot(t(R_E_moa_MJ.m2.hr1), type = "l", lwd = 2, col = cols, xaxt = "n",
#'     xlab = "Geographical latitude (degrees)",
#'     ylab = "Mean hourly solar irradiance (MJ m-2 hr-1)")
#' axis(1, at = seq(1, ncol(R_E_moa_MJ.m2.hr1)), labels = colnames(R_E_moa_MJ.m2.hr1))
#' legend(1, 2, legend = rownames(R_E_moa_MJ.m2.hr1), col = cols, lty = 1 : 2, lwd = 2, xpd = TRUE)
#'
#' \donttest{
#' # Daylength at latitude 75N in the year 2000, according to the two different radiation models
#' DL_hr <- matrix(nrow = 2, ncol = 366, dimnames = list(c("Solar123", "SPLASH"), seq(1, 366)))
#' DL_hr[1, ] <- cliAvgDlySolIrrPoints(75., 0., daily = TRUE, mlyOpVar = "DL")$dly$DL_hr
#' DL_hr[2, ] <- cliAvgDlySolIrrPoints(75., 0., aprchSIM = "SPLASH", daily = TRUE,
#'     mlyOpVar = "DL")$dly$DL_hr
#' cols <- c("black", "green")
#' matplot(t(DL_hr), type = "l", lwd = 2, col = cols, xaxt = "n", xlab = "Day number in the year",
#'     ylab = "Daylength (hr)")
#' axis(1, at = seq(1, ncol(DL_hr)), labels = colnames(DL_hr))
#' legend(1, 20, legend = rownames(DL_hr), col = cols, lty = 1 : 2, lwd = 2, xpd = TRUE)
#' }
#'
#' @importFrom stats complete.cases setNames
#' @importFrom strex match_arg
#' @import palinsol
#'
#' @export
#'
cliAvgDlySolIrrPoints <- function(lat, elv = NULL, year = 2000, aprchSIM = c("Solar123", "SPLASH"),
                                  aprchTR = c("daily", "hourly"), daily = FALSE, mlyOpVar = c("R_E")) {

  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  aprchSIM <- strex::match_arg(aprchSIM)
  aprchTR <- strex::match_arg(aprchTR)

  # Vectorization of scalar variables
  cv.scl <- c("elv", "year")
  if (any(sapply(cv.scl, function(x) { (length(get(x)) == 1 & is.numeric(get(x))) | identical(get(x), NA) }))) {
    lgth <- errorHandling(lat = lat)$lgth
    list2env(sapply(cv.scl, function(x) {
      if ((length(get(x)) == 1 & is.numeric(get(x))) | identical(get(x), NA)) {
        assign(x, rep(get(x), lgth)) } else { assign(x, get(x)) } },
      simplify = FALSE), envir = environment())
  }
  err_han <- errorHandling(lat = lat, elv = elv, year = year)
  list2env(Filter(Negate(is.null), err_han), envir = environment())

  mlyOpVar <- strex::match_arg(mlyOpVar, c("R_E", "DL"), several_ok = TRUE)


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
  cal <- setNames(lapply(vector("list", 2), function(x) vector()), c("nrml", "leap"))
  vld <- Reduce("&", list(do.call(complete.cases, list(lat, elv, year))))
  cal$nrml <- as.numeric(which(Reduce("&", list(vld, sapply(n_doy, function(x) { x == N_doy[1] & !is.na(x) })))))
  cal$leap <- as.numeric(which(Reduce("&", list(vld, sapply(n_doy, function(x) { x == N_doy[2] & !is.na(x) })))))

  applCal <- match(sort(unique(n_doy[!is.na(n_doy)])), N_doy)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 02. Prepare the output objects
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (daily) {
    dlyOpVarN <- if (any(mlyOpVar %in% "DL")) { c("DL_hr") } else { vector("character") }
    if (any(mlyOpVar %in% "R_E")) {
      dlyOpVarN <- if (aprchTR == "hourly") { c(dlyOpVarN, "R_E_MJ.m2.hr1") } else { c(dlyOpVarN, "R_E_MJ.m2.dy1") }
    }
    dlyVal <- setNames(replicate(length(dlyOpVarN),
                                 matrix(NA, nrow = lgth,
                                        ncol = ifelse(!identical(integer(0), applCal), max(N_doy[applCal]), 365),
                                        dimnames = list(NULL, NULL)),
                                 simplify = FALSE),
                       dlyOpVarN)
  }

  mlyOpVarN <- if (any(mlyOpVar %in% "DL")) { c("DL_moa_hr") } else { vector("character") }
  if (any(mlyOpVar %in% "R_E")) {
    mlyOpVarN <- if (aprchTR == "hourly") {
      c(mlyOpVarN, "R_E_moa_MJ.m2.hr1")
    } else {
      c(mlyOpVarN, "R_E_moa_MJ.m2.dy1")
    }
  }
  mlyVal <- setNames(replicate(length(mlyOpVarN),
                               matrix(NA, nrow = lgth, ncol = n_moy, dimnames = list(NULL, month.abb)),
                               simplify = FALSE),
                     mlyOpVarN)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 03. Run the algorithm for each calendar type
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (i_cal in applCal) {

    slctd <- cal[[c("nrml", "leap")[i_cal]]]
    tstp <- seq(1, N_doy[i_cal])

    i_doy <- 0
    for (i_moy in 1 : n_moy) { # monthly

      for (i_dom in 1 : MvIDNM[c("nrml", "leap")[i_cal], i_moy]) { # daily

        i_doy <- MvCDNM[c("nrml", "leap")[i_cal], i_moy] + i_dom

        # Calculate daily solar irradiance, and declination for a given day
        if (aprchSIM == "SPLASH") {
          dlyQty <- calcDlySolRad(year[slctd], rep(i_doy, length(slctd)), lat[slctd], argCkd = T)
          # ref: Eq 2 in Allen (1996)
          if (is.null(elv[slctd])) {
            R_E_MJ.m2.dy1 <- dlyQty$H_0_j.m2 * 10. ^ -6
          } else {
            R_E_MJ.m2.dy1 <- ((c$c + c$d) * (1. + (2.67e-5) * elv[slctd])) * dlyQty$H_0_j.m2 * 10. ^ -6
          }
          # ref: Eq 1.6.11 in Duffie and Beckman (1991)
          if (any(mlyOpVar %in% "DL") | aprchTR == "hourly") { DL_hr <- (2. / 15.) * dlyQty$h_s_deg }
          if (aprchTR == "hourly") { R_E_MJ.m2.hr1 <- ifelse(DL_hr == 0., 0., R_E_MJ.m2.dy1 / DL_hr) }
        } else {
          dlyQty <- calcDlySolIrr(year[slctd], rep(i_doy, length(slctd)), lat[slctd], elv[slctd], argCkd = T)
          R_E_MJ.m2.hr1 <- dlyQty$R_E_MJ.m2.hr1
          if (any(mlyOpVar %in% "DL") | aprchTR == "daily") { DL_hr <- dlyQty$DL_hr }
          # ref: Eq 3.2 in Yin (1997b)
          if (aprchTR == "daily") { R_E_MJ.m2.dy1 <- R_E_MJ.m2.hr1 * DL_hr }
        }

        # Save the daily values (if necessary)
        if (daily) {
          for (i_var in 1 : (length(dlyOpVarN))) {
            dlyVal[[dlyOpVarN[i_var]]][slctd, i_doy] <- get(dlyOpVarN[i_var])
          }
        }

        # Update the monthly means
        if (any(mlyOpVar %in% "R_E")) {
          if (aprchTR == "hourly") {
            mlyVal$R_E_moa_MJ.m2.hr1[slctd, i_moy] <- sapply(1 : length(slctd), function(i) {
              plus(c(mlyVal$R_E_moa_MJ.m2.hr1[slctd[i], i_moy],
                     R_E_MJ.m2.hr1[i] / MvIDNM[c("nrml", "leap")[i_cal], i_moy]))
            })
          } else {
            mlyVal$R_E_moa_MJ.m2.dy1[slctd, i_moy] <- sapply(1 : length(slctd), function(i) {
              plus(c(mlyVal$R_E_moa_MJ.m2.dy1[slctd[i], i_moy],
                     R_E_MJ.m2.dy1[i] / MvIDNM[c("nrml", "leap")[i_cal], i_moy]))
            })
          }
        }

        if (any(mlyOpVar %in% "DL")) {
          mlyVal$DL_moa_hr[slctd, i_moy] <- sapply(1 : length(slctd), function(i) {
            plus(c(mlyVal$DL_moa_hr[slctd[i], i_moy], DL_hr[i] / MvIDNM[c("nrml", "leap")[i_cal], i_moy]))
          })
        }

      } # daily

    } # monthly

  }

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  if (daily) {
    rslt <- list(dly = dlyVal, mly = mlyVal)
  } else {
    rslt <- mlyVal
  }
  return(rslt)

}


plus <- function(x) { if(all(is.na(x))) { c(x[0], NA)} else { sum(x, na.rm = TRUE) } }
