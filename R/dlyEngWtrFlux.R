#' Estimator for Daily Amounts of Energy and Water Fluxes
#'
#' @description Estimates the daily amounts of energy and water fluxes and the associated monthly bioclimatic
#'     variables, by using the SPLASH algorithm described by Davis et al. (2017). This version of the algorithm is
#'     directly suitable for paleoclimate applications because it takes into account the time variability of the
#'     Earth's orbital elements, and thus changes in the seasonal cycle of insolation.
#'
#' @param TEMP 'numeric' R object with one-year time series of daily mean air temperature (in °C)
#' @param PREC 'numeric' R object with one-year time series of daily precipitation sum (in mm)
#' @param BSDF 'numeric' R object with one-year time series of daily fractional sunshine duration (dimensionless)
#' @param lat 'numeric' vector with the latitude coordinates (in decimal degrees)
#' @param elv 'numeric' vector with the elevation values (in meters above sea level)
#' @param year 'numeric' vector with values of the year (using astronomical year numbering)
#' @param MSMC 'numeric' vector with values of the maximum soil moisture capacity (aka 'bucket size') (in mm)
#' @param daily 'logical' scalar that indicates whether or not daily values should also be computed.
#' @param mlyOpVar 'character' vector of at least one length that indicates the bioclimatic variable(s) for which
#'     monthly time series are to be calculated. Valid values are as follows: \cr
#'     (a) \code{'EET'} - monthly amounts of the equilibrium evapotranspiration (in mm); \cr
#'     (b) \code{'PET'} - monthly amounts of the potential evapotranspiration (in mm); \cr
#'     (c) \code{'AET'} - monthly amounts of the actual evapotranspiration (in mm); \cr
#'     (d) \code{'PTC'} - monthly values of the Priestley–Taylor coefficient (dimensionless); \cr
#'     (e) \code{'CWD'} - monthly values of the climatic water deficit (in mm).
#'
#' @details To estimate the daily radiation, evapotranspiration and soil moisture for an equilibrium year, the SPLASH
#'     algorithm described by Davis et al. (2017) is implemented with two slight amendments. In accordance with Davis
#'     et al. (2017), daily insolation (incoming solar radiation at the top of the atmosphere) is estimated by using
#'     Eq 1.10.3 in Duffie and Beckman (1991), with the remark that orbital parameters of the Earth are not assumed
#'     to be constant. Temporal variability of orbital parameters is considered through the implementation of the
#'     procedure as proposed by Berger and Loutre (1991). To simulate seasonal changes in the climatic water balance,
#'     the simple 'bucket model' proposed by Cramer and Prentice (1988) is applied in accordance with the SPLASH
#'     v.1.0 model. In this model, the daily value of actual evapotranspiration is estimated as an the analytical
#'     integral of the minimum of the instantaneous evaporative supply and demand rates over a single day (see Eq 27
#'     in Davis et  al. (2017)). The SPLASH algorithm is modified in a further aspect: in the 'bucket model', the
#'     'bucket size' is freely changeable, i.e., it can be specified regionally. Its value is set to 150 mm by
#'     default, in accordance with Cramer and Prentice (1988). \cr
#'     The function provides daily estimates for the following key quantities: daily insolation
#'     (\code{'H_0_J.m2.dy1'}), daily net surface radiation (\code{'H_np_J.m2.dy1'}, and \code{'H_nn_J.m2.dy1'});
#'     photosynthetic photon flux density (\code{'PPFD_mol.m2.dy1'}); daily condensation, soil moisture and runoff
#'     (\code{'CN_mm.dy1'}, \code{'SM_mm.dy1'}, and \code{'RO_mm.dy1'}, respectively); and daily equilibrium,
#'     potential, and actual evapotranspiration (\code{'EET_mm.dy1'}, \code{'PET_mm.dy1'}, and \code{'AET_mm.dy1'}).
#'     It also integrates daily data for bioclimatic variables relevant to ecoclimatological studies at a monthly
#'     timescale: monthly equilibrium, potential and actual evapotranspiration (\code{'EET_mo_mm.mo1'},
#'     \code{'PET _mo_mm.mo1'}, and \code{'AET_mo_mm.mo1'}), monthly Priestley–Taylor coefficient (\code{'PTC_mo'}),
#'     monthly climatic water deficit (\code{'CWD_mo_mm.mo1'}).
#'
#' @return If daily values are also requested (\code{daily = TRUE}), the function returns a list of lists with daily
#'     and monthly data. If \code{daily = FALSE}, the return object is a list with the monthly values. The character
#'     vector \code{'mlyOpVar'} determines for which variables are integrated at monthly scale (for explanations see
#'     'Details'). Daily data is available for the following quantities:
#'
#'     \itemize{
#'       \item{\code{H_0_J.m2.dy1}: daily solar irradiation (in J m-2)}
#'       \item{\code{H_np_J.m2.dy1}: daily positive (daytime) net surface radiation (in J m-2)}
#'       \item{\code{H_nn_J.m2.dy1}: daily negative (nighttime) net surface radiation (in J m-2)}
#'       \item{\code{PPFD_mol.m2.dy1}: daily photosynthetically active radiation (in mol m-2)}
#'       \item{\code{CN_mm.dy1}: daily condensation (in mm)}
#'       \item{\code{SM_mm.dy1}: daily soil moisture (in mm)}
#'       \item{\code{RO_mm.dy1}: daily runoff (in mm)}
#'       \item{\code{EET_mm.dy1}: daily equilibrium evapotranspiration (in mm)}
#'       \item{\code{PET_mm.dy1}: daily potential evapotranspiration (in mm)}
#'       \item{\code{AET_mm.dy1}: daily actual evapotranspiration (in mm)}
#'     }
#'
#'     Each matrix in the list of daily data consists of 365 or 366 columns, while monthly data are available, of
#'     course, as 12-column matrices. The former are accessible in the list \code{'dly'}, while the latter can be
#'     found in the list labelled as \code{'mly'}.
#'
#' @note As with any function with a point mode, a set of basic input data is defined here. In this case, they are as
#'    follows: \code{'TEMP'} (one-year time series of daily mean air temperature), \code{'PREC'} (one-year time series
#'    of daily precipitation sum), and \code{'BSDF'} (one-year time series of daily mean relative sunshine
#'    duration). The objects \code{'TEMP'}, \code{'PREC'} and \code{'BSDF'} must be either vectors of length 365 (or
#'    366) or 365-column (or 366-column) matrices. The first dimensions of these matrices have to be the same length.
#'    The function automatically converts vectors into single-row matrices during the error handling, and then uses
#'    these matrices. The first dimensions of these matrices determines the number of rows in the result matrices. In
#'    the case of arguments that do not affect the course of the calculation procedure or the structure of the
#'    return object, scalar values (i.e., 'numeric' vector of length 1) may also be allowed. In this case, they are as
#'    follows: \code{'lat'} (latitude coordinates in decimal degrees), \code{'elv'} (elevation in meters above sea
#'    level), \code{'year'} (year using astronomical year numbering), and \code{'MSMC'} ('bucket size' in mm). These
#'    scalars are converted to vectors by the function during the error handling, and these vectors are applied in
#'    the further calculations. If these data are stored in vectors of length at least 2, their length must be the
#'    same size of first dimension of the matrices containing the basic data.
#'
#' @references
#'
#' \cite{Berger A, Loutre MF (1991) Insolation values for the climate of the last 10 million years. Quat Sci Rev
#'     10(4):297-317. \doi{10.1016/0277-3791(91)90033-Q}}
#'
#' \cite{Cramer W, Prentice IC (1988) Simulation of regional soil moisture deficits on a European scale. Nor J Geogr
#'     42(2-3):149–151. \doi{10.1080/00291958808552193}}
#'
#' \cite{Davis TW, Prentice IC, Stocker BD, Thomas RT, Whitley RJ, Wang H, Evans BJ, Gallego-Sala AV, Sykes MT,
#'     Cramer W (2017) Simple process-led algorithms for simulating habitats (SPLASH v.1.0): robust indices of
#'     radiation, evapotranspiration and plant-available moisture. Geosci Model Dev 10(2):689–708.
#'     \doi{10.5194/gmd-10-689-2017}}
#'
#' \cite{Duffie JA, Beckman WA (1991) Solar Engineering of Thermal Processes. Second Edition. Wiley-Interscience,
#'     New York, NY}
#'
#' @examples
#' \donttest{
#' library(graphics)
#'
#' # Loading mandatory data for the Example 'Points'
#' data(inp_exPoints)
#'
#' with(inp_exPoints, {
#' # Estimates the daily amounts of energy and water fluxes with default settings,
#' # at a grid cell near Szeged, Hungary (46.3N, 20.2E) (for the normal period 1981-2010)
#' year <- trunc(mean(seq(1981, 2010)))
#' wea <- dlyWeaGenPoints(colMeans(temp), colMeans(prec), colMeans(bsdf), year = year)
#' ewf <- dlyEngWtrFluxPoints(wea$TEMP, wea$PREC, wea$BSDF, lat, lon, elv, year = year)
#'
#' # Check daily energy and water fluxes
#' opar <- par(no.readonly = TRUE)
#' par(mfrow = c(4, 1))
#' var <- list(t(ewf$dly$H_np_J.m2.dy1) * 1e-6, t(ewf$dly$SM_mm.dy1), t(wea$PREC))
#' lbl <- list(expression(italic(H[N])~(MJ~m^{-2})), expression(italic(SM[n])~(mm)),
#'     expression(italic(P[n])~(mm)))
#' at <- list(seq(0, 16, 4), seq(0, 80, 20), seq(0, 4))
#' txt <- list("(a)", "(b)", "(c)")
#' for (i in 1 : length(var)) {
#'   par(mar = c(1, 5, 1, 1))
#'   plot(var[[i]], type = "l", lwd = 2, xlab = NA, ylab = NA, axes = FALSE)
#'   axis(side = 1, las = 1, tck = -0.03, labels = NA, at = seq(-60, 720, 30))
#'   axis(side = 2, las = 1, tck = -0.03, labels = NA, at = at[[i]])
#'   axis(side = 2, las = 1, lwd = 0, line = -0.4, cex.axis = 1.6, at = at[[i]])
#'   mtext(side = 2, lbl[[i]], line = 3, cex = 1.1)
#'   text(-12, max(at[[i]]) / 4, txt[[i]], pos = 4, cex = 1.7)
#' }
#' par(mar = c(2, 5, 1, 1))
#' plot(t(ewf$dly$PET_mm.dy1), type = "l", lwd = 2, xlab = NA, ylab = NA, axes = FALSE,
#'   ylim = c(0, max(t(ewf$dly$PET_mm.dy1))))
#' lines(t(ewf$dly$AET_mm.dy1), lty = 2, lwd = 2, col = "green")
#' axis(side = 1, las = 1, tck = -0.03, labels = NA, at = seq(-60, 720, 30))
#' axis(side = 1, las = 1, lwd = 0, line = -0.4, at = seq(-60, 720, 30), cex.axis = 1.6)
#' axis(side = 2, las = 1, tck = -0.03, labels = NA, at = seq(-1, 6, 1))
#' axis(side = 2, las = 1, lwd = 0, line = -0.4, cex.axis = 1.6, at = seq(-1, 6, 1))
#' legend("topright", legend = c(expression(italic(E[n]^{q})), expression(italic(E[n]^{a}))),
#'     col = c("black", "green"), lty = c(1, 2), cex = 1.6, inset = 0.02,
#'     adj = c(0.5, 0.5), lwd = c(2, 2), horiz = TRUE, bty = "n", seg.len = 1)
#' mtext(side = 2, expression(italic(E[n])~(mm)), line = 3, cex = 1.1)
#' text(-12, 1.5, "(d)", pos = 4, cex = 1.7)
#' par(opar)
#'
#' # Check monthly water balance quantities
#' plot(t(ewf$mly$PET_mo_mm.mo1), type = "l", lwd = 2, ylim = c(0, 1.1 * max(ewf$mly$PET_mo_mm.mo1)),
#'     xlab = NA, ylab = NA, axes = FALSE)
#' lines(t(ewf$mly$EET_mo_mm.mo1), lty = 1, lwd = 2, col = "green")
#' lines(t(ewf$mly$AET_mo_mm.mo1), lty = 2, lwd = 2, col = "blue")
#' box(lwd = 2)
#' axis(side = 1, las = 1, tck = -0.02, labels = NA, at = seq(1, 12))
#' axis(side = 1, las = 1, lwd = 0, line = -0.4, labels = month.abb, at = seq(1, 12), cex.axis = 1.2)
#' axis(side = 2, las = 1, tck = -0.02, labels = NA, at = seq(-20, 200, 20))
#' axis(side = 2, las = 1, lwd = 0, line = -0.4, at = seq(-20, 200, 20), cex.axis = 1.2)
#' mtext(side = 2, expression(list(Evapotranspiration, mm~month^{-1})), line = 2, cex = 1.2)
#' legend("top", legend = c("Potential", "Equilibrium", "Actual"), col = c("black", "green", "blue"),
#'     lty = c(1, 1, 2), lwd = c(2, 2, 2), inset = 0.01, x.intersp = 1.1, y.intersp = 2.0,
#'     horiz = TRUE, bty = "n", cex = 1.2)
#' })
#' }
#'
#' @importFrom stats complete.cases setNames
#' @importFrom strex match_arg
#'
#' @export
#'
dlyEngWtrFluxPoints <- function(TEMP, PREC, BSDF, lat, elv, year = 2000, MSMC = 150., daily = TRUE,
                                mlyOpVar = c("EET", "PET", "AET")) {

  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  errorChecking(year = year, MSMC = MSMC)

  # Vectorization of scalar variables
  cv.scl <- c("lat", "elv", "year", "MSMC")
  if (any(sapply(cv.scl, function(x) { (length(get(x)) == 1 & is.numeric(get(x))) | identical(get(x), NA) }))) {
    lgth <- errorHandling(TEMP = TEMP, PREC = PREC, BSDF = BSDF)$lgth
    list2env(sapply(cv.scl, function(x) {
      if ((length(get(x)) == 1 & is.numeric(get(x))) | identical(get(x), NA)) {
        assign(x, rep(get(x), lgth)) } else { assign(x, get(x)) } },
      simplify = FALSE), envir = environment())
  }
  err_han <- errorHandling(TEMP = TEMP, PREC = PREC, BSDF = BSDF, lat = lat, elv = elv, year = year, MSMC = MSMC)
  list2env(Filter(Negate(is.null), err_han), envir = environment())

  cv.arg <- c("BSDF", "TEMP", "PREC", "lat", "elv")
  for (i in 1 : length(cv.arg)) {
    if (is.null(get(cv.arg[i]))) { stop("Invalid argument: '", cv.arg[i], "' is missing, with no default.") }
  }

  mlyOpVar <- strex::match_arg(mlyOpVar, c("EET", "PET", "AET", "PTC", "CWD"), several_ok = TRUE)


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
  if (length(unique(lapply(list(TEMP, PREC, BSDF), function(VAR) { dim(VAR)[2] }))) == 1) {

    # Check whether or not the first n_doy values of each storage vector are valid (notNA)
    vldIndivIn <- sapply(list(TEMP = TEMP, PREC = PREC, BSDF = BSDF), function(VAR) {
      mapply(function(x, y) { sum(is.na(x[1 : y])) }, split(VAR, row(VAR)), replace(n_doy, is.na(n_doy), 1)) == 0
    })
    vldIn <- apply(matrix(vldIndivIn, ncol = 3), 1, function(x) { Reduce("&", x) })

    # Check whether or not values after the first n_doy values of each storage vector are valid (NA)
    vldIndivOut <- sapply(list(TEMP = TEMP, PREC = PREC, BSDF = BSDF), function(VAR) {
      sapply(mapply(function(x, y) { setdiff(x, x[1 : y]) }, split(TEMP, row(TEMP)), replace(n_doy, is.na(n_doy), 1)),
             function(x) { any(identical(x, numeric(0)), is.na(x)) })
    })
    vldOut <- apply(matrix(vldIndivOut, ncol = 3), 1, function(x) { Reduce("&", x) })

    vld <- Reduce("&", list(vldIn, vldOut))

    vld <- Reduce("&", list(do.call(complete.cases, list(year, lat, elv, MSMC)), vld))
    cal$nrml <- as.numeric(which(Reduce("&", list(vld, sapply(n_doy, function(x) { x == N_doy[1] & !is.na(x) })))))
    cal$leap <- as.numeric(which(Reduce("&", list(vld, sapply(n_doy, function(x) { x == N_doy[2] & !is.na(x) })))))

  } else {
    stop("The second dimension of objects transformed to matrices that contain meteorological variables ",
         "('TEMP', 'PREC' and 'BSDF') is not the same.")
  }
  applCal <- match(sort(unique(n_doy[!is.na(n_doy)])), N_doy)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 02. Prepare the output objects
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (daily) {
    dlyOpVarN <- opVarChoices$Variable[opVarChoices$Timescale == "daily"]
  } else {
    dlyOpVarN <- c("SM_mm.dy1")
  }
  dlyVal <- setNames(replicate(length(dlyOpVarN),
                               matrix(NA, nrow = lgth,
                                      ncol = ifelse(!identical(integer(0), applCal), max(N_doy[applCal]), 365),
                                      dimnames = list(NULL, NULL)),
                               simplify = FALSE),
                     dlyOpVarN)

  mlyOpVarN <- opVarChoices$Variable[opVarChoices$Timescale == "monthly"]
  mlyVal <- setNames(replicate(length(mlyOpVarN),
                               matrix(0., nrow = lgth, ncol = n_moy, dimnames = list(NULL, month.abb)),
                               simplify = FALSE),
                     mlyOpVarN)

  if (all(sapply(cal, function(x) { identical(x, numeric(0)) })) | all(!vld)) {
    slctdMlyOpVar <- replace(paste0(mlyOpVar, "_mo_mm.mo1"),
                             paste0(mlyOpVar, "_mo_mm.mo1") == "PTC_mo_mm.mo1", "PTC_mo")
    if (daily) {
      rslt <- list(dly = dlyVal, mly = mlyVal[slctdMlyOpVar])
    } else {
      rslt <- mlyVal[slctdMlyOpVar]
    }
    return(rslt)
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 03. Run the algorithm for each calendar type
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (i_cal in applCal) {

    slctd <- cal[[c("nrml", "leap")[i_cal]]]
    tstp <- seq(1, N_doy[i_cal])

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 03.01 Spin up the soil moisture content to the model climatology
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dlyVal$SM_mm.dy1[slctd, tstp] <-
      spin_up(year[slctd], lat[slctd], elv[slctd], BSDF[slctd, tstp, drop = FALSE], TEMP[slctd, tstp, drop = FALSE],
              PREC[slctd, tstp, drop = FALSE], MSMC[slctd], argCkd = TRUE)$SM_mm.dy1

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 03.02 Run the SPLASH algorithm for a full year
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    i_doy <- 0
    for (i_moy in 1 : n_moy) { # monthly

      for (i_dom in 1 : MvIDNM[c("nrml", "leap")[i_cal], i_moy]) { # daily

        i_doy <- MvCDNM[c("nrml", "leap")[i_cal], i_moy] + i_dom

        # Calculate daily radiation, condensation, and evaporation fluxes for a given day
        W_n <- dlyVal$SM_mm.dy1[slctd, ifelse((i_doy - 1) < 1, N_doy[i_cal], (i_doy - 1))]
        dlyQty <- run_one_day(year[slctd], rep(i_doy, length(slctd)), lat[slctd], elv[slctd],
                              BSDF[slctd, , drop = FALSE][, i_doy],
                              TEMP[slctd, , drop = FALSE][, i_doy],
                              PREC[slctd, , drop = FALSE][, i_doy],
                              W_n, MSMC[slctd], argCkd = TRUE)

        # Update the daily soil moisture and save the other daily values (if necessary)
        if (daily) {
          for (i_var in 1 : (length(dlyOpVarN))) {
            dlyVal[[dlyOpVarN[i_var]]][slctd, i_doy] <- dlyQty[[dlyOpVarN[i_var]]]
          }
        } else {
          dlyVal$SM_mm.dy1[slctd, i_doy] <- dlyQty$SM_mm.dy1
        }

        # Update the monthly totals of evapotranspiration
        mlyVal$EET_mo_mm.mo1[slctd, i_moy] <- mlyVal$EET_mo_mm.mo1[slctd, i_moy] + dlyQty$EET_mm.dy1
        mlyVal$PET_mo_mm.mo1[slctd, i_moy] <- mlyVal$PET_mo_mm.mo1[slctd, i_moy] + dlyQty$PET_mm.dy1
        mlyVal$AET_mo_mm.mo1[slctd, i_moy] <- mlyVal$AET_mo_mm.mo1[slctd, i_moy] + dlyQty$AET_mm.dy1

      } # daily

      # Calculate the monthly water balance indices
      # ref: Eqs 34 and 33 in Davis et al. (2017)
      mlyVal$PTC_mo[slctd, i_moy] <- mlyVal$AET_mo_mm.mo1[slctd, i_moy] / mlyVal$EET_mo_mm.mo1[slctd, i_moy]
      mlyVal$CWD_mo_mm.mo1[slctd, i_moy] <- mlyVal$PET_mo_mm.mo1[slctd, i_moy] - mlyVal$AET_mo_mm.mo1[slctd, i_moy]

    } # monthly

  }

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  slctdMlyOpVar <- replace(paste0(mlyOpVar, "_mo_mm.mo1"),
                           paste0(mlyOpVar, "_mo_mm.mo1") == "PTC_mo_mm.mo1", "PTC_mo")
  if (daily) {
    rslt <- list(dly = dlyVal, mly = mlyVal[slctdMlyOpVar])
  } else {
    rslt <- mlyVal[slctdMlyOpVar]
  }
  return(rslt)

}


# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script contains functions for running the SPLASH algorithm for point-based data, i.e.:
#   spin_up(year, lat, elv, BSDF, TEMP, PREC, MSMC, argCkd = FALSE)
#   run_one_day(year, n, lat, elv, S_f, T_a, P_n, W_n, W_max, quick = TRUE, argCkd = FALSE)
#
#### DEFINE FUNCTIONS ################################################################################################

# ********************************************************************************************************************
# Name:     spin_up
# Inputs:   - double, year (using astronomical year numbering) (year)
#           - double, latitude, deg (lat)
#           - double, elevation, m (elv)
#           - double, one-year time series of daily fractional sunshine duration, unitless (BSDF)
#           - double, one-year time series of daily mean air temperature, deg C (TEMP)
#           - double, one-year time series of daily precipitation sum, mm (PREC)
#           - double, maximum soil moisture capacity /'bucket size'/, mm (MSMC)
#           - logical, indicates whether or not the checking and correction of arguments should be omitted (argCkd)
# Returns:  - list object (rslt)
#             $SM_mm.dy1 ... daily soil moisture, mm
#             $spinCount ... number of steps required to balance
# Features: This function updates the daily soil moisture until it reaches equilibrium.
# Depends:  - run_one_day() ... daily radiation, condensation, and evaporation fluxes
# ********************************************************************************************************************
spin_up <- function(year, lat, elv, BSDF, TEMP, PREC, MSMC, argCkd = FALSE) {

  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  err_han <- errorHandling(year = year, lat = lat, elv = elv, BSDF = BSDF, TEMP = TEMP, PREC = PREC, MSMC = MSMC,
                           onlyLgthChecking = argCkd)
  list2env(Filter(Negate(is.null), err_han), envir = environment())

  cv.arg <- c("year", "lat", "elv", "BSDF", "TEMP", "PREC", "MSMC")
  for (i in 1 : length(cv.arg)) {
    if (is.null(get(cv.arg[i]))) { stop("Invalid argument: '", cv.arg[i], "' is missing, with no default.") }
  }


  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  rslt <- list()

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 01. Select the calendar(s) applied
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  n_doy <- rep(NA, length(year))
  n_doy[!is.na(year)] <- ifelse(is.leap.year(year[!is.na(year)]), 366, 365)

  N_doy <- c(365, 366)
  cal <- setNames(lapply(vector("list", 2), function(x) vector()), c("nrml", "leap"))
  if (length(unique(lapply(list(TEMP, PREC, BSDF), function(VAR) { dim(VAR)[2] }))) == 1) {

    # Check whether or not the first n_doy values of each storage vector are valid (notNA)
    vldIndivIn <- sapply(list(TEMP = TEMP, PREC = PREC, BSDF = BSDF), function(VAR) {
      mapply(function(x, y) { sum(is.na(x[1 : y])) }, split(VAR, row(VAR)), replace(n_doy, is.na(n_doy), 1)) == 0
    })
    vldIn <- apply(matrix(vldIndivIn, ncol = 3), 1, function(x) { Reduce("&", x) })

    # Check whether or not values after the first n_doy values of each storage vector are valid (NA)
    vldIndivOut <- sapply(list(TEMP = TEMP, PREC = PREC, BSDF = BSDF), function(VAR) {
      sapply(mapply(function(x, y) { setdiff(x, x[1 : y]) }, split(TEMP, row(TEMP)), replace(n_doy, is.na(n_doy), 1)),
             function(x) { any(identical(x, numeric(0)), is.na(x)) })
    })
    vldOut <- apply(matrix(vldIndivOut, ncol = 3), 1, function(x) { Reduce("&", x) })

    vld <- Reduce("&", list(vldIn, vldOut))

    vld <- Reduce("&", list(do.call(complete.cases, list(year, lat, elv, MSMC)), vld))
    cal$nrml <- as.numeric(which(Reduce("&", list(vld, sapply(n_doy, function(x) { x == N_doy[1] & !is.na(x) })))))
    cal$leap <- as.numeric(which(Reduce("&", list(vld, sapply(n_doy, function(x) { x == N_doy[2] & !is.na(x) })))))

  } else {
    stop("The second dimension of objects transformed to matrices that contain meteorological variables ",
         "('TEMP', 'PREC' and 'BSDF') is not the same.")
  }
  applCal <- match(sort(unique(n_doy[!is.na(n_doy)])), N_doy)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 02. Set the auxiliary vectors and result objects
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rslt$SM_mm.dy1 <- matrix(nrow = lgth, ncol = ifelse(!identical(integer(0), applCal), max(N_doy[applCal]), 365))
  rslt$spinCount <- rep(NA, lgth)
  start_SM <- end_SM <- diff_SM <- rep(NA, lgth)
  if (identical(integer(0), applCal) | all(!vld)) { return(rslt) }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 03. Run the algorithm for each calendar type
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (i_cal in applCal) {

    slctd <- cal[[c("nrml", "leap")[i_cal]]]
    tstp <- seq(1, N_doy[i_cal])

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 03.01 Set an empty bucket for each day (the daily soil moisture equal to zero)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rslt$SM_mm.dy1[slctd, tstp] <- matrix(0., nrow = length(slctd), ncol = N_doy[i_cal])

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 03.02 Run the model for one whole year
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (i_doy in tstp) {
      rslt$SM_mm.dy1[slctd, i_doy] <-
        run_one_day(year[slctd], rep(i_doy, length(slctd)), lat[slctd], elv[slctd],
                    BSDF[slctd, , drop = FALSE][, i_doy],
                    TEMP[slctd, , drop = FALSE][, i_doy],
                    PREC[slctd, , drop = FALSE][, i_doy],
                    rslt$SM_mm.dy1[slctd, , drop = FALSE][, ifelse(i_doy == 1, N_doy[i_cal], i_doy - 1)],
                    MSMC[slctd], quick = TRUE, argCkd = TRUE)$SM_mm.dy1
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 03.03 Calculate the difference in soil moisture between the initial and the final conditions
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    start_SM[slctd] <- rslt$SM_mm.dy1[slctd, 1]
    end_SM[slctd] <- run_one_day(year[slctd], rep(1, length(slctd)), lat[slctd], elv[slctd],
                                 BSDF[slctd, , drop = FALSE][, 1],
                                 TEMP[slctd, , drop = FALSE][, 1],
                                 PREC[slctd, , drop = FALSE][, 1],
                                 rslt$SM_mm.dy1[slctd, , drop = FALSE][, N_doy[i_cal]],
                                 MSMC[slctd], quick = TRUE, argCkd = TRUE)$SM_mm.dy1
    diff_SM[slctd] <- abs(end_SM[slctd] - start_SM[slctd])

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 03.04 Equilibrate this difference
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rslt$spinCount[slctd] <- 1
    slctd <- which(diff_SM > 0.1)

    while (!identical(integer(0), slctd)) {

      # Run the model again
      for (i_doy in tstp) {
        rslt$SM_mm.dy1[slctd, i_doy] <-
          run_one_day(year[slctd], rep(i_doy, length(slctd)), lat[slctd], elv[slctd],
                      BSDF[slctd, , drop = FALSE][, i_doy],
                      TEMP[slctd, , drop = FALSE][, i_doy],
                      PREC[slctd, , drop = FALSE][, i_doy],
                      rslt$SM_mm.dy1[slctd, , drop = FALSE][, ifelse(i_doy == 1, N_doy[i_cal], i_doy - 1)],
                      MSMC[slctd], quick = TRUE, argCkd = TRUE)$SM_mm.dy1
      }

      # Calculate the difference again
      start_SM[slctd] <- rslt$SM_mm[slctd, 1]
      end_SM[slctd] <- run_one_day(year[slctd], rep(1, length(slctd)), lat[slctd], elv[slctd],
                                   BSDF[slctd, , drop = FALSE][, 1],
                                   TEMP[slctd, , drop = FALSE][, 1],
                                   PREC[slctd, , drop = FALSE][, 1],
                                   rslt$SM_mm.dy1[slctd, , drop = FALSE][, N_doy[i_cal]],
                                   MSMC[slctd], quick = TRUE, argCkd = TRUE)$SM_mm.dy1
      diff_SM[slctd] <- abs(end_SM[slctd] - start_SM[slctd])

      rslt$spinCount[slctd] <- rslt$spinCount[slctd] + 1
      slctd <- which(diff_SM > 0.1)
    }

  }

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  return(rslt)

}


# ********************************************************************************************************************
# Name:     run_one_day
# Inputs:   - double, year (using astronomical year numbering) (year)
#           - double, day of the year (n)
#           - double, latitude, deg (lat)
#           - double, elevation, m (elv)
#           - double, daily fractional sunshine duration, unitless (S_f)
#           - double, daily mean air temperature, deg C (T_a)
#           - double, daily precipitation sum, mm (P_n)
#           - double, previous day's soil moisture, mm (W_n)
#           - double, maximum soil moisture capacity /'bucket size'/, mm (W_max)
#           - logical, indicates whether or not daily values should be included the results (quick)
#           - logical, indicates whether or not the checking and correction of arguments should be omitted (argCkd)
# Returns:  - list object (rslt)
#             $H_0_J.m2.dy1 ...... daily solar irradiation, J m-2
#             $H_np_J.m2.dy1 ..... daily positive (daytime) net surface radiation, J m-2
#             $H_nn_J.m2.dy1 ..... daily negative (nightime) net surface radiation, J m-2
#             $PPFD_mol.m2.dy1 ... daily photosynthetically active radiation, mol m-2
#             $CN_mm.dy1 ......... daily condensation, mm
#             $EET_mm.dy1 ........ daily equilibrium evapotranspiration, mm
#             $PET_mm.dy1 ........ daily potential evapotranspiration, mm
#             $AET_mm.dy1 ........ daily actual evapotranspiration, mm
#             $SM_mm.dy1 ......... daily soil moisture, mm
#             $RO_mm.dy1 ......... daily runoff, mm
# Features: This function runs the SPLASH algorithm at a single location for a single day.
# Depends:  - c$S_c ............. supply rate constant, mm hr-1
#           - calcDlyEvapot() ... daily radiation, condensation, and evaporation fluxes
# Ref:      - Davis TW, Prentice IC, Stocker BD, Thomas RT, Whitley RJ, Wang H, Evans BJ, Gallego-Sala AV, Sykes MT,
#             Cramer W (2017) Simple process-led algorithms for simulating habitats (SPLASH v.1.0): robust indices of
#             radiation, evapotranspiration and plant-available moisture. Geosci Model Dev 10(2):689–708.
#             DOI: 10.5194/gmd-10-689-2017
# ********************************************************************************************************************
run_one_day <- function(year, n, lat, elv, S_f, T_a, P_n, W_n, W_max, quick = FALSE, argCkd = FALSE) {

  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  err_han <- errorHandling(year = year, n = n, lat = lat, elv = elv, S_f = S_f, T_a = T_a, P_n = P_n, W_n = W_n,
                           W_max = W_max, onlyLgthChecking = argCkd)
  list2env(Filter(Negate(is.null), err_han), envir = environment())

  cv.arg <- c("year", "n", "lat", "elv", "S_f", "T_a", "P_n", "W_n", "W_max")
  for (i in 1 : length(cv.arg)) {
    if (is.null(get(cv.arg[i]))) { stop("Invalid argument: '", cv.arg[i], "' is missing, with no default.") }
  }

  cv.opv <- c("SM_mm.dy1", "RO_mm.dy1")
  if (!quick) {
    cv.opv <- c(cv.opv, c("H_0_J.m2.dy1", "H_np_J.m2.dy1", "H_nn_J.m2.dy1", "PPFD_mol.m2.dy1", "CN_mm.dy1",
                          "EET_mm.dy1", "PET_mm.dy1", "AET_mm.dy1"))
  }

  incl <- which(complete.cases(year, n, lat, elv, S_f, T_a, P_n, W_n, W_max))


  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  rslt <- setNames(lapply(vector("list", length(cv.opv)), function(x) rep(NA, lgth)), cv.opv)
  if (identical(integer(0), incl)) {
    return(rslt)
  } else {
    for (i in 1 : length(cv.arg)) { assign(cv.arg[i], get(cv.arg[i])[incl]) }
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 01. Calculate the evaporative supply rate, mm hr-1
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq 21 in Davis et al. (2017)
  S_w <- c$S_c * W_n / W_max

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 02. Compute daily radiation, condensation, and evaporation fluxes
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ET <- calcDlyEvapot(year, n, lat, elv, S_f, T_a, S_w, argCkd = TRUE)
  if (!quick) {
    rslt$H_0_J.m2.dy1[incl] <- ET$H_0_j.m2
    rslt$H_np_J.m2.dy1[incl] <- ET$H_np_j.m2
    rslt$H_nn_J.m2.dy1[incl] <- ET$H_nn_j.m2
    rslt$PPFD_mol.m2.dy1[incl] <- ET$Q_n_mol.m2
    rslt$CN_mm.dy1[incl] <- ET$Cond_mm
    rslt$EET_mm.dy1[incl] <- ET$EET_mm
    rslt$PET_mm.dy1[incl] <- ET$PET_mm
    rslt$AET_mm.dy1[incl] <- ET$AET_mm
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 03. Update the daily soil moisture, mm
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq 1 i Davis et al. (2017)
  SM <- W_n + P_n + ET$Cond_mm - ET$AET_mm

  SM <- ifelse(SM > W_max, W_max,   # The bucket is full. Set the soil moisture to the maximum capacity.
               ifelse(SM < 0., 0.,   # The bucket is empty. Set the soil moisture equal to zero
                      SM))

  RO <- ifelse(SM > W_max, SM - W_max,   # The bucket is full. Add the remaining water to the runoff
               0.)   # The bucket is empty. Set the runoff equal to zero

  if (!quick) {
    rslt$AET_mm.dy1[incl] <- ifelse(SM < 0., rslt$AET_mm.dy1[incl] + SM, rslt$AET_mm.dy1[incl])
  }

  rslt$SM_mm.dy1[incl] <- SM
  rslt$RO_mm.dy1[incl] <- RO

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  return(rslt)

}
