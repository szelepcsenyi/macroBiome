# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script contains functions to calculate daily irradiance, i.e.:
#   calcDtMOptAM(lat, DCL, Z_0, Z_n, t_1, t_80, argCkd = FALSE)
#   calcDlySolIrr(year, n, lat, elv = NULL, aprchI_0 = c("Yin1999", "Iqbal1983"), aprchDCL = c("Brock1981", "CBM"),
#      argCkd = FALSE)
#
#### DEFINE FUNCTIONS ################################################################################################

# ********************************************************************************************************************
# Name:     calcDlySolIrr
# Inputs:   - double, year (using astronomical year numbering) (year)
#           - double, day of the year (n)
#           - double, latitude, deg (lat)
#           - double, elevation, m (elv)
#           - character, name of the approach for calculating values of the solar constant (aprchI_0)
#           - character, name of the approach for calculating values of the solar declination (aprchDCL)
#           - logical, indicates whether or not the checking and correction of arguments should be omitted (argCkd)
# Returns:  - list object (l.SolIrr)
#             $I_0_MJ.m2.hr1 ... solar constant, MJ m-2 hr-1
#             $DCL_deg ......... solar declination, deg
#             $Z_0_deg ......... solar zenith angle at solar noon, deg
#             $Z_n_deg ......... solar zenith angle at sunset (or at solar midnight when the sun does not rise), deg
#             $t_1_hr .......... sunrise/sunset hour (solar hour at sunrise and sunset), hr
#             $t_80_hr ......... solar hour corresponding to the 80° zenith angle, hr
#             $m_a ............. daytime mean optical air mass, unitless
#             $tau ............. transmission coefficient, unitless
#             $f_m ............. an exponent depending on optical air mass, unitless
#             $cosZ_a .......... diurnal average of the cosine of the solar zenith angle, deg
#             $R_E_MJ.m2.hr1 ... mean hourly solar irradiance under cloudless-sky conditions, MJ m-2 hr-1
#             $DL_hr ........... daylength, hr
# Features: This function calculates daily solar irradiance data.
# Depends:  - calcDtMOptAM() ... daytime mean optical air mass
# Ref:      - Brock TD (1981) Calculating solar radiation for ecological studies. Ecol Model 14(1–2):1-19.
#             DOI: 10.1016/0304-3800(81)90011-9
#           - Duffie JA, Beckman WA (1980) Solar Engineering of Thermal Processes. Wiley-Interscience, New York, NY
#           - Forsythe WC, Rykiel Jr. EJ, Stahl RS, Wu H, Schoolfield RM (1995) A model comparison for daylength as a
#             function of latitude and day of year. Ecol Model 80(1):87–95. DOI: 10.1016/0304-3800(94)00034-F
#           - Iqbal M (1983) An Introduction to Solar Radiation. Academic Press, London, UK.
#             DOI: 10.1016/B978-0-12-373750-2.X5001-0
#           - Linacre E (1992) Climate, Data and Resources. A Reference and Guide. Routledge, New York, NY
#           - Yin X (1997a) Calculating daytime mean relative air mass. Agric For Meteorol 87(2-3):85-90.
#             DOI: 10.1016/S0168-1923(97)00029-4
#           - Yin X (1997b) Optical Air Mass: Daily Integration and its Applications. Meteorol Atmos Phys
#             63(3-4):227-233. DOI: 10.1007/BF01027387
#           - Yin X (1999) Bright Sunshine Duration in Relation to Precipitation, Air Temperature and Geographic
#             Location. Theor Appl Climatol 64(1–2):61–68. DOI: 10.1007/s007040050111
# ********************************************************************************************************************
calcDlySolIrr <- function(year, n, lat, elv = NULL, aprchI_0 = c("Yin1999", "Iqbal1983"),
                          aprchDCL = c("Brock1981", "CBM"), argCkd = FALSE) {

  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  aprchI_0 <- strex::match_arg(aprchI_0)
  aprchDCL <- strex::match_arg(aprchDCL)

  err_han <- errorHandling(year = year, n = n, lat = lat, elv = elv, onlyLgthChecking = argCkd)
  list2env(Filter(Negate(is.null), err_han), envir = environment())

  cv.arg <- c("year", "n", "lat")
  for (i in 1 : length(cv.arg)) {
    if (is.null(get(cv.arg[i]))) { stop("Invalid argument: '", cv.arg[i], "' is missing, with no default.") }
  }

  cv.opv <- c("I_0_MJ.m2.hr1", "DCL_deg", "Z_0_deg", "Z_n_deg", "t_1_hr", "t_80_hr", "m_a", "tau", "f_m", "cosZ_a",
              "R_E_MJ.m2.hr1", "DL_hr")

  if (missing(elv)) {
    incl <- which(stats::complete.cases(year = year, n = n, lat = lat))
    elv <- NULL
  } else {
    incl <- which(stats::complete.cases(year = year, n = n, lat = lat, elv = elv))
    cv.arg <- c(cv.arg, c("elv"))
  }


  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  l.SolIrr <- stats::setNames(lapply(vector("list", length(cv.opv)), function(x) rep(NA, lgth)), cv.opv)
  if (identical(integer(0), incl)) {
    return(l.SolIrr)
  } else {
    for (i in 1 : length(cv.arg)) { assign(cv.arg[i], get(cv.arg[i])[incl]) }
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 01. Calculate the solar constant (I_0), MJ m-2 hr-1
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (aprchI_0 == "Iqbal1983") {    # ref: Iqbal (1983)

    # 'Average' solar constant, I_SC (in MJ m-2 hr-1)
    I_SC <- 4.871    # ref: 'text' in Iqbal (1983)

    # Correction factor for the eccentricity of Earth’s orbit (unitless)
    # ref: Eq 1.2.3 in Iqbal (1983); Eq 1.4.1 in Duffie and Beckman (1980)
    E_0 <- 1. + (0.033 * cos(360. * n / 365. * pi /180.))

    # Solar constant, I_0 (in MJ m-2 hr-1)
    I_0 <- I_SC * E_0    # ref: Eq 4.2.1 in Iqbal (1983)

  } else {                          # ref: Yin (1999); Brock (1981)

    # 'Average' solar constant, I_SC (in MJ m-2 hr-1)
    I_SC <- 4.9212    # ref: 'text' in Yin (1997b); Yin (1999)

    # Radius vector of the Earth, R_1 (unitless)
    # ref: Eq 2 in Brock (1981)
    R_1 <- 1. / (1. + (0.033 * cos(360. * (n - 1.) / 365. * pi / 180.)))**0.5

    # Solar constant, I_0 (in MJ m-2 hr-1)
    I_0 <- I_SC /R_1 ** 2.    # ref: 'text' in Brock (1981)

  }
  l.SolIrr$I_0_MJ.m2.hr1[incl] <- I_0

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 02. Calculate the solar declination (DCL), deg
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (aprchDCL == "CBM") {    # ref: Eqs 1 and 2 in Forsythe et al. (1995)
    DCL <- rad2deg(asin(0.39795 * cos(0.2163108 + 2. * atan(0.9671396 * tan(0.00860 * (n - 186.))))))
  } else {                    # ref: Eq A7 in Yin (1997a); Eq 1 in Brock (1981); Eq 4 in Forsythe et al. (1995)
    DCL <- 23.45 * sin(deg2rad(360. / 365. * (n + 283.)))
  }
  l.SolIrr$DCL_deg[incl] <- DCL

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 03. Calculate the solar zenith angle at solar noon (Z_0), deg
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq 2.10 in Yin (1997b)
  Z_0 <- ifelse(abs(lat - DCL) < 90., lat - DCL, 90.)
  l.SolIrr$Z_0_deg[incl] <- Z_0

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 04. Calculate the solar zenith angle at sunset (or at solar midnight when the sun does not rise) (Z_n), deg
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq 2.9 in Yin (1997b)
  Z_n <- ifelse(abs(lat + DCL) < 90., 90.,
                rad2deg(acos(sin(deg2rad(DCL)) * sin(deg2rad(lat)) - cos(deg2rad(DCL)) * cos(deg2rad(lat)))))
  l.SolIrr$Z_n_deg[incl] <- Z_n

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 05. Calculate the daytime mean optical air mass (and the value of an associated experiential function)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Coefficents for calculating optical air mass
  # ref: Eqs 2.2 and 2.3 in Yin (1997b)
  c_LE80 <- matrix(c(0.008307, 1.021, -0.01259), nrow = 1, dimnames = list("c_LE80", c("c1", "c2", "c3")))
  c_GT80 <- matrix(c(0.03716, 1.538, -1.689), nrow = 1, dimnames = list("c_GT80", c("c1", "c2", "c3")))
  coef_cho <- rbind(c_LE80, c_GT80)

  # Optical air mass at solar noon (corresponding to cosine zenith angle at its maximum), m_o (unitless)
  # ref: Eqs 2.2 and 2.3 in Yin (1997b)
  coef_val <- matrix(coef_cho[ifelse(Z_0 <= 80., "c_LE80", "c_GT80"), ], nrow = length(incl), ncol = ncol(coef_cho),
                     dimnames = list(NULL, colnames(coef_cho)))
  # ref: Eq 2.1 in Yin (1997b)
  m_o <- as.numeric(coef_val[, "c2"] / (coef_val[, "c1"] + cos(deg2rad(Z_0))) + coef_val[, "c3"])

  # Optical air mass corresponding to cosine zenith angle at its medium, m_c (unitless)
  # ref: definitions for Eq 3.3 in Yin (1997b)
  Z_c <- rad2deg(acos((cos(deg2rad(Z_0)) + cos(deg2rad(Z_n))) / 2.))
  # ref: Eqs 2.2 and 2.3 in Yin (1997b)
  coef_val <- matrix(coef_cho[ifelse(Z_c <= 80., "c_LE80", "c_GT80"), ], nrow = length(incl), ncol = ncol(coef_cho),
                     dimnames = list(NULL, colnames(coef_cho)))
  # ref: Eq 2.1 in Yin (1997b)
  m_c <- as.numeric(coef_val[, "c2"] / (coef_val[, "c1"] + cos(deg2rad(Z_c))) + coef_val[, "c3"])

  # Optical air mass corresponding to cosine zenith angle at its bottom-quarter range point, m_l (unitless)
  # ref: definitions for Eq 3.3 in Yin (1997b)
  Z_l <- rad2deg(acos((cos(deg2rad(Z_0)) + 3. * cos(deg2rad(Z_n))) / 4.))
  # ref: Eqs 2.2 and 2.3 in Yin (1997b)
  coef_val <- matrix(coef_cho[ifelse(Z_l <= 80., "c_LE80", "c_GT80"), ], nrow = length(incl), ncol = ncol(coef_cho),
                     dimnames = list(NULL, colnames(coef_cho)))
  # ref: Eq 2.1 in Yin (1997b)
  m_l <- as.numeric(coef_val[, "c2"] / (coef_val[, "c1"] + cos(deg2rad(Z_l))) + coef_val[, "c3"])

  # Sunrise/sunset hour (solar hour at sunrise and sunset), t_1 (hr)
  # ref: Eq 2.5 in Yin (1997b)
  # t_1 <- ifelse(abs(lat - DCL) >= 90., 0.,
  #               ifelse(abs(lat + DCL) >= 90., 12.,
  #                      (12. / pi) * acos(-1. * tan(deg2rad(lat)) * tan(deg2rad(DCL)))))
  t_1 <- vector(length = length(incl))
  t_1[abs(lat - DCL) >= 90. & abs(lat + DCL) < 90.] <- 0.
  t_1[abs(lat - DCL) < 90. & abs(lat + DCL) >= 90.] <- 12.
  sls <- abs(lat - DCL) < 90. & abs(lat + DCL) < 90.
  if (!identical(sls, FALSE)) { t_1[sls] <- (12. / pi) * acos(-1. * tan(deg2rad(lat[sls])) * tan(deg2rad(DCL[sls]))) }
  l.SolIrr$t_1_hr[incl] <- t_1

  # Solar hour corresponding to the 80° zenith angle, t_80 (hr)
  # ref: Eq 2.8 in Yin (1997b)
  ratio <- (cos(deg2rad(80.)) - (sin(deg2rad(lat)) * sin(deg2rad(DCL)))) / cos(deg2rad(lat)) * cos(deg2rad(DCL))
  t_80 <- vector(length = length(incl))
  t_80[ratio >= 1.] <- 0.
  t_80[ratio <= -1.] <- 12.
  t_80[ratio < 1. & ratio > -1.] <- (1. / 15.) * rad2deg(acos(ratio[ratio < 1. & ratio > -1.]))
  l.SolIrr$t_80_hr[incl] <- t_80

  # Daytime mean optical air mass, m_a (unitless)
  # ref: 2.7 Eq in Yin (1997b)
  # print(c("lat", "DCL", "Z_0", "Z_n", "t_1", "t_80"))
  # print(c(lat, DCL, Z_0, Z_n, t_1, t_80))
  m_a <- calcDtMOptAM(lat, DCL, Z_0, Z_n, t_1, t_80, argCkd = TRUE)

  # Checking whether use of elevation-correction for the optical air mass values
  if (!is.null(elv)) {
    # Station atmospheric pressure, sap (unitless)
    # ref: Eq 2.11 in Yin (1997b); Linacre (1992)
    sap <- exp(-elv / 8000.)

    m_o <- m_o * sap
    m_c <- m_c * sap
    m_l <- m_l * sap
    m_a <- m_a * sap
  }
  l.SolIrr$m_a[incl] <- m_a

  # Transmission coefficient, tau (unitless)
  # ref: 'text' Yin (1997b), Linacre (1992)
  tau <- 0.75
  l.SolIrr$tau[incl] <- tau

  # An exponent depending on optical air mass, f_m (unitless)
  # ref: Eq 3.3 in Yin (1997b)
  f_m <- 0.2323 + 0.01452 * (m_a + m_l) * exp(1.403 * tau) - 0.1528 * m_o + m_c + 0.487 * (m_c - m_l)
  l.SolIrr$f_m[incl] <- f_m

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 06. Calculate the diurnal average of the cosine of the solar zenith angle
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Half-length of the illuminated day, Hr (in radians)
  # ref: Eq A2 in Yin (1997a)
  Hr <- vector(length = length(incl))
  Hr[abs(lat - DCL) >= 90. & abs(lat + DCL) < 90.] <- 0.
  Hr[abs(lat - DCL) < 90. & abs(lat + DCL) >= 90.] <- pi
  sls <- abs(lat - DCL) < 90. & abs(lat + DCL) < 90.
  if (!identical(sls, FALSE)) { Hr[sls] <- acos(-1. * tan(deg2rad(lat[sls])) * tan(deg2rad(DCL[sls]))) }

  # Diurnal average of the cosine of the solar zenith angle, cosZ_a (unitless)
  # ref: Eq A5 in Yin (1997a)
  cosZ_a <- vector(length = length(incl))
  cosZ_a[abs(lat - DCL) >= 90.] <- 0.
  sls <- abs(lat - DCL) < 90.
  if (!identical(sls, FALSE)) { cosZ_a[sls] <- sin(deg2rad(lat[sls])) * sin(deg2rad(DCL[sls])) +
    (1. / Hr[sls]) * cos(deg2rad(lat[sls])) * cos(deg2rad(DCL[sls])) * sin(Hr[sls]) }
  l.SolIrr$cosZ_a[incl] <- cosZ_a

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 07. Calculate the mean hourly solar irradiance under cloudless-sky conditions (R_E), MJ m-2 hr-1
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: 'text' in Yin (1999); Eq 3.2 in Yin (1997b)
  R_E <- I_0 * cosZ_a * (tau ** f_m)
  l.SolIrr$R_E_MJ.m2.hr1[incl] <- R_E

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 08. Calculate the daylength (DL), hr
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  DL <- vector(length = length(incl))
  DL[abs(lat - DCL) >= 90. & abs(lat + DCL) < 90.] <- 0.
  DL[abs(lat - DCL) < 90. & abs(lat + DCL) >= 90.] <- 24.
  sls <- abs(lat - DCL) < 90. & abs(lat + DCL) < 90.
  if (aprchDCL == "CBM") {  # ref: Eqs 3 in Forsythe et al. (1995)
    DL[sls] <- 24. - (24. / pi) * acos((sin(deg2rad(0.8333)) + sin(deg2rad(lat[sls])) * sin(deg2rad(DCL[sls]))) /
                                         (cos(deg2rad(lat[sls])) * cos(deg2rad(DCL[sls]))))
  } else {                  # ref: 'text' in Yin (1997a); Eq 4 in Brock (1981); Eqs 5 and 6 in Forsythe et al. (1995)
    DL[sls] <- 24. / 180. * rad2deg(acos(-1. * tan(deg2rad(lat[sls])) * tan(deg2rad(DCL[sls]))))
  }
  l.SolIrr$DL_hr[incl] <- DL

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  return(l.SolIrr)

}


# ********************************************************************************************************************
# Name:     calcDtMOptAM
# Inputs:   - double, latitude, deg (lat)
#           - double, solar declination, deg (DCL)
#           - double, solar zenith angle at solar noon, deg (Z_0)
#           - double, solar zenith angle at sunset (or at solar midnight when the sun does not rise), deg (Z_n)
#           - double, sunrise/sunset hour (solar hour at sunrise and sunset), hr (t_1)
#           - double, solar hour corresponding to the 80° zenith angle, hr (t_80)
#           - logical, indicates whether or not the checking and correction of arguments should be omitted (argCkd)
# Returns:  - double, daytime mean optical air mass, unitless (m_a)
# Features: This function returns the daytime mean optical air mass by using the Equation 2.7 in Yin (1997b).
# Ref:      - Yin X (1997b) Optical Air Mass: Daily Integration and its Applications. Meteorol Atmos Phys
#             63(3-4):227-233. DOI: 10.1007/BF01027387
# ********************************************************************************************************************
calcDtMOptAM <- function(lat, DCL, Z_0, Z_n, t_1, t_80, argCkd = FALSE) {

  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

  cv.arg <- c("lat", "DCL", "Z_0", "Z_n", "t_1", "t_80")
  for (i in 1 : length(cv.arg)) {
    if (is.null(get(cv.arg[i]))) { stop("Invalid argument: '", cv.arg[i], "' is missing, with no default.") }
  }

  if (!argCkd) {

    maVarFeatures <- data.frame(lgth = rep(1, 6), lwst_val = c(rep(-90, 2), rep(0, 4)),
                                hgst_val = c(rep(90, 4), rep(24, 2)),
                                rng_wrd = c(rep("values from the interval [-90, 90]", 2),
                                            rep("values from the interval [0, 90]", 2),
                                            rep("values from the interval [0, 24]", 2)))
    rownames(maVarFeatures) <- cv.arg

    # Checking the class and content of objects
    for (i in 1 : length(cv.arg)) {
      arg <- get(cv.arg[i])
      if (!is.list(arg) & length(arg) >= 1L &
          length(which(!is.na(match(class(arg), c("numeric", "integer"))))) >= 1L) {
        if (!(is.numeric(arg) &
            all(arg[!is.na(arg)] >= maVarFeatures[cv.arg[i], "lwst_val"] &
                arg[!is.na(arg)] <= maVarFeatures[cv.arg[i], "hgst_val"]))) {
          stop("Invalid argument: '", cv.arg[i], "' has to be a vector with ", maVarFeatures[cv.arg[i], "rng_wrd"],
               ".")
        }
      } else {
        stop("Invalid argument: '", cv.arg[i], "' has to be a vector with ", maVarFeatures[cv.arg[i], "rng_wrd"],
             ".")
      }
    }

    ls.arg <- stats::setNames(list(lat, DCL, Z_0, Z_n, t_1, t_80), cv.arg)
    if (!all(sapply(ls.arg, length) == length(ls.arg[[1]]))) {
      stop("Invalid argument: The vectors ",
           paste0(paste(cv.arg[1 : (length(cv.arg) - 1)], sep = "", collapse = ", "), " and ",
                  cv.arg[length(cv.arg)]),
           " have not the same length.")
    }

  }



  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

  # Coefficents for calculating optical air mass
  # ref: Eqs 2.2 and 2.3 in Yin (1997b)
  c_LE80 <- matrix(c(0.008307, 1.021, -0.01259), nrow = 1, dimnames = list("c_LE80", c("c1", "c2", "c3")))
  c_GT80 <- matrix(c(0.03716, 1.538, -1.689), nrow = 1, dimnames = list("c_GT80", c("c1", "c2", "c3")))

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 01. Select the equation(s) applied
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  eq <- stats::setNames(lapply(vector("list", 4), function(x) vector()), c("a", "b", "c", "d"))
  eq$a <- which(t_1 == 0.)
  eq$b <- which(t_1 != 0. & Z_n <= 80.)
  eq$c <- which(t_1 != 0. & Z_0 >= 80.)
  eq$d <- which(t_1 != 0. & Z_n > 80. & Z_0 < 80.)

  applEq <- which(lapply(eq, function(x) { !identical(x, integer(0)) } ) == T)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 02. Set the result vector
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  m_a <- rep(NA, length(lat))

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 03. Run the algorithm for each equation type
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (i_eq in applEq) {

    slctd <- eq[[c("a", "b", "c", "d")[i_eq]]]

    sn <- 0
    while(sn <= 2) {

      if (c("a", "b", "c", "d")[i_eq] == "a") {
        m_a[slctd] <- 39.7
        break
      }

      if (c("a", "b", "c", "d")[i_eq] == "b") {
        t_1g <- t_1[slctd]
        cg <- c_LE80
      }

      if (c("a", "b", "c", "d")[i_eq] == "c") {
        t_1g <- t_1[slctd]
        cg <- c_GT80
      }

      if (c("a", "b", "c", "d")[i_eq] == "d") {
        sn <- sn + 1

        if (sn == 1) t_1g <- t_80[slctd]
        if (sn == 2) t_1g <- t_1[slctd]
        if (sn == 3) t_1g <- t_80[slctd]

        if (sn == 1) cg <- c_LE80
        if (sn == 2) cg <- c_GT80
        if (sn == 3) cg <- c_GT80
      }

      a_ar <- cg[1] + sin(deg2rad(lat[slctd])) * sin(deg2rad(DCL[slctd]))
      if (sn < 2) {
        b_ar <- cos(deg2rad(lat[slctd])) * cos(deg2rad(DCL[slctd]))
      }

      mbr2 <- mbr3 <- mbr4 <- rep(NA, length(slctd))

      sls <- a_ar > b_ar
      if (!identical(sls, integer(0))) {
        mbr2[sls] <- cg[2] / sqrt(a_ar[sls] ** (2.) - b_ar[sls] ** (2.))
        numer <- b_ar[sls] + a_ar[sls] * cos(deg2rad(15. * t_1g[sls]))
        denom <- a_ar[sls] + b_ar[sls] * cos(deg2rad(15. * t_1g[sls]))
        mbr3[sls] <- acos(numer / denom)
      }

      sls <- a_ar < b_ar
      if (!identical(sls, integer(0))) {
        mbr2[sls] <- cg[2] / sqrt(b_ar[sls] ** (2.) - a_ar[sls] ** (2.))
        sqr1 <- sqrt((b_ar[sls] + a_ar[sls]) * (1. + cos(deg2rad(15. * t_1g[sls]))))
        sqr2 <- sqrt((b_ar[sls] - a_ar[sls]) * (1. - cos(deg2rad(15. * t_1g[sls]))))
        mbr3[sls] <- log((sqr1 + sqr2) / (sqr1 - sqr2))
      }

      sls <- a_ar == b_ar
      if (!identical(sls, integer(0))) {
        mbr2[sls] <- cg[2] / a_ar[sls]
        mbr3[sls] <- tan(deg2rad(15. * t_1g[sls] / 2.))
      }

      mbr4 <- cg[3] * t_1g

      F_eq6 <- 1. / deg2rad(15.) * mbr2 * mbr3 + mbr4

      if (sn == 0) { F0_eq6 <- F_eq6 }
      if (sn == 1) { F1_eq6 <- F_eq6 }
      if (sn == 2) { F2_eq6 <- F_eq6 }
      if (sn == 3) { F3_eq6 <- F_eq6 }

      if (sn == 1 | sn == 2) { next }

      if (sn == 0) {
        m_a[slctd] <- 1. / t_1g * F0_eq6
        break
      }

      if (sn == 3) {
        m_a[slctd] <- 1. / t_1[slctd] * (F1_eq6 + F2_eq6 - F3_eq6)
        break
      }

    }

  }

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  return(m_a)

}


deg2rad <- function(deg) { deg * (pi / 180.) }
rad2deg <- function(rad) { rad * (180. / pi) }
