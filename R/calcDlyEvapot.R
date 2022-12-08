# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script contains functions to calculate daily radiation, condensation, and evapotranspiration, i.e.:
#   calcDlyEvapot(year, n, lat, elv, S_f, T_a, S_w = NULL, argCkd = FALSE)
#   elv2prs(elv)
#   calcSlpSatPrs(T_a)
#   calcLtHtVap(T_a)
#   calcDenWtr(T_a, P_atm, argCkd = FALSE)
#   calcPsychroConst(T_a, P_atm, argCkd = FALSE)
#   calcSpHtCap(T_a, argCkd = FALSE)
#
#### DEFINE FUNCTIONS ################################################################################################

# ********************************************************************************************************************
# Name:     calcDlyEvapot
# Inputs:   - double, year (using astronomical year numbering) (year)
#           - double, day of the year (n)
#           - double, latitude, deg (lat)
#           - double, elevation, m (elv)
#           - double, daily fractional sunshine duration, unitless (S_f)
#           - double, daily mean air temperature, deg C (T_a)
#           - double, evaporative supply rate, mm hr-1 (S_w)
#           - logical, indicates whether or not the checking and correction of arguments should be omitted (argCkd)
# Returns:  - list object (l.Evapot)
#             $nu_deg ........ true anomaly, deg
#             $lam_deg ....... true longitude, deg
#             $d_r ........... distance factor, unitless
#             $dcl_deg ....... declination angle, deg
#             $h_s_deg ....... sunset hour angle, deg
#             $H_0_j.m2 ...... daily top-of-the-atmosphere solar radiation, J m-2
#             $tau ........... atmospheric transmittivity, unitless
#             $Q_n_mol.m2 .... daily photosynthetic photon flux density, mol m-2
#             $I_LW_w.m2 ..... net long-wave radiation flux, W m-2
#             $h_n_deg ....... net radiation cross-over angle, deg
#             $H_np_j.m2 ..... daily daytime net radiation, J m-2
#             $H_nn_j.m2 ..... daily nighttime net radiation, J m-2
#             $E_con_m3.j .... water-to-energy conversion factor, m3 J-1
#             &P_atm_pa ...... atmospheric pressure at a given elevation, Pa
#             &s_pa.k......... slope of saturation vapor pressure-temperature curve, Pa K-1
#             &L_v_j.kg ...... latent heat of vaporization of water, J kg-1
#             &rho_w_kg.m3 ... density of water, kg m-3
#             &gam_pa.k ...... psychrometric constant, Pa K-1
#             $Cond_mm ....... daily condensation, mm
#             $EET_mm ........ daily equilibrium evapotranspiration, mm
#             $PET_mm ........ daily potential evapotranspiration, mm
#             $AET_mm ........ daily actual evapotranspiration, mm
# Features: This function calculates daily radiation, condensation, and evaporation fluxes.
# Depends:  - c$ome ................ entrainment factor, unitless
#           - calcDlySolRad() ...... daily radiation fluxes
#           - elv2prs() ............ elevation-dependent atmospheric pressure
#           - calcSlpSatPrs() ...... slope of the saturation pressure temperature curve
#           - calcLtHtVap() ........ latent heat of vaporization
#           - calcDenWtr() ......... density of water
#           - calcPsychroConst() ... psychrometric constant
# ********************************************************************************************************************
calcDlyEvapot <- function(year, n, lat, elv, S_f, T_a, S_w = NULL, argCkd = FALSE) {
  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  err_han <- errorHandling(year = year, n = n, lat = lat, elv = elv, S_f = S_f, T_a = T_a, S_w = S_w,
                           onlyLgthChecking = argCkd)
  list2env(Filter(Negate(is.null), err_han), envir = environment())

  cv.arg <- c("year", "n", "lat", "elv", "S_f", "T_a")
  for (i in 1 : length(cv.arg)) {
    if (is.null(get(cv.arg[i]))) { stop("Invalid argument: '", cv.arg[i], "' is missing, with no default.") }
  }

  cv.opv <- c("H_0_j.m2", "H_np_j.m2", "H_nn_j.m2", "Q_n_mol.m2", "P_atm_pa", "s_pa.k", "L_v_j.kg", "rho_w_kg.m3",
              "gam_pa.k", "E_con_m3.j", "Cond_mm", "EET_mm", "PET_mm")

  if (missing(S_w)) {
    incl <- which(stats::complete.cases(year = year, n = n, lat = lat, elv = elv, S_f = S_f, T_a = T_a))
    S_w <- NULL
  } else {
    cv.opv <- c(cv.opv, "AET_mm")
    incl <- which(stats::complete.cases(year = year, n = n, lat = lat, elv = elv, S_f = S_f, T_a = T_a, S_w = S_w))
    cv.arg <- c(cv.arg, c("S_w"))
  }


  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  l.Evapot <- stats::setNames(lapply(vector("list", length(cv.opv)), function(x) rep(NA, lgth)), cv.opv)
  if (identical(integer(0), incl)) {
    return(l.Evapot)
  } else {
    for (i in 1 : length(cv.arg)) { assign(cv.arg[i], get(cv.arg[i])[incl]) }
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 01. Calculate daily radiation fluxes
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  l.SolRad <- calcDlySolRad(year, n, lat, elv, S_f, T_a, argCkd = TRUE)
  if (!missing(S_w)) {
    r_u <- l.SolRad$r_u
    r_v <- l.SolRad$r_v
    r_w <- l.SolRad$r_w_W.m2
    I_LW <- l.SolRad$I_LW_w.m2
    h_n <- l.SolRad$h_n_deg
  }
  H_np <- l.SolRad$H_np_j.m2
  H_nn <- l.SolRad$H_nn_j.m2
  l.Evapot$H_0_j.m2[incl]  <- l.SolRad$H_0_j.m2
  l.Evapot$H_np_j.m2[incl] <- l.SolRad$H_np_j.m2
  l.Evapot$H_nn_j.m2[incl] <- l.SolRad$H_nn_j.m2
  l.Evapot$Q_n_mol.m2[incl] <- l.SolRad$Q_n_mol.m2

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 02. Calculate the water-to-energy conversion factor (E_con), m3 J-1
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Atmospheric pressure at a given elevation, Pa
  P_atm <- elv2prs(elv)
  l.Evapot$P_atm_pa[incl] <- P_atm

  # Slope of saturation vapor pressure-temperature curve, Pa K-1
  s <- calcSlpSatPrs(T_a)
  l.Evapot$s_pa.k[incl] <- s

  # Latent heat of vaporization of water, J kg-1
  L_v <- calcLtHtVap(T_a)
  l.Evapot$L_v_j.kg[incl] <- L_v

  # Density of water, kg m-3
  rho_w <- calcDenWtr(T_a, P_atm)
  l.Evapot$rho_w_kg.m3[incl] <- rho_w

  # Psychrometric constant, Pa K-1
  gam <- calcPsychroConst(T_a, P_atm, argCkd = TRUE)
  l.Evapot$gam_pa.k[incl] <- gam

  # Water-to-energy conversion factor, m3 J-1
  # ref: Eq 19 in Davis et al. (2017)
  E_con <- s / (L_v * rho_w * (s + gam))
  l.Evapot$E_con_m3.j[incl] <- E_con

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 03. Calculate the daily condensation (C_n), mm
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq 18 in Davis et al. (2017)
  C_n <- (1e3) * E_con * abs(H_nn)
  l.Evapot$Cond_mm[incl] <- C_n

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 04. Estimate the daily equilibrium evapotranspiration (EET_d), mm
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq 24 in Davis et al. (2017)
  EET_d <- (1e3) * E_con * H_np
  l.Evapot$EET_mm[incl] <- EET_d

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 05. Estimate the daily potential evapotranspiration (PET_d), mm
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq 25 in Davis et al. (2017)
  PET_d <- (1. + c$ome) * EET_d
  l.Evapot$PET_mm[incl] <- PET_d

  if (!missing(S_w)) {
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 06. Calculate an intermediate variable (r_x), mm m2 W-1 hr-1
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ref: 'text' in Davis et al. (2017)
    r_x <- (3.6e6) * (1. + c$ome) * E_con

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 07. Calculate the intersection hour angle (h_i), deg
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ref: Eq 28 in Davis et al. (2017)
    cos_h_i <- S_w / (r_x * r_w * r_v) + I_LW / (r_w * r_v) - r_u / r_v
    h_i <- vector(length = length(cos_h_i))
    h_i[cos_h_i >= 1.] <- 0.   # Supply is in excess of demand during the entire day
    h_i[cos_h_i <= -1.] <- 180.   # Supply limits demand during the entire day
    h_i[cos_h_i < 1. & cos_h_i > -1.] <- rad2deg(acos(cos_h_i[cos_h_i < 1. & cos_h_i > -1.]))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 08. Estimate the daily actual evapotranspiration (AET_d), mm
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ref: Eq 27b in Davis et al. (2017)
    AET_d <- (24. / pi) * (S_w * deg2rad(h_i) + r_x * r_w * r_v * (sin(deg2rad(h_n)) - sin(deg2rad(h_i))) +
                             (r_x * r_w * r_u - r_x * I_LW) * deg2rad(h_n - h_i))
    l.Evapot$AET_mm[incl] <- AET_d
  }

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  return(l.Evapot)
}


# ********************************************************************************************************************
# Name:     elv2prs
# Inputs:   - double, elevation, m (elv)
# Returns:  - double, atmospheric pressure, Pa
# Features: This function calculates the atmospheric pressure at a given elevation.
# Depends:  - c$g ..... standard gravity, m s-2
#           - c$L ..... temperature lapse rate, K m-1
#           - c$M_a ... molecular weight of dry air, kg mol-1
#           - c$P_0 ... standard sea-level pressure, Pa
#           - c$R ..... universal gas constant, J mol-1 K-1
#           - c$T_0 ... base temperature, deg C
# Ref:      - Allen RG, Pereira LS, Raes D, Smith M (1998) FAO Irrigation and Drainage Paper No. 56. Technical Report.
#             Food and Agriculture Organization of the United Nations, Rome, Italy
# ********************************************************************************************************************
elv2prs <- function(elv) {
  # ref: Eq 7 in Allen et al. (1998), Eq 20 in Davis et al. (2017)
  c$P_0 * (1. - c$L * elv / c$T_0) ^ (c$g * c$M_a /(c$R * c$L))
}


# ********************************************************************************************************************
# Name:     calcSlpSatPrs
# Inputs:   - double, air temperature, deg C (T_a)
# Returns:  - double, rate of increase of saturated vapor pressure with temperature, Pa K-1
# Features: This function calculates the slope of saturation vapor pressure-temperature curve.
# Ref:      - Prentice IC, Sykes MT, Cramer W (1993) A simulation model for the transient effects of climate change on
#             forest landscapes. Ecol Model 65(1-2):51-70. DOI: 10.1016/0304-3800(93)90126-D
#           - Allen RG, Pereira LS, Raes D, Smith M (1998) FAO Irrigation and Drainage Paper No. 56. Technical Report.
#             Food and Agriculture Organization of the United Nations, Rome, Italy
# ********************************************************************************************************************
calcSlpSatPrs <- function(T_a) {
  # ref: Eq 6 in Prentice et al. (1993); Eq 13 in Allen et al. (1998); Eq B1 in Davis et al. (2017)
  (17.269) * (237.3) * (610.78) * exp(17.269 * T_a / (237.3 + T_a)) / (237.3 + T_a) ^ 2.
}


# ********************************************************************************************************************
# Name:     calcLtHtVap
# Inputs:   - double, air temperature, deg C (T_a)
# Returns:  - double, temperature-dependent enthalpy of vaporization, J kg-1
# Features: This function calculates the latent heat of vaporization of water.
# Ref:      - Henderson‐Sellers B (1984) A new formula for latent heat of vaporization of water as a function of
#             temperature. Q J R Meteorol Soc 11(466):1186-1190. DOI: 10.1002/qj.49711046626
# ********************************************************************************************************************
calcLtHtVap <- function(T_a) {
  # ref: Eq 8 in Henderson-Sellers (1984); Eq B2 in Davis et al. (2017)
  1.91846e6 * ((T_a + 273.15) / (T_a + 273.15 - 33.91)) ^ 2.
}


# ********************************************************************************************************************
# Name:     calcDenWtr
# Inputs:   - double, air temperature, deg C (T_a)
#           - double, atmospheric pressure, Pa (P_atm)
#           - logical, indicates whether or not the checking and correction of arguments should be omitted (argCkd)
# Returns:  - double, temperature- and pressure-dependent density of pure water, kg m-3
# Features: This function calculates the temperature- and pressure-dependent density of pure water.
# Ref:      - Chen C‐T, Fine RA, Millero FJ (1977) The equation of state of pure water determined from sound speeds.
#             J Chem Phys 66(5):2142-2144. DOI: 10.1063/1.434179
# ********************************************************************************************************************
calcDenWtr <- function(T_a, P_atm, argCkd = FALSE) {
  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  if (!argCkd) {
    err_han <- errorHandling(T_a = T_a, P_atm = P_atm)
    list2env(Filter(Negate(is.null), err_han), envir = environment())
  }

  cv.arg <- c("T_a", "P_atm")
  for (i in 1 : length(cv.arg)) {
    if (is.null(get(cv.arg[i]))) { stop("Invalid argument: '", cv.arg[i], "' is missing, with no default.") }
  }


  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 01. Calculate the density of water at 1 atm, g cm-3
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq 2 in Chen et al. (1977); Eq B4 in Davis et al. (2017)
  rho_0 <- 0.99983952 + (6.788260e-5) * T_a + (-9.08659e-6) * T_a ^ 2. + (1.022130e-7) * T_a ^ 3. +
    (-1.35439e-9) * T_a ^ 4. + (1.471150e-11) * T_a ^ 5. + (-1.11663e-13) * T_a ^ 6. + (5.044070e-16) * T_a ^ 7. +
    (-1.00659e-18) * T_a ^ 8.

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 02. Calculate the bulk modulus of water at 1 atm, bar
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq 3 in Chen et al. (1977); Eq B5 in Davis et al. (2017)
  K_0 <- 19652.17 + 148.1830 * T_a + (-2.29995) * T_a ^ 2. + 0.01281 * T_a  ^ 3. + (-4.91564e-5) * T_a  ^ 4. +
    (1.035530e-7) * T_a  ^ 5.

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 03. Calculate temperature-dependent coefficients, unitless and bar-1
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eqs 7 and 8 in Chen et al. (1977); Eqs B6 and B7 in Davis et al. (2017)
  C_A <- 3.26138 + (5.223e-4) * T_a + (1.324e-4) * T_a ^ 2. + (-7.655e-7) * T_a ^ 3. + (8.584e-10) * T_a ^ 4.

  C_B <- (7.2061e-5) + (-5.8948e-6) * T_a + (8.69900e-8) * T_a ^ 2. + (-1.0100e-9) * T_a ^ 3. +
    (4.3220e-12) * T_a ^ 4.

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 04. Convert the pressure value from pascals to bars (1 bar = 100 000 Pa)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  P_atm_bar <- (1e-5) * P_atm

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 05. Calculate the density of water based on the temperature and pressure, kg m-3
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq B3 in Davis et al. (2017)
  rho_w <- (1e3) * rho_0 * (K_0 + C_A * P_atm_bar + C_B * P_atm_bar ^ 2.) /
    (K_0 + C_A * P_atm_bar + C_B * P_atm_bar ^ 2. - P_atm_bar)

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  return(rho_w)
}


# ********************************************************************************************************************
# Name:     calcPsychroConst
# Inputs:   - double, air temperature, deg C (T_a)
#           - double, atmospheric pressure, Pa (P_atm)
#           - logical, indicates whether or not the checking and correction of arguments should be omitted (argCkd)
# Returns:  - double, temperature- and pressure-dependent psychrometric constant, Pa K-1
# Features: This function calculates the temperature- and pressure-dependent psychrometric constant.
# Depends:  - c$M_a ... molecular weight of dry air, kg mol-1
#           - c$M_v ... molecular weight of water vapor, kg mol-1
#           - calcLtHtVap() ... latent heat of vaporization
#           - calcSpHtCap() ... specific heat capacity of humid air
# Ref:      - Allen RG, Pereira LS, Raes D, Smith M (1998) FAO Irrigation and Drainage Paper No. 56. Technical Report.
#             Food and Agriculture Organization of the United Nations, Rome, Italy
# ********************************************************************************************************************
calcPsychroConst <- function(T_a, P_atm, argCkd = FALSE) {
  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  if (!argCkd) {
    err_han <- errorHandling(T_a = T_a, P_atm = P_atm)
    list2env(Filter(Negate(is.null), err_han), envir = environment())
  }

  cv.arg <- c("T_a", "P_atm")
  for (i in 1 : length(cv.arg)) {
    if (is.null(get(cv.arg[i]))) { stop("Invalid argument: '", cv.arg[i], "' is missing, with no default.") }
  }


  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

  # Calculate the specific heat capacity of humid air, J kg-1 K-1
  C_p <- calcSpHtCap(T_a, argCkd = TRUE)

  # Calculate the latent heat of vaporization of water, J kg-1
  L_v <- calcLtHtVap(T_a)

  # Calculate the psychrometric constant, Pa K-1
  # ref: Eq 3.10 in Allen et al. (1998); Eq B8 in Davis et al. (2017)
  gam <- C_p * c$M_a * P_atm / (c$M_v * L_v)

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  return(gam)
}


# ********************************************************************************************************************
# Name:     calcSpHtCap
# Inputs:   - double, air temperature, deg C (T_a)
#           - logical, indicates whether or not the checking and correction of the argument should be omitted (argCkd)
# Returns:  - double, specific heat capacity of humid air, J kg-1 K-1
# Features: This function calculates the specific heat capacity of humid air.
# Ref:      - Tsilingiris PT (2008) Thermophysical and transport properties of humid air at temperature range between
#             0 and 100 °C. Energy Conversion Manage 49(5):1098-1110. DOI: 10.1016/j.enconman.2007.09.015
# ********************************************************************************************************************
calcSpHtCap <- function(T_a, argCkd = FALSE) {
  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  if (!argCkd) {
    err_han <- errorHandling(T_a = T_a)
    list2env(Filter(Negate(is.null), err_han), envir = environment())
  }

  cv.arg <- c("T_a")
  for (i in 1 : length(cv.arg)) {
    if (is.null(get(cv.arg[i]))) { stop("Invalid argument: '", cv.arg[i], "' is missing, with no default.") }
  }


  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  T_a <- ifelse(T_a < 0., 0., ifelse(T_a > 100., 100., T_a))

  # ref: Eq 47 in Tsilingiris (2008);  Eq B9 in Davis et al. (2017)
  C_p <- 1.0045714270 + (2.050632750e-3) * T_a + (-1.631537093e-4) * T_a ^ 2. + (6.212300300e-6) * T_a ^ 3. +
    (-8.830478888e-8) * T_a ^ 4. + (5.071307038e-10) * T_a  ^ 5.
  C_p <- (1e3) * C_p

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  return(C_p)
}
