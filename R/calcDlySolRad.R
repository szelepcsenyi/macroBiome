# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script contains functions to calculate daily radiation, i.e.:
#   calcOrbPos(year, n, argCkd = FALSE)
#   calcDlySolRad(year, n, lat, elv = NULL, S_f = NULL, T_a = NULL, argCkd = FALSE)
#
#### DEFINE FUNCTIONS ################################################################################################

# ********************************************************************************************************************
# Name:     calcDlySolRad
# Inputs:   - double, year (using astronomical year numbering) (year)
#           - double, day of the year (n)
#           - double, latitude, deg (lat)
#           - double, elevation, m (elv)
#           - double, daily fractional sunshine duration, unitless (S_f)
#           - double, daily mean air temperature, deg C (T_a)
#           - logical, indicates whether or not the checking and correction of arguments should be omitted (argCkd)
# Returns:  - list object (l.SolRad)
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
# Features: This function calculates daily radiation fluxes.
# Depends:  - c$A ............ empirical constant for the net long-wave radiation flux, unitless
#           - c$bet_sw ....... shortwave albedo, unitless
#           - c$bet_vis ...... visible light albedo, unitless
#           - c$b ............ empirical constant for the net long-wave radiation flux, unitless
#           - c$c ............ empirical constant for the net short-wave radiation flux, unitless
#           - c$d ............ empirical constant for the net short-wave radiation flux, unitless
#           - c$fFEC ......... flux-to-energy conversion factor, umol J-1
#           - c$Gsc .......... solar constant, W m-2
#           - calcOrbPos() ... orbital position of the Earth
# Ref:      - Allen RG (1996) Assessing Integrity of Weather Data for Reference Evapotranspiration Estimation. J Irrig
#             Drain Eng 122(2):97-106. DOI: 10.1061/(ASCE)0733-9437(1996)122:2(97)
#           - Berger AL (1978) Long-term variations of daily insolation and Quaternary climatic changes. J Atmos Sci
#             35(12):2361-2367. DOI: 10.1175/1520-0469(1978)035<2362:ltvodi>2.0.co;2
#           - Berger A, Loutre MF (1991) Insolation values for the climate of the last 10 million years. Quat Sci Rev
#             10(4):297-317. DOI: 10.1016/0277-3791(91)90033-Q
#           - Berger A, Loutre M‐F, Tricot C (1993) Insolation and Earth's orbital period. J Geophys Res Atmos
#             98(D6):10341-10362. DOI: 10.1029/93JD00222
#           - Duffie JA, Beckman WA (1991) Solar Engineering of Thermal Processes. Second Edition.
#             Wiley-Interscience, New York, NY
#           - Linacre ET (1968) Estimating the net-radiation flux. Agric Meteorol 5(1):49-63.
#             DOI: 10.1016/0002-1571(68)90022-8
#           - Woolf HM (1968) On the computation of solar evaluation angles and the determination of sunrise and
#             sunset times, Technical Report. NASA-TM-X-164, National Aeronautics and Space Administration,
#             Washington, DC
# ********************************************************************************************************************
calcDlySolRad <- function(year, n, lat, elv = NULL, S_f = NULL, T_a = NULL, argCkd = FALSE) {
  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  err_han <- errorHandling(year = year, n = n, lat = lat, elv = elv, S_f = S_f, T_a = T_a, onlyLgthChecking = argCkd)
  list2env(Filter(Negate(is.null), err_han), envir = environment())

  cv.arg <- c("year", "n", "lat")
  for (i in 1 : length(cv.arg)) {
    if (is.null(get(cv.arg[i]))) { stop("Invalid argument: '", cv.arg[i], "' is missing, with no default.") }
  }

  cv.opv <- c("nu_deg", "lam_deg", "d_r", "dcl_deg", "r_u", "r_v", "h_s_deg", "H_0_j.m2")

  if (any(missing(elv), missing(S_f), missing(T_a))) {
    incl <- which(stats::complete.cases(year = year, n = n, lat = lat))
    elv <- S_f <- T_a <- NULL
  } else {
    cv.opv <- c(cv.opv, c("tau", "Q_n_mol.m2", "I_LW_w.m2", "r_w_W.m2", "h_n_deg", "H_np_j.m2", "H_nn_j.m2"))
    incl <- which(stats::complete.cases(year = year, n = n, lat = lat, elv = elv, S_f = S_f, T_a = T_a))
    cv.arg <- c(cv.arg, c("elv", "S_f", "T_a"))
  }


  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  l.SolRad <- stats::setNames(lapply(vector("list", length(cv.opv)), function(x) rep(NA, lgth)), cv.opv)
  if (identical(integer(0), incl)) {
    return(l.SolRad)
  } else {
    for (i in 1 : length(cv.arg)) { assign(cv.arg[i], get(cv.arg[i])[incl]) }
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 01. Determine the orbital position of the Earth (at a paleoclimatic scale)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Berger and Loutre (1991), Berger (1978)
  orbPos <- calcOrbPos(year, n, argCkd = FALSE)
  e_pc   <- orbPos$e_pc
  eps_pc <- orbPos$eps_pc
  nu     <- orbPos$nu
  lam    <- orbPos$lam
  l.SolRad$nu_deg[incl]  <- nu
  l.SolRad$lam_deg[incl] <- lam

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 02. Calculate the distance factor (d_r), unitless
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq 3 in Berger et al. (1993); Eq 3 in Davis et al. (2017)
  d_r  <- ((1. + e_pc * cos(deg2rad(nu))) / (1. - e_pc ^ 2.)) ^ 2.
  l.SolRad$d_r[incl] <- d_r

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 03. Calculate the declination angle (dcl), deg
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq 1.2 in Woolf (1968); Eq 5 in Davis et al. (2017)
  dcl <- rad2deg(asin(sin(deg2rad(lam)) * sin(deg2rad(eps_pc))))
  l.SolRad$dcl_deg[incl] <- dcl

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 04. Calculate intermediate variables (r_u and r_v), unitless
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: 'text' in Davis et al. (2017)
  r_u <- sin(deg2rad(dcl)) * sin(deg2rad(lat))
  r_v <- cos(deg2rad(dcl)) * cos(deg2rad(lat))
  l.SolRad$r_u[incl] <- r_u
  l.SolRad$r_v[incl] <- r_v

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 05. Calculate the sunset hour angle (h_s), deg
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Note: (r_u / r_v) == (tan(dcl) * tan(lat))
  # ref: Eq 8 in Davis et al. (2017)
  cos_h_s <- -1. * r_u / r_v
  h_s <- vector(length = length(cos_h_s))
  h_s[cos_h_s >= 1.] <- 0.   # Polar night (no sunrise)
  h_s[cos_h_s <= -1.] <- 180.   # Polar day (no sunset)
  h_s[cos_h_s < 1. & cos_h_s > -1.] <- rad2deg(acos(cos_h_s[cos_h_s < 1. & cos_h_s > -1.]))
  l.SolRad$h_s_deg[incl] <- h_s

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 06. Calculate the integrated daily extraterrestrial radiation on a horizontal surface (H_0), J m-2
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq 1.10.3 in Duffie and Beckman (1991); Eq 7 in Davis et al. (2017)
  H_0 <- (86400. / pi) * c$Gsc * d_r * (r_u * deg2rad(h_s) + r_v * sin(deg2rad(h_s)))
  l.SolRad$H_0_j.m2[incl] <- H_0

  if (!any(missing(elv), missing(S_f), missing(T_a))) {
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 07. Calculate the atmospheric transmittivity (tau), unitless
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ref:  Eq 11 in Linacre (1968); Eq 2 in Allen (1996); Eqs 11 and 12 in Davis et al. (2017)
    tau <- (c$c + c$d * S_f) * (1. + (2.67e-5) * elv)
    l.SolRad$tau[incl] <- tau

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 08. Calculate the daily photosynthetically active radiation (Q_n), mol m-2
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ref: Eq 17 in Davis et al. (2017)
    Q_n <- (1e-6) * c$fFEC * (1. - c$bet_vis) * tau * H_0
    l.SolRad$Q_n_mol.m2[incl] <- Q_n

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 09. Estimate the net long-wave radiation flux (I_LW), W m-2
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ref: Eq 13 in Davis et al. (2017)
    I_LW <- (c$b + (1. - c$b) * S_f) * (c$A - T_a)
    l.SolRad$I_LW_w.m2[incl] <- I_LW

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 10. Calculate an intermediate variable (r_w), W m-2
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ref: 'text' in Davis et al. (2017)
    r_w <- (1. - c$bet_sw) * tau * c$Gsc * d_r
    l.SolRad$r_w_W.m2[incl] <- r_w

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 11. Calculate the net radiation cross-over angle (h_n), deg
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ref: Eq 15 in Davis et al. (2017)
    cos_h_n <- (I_LW - r_w * r_u) / (r_w * r_v)
    h_n <- vector(length = length(cos_h_n))
    h_n[cos_h_n >= 1.] <- 0.   # Net radiation is negative all day
    h_n[cos_h_n <= -1.] <- 180.   # Net radiation is positive all day
    h_n[cos_h_n < 1. & cos_h_n > -1.] <- rad2deg(acos(cos_h_n[cos_h_n < 1. & cos_h_n > -1.]))
    l.SolRad$h_n_deg[incl] <- h_n

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 12. Calculate the positive (daytime) net surface radiation (H_np), J m-2
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ref: Eq 14 in Davis et al. (2017)
    H_np <- (86400. / pi) * (deg2rad(h_n) * (r_w * r_u - I_LW) + r_w * r_v * sin(deg2rad(h_n)))
    l.SolRad$H_np_j.m2[incl] <- H_np

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 13. Calculate the negative (nighttime) net surface radiation (H_nn), J m-2
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ref: Eq 16 in Davis et al. (2017)
    H_nn <- (86400./ pi) * (r_w * r_v * (sin(deg2rad(h_s)) - sin(deg2rad(h_n))) + r_w * r_u * deg2rad(h_s - h_n) -
                              I_LW * (pi - deg2rad(h_n)))
    l.SolRad$H_nn_j.m2[incl] <- H_nn

  }

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  return(l.SolRad)
}


# ********************************************************************************************************************
# Name:     calcOrbPos
# Inputs:   - double, year (using astronomical year numbering) (year)
#           - double, day of the year (n)
#           - logical, indicates whether or not the checking and correction of arguments should be omitted (argCkd)
# Returns:  - list object (l.OrbPos)
#             $nu ......... true anomaly, deg
#             $lam ........ true longitude, deg
#             $e_pc ....... eccentricity of Earth’s elliptical orbit, unitless
#             $eps_pc ..... obliquity of Earth's elliptical orbit, deg
#             $w_ome_pc ... longitude of perihelion, deg
# Features: This function returns the heliocentric coordinates (true anomaly and true longitude) of the Earth for a
#           given day, and the orbital parameters of the Earth for a given year/epoch.
# Ref:      - Berger AL (1978) Long-term variations of daily insolation and Quaternary climatic changes. J Atmos Sci
#             35(12):2361-2367. DOI: 10.1175/1520-0469(1978)035<2362:ltvodi>2.0.co;2
#           - Berger A, Loutre MF (1991) Insolation values for the climate of the last 10 million years. Quat Sci Rev
#             10(4):297-317. DOI: 10.1016/0277-3791(91)90033-Q
# ********************************************************************************************************************
calcOrbPos <- function(year, n, argCkd = FALSE) {
  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  err_han <- errorHandling(year = year, n = n, onlyLgthChecking = argCkd)
  list2env(Filter(Negate(is.null), err_han), envir = environment())

  cv.arg <- c("year", "n")
  for (i in 1 : length(cv.arg)) {
    if (is.null(get(cv.arg[i]))) { stop("Invalid argument: '", cv.arg[i], "' is missing, with no default.") }
  }

  cv.opv <- c("e_pc", "eps_pc", "w_ome_pc", "nu", "lam")

  incl <- which(stats::complete.cases(year, n))

  for (i in 1 : length(cv.arg)) { assign(cv.arg[i], get(cv.arg[i])[incl]) }


  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  l.OrbPos <- stats::setNames(lapply(vector("list", length(cv.opv)), function(x) rep(NA, lgth)), cv.opv)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 01. Calculate the number of days in a given year (N), days
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  N <- ifelse(is.leap.year(year), 366, 365)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 02. Calculate orbital parameters (at a paleoclimatic scale)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Berger and Loutre (1991)
  orbParam <- do.call(rbind, lapply(year - 1950, ber90, degree = TRUE))
  e_pc <- as.numeric(orbParam[, "ecc"])
  eps_pc <- as.numeric(orbParam[, "eps"])
  w_ome_pc <- as.numeric(orbParam[, "varpi"])
  l.OrbPos$e_pc[incl] <- e_pc
  l.OrbPos$eps_pc[incl] <- eps_pc
  l.OrbPos$w_ome_pc[incl] <- w_ome_pc

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 03. Calculate intermediate variables
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  e_pc2 <- e_pc ^ 2.
  e_pc3 <- e_pc ^ 3.
  bet <- sqrt(1. - e_pc ^ 2.)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 04. Calculate the mean longitude of the vernal equinox, deg
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq A1 in Davis et al. (2017)
  lam_m0 <- rad2deg(2. * ((e_pc / 2. + e_pc3 / 8.) * (1. + bet) * sin(deg2rad(w_ome_pc)) -
                            e_pc2/4. * (0.5 + bet) * sin(2. * deg2rad(w_ome_pc)) +
                            e_pc3 / 8. * (1. / 3. + bet) * sin(3. * deg2rad(w_ome_pc))))

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 05. Calculate the mean longitude for a day of the year, deg
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq A2 in Davis et al. (2017)
  lam_m <- lam_m0 + (n - 80.) * (360. / N)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 06. Calculate the mean anomaly, deg
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq A3 in Davis et al. (2017)
  nu_m <- lam_m - w_ome_pc

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 07. Calculate the true anomaly (uncorrected), deg
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq A4 in Davis et al. (2017)
  nu <- rad2deg(deg2rad(nu_m) + (2. * e_pc - e_pc3 / 4.) * sin(deg2rad(nu_m)) +
                  5. / 4. * e_pc2 * sin(2. * deg2rad(nu_m)) +
                  13. / 12.* e_pc3 * sin(3. * deg2rad(nu_m)))

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 08. Calculate the true longitude, deg
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq A5 in Davis et al. (2017)
  lam <- nu + w_ome_pc
  lam <- ifelse(lam < 0., lam + 360., ifelse(lam > 360., lam - 360., lam))

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 09. Correct the true anomaly, deg
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq A5 in Davis et al. (2017)
  nu <- lam - w_ome_pc
  nu <- ifelse(nu < 0., nu + 360, nu)

  l.OrbPos$nu[incl] <- nu
  l.OrbPos$lam[incl] <- lam

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  return(l.OrbPos)
}


deg2rad <- function(deg) { deg * (pi / 180.) }
rad2deg <- function(rad) { rad * (180. / pi) }
is.leap.year <- function(year) { year%%4 == 0 & (year%%100 != 0 | year%%400 == 0) }
