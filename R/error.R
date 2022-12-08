# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script contains functions for error checking and handling, i.e.:
#   errorHandling(temp = NULL, prec = NULL, pint = NULL, bsdf = NULL, TEMP = NULL, PREC = NULL, BSDF = NULL,
#      lat = NULL, lon = NULL, elv = NULL, year = NULL, n = NULL, S_f = NULL, T_a = NULL, S_w = NULL, P_atm = NULL,
#      P_n = NULL, W_n = NULL, W_max = NULL, MSMC = NULL, onlyLgthChecking = FALSE)
#   errorChecking(year = 2000, MSMC = 150., dvVAR = rep(0.7, 12), dvTEMP = rep(0.7, 12), dvPREC = rep(0.7, 12))
#   errorCheckingGrid(rs.temp = NULL, rs.prec = NULL, rs.bsdf = NULL, rl.elv = NULL, rl.MSMC = NULL)
#
#### DEFINE FUNCTIONS ################################################################################################

# ********************************************************************************************************************
# Name:     errorHandling
# Inputs:   - double, one-year time series of monthly mean air temperature, deg C (temp)
#           - double, one-year time series of monthly precipitation sum, mm (prec)
#           - double, one-year time series of monthly mean precipitation intensity, mm dy-1 (pint)
#           - double, one-year time series of monthly mean relative sunshine duration, unitless (bsdf)
#           - double, one-year time series of daily mean air temperature, deg C (TEMP)
#           - double, one-year time series of daily precipitation sum, mm (PREC)
#           - double, one-year time series of daily fractional sunshine duration, unitless (BSDF)
#           - double, latitude, deg (lat)
#           - double, longitude, deg (lon)
#           - double, elevation, m (elv)
#           - double, year (using astronomical year numbering) (year)
#           - double, day of the year (n)
#           - double, daily fractional sunshine duration, unitless (S_f)
#           - double, daily mean air temperature, deg C (T_a)
#           - double, evaporative supply rate, mm hr-1 (S_w)
#           - double, atmospheric pressure, Pa (P_atm)
#           - double, daily precipitation sum, mm (P_n)
#           - double, previous day's soil moisture, mm (W_n)
#           - double, maximum soil moisture capacity /'bucket size'/, mm (W_max / MSMC)
#           - logical, indicates that the function only need to check content of non-basic input data and to
#             determine the length dimension of objects (onlyLgthChecking)
# Returns:  - warning message if any input arguments do not meet requirements
# Features: This function checks and corrects the input arguments to each function of the package meet the
#           requirements.
# ********************************************************************************************************************
errorHandling <- function(temp = NULL, prec = NULL, pint = NULL, bsdf = NULL, TEMP = NULL, PREC = NULL, BSDF = NULL,
                          lat = NULL, lon = NULL, elv = NULL, year = NULL, n = NULL, S_f = NULL, T_a = NULL,
                          S_w = NULL, P_atm = NULL, P_n = NULL, W_n = NULL, W_max = NULL, MSMC = NULL,
                          onlyLgthChecking = FALSE) {

 ## Set auxiliary objects required for error handling
  # Labels for objects containing monthly meteorological variables
  cv.mly_var <- c("temp", "prec", "pint", "bsdf")

  # Labels for objects containing daily meteorological variables
  cv.dly_var <- c("TEMP", "PREC", "BSDF")

  # Labels for objects containing meteorological variables
  cv.met_var <- c(cv.mly_var, cv.dly_var)

  # Labels for objects containing geographical data
  cv.geo_dta <- c("lat", "lon", "elv")

  # Labels for objects containing temporal data
  cv.tle_dta <- c("year", "n")

  # Labels for objects containing location data
  cv.loc_dta <- c(cv.geo_dta, cv.tle_dta)

  # Labels for objects containing physical quantities
  cv.phy_qty <- c("S_f", "T_a", "S_w", "P_atm", "P_n", "W_n", "W_max", "MSMC")

  # Labels for valid classes of objects containing meteorological variables
  cv.vld_cls <- c("numeric", "integer", "ts", "xts", "zoo", "matrix", "data.frame", "logical")

  # Labels for time series classes of objects containing meteorological variables
  cv.1ts_cls <- c("ts", "xts", "zoo")

  # Labels for valid classes of objects containing location data
  cv.0ts_cls <- setdiff(cv.vld_cls, cv.1ts_cls)

  # An auxiliary matrix for checking the dimensions of point data
  cv.pt_obj <- c(cv.met_var, cv.loc_dta, cv.phy_qty)
  df.pt_len <- data.frame("object" = cv.pt_obj, "length" = as.numeric(NA))

  # Labels for objects containing real numbers
  cv.yre_num <- setdiff(cv.pt_obj, cv.tle_dta)

  # Labels for objects containing non-basic input data
  cv.nba_dta <- c("n", "S_w", "P_atm", "W_n")

  if (!onlyLgthChecking) {

   ## Correcting missing values of objects in order to facilitate content-checking
    for (i_i in 1 : length(cv.tle_dta)) {
      tle_dta <- get(cv.tle_dta[i_i])
      if (!is.null(tle_dta)) {
        assign(cv.tle_dta[i_i], replace(tle_dta, is.na(tle_dta), NA_integer_))
      }
    }

    for (i_r in 1 : length(cv.yre_num)) {
      yre_num <- get(cv.yre_num[i_r])
      if (!is.null(yre_num)) {
        assign(cv.yre_num[i_r], replace(yre_num, is.na(yre_num), NA_real_))
      }
    }


   ## Checking and correcting objects
    # Checking the class of objects containing meteorological variables
    for (i_mv in 1 : length(cv.met_var)) {
      met_var <- get(cv.met_var[i_mv])
      if (!is.null(met_var)) {
        if (length(which(!is.na(match(class(met_var), cv.vld_cls)))) < 1L) {
          stop("Invalid argument: 'class(", cv.met_var[i_mv],
               ")' must be in c('numeric', 'integer', 'ts', 'xts', 'zoo', 'matrix', 'data.frame', 'logical').")
        }
      }
    }

    # Convert a data frame-class object containing a meteorological variable
    # to a matrix in order to check its content
    for (i_mv in 1 : length(cv.met_var)) {
      met_var <- get(cv.met_var[i_mv])
      if (!is.null(met_var)) {
        if (any("data.frame" %in% class(met_var))) {
          assign(cv.met_var[i_mv], as.matrix(met_var))
        }
      }
    }

    # Checking and correcting the data structure and content of objects containing monthly meteorological variables
    for (i_m in 1 : length(cv.mly_var)) {
      mly_var <- get(cv.mly_var[i_m])
      if (!is.null(mly_var)) {
        if (is.null(dim(mly_var))) {
          if (!is.list(mly_var) & length(mly_var) == 12L) {
            if (is.numeric(mly_var) &
                all(mly_var[!is.na(mly_var)] >= varFeatures[cv.mly_var[i_m], "lwst_val"] &
                    mly_var[!is.na(mly_var)] <= varFeatures[cv.mly_var[i_m], "hgst_val"])) {
              row_ttl <- if (length(which(!is.na(match(class(mly_var), cv.1ts_cls)))) >= 1)
                stats::time(mly_var) else NULL
              col_ttl <- month.abb
              assign(cv.mly_var[i_m], matrix(mly_var, ncol = 12L, dimnames = list(row_ttl, col_ttl), byrow = TRUE))
            } else {
              stop("Invalid argument: '", cv.mly_var[i_m], "' has to be a vector of length 12 or ",
                   "a 12-column matrix-like object with ", varFeatures[cv.mly_var[i_m], "rng_wrd"], ".")
            }
          } else {
            stop("Invalid argument: '", cv.mly_var[i_m], "' has to be a vector of length 12 or ",
                 "a 12-column matrix-like object with ", varFeatures[cv.mly_var[i_m], "rng_wrd"], ".")
          }
        } else if (length(dim(mly_var)) == 2L & dim(mly_var)[1L] >= 1L & dim(mly_var)[2L] == 12L) {
          if (is.numeric(mly_var) &
              all(mly_var[!is.na(mly_var)] >= varFeatures[cv.mly_var[i_m], "lwst_val"] &
                  mly_var[!is.na(mly_var)] <= varFeatures[cv.mly_var[i_m], "hgst_val"])) {
            row_ttl <- if (length(which(!is.na(match(class(mly_var), cv.1ts_cls)))) >= 1)
              stats::time(mly_var) else dimnames(mly_var)[[1L]]
            col_ttl <- if (is.null(dimnames(mly_var)[[2L]])) month.abb else dimnames(mly_var)[[2L]]
            assign(cv.mly_var[i_m], matrix(mly_var, ncol = 12L, dimnames = list(row_ttl, col_ttl)))
          } else {
            stop("Invalid argument: '", cv.mly_var[i_m], "' has to be a vector of length 12 or ",
                 "a 12-column matrix-like object with ", varFeatures[cv.mly_var[i_m], "rng_wrd"], ".")
          }
        } else {
          stop("Invalid argument: '", cv.mly_var[i_m], "' has to be a vector of length 12 or ",
               "a 12-column matrix-like object with ", varFeatures[cv.mly_var[i_m], "rng_wrd"], ".")
        }
      }
    }

    # Checking and correcting the data structure and content of objects containing daily meteorological variables
    for (i_d in 1 : length(cv.dly_var)) {
      dly_var <- get(cv.dly_var[i_d])
      if (!is.null(dly_var)) {
        if (is.null(dim(dly_var))) {
          if (!is.list(dly_var) & any(length(dly_var) == c(365L, 366L))) {
            if (is.numeric(dly_var) &
                all(dly_var[!is.na(dly_var)] >= varFeatures[cv.dly_var[i_d], "lwst_val"] &
                    dly_var[!is.na(dly_var)] <= varFeatures[cv.dly_var[i_d], "hgst_val"])) {
              row_ttl <- if (length(which(!is.na(match(class(dly_var), cv.1ts_cls)))) >= 1)
                stats::time(dly_var) else NULL
              col_ttl <- seq(1, length(dly_var))
              assign(cv.dly_var[i_d], matrix(dly_var, ncol = length(dly_var), dimnames = list(row_ttl, col_ttl),
                                             byrow = TRUE))
            } else {
              stop("Invalid argument: '", cv.dly_var[i_d], "' has to be a vector of length 365 (or 366) or ",
                   "a 365-column (or 366-column) matrix-like object with ", varFeatures[cv.dly_var[i_d], "rng_wrd"],
                   ".")
            }
          } else {
            stop("Invalid argument: '", cv.dly_var[i_d], "' has to be a vector of length 365 (or 366) or ",
                 "a 365-column (or 366-column) matrix-like object with ", varFeatures[cv.dly_var[i_d], "rng_wrd"],
                 ".")
          }
        } else if (length(dim(dly_var)) == 2L & dim(dly_var)[1L] >= 1L & any(dim(dly_var)[2L] == c(365L, 366L))) {
          if (is.numeric(dly_var) &
              all(dly_var[!is.na(dly_var)] >= varFeatures[cv.dly_var[i_d], "lwst_val"] &
                  dly_var[!is.na(dly_var)] <= varFeatures[cv.dly_var[i_d], "hgst_val"])) {
            row_ttl <- if (length(which(!is.na(match(class(dly_var), cv.1ts_cls)))) >= 1)
              stats::time(dly_var) else dimnames(dly_var)[[1L]]
            col_ttl <- if (is.null(dimnames(dly_var)[[2L]])) seq(1, dim(dly_var)[2L]) else dimnames(dly_var)[[2L]]
            assign(cv.dly_var[i_d], matrix(dly_var, ncol = dim(dly_var)[2L], dimnames = list(row_ttl, col_ttl)))
          } else {
            stop("Invalid argument: '", cv.dly_var[i_d], "' has to be a vector of length 365 (or 366) or ",
                 "a 365-column (or 366-column) matrix-like object with ", varFeatures[cv.dly_var[i_d], "rng_wrd"],
                 ".")
          }
        } else {
          stop("Invalid argument: '", cv.dly_var[i_d], "' has to be a vector of length 365 (or 366) or ",
               "a 365-column (or 366-column) matrix-like object with ", varFeatures[cv.dly_var[i_d], "rng_wrd"],
               ".")
        }
      }
    }


    # Checking the class of objects containing location data and physical quantities
    for (i_ld in 1 : length(c(cv.loc_dta, cv.phy_qty))) {
      loc_dta <- get(c(cv.loc_dta, cv.phy_qty)[i_ld])
      if (!is.null(loc_dta)) {
        if (length(which(!is.na(match(class(loc_dta), cv.0ts_cls)))) < 1L) {
          stop("Invalid argument: 'class(", c(cv.loc_dta, cv.phy_qty)[i_ld],
               ")' must be in c('numeric', 'integer', 'matrix', 'data.frame', 'logical').")
        }
      }
    }

    # Convert a data frame-class object containing a location data or physical quantity
    # to a matrix in order to check its content
    for (i_ld in 1 : length(c(cv.loc_dta, cv.phy_qty))) {
      loc_dta <- get(c(cv.loc_dta, cv.phy_qty)[i_ld])
      if (!is.null(loc_dta)) {
        if (any("data.frame" %in% class(loc_dta))) {
          assign(c(cv.loc_dta, cv.phy_qty)[i_ld], as.matrix(loc_dta))
        }
      }
    }

    # Checking and correcting the data structure and content of objects containing geographical data
    for (i_gd in 1 : length(cv.geo_dta)) {
      geo_dta <- get(cv.geo_dta[i_gd])
      if (!is.null(geo_dta)) {
        if (is.null(dim(geo_dta))) {
          if (!is.list(geo_dta) & length(geo_dta) >= 1L) {
            if (!is.numeric(geo_dta) |
                any(geo_dta[!is.na(geo_dta)] < varFeatures[cv.geo_dta[i_gd], "lwst_val"] |
                    geo_dta[!is.na(geo_dta)] > varFeatures[cv.geo_dta[i_gd], "hgst_val"])) {
              stop("Invalid argument: '", cv.geo_dta[i_gd], "' has to be a vector or ",
                   "a one-column matrix-like object with ", varFeatures[cv.geo_dta[i_gd], "rng_wrd"], ".")
            }
          } else {
            stop("Invalid argument: '", cv.geo_dta[i_gd], "' has to be a vector or ",
                 "a one-column matrix-like object with ", varFeatures[cv.geo_dta[i_gd], "rng_wrd"], ".")
          }
        } else if (length(dim(geo_dta)) == 2L & dim(geo_dta)[1L] >= 1L) {
          if (dim(geo_dta)[2L] == 1L) {
            if (is.numeric(geo_dta) &
                all(geo_dta[!is.na(geo_dta)] >= varFeatures[cv.geo_dta[i_gd], "lwst_val"] &
                    geo_dta[!is.na(geo_dta)] <= varFeatures[cv.geo_dta[i_gd], "hgst_val"])) {
              assign(cv.geo_dta[i_gd], as.vector(geo_dta))
            } else {
              stop("Invalid argument: '", cv.geo_dta[i_gd], "' has to be a vector or ",
                   "a one-column matrix-like object with ", varFeatures[cv.geo_dta[i_gd], "rng_wrd"], ".")
            }
          } else {
            stop("Invalid argument: '", cv.geo_dta[i_gd], "' has to be a vector or ",
                 "a one-column matrix-like object with ", varFeatures[cv.geo_dta[i_gd], "rng_wrd"], ".")
          }
        } else {
          stop("Invalid argument: '", cv.geo_dta[i_gd], "' has to be a vector or ",
               "a one-column matrix-like object with ", varFeatures[cv.geo_dta[i_gd], "rng_wrd"], ".")
        }
      }
    }

    # Checking and correcting the data structure and content of objects containing temporal data
    for (i_td in 1 : length(cv.tle_dta)) {
      tle_dta <- get(cv.tle_dta[i_td])
      if (!is.null(tle_dta)) {
        if (is.null(dim(tle_dta))) {
          if (!is.list(tle_dta) & length(tle_dta) >= 1L) {
            if (!is.numeric(tle_dta) |
                any(tle_dta[!is.na(tle_dta)] < varFeatures[cv.tle_dta[i_td], "lwst_val"] |
                    tle_dta[!is.na(tle_dta)] > varFeatures[cv.tle_dta[i_td], "hgst_val"]) |
                any(tle_dta[!is.na(tle_dta)] %% 1 != 0)) {
              stop("Invalid argument: '", cv.tle_dta[i_td], "' has to be a vector or ",
                   "a one-column matrix-like object with ", varFeatures[cv.tle_dta[i_td], "rng_wrd"], ".")
            }
          } else {
            stop("Invalid argument: '", cv.tle_dta[i_td], "' has to be a vector or ",
                 "a one-column matrix-like object with ", varFeatures[cv.tle_dta[i_td], "rng_wrd"], ".")
          }
        } else if (length(dim(tle_dta)) == 2L & dim(tle_dta)[1L] >= 1L) {
          if (dim(tle_dta)[2L] == 1L) {
            if (is.numeric(tle_dta) &
                all(tle_dta[!is.na(tle_dta)] >= varFeatures[cv.tle_dta[i_td], "lwst_val"] &
                    tle_dta[!is.na(tle_dta)] <= varFeatures[cv.tle_dta[i_td], "hgst_val"]) &
                all(tle_dta[!is.na(tle_dta)] %% 1 == 0)) {
              assign(cv.tle_dta[i_td], as.vector(tle_dta))
            } else {
              stop("Invalid argument: '", cv.tle_dta[i_td], "' has to be a vector or ",
                   "a one-column matrix-like object with ", varFeatures[cv.tle_dta[i_td], "rng_wrd"], ".")
            }
          } else {
            stop("Invalid argument: '", cv.tle_dta[i_td], "' has to be a vector or ",
                 "a one-column matrix-like object with ", varFeatures[cv.tle_dta[i_td], "rng_wrd"], ".")
          }
        } else {
          stop("Invalid argument: '", cv.tle_dta[i_td], "' has to be a vector or ",
               "a one-column matrix-like object with ", varFeatures[cv.tle_dta[i_td], "rng_wrd"], ".")
        }
      }
    }

    # Checking and correcting the data structure and content of objects containing physical quantities
    for (i_pq in 1 : length(cv.phy_qty)) {
      phy_qty <- get(cv.phy_qty[i_pq])
      if (!is.null(phy_qty)) {
        if (is.null(dim(phy_qty))) {
          if (!is.list(phy_qty) & length(phy_qty) >= 1L) {
            if (!is.numeric(phy_qty) |
                any(phy_qty[!is.na(phy_qty)] < varFeatures[cv.phy_qty[i_pq], "lwst_val"] |
                    phy_qty[!is.na(phy_qty)] > varFeatures[cv.phy_qty[i_pq], "hgst_val"])) {
              stop("Invalid argument: '", cv.phy_qty[i_pq], "' has to be a vector or ",
                   "a one-column matrix-like object with ", varFeatures[cv.phy_qty[i_pq], "rng_wrd"], ".")
            }
          } else {
            stop("Invalid argument: '", cv.phy_qty[i_pq], "' has to be a vector or ",
                 "a one-column matrix-like object with ", varFeatures[cv.phy_qty[i_pq], "rng_wrd"], ".")
          }
        } else if (length(dim(phy_qty)) == 2L & dim(phy_qty)[1L] >= 1L) {
          if (dim(phy_qty)[2L] == 1L) {
            if (is.numeric(phy_qty) &
                all(phy_qty[!is.na(phy_qty)] >= varFeatures[cv.phy_qty[i_pq], "lwst_val"] &
                    phy_qty[!is.na(phy_qty)] <= varFeatures[cv.phy_qty[i_pq], "hgst_val"])) {
              assign(cv.phy_qty[i_pq], as.vector(phy_qty))
            } else {
              stop("Invalid argument: '", cv.phy_qty[i_pq], "' has to be a vector or ",
                   "a one-column matrix-like object with ", varFeatures[cv.phy_qty[i_pq], "rng_wrd"], ".")
            }
          } else {
            stop("Invalid argument: '", cv.phy_qty[i_pq], "' has to be a vector or ",
                 "a one-column matrix-like object with ", varFeatures[cv.phy_qty[i_pq], "rng_wrd"], ".")
          }
        } else {
          stop("Invalid argument: '", cv.phy_qty[i_pq], "' has to be a vector or ",
               "a one-column matrix-like object with ", varFeatures[cv.phy_qty[i_pq], "rng_wrd"], ".")
        }
      }
    }

  }


  # Checking the content of the object containing non-basic data
  for (i_nd in 1 : length(cv.nba_dta)) {
    nba_dta <- get(cv.nba_dta[i_nd])
    if (!is.null(nba_dta)) {
      if (!is.numeric(nba_dta) |
          any(nba_dta[!is.na(nba_dta)] < varFeatures[cv.nba_dta[i_nd], "lwst_val"] |
              nba_dta[!is.na(nba_dta)] > varFeatures[cv.nba_dta[i_nd], "hgst_val"])) {
        stop("Invalid argument: '", cv.nba_dta[i_nd], "' has to be a vector or ",
             "a one-column matrix-like object with ", varFeatures[cv.nba_dta[i_nd], "rng_wrd"], ".")
      }
    }
  }


  # Checking whether the vectors and the first dimensions of the matrices have the same length
  for (i_pt in 1 : length(cv.pt_obj)) {
    pt_obj <- get(cv.pt_obj[i_pt])
    if (!is.null(pt_obj)) {
      df.pt_len[i_pt, "length"] <- ifelse(is.null(dim(pt_obj)), length(pt_obj), dim(pt_obj)[1L])
    }
  }
  df.dc_pt_len <- df.pt_len[!is.na(df.pt_len[, "length"]), ]
  if (!(diff(range(df.dc_pt_len[, "length"])) < .Machine$double.eps ^ 0.5)) {
    cv.iv_obj <- paste0("'", df.dc_pt_len$object, "'")
    if (all(df.dc_pt_len$object %in% cv.met_var)) {
      stop("Invalid argument: The first dimensions of the objects ",
           paste0(paste(cv.iv_obj[1 : (length(cv.iv_obj) - 1)], sep = "", collapse = ", "), " and ",
                  cv.iv_obj[length(cv.iv_obj)]),
           " have not the same length.")
    } else if (all(df.dc_pt_len$object %in% c(cv.loc_dta, cv.phy_qty))) {
      stop("Invalid argument: The vectors ",
           paste0(paste(cv.iv_obj[1 : (length(cv.iv_obj) - 1)], sep = "", collapse = ", "), " and ",
                  cv.iv_obj[length(cv.iv_obj)]),
           " have not the same length.")
    } else {
      cv.iv_mvo <- paste0("'", df.dc_pt_len$object[df.dc_pt_len$object %in% cv.met_var], "'")
      cv.iv_ldo <- paste0("'", df.dc_pt_len$object[df.dc_pt_len$object %in% c(cv.loc_dta, cv.phy_qty)], "'")
      stop("Invalid argument: The first dimension(s) of the object(s) ",
           ifelse(length(cv.iv_mvo) == 1, cv.iv_mvo,
                  paste0(paste(cv.iv_mvo[1 : (length(cv.iv_mvo) - 1)], sep = "", collapse = ", "), " and ",
                         cv.iv_mvo[length(cv.iv_mvo)])),
           " and the vector(s) ",
           ifelse(length(cv.iv_ldo) == 1, cv.iv_ldo,
                  paste0(paste(cv.iv_ldo[1 : (length(cv.iv_ldo) - 1)], sep = "", collapse = ", "), " and ",
                         cv.iv_ldo[length(cv.iv_ldo)])),
           " have not the same length.")
    }
  }

  err_han <- mget(cv.pt_obj)
  err_han$lgth <- unique(df.dc_pt_len$length)

  return(err_han)

}


# ********************************************************************************************************************
# Name:     errorChecking
# Inputs:   - double, year (using astronomical year numbering) (year)
#           - double, maximum soil moisture capacity /'bucket size'/, mm (MSMC)
#           - double, monthly time series of the damping variable for the given climate variable (dvVAR)
#           - double, monthly time series of the damping variable for the air temperature data (dvTEMP)
#           - double, monthly time series of the damping variable for the precipitation data (dvPREC)
# Returns:  - warning message if any input arguments do not meet requirements
# Features: This function checks the input arguments to each function of the package meet the requirements.
# ********************************************************************************************************************
errorChecking <- function(year = 2000, MSMC = 150., dvVAR = rep(0.7, 12), dvTEMP = rep(0.7, 12),
                          dvPREC = rep(0.7, 12)) {

  # Labels for objects to be checked
  cv.ckd_obj <- c("year", "MSMC", "dvVAR", "dvTEMP", "dvPREC")

  # Labels for objects containing damping variables that are used in an iterative interpolation technique
  cv.dmp_var <- c("dvVAR", "dvTEMP", "dvPREC")

  # Labels for valid classes of objects to be checked
  cv.0ts_cls <- c("numeric", "integer", "matrix", "data.frame", "logical")

  # Checking the class of objects
  for (i_co in 1 : length(cv.ckd_obj)) {
    ckd_obj <- get(cv.ckd_obj[i_co])
    if (length(which(!is.na(match(class(ckd_obj), cv.0ts_cls)))) < 1L) {
      stop("Invalid argument: 'class(", cv.ckd_obj[i_co],
           ")' must be in c('numeric', 'integer', 'matrix', 'data.frame', 'logical').")
    }
  }

  # Checking the content of the object containing damping variables that are used
  # in an iterative interpolation technique
  for (i in 1 : length(cv.dmp_var)) {
    dmp_var <- get(cv.dmp_var[i])
    if (is.null(dim(dmp_var)) & !is.list(dmp_var) & length(dmp_var) == 12L) {
      if (!is.numeric(dmp_var) | any(dmp_var < 0.) | any(dmp_var > 1.)) {
        stop("Invalid argument: '", cv.dmp_var[i],
             "' has to be a vector of length 12 with values from the interval [0, 1].")
      }
    } else {
      stop("Invalid argument: '", cv.dmp_var[i],
           "' has to be a vector of length 12 with values from the interval [0, 1].")
    }
  }

}


# ********************************************************************************************************************
# Name:     errorCheckingGrid
# Inputs:   - S4, one-year time series of monthly mean air temperature, deg C (rs.temp)
#           - S4, one-year time series of monthly precipitation sum, mm (rs.prec)
#           - S4, one-year time series of monthly mean relative sunshine duration, unitless (rs.bsdf)
#           - S4, elevation, m (rl.elv)
#           - S4, maximum soil moisture capacity /'bucket size'/, mm (rl.MSMC)
# Returns:  - warning message if any input arguments do not meet requirements
# Features: This function checks the input arguments to each function of the package meet the requirements.
# ********************************************************************************************************************
errorCheckingGrid <- function(rs.temp = NULL, rs.prec = NULL, rs.bsdf = NULL, rl.elv = NULL, rl.MSMC = NULL) {

  # Labels for objects to be checked
  cv.ckd_obj <- c("rs.temp", "rs.prec", "rs.bsdf", "rl.elv", "rl.MSMC")

  # Labels for valid classes of objects to be checked
  cv.rstr_cls <- c("RasterLayer", "RasterBrick", "RasterStack")

  # Number of layers in each Raster* object
  nv.num_lyr <- c(12, 12, 12, 1, 1)

  # Checking the class of objects
  for (i_co in 1 : length(cv.ckd_obj)) {
    ckd_obj <- get(cv.ckd_obj[i_co])
    if (!is.null(ckd_obj)) {
      if (length(which(!is.na(match(class(ckd_obj), cv.rstr_cls)))) < 1L) {
        stop("Invalid argument: '", cv.ckd_obj[i_co], "' has to be a RasterLayer, RasterBrick, or a RasterStack.")
      }
    }
  }

  # Checking whether Raster* objects contain enough layers
  for (i_co in 1 : length(cv.ckd_obj)) {
    ckd_obj <- get(cv.ckd_obj[i_co])
    num_lyr <- nv.num_lyr[i_co]
    if (!is.null(ckd_obj)) {
      if (raster::nlayers(ckd_obj) != num_lyr) {
        stop("Invalid argument: '", cv.ckd_obj[i_co], "' has to be a Raster* object with ", num_lyr,
             ifelse(num_lyr == 1, " layer.", " layers."))
      }
    }
  }

  # Checking whether Raster* objects have the same extent, number of rows and columns, projection, and resolution
  raster::rasterOptions(tolerance = 0.5)
  if (!raster::compareRaster(Filter(Negate(is.null), as.list(rs.temp, rs.prec, rs.bsdf, rl.elv, rl.MSMC)),
                             res = TRUE, rotation = FALSE, stopiffalse = FALSE, showwarning = TRUE)) {
    cv.gr_obj <- cv.ckd_obj
    cv.gr_obj <- cv.gr_obj[sapply(cv.gr_obj, function(x) { !is.null(get(x)) })]
    stop("The Raster* objects ",
         paste0(paste(cv.gr_obj[1 : (length(cv.gr_obj) - 1)], sep = "", collapse = ", "), " and ",
                cv.gr_obj[length(cv.gr_obj)]),
         " must have the same extent, number of rows and columns, projection, and resolution.")
  }

}
