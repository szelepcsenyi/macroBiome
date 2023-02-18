#' Vegetation Classifier Using the KGC System
#'
#' @description Calculates the values of bioclimatic indices used in the Köppen-Geiger classification (KGC) system
#'     (Köppen 1936), and designates the KGC type using these values, by using the monthly time series of temperature
#'     and precipitation. The classification scheme is based on the procedure described by Köppen (1936) and follows
#'     the modifications described by Peel et al. (2007).
#'
#' @param temp 'numeric' R object with one-year time series of monthly mean air temperature (in °C)
#' @param prec 'numeric' R object with one-year time series of monthly precipitation sum (in mm)
#' @param verbose 'logical' scalar that indicates whether or not values of the bioclimatic indices used should be
#'     added to the output.
#'
#' @details To classify vegetation, the KGC system developed by Köppen (1936) and fine-tuned by Peel. et al. (2007)
#'     uses the values of the following 13 bioclimatic indices:
#'
#'     \itemize{
#'       \item{\code{tap}: Total Annual Precipitation (in mm)}
#'       \item{\code{mat}: Mean Annual Temperature (in °C)}
#'       \item{\code{tw}: Mean Temperature of the Warmest Month (in °C)}
#'       \item{\code{tc}: Mean Temperature of the Coldest Month (in °C)}
#'       \item{\code{tm10}: Number of Months with Mean Temperature above 10°C (dimensionless)}
#'       \item{\code{pdry}: Precipitation Sum of the Driest Month (in mm)}
#'       \item{\code{psdry}: Precipitation Sum of the Driest Month in the Summer Half-Year (in mm)}
#'       \item{\code{pwdry}: Precipitation Sum of the Driest Month in the Winter Half-Year (in mm)}
#'       \item{\code{pswet}: Precipitation Sum of the Wettest Month in the Summer Half-Year (in mm)}
#'       \item{\code{pwwet}: Precipitation Sum of the Wettest Month in the Winter Half-Year (in mm)}
#'       \item{\code{ps}: Precipitation Sum of the Summer Half-Year (in mm)}
#'       \item{\code{pw}: Precipitation Sum of the Winter Half-Year (in mm)}
#'       \item{\code{pth}: Dryness Threshold (in mm)}
#'     }
#'
#'     For details about calculating bioclimatic indices, see the function
#'     \code{\link[macroBiome]{cliBioCliIdxPoints}}. Since \code{pth} is more of a technical measure, it is not
#'     calculated by the function \code{\link[macroBiome]{cliBioCliIdxPoints}}. The value of \code{pth} depends on
#'     mean annual temperature and annual cycle of precipitation: \code{pth = 2 * mat} if >70% of precipitation falls
#'     in winter half-year, (b) \code{pth = 2 * mat + 28} if >70% of precipitation falls in summer half-year, and
#'     otherwise \code{pth = 2 * mat + 14}. For this index, the same definitions are used for seasons as in the
#'     function \code{\link[macroBiome]{cliBioCliIdxPoints}}, i.e., summer (winter) half-year is defined as the
#'     warmer (cooler) six month period of AMJJAS (from April to September) and ONDJFM (from October to March). \cr
#'     Numerous variants of the Köppen classification system are known (e.g., Köppen-Geiger classification: Köppen
#'     1936; Köppen-Trewartha classification: Trewartha and Horn 1980). Here, one of the most widely used versions of
#'     the Köppen-Geiger classification system is implemented, in accordance with works of Peel et al. (2007) and
#'     Beck et al. (2018). This classification system is the same as that presented by Köppen (1936) with three
#'     differences. First, classes \code{'C'} (temperate) and \code{'D'} (cold) are distinguished using a 0°C
#'     threshold instead of a -3°C threshold, following Russell (1931). Second, the sub-classes of the class
#'     \code{'B'} (arid) are identified depending on whether 70% of precipitation falls in the summer or winter
#'     half-year. Third, the sub-classes \code{'s'} (dry summer) and \code{'w'} (dry winter) within the classes
#'     \code{'C'} and \code{'D'} are made mutually exclusive by assigning \code{'s'} when more precipitation falls in
#'     winter than in summer and assigning \code{'w'} otherwise. In this version, a total of 30 KGC types are
#'     distinguished (see \code{\link[macroBiome]{vegClsNumCodes}}).
#'
#' @return Depending on the setting, a data frame with one or more columns where the KGC types are stored in the
#'     last (character) column, while the additional columns contain the values of bioclimatic indices used. The
#'     abbreviations of KGC types can be found in the data frame \code{\link[macroBiome]{vegClsNumCodes}}. If
#'     \code{verbose = FALSE}, the return object is a one-column data frame with the KGC types.
#'
#' @note As with any function with a point mode, a set of basic input data is defined here. In this case, they are as
#'     follows: \code{'temp'} (one-year time series of monthly mean air temperature), and \code{'prec'} (one-year
#'     time series of monthly precipitation sum). The objects \code{'temp'} and \code{'pre'} must be either vectors
#'     of length 12 or 12-column matrices. The first dimensions of these matrices have to be the same length. The
#'     function automatically converts vectors into single-row matrices during the error handling, and then uses
#'     these matrices. The first dimensions of these matrices determines the number of rows in the result matrix.
#'
#' @references
#'
#' \cite{Beck HE, Zimmermann NE, McVicar TR, Vergopolan N, Berg A, Wood EF (2018) Present and future Köppen-Geiger
#'     climate classification maps at 1-km resolution. Sci Data 5:180214. \doi{10.1038/sdata.2018.214}}
#'
#' \cite{Köppen W (1936) Das geographische System der Klimate. In: Köppen W, Geiger R (eds) Handbuch der
#'     Klimatologie. Verlag von Gebrüder Borntraeger, Berlin, Germany, pp 1–44}
#'
#' \cite{Peel MC, Finlayson BL, McMahon TA (2007) Updated world map of the Köppen-Geiger climate classification.
#'     Hydrol Earth Syst Sci 11(5):1633–1644. \doi{10.5194/hess-11-1633-2007}}
#'
#' \cite{Russell RJ (1931) Dry Climates of the United States: I. Climatic Map. University of California, Publications
#'     in Geography 5:1–41}
#'
#' \cite{Trewartha GT, Horn LH (1980) An Introduction to Climate. Fifth Edition. McGraw-Hill, New York, NY}
#'
#' @examples
#' # Loading mandatory data for the Example 'Points'
#' data(inp_exPoints)
#'
#' # Designate the KGC type (using the related bioclimatic indices),
#' # at a grid cell near Szeged, Hungary (46.3N, 20.2E) (for the normal period 1981-2010)
#' with(inp_exPoints, {
#' KGC <- cliKoppenPoints(colMeans(temp), colMeans(prec), verbose = TRUE)
#' numCode <- which(sapply(vegClsNumCodes$Code.KGC, identical, KGC[, "vegCls"]))
#' cbind(KGC[,-c(14)], vegClsNumCodes[numCode, c("Name.KGC", "Code.KGC")])
#' })
#'
#' @importFrom stats setNames
#'
#' @export
#'
cliKoppenPoints <- function(temp, prec, verbose = FALSE) {

  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  err_han <- errorHandling(temp = temp, prec = prec)
  list2env(Filter(Negate(is.null), err_han), envir = environment())

  cv.arg <- c("temp", "prec")
  for (i in 1 : length(cv.arg)) {
    if (is.null(get(cv.arg[i]))) { stop("Invalid argument: '", cv.arg[i], "' is missing, with no default.") }
  }

  cv.bci <- c("tap", "mat", "tw", "tc", "tm10", "pdry", "psdry", "pwdry", "pswet", "pwwet", "ps", "pw")

  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 01. Calculate values of each bioclimatic index required to classify vegetation
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  bioCliIdx <- cliBioCliIdxPoints(temp, prec, bciOpVar = cv.bci, argCkd = T)
  list2env(unclass(as.data.frame(bioCliIdx)), envir = environment())

  # Dryness Threshold, mm (pth) | ref: table caption of Table 1 in Peel et al. (2007)
  pth <- ifelse(0.7 * tap < pw, 2. * mat, ifelse(0.7 * tap >= ps, 2. * mat + 28., 2. * mat + 14.))

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 02. Determine the vegetation type by using the Köppen-Geiger classification system
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Set a vector containing all types used by the KGC system
  cls <- vegCkgcTypes$CODE.KGC

  # Classify the vegetation type under certain boundary conditions
  # Arid
  B   <- tap < 10. * pth
  BW  <- B & tap < 5. * pth
  BWh <- BW & mat >= 18.
  BWk <- BW & mat < 18.
  BS  <- B & tap >= 5. * pth
  BSh <- BS & mat >= 18.
  BSk <- BS & mat < 18.

  # Tropical
  A  <- tc >= 18 & !B
  Af <- A & pdry >= 60.
  Am <- A & !Af & pdry >= 100. - tap / 25.
  Aw <- A & !Af & pdry < 100. - tap / 25.

  # Temperate
  C <- tw > 10. & tc > 0. & tc < 18. & !B
  Cs  <- C & psdry < 40. & psdry < (pwwet / 3.)
  Cw  <- C & pwdry < pswet / 10.
  Cs[Cs & Cw & ps > pw] <- FALSE
  Cw[Cs & Cw & ps <= pw] <- FALSE
  Csa <- Cs & tw >= 22.
  Csb <- Cs & !Csa & tm10 >= 4
  Csc <- Cs & !Csa & !Csb & tm10 >= 1 & tm10 < 4
  Cwa <- Cw & tw >= 22.
  Cwb <- Cw & !Cwa & tm10 >= 4
  Cwc <- Cw & !Cwa & !Cwb & tm10 >= 1 & tm10 < 4
  Cf  <- C & !Cs & !Cw
  Cfa <- Cf & tw >= 22.
  Cfb <- Cf & !Cfa & tm10 >= 4
  Cfc <- Cf & !Cfa & !Cfb & tm10 >= 1 & tm10 < 4

  # Cold
  D <- tw > 10. & tc <= 0. & !B
  Ds  <- D & psdry < 40. & psdry < (pwwet / 3.)
  Dw  <- D & pwdry < pswet / 10.
  Ds[Ds & Dw & ps > pw] <- FALSE
  Dw[Ds & Dw & ps <= pw] <- FALSE
  Dsa <- Ds & tw >= 22.
  Dsb <- Ds & !Dsa & tm10 >= 4
  Dsd <- Ds & !Dsa & !Dsb & tc < -38.
  Dsc <- Ds & !Dsa & !Dsb & !Dsd
  Dwa <- Dw & tw >= 22.
  Dwb <- Dw & !Dwa & tm10 >= 4
  Dwd <- Dw & !Dwa & !Dwb & tc < -38.
  Dwc <- Dw & !Dwa & !Dwb & !Dwd
  Df  <- D & !Ds & !Dw
  Dfa <- Df & tw >= 22.
  Dfb <- Df & !Dfa & tm10 >= 4
  Dfd <- Df & !Dfa & !Dfb & tc < -38.
  Dfc <- Df & !Dfa & !Dfb & !Dfd

  # Polar
  E <- tw <= 10. & !B
  ET <- E & tw > 0.
  EF <- E & tw <= 0.

  psblVegCls <- do.call(cbind, as.list(setNames(mget(cls), cls)))
  vegCls <- sapply(1 : nrow(psblVegCls), function(i) { ifelse(length(which(psblVegCls[i, ])) == 0, NA,
                                                              names(which(psblVegCls[i, ]))) })

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  if (verbose) {
    rslt <- data.frame(cbind(bioCliIdx, pth), vegCls = vegCls)
  } else {
    rslt <- data.frame(vegCls = vegCls)
  }
  return(rslt)

}


#' Vegetation Classifier Using the KGC System
#'
#' @description Calculates the values of bioclimatic indices used in the Köppen-Geiger classification (KGC) system
#'     (Köppen 1936), and designates the KGC type using these values, by using the monthly time series of temperature
#'     and precipitation. The classification scheme is based on the procedure described by Köppen (1936) and follows
#'     the modifications described by Peel et al. (2007).
#'
#' @param rs.temp multi-layer Raster* object with one-year time series of monthly mean air temperature (in °C)
#' @param rs.prec multi-layer Raster* object with one-year time series of monthly precipitation sum (in mm)
#' @param verbose 'logical' scalar that indicates whether or not values of the bioclimatic indices used should be
#'     added to the output.
#' @param filename output filename
#' @param ... additional arguments passed on to \code{\link[raster]{writeRaster}}
#'
#' @details See \code{\link[macroBiome]{cliKoppenPoints}}.
#'
#' @return Depending on the setting, a RasterStack with one or more layers where the numeric integers encoding the
#'     KGC type are stored at the last layer, while the additional layers contain the values of bioclimatic indices
#'     used. The meaning of integers is given in the data frame \code{\link[macroBiome]{vegClsNumCodes}}. If
#'     \code{verbose = FALSE}, the return object is a single-layer RasterStack with numeric integers encoding the KGC
#'     type.
#'
#' @note The objects \code{'rs.temp'} and \code{'rs.prec'} must be 12-layer Raster* objects. These Raster* objects
#'     must have the same bounding box, projection, and resolution.
#'
#' @references
#'
#' \cite{Köppen W (1936) Das geographische System der Klimate. In: Köppen W, Geiger R (eds) Handbuch der
#'     Klimatologie. Verlag von Gebrüder Borntraeger, Berlin, Germany, pp 1–44}
#'
#' \cite{Peel MC, Finlayson BL, McMahon TA (2007) Updated world map of the Köppen-Geiger climate classification.
#'     Hydrol Earth Syst Sci 11(5):1633–1644. \doi{10.5194/hess-11-1633-2007}}
#'
#' @examples
#' # Loading mandatory data for the Example 'Climate Normal Grid'
#' data(inp_exClnrGrid)
#'
#' # Designate the KGC types (using the related bioclimatic indices)
#' # for Csongrad-Csanad County (for the normal period 1981-2010)
#' with(inp_exClnrGrid, {
#' rs.KGC <- cliKoppenGrid(temp, prec, verbose = TRUE)
#' rs.KGC
#' })
#'
#' @import raster
#'
#' @export
#'
cliKoppenGrid <- function(rs.temp, rs.prec, verbose = FALSE, filename = "", ...) {

  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  cv.arg <- c("rs.temp", "rs.prec")
  for (i in 1 : length(cv.arg)) {
    if (is.null(get(cv.arg[i]))) { stop("Invalid argument: '", cv.arg[i], "' is missing, with no default.") }
  }

  errorCheckingGrid(rs.temp = rs.temp, rs.prec = rs.prec)


  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

  n_lyr <- ifelse(verbose, 14, 1)

  rs.rslt <- brick(rs.temp, nl = n_lyr)

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
      cv.arg <- cv.mly_var
      for (i_arg in 1 : length(cv.arg)) {
        assign(substring(cv.arg[i_arg], 4), getValues(get(cv.arg[i_arg]), row = bs$row[i], nrows = bs$nrows[i]))
      }
      df.rslt <- cliKoppenPoints(temp, prec, verbose = verbose)
      numCode <- as.numeric(rownames(vegCkgcTypes))[match(df.rslt[["vegCls"]], vegCkgcTypes$CODE.KGC)]
      df.rslt[["vegCls"]] <- numCode

      rs.rslt <- writeValues(rs.rslt, as.matrix(df.rslt), bs$row[i])
      pbStep(pb, i)
    }
    rs.rslt <- writeStop(rs.rslt)
  } else {
    for (i in 1 : bs$n) {
      cv.mly_var <- c("rs.temp", "rs.prec")
      cv.arg <- cv.mly_var
      for (i_arg in 1 : length(cv.arg)) {
        assign(substring(cv.arg[i_arg], 4), getValues(get(cv.arg[i_arg]), row = bs$row[i], nrows = bs$nrows[i]))
      }
      df.rslt <- cliKoppenPoints(temp, prec, verbose = verbose)
      numCode <- as.numeric(rownames(vegCkgcTypes))[match(df.rslt[["vegCls"]], vegCkgcTypes$CODE.KGC)]
      df.rslt[["vegCls"]] <- numCode

      cols <- bs$row[i] : (bs$row[i] + bs$nrows[i] - 1)
      arr[, cols, ] <- array(as.matrix(df.rslt), dim = c(bs$nrows[i], ncol(rs.rslt), nlayers(rs.rslt)))
      pbStep(pb, i)
    }
    for (lyr in 1 : nlayers(rs.rslt)) {
      rs.rslt <- setValues(rs.rslt, as.vector(arr[, , lyr]), layer = lyr)
    }
  }

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  names(rs.rslt) <- colnames(df.rslt)
  return(stack(rs.rslt))

}
