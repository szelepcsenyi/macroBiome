#' Vegetation Classifier Using the HLZ System
#'
#' @description Calculates the values of bioclimatic indices used in the Holdridge life zone (HLZ) system (Holdridge
#'     1947, 1967), and designates the HLZ type using these values, by using the monthly time series of temperature
#'     and precipitation.
#'
#' @param temp 'numeric' R object with one-year time series of monthly mean air temperature (in °C)
#' @param prec 'numeric' R object with one-year time series of monthly precipitation sum (in mm)
#' @param verbose 'logical' scalar that indicates whether or not values of the bioclimatic indices used should be
#'     added to the output.
#'
#' @details To classify vegetation, the HLZ system developed by Holdridge (1947, 1967) uses the values of the
#'     following 3 bioclimatic indices:
#'
#'     \itemize{
#'       \item{\code{abt}: Mean Annual Biotemperature (Eq 1 in Szelepcsényi et al. (2014); in °C)}
#'       \item{\code{tap}: Total Annual Precipitation (in mm)}
#'       \item{\code{per}: Potential Evapotranspiration Ratio (Eq 4 in Szelepcsényi et al. (2014); dimensionless)}
#'     }
#'
#'     For details about calculating bioclimatic indices, see the function
#'     \code{\link[macroBiome]{cliBioCliIdxPoints}}. \cr
#'     The HLZ system classifies the vegetation type based on the distance from the ideal (theoretical) point in the
#'     3-dimensional space of bioclimatic indices. Numerous variants of the HLZ system are known (e.g.,
#'     Henderson-Sellers 1994; Yates et al. 2000). Here, one of its most widely used versions ('version with no
#'     altitudinal belts') is implemented, in accordance with works of Szelepcsényi et al. (2014, 2018). In this
#'     version, a total of 39 HLZ types are distinguished (see \code{\link[macroBiome]{vegClsNumCodes}}).
#'
#' @return Depending on the setting, a data frame with one or more columns where the HLZ types are stored in the last
#'     (character) column, while the additional columns contain the values of bioclimatic indices used. The
#'     abbreviations of HLZ types can be found in the data frame \code{\link[macroBiome]{vegClsNumCodes}}. If
#'     \code{verbose = FALSE}, the return object is a one-column data frame with the HLZ types.
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
#' \cite{Henderson-Sellers A (1994) Global terrestrial vegetation ‘prediction’: the use and abuse of climate and
#'     application models. Prog Phys Geogr 18(2):209–246. \doi{10.1177/030913339401800203}}
#'
#' \cite{Holdridge LR (1947) Determination of World Plant Formations From Simple Climatic Data. Science
#'     105(2727):367–368. \doi{10.1126/science.105.2727.367}}
#'
#' \cite{Holdridge LR (1967) Life zone ecology. Tropical Science Center, San Jose, Costa Rica}
#'
#' \cite{Szelepcsényi Z, Breuer H, Sümegi P (2014) The climate of Carpathian Region in the 20th century based on the
#'     original and modified Holdridge life zone system. Cent Eur J Geosci 6(3):293–307.
#'     \doi{10.2478/s13533-012-0189-5}}
#'
#' \cite{Szelepcsényi Z, Breuer H, Kis A, Pongrácz R, Sümegi P (2018) Assessment of projected climate change in the
#'     Carpathian Region using the Holdridge life zone system. Theor Appl Climatol 131(1–2):593–610.
#'     \doi{10.1007/s00704-016-1987-3}}
#'
#' \cite{Yates DN, Kittel TGF, Cannon RF (2000) Comparing the Correlative Holdridge Model to Mechanistic
#'     Biogeographical Models for Assessing Vegetation Distribution Response to Climatic Change. Clim Chang
#'     44(1–2):59–87. \doi{10.1023/A:1005495908758}}
#'
#' @examples
#' # Loading mandatory data for the Example 'Points'
#' data(inp_exPoints)
#'
#' # Designate the HLZ type (using the related bioclimatic indices),
#' # at a grid cell near Szeged, Hungary (46.3N, 20.2E) (for the normal period 1981-2010)
#' with(inp_exPoints, {
#' HLZ <- cliHoldridgePoints(colMeans(temp), colMeans(prec), verbose = TRUE)
#' numCode <- which(sapply(vegClsNumCodes$Code.HLZ, identical, HLZ[, "vegCls"]))
#' cbind(HLZ[,-c(4)], vegClsNumCodes[numCode, c("Name.HLZ", "Code.HLZ")])
#' })
#'
#' @importFrom stats setNames
#'
#' @export
#'
cliHoldridgePoints <- function(temp, prec, verbose = FALSE) {

  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  err_han <- errorHandling(temp = temp, prec = prec)
  list2env(Filter(Negate(is.null), err_han), envir = environment())

  cv.arg <- c("temp", "prec")
  for (i in 1 : length(cv.arg)) {
    if (is.null(get(cv.arg[i]))) { stop("Invalid argument: '", cv.arg[i], "' is missing, with no default.") }
  }

  cv.bci <- c("abt", "tap", "per")

  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 01. Calculate values of each bioclimatic index required to classify vegetation
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  bioCliIdx <- cliBioCliIdxPoints(temp, prec, bciOpVar = cv.bci, argCkd = T)
  list2env(setNames(split(bioCliIdx, col(bioCliIdx)), colnames(bioCliIdx)), envir = environment())

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 02. Set the result object containing vegetation class
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  vegCls <- psblVegCls <- rep(NA_character_, length = lgth)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 03. Determine the vegetation class by using the Holdridge life zone system
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Set a board containing the distance from the optimal value for each pair of indices and classes
  distBds <- array(NA, dim = c(nrow(hlzDefSubset), length(cv.bci), lgth),
                   dimnames = list(hlzDefSubset$ID, cv.bci, NULL))
  for (i_bci in 1 : length(cv.bci)) {
    bciLbl <- cv.bci[i_bci]
    bci <- get(bciLbl)
    for (i_vcl in 1 : nrow(hlzDefSubset)) {
      optVal <- log2(hlzDefSubset[[bciLbl]][i_vcl]) + 0.5
      distBds[i_vcl, i_bci, ] <- ifelse(lapply(bci, function(x) {x == 0}), NA, (log2(bci) - optVal) ** 2.)
    }
  }
  vld <- as.logical(apply(matrix(mapply(function(x) { !is.na(x) },
                                        data.frame(abt = get("abt"), tap = get("tap"), per = get("per"))),
                                 nrow = lgth), 1, prod))
  psblVegCls[vld] <- lapply(seq(1, lgth)[vld], function(i) {
    which(sqrt(rowSums(distBds[, , i])) == min(sqrt(rowSums(distBds[, , i]))))
  })

  # Set some magic numbers
  # Polar temperature line
  plTL <- 1.5

  # Minimum value of the total annual precipitation
  mn_tap <- 62.5

  # Minimum value of the potential evapotranspiration ratio
  mn_per <- 0.125

  # Frost or critical temperature line
  frTL <- 2. ** (log2(12.) + 0.5)

  # Vegetation classes in warm temperate and subtropical regions
  WtCls <- seq(17, 23)
  names(WtCls) <- c("WtD", "WtDs", "WtTs", "WtDf", "WtMf", "WtWf", "WtRf")
  StCls <- seq(24, 30)
  names(StCls) <- c("StD", "StDs", "StTw", "StDf", "StMf", "StWf", "StRf")

  # Classify the vegetation type under certain boundary conditions
  if (!identical(vld, integer(0))) {

    gr <- which(vld)
    grI       <- gr[which(get("tap")[gr] < mn_tap | get("per")[gr] < mn_per)]
    grII      <- gr[which(get("tap")[gr] >= mn_tap & get("per")[gr] >= mn_per)]
    if (!identical(grII, integer(0))) {
      grIIA     <- grII[which(get("abt")[grII] < plTL)]
      grIIB     <- grII[which(get("abt")[grII] >= plTL)]
    } else {
      grIIA <- grIIB <- integer(0)
    }
    if (!identical(grIIB, integer(0))) {
      grIIB1    <- grIIB[which(sapply(psblVegCls[grIIB], function(x) { length(x) == 1 }))]
      grIIB2    <- grIIB[which(sapply(psblVegCls[grIIB], function(x) { length(x) > 1 }))]
    } else {
      grIIB1 <- grIIB2 <- integer(0)
    }
    if (!identical(grIIB2, integer(0))) {
      grIIB2a   <- grIIB2[which(sapply(psblVegCls[grIIB2], function(x) { !any(is.element(x, c(WtCls, StCls))) }))]
      grIIB2b   <- grIIB2[which(sapply(psblVegCls[grIIB2], function(x) { any(is.element(x, c(WtCls, StCls))) }))]
    } else {
      grIIB2a <- grIIB2b <- integer(0)
    }
    if (!identical(grIIB2b, integer(0))) {
      grIIB2bwt <- grIIB2b[which(get("abt")[grIIB2b] < frTL)]
      grIIB2bst <- grIIB2b[which(get("abt")[grIIB2b] >= frTL)]
    } else {
      grIIB2bwt <- grIIB2bst <- integer(0)
    }

    # I. Bare soil and no vegetation
    if (!identical(grI, integer(0))) {
      vegCls[grI]       <- "BaSl"
    }
    # II.A Polar desert
    if (!identical(grIIA, integer(0))) {
      vegCls[grIIA]     <- "PD"
    }
    # II.B.1 One possible valid vegetation type
    if (!identical(grIIB1, integer(0))) {
      vegCls[grIIB1]    <- sapply(psblVegCls[grIIB1], names)
    }
    # II.B.2.a Two or more possible valid vegetation types (without warm temperate and subtropical regions)
    if (!identical(grIIB2a, integer(0))) {
      vegCls[grIIB2a]   <- sapply(psblVegCls[grIIB2a], function(x) { names(x)[which(x == min(x))] })
    }
    # II.B.2.bwt Two or more possible valid vegetation types (with warm temperate region)
    if (!identical(grIIB2bwt, integer(0))) {
      vegCls[grIIB2bwt] <- sapply(psblVegCls[grIIB2bwt], function(x) {
        names(x[!is.element(x, StCls)])[which(x[!is.element(x, StCls)] == min(x[!is.element(x, StCls)]))]
      })
    }
    # II.B.2.bst Two or more possible valid vegetation types (with subtropical region)
    if (!identical(grIIB2bst, integer(0))) {
      vegCls[grIIB2bst] <- sapply(psblVegCls[grIIB2bst], function(x) {
        names(x[!is.element(x, WtCls)])[which(x[!is.element(x, WtCls)] == min(x[!is.element(x, WtCls)]))]
      })
    }

  }

  # ~~~~ RETURN VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  if (verbose) {
    rslt <- data.frame(do.call(cbind, mget(cv.bci)), vegCls = vegCls)
  } else {
    rslt <- data.frame(vegCls = vegCls)
  }
  return(rslt)

}


#' Vegetation Classifier Using the HLZ System
#'
#' @description Calculates the values of bioclimatic indices used in the Holdridge life zone (HLZ) system
#'     (Holdridge 1947, 1967), and designates the HLZ type using these values, for a given region, by using the
#'     monthly time series of temperature and precipitation.
#'
#' @param rs.temp multi-layer Raster* object with one-year time series of monthly mean air temperature (in °C)
#' @param rs.prec multi-layer Raster* object with one-year time series of monthly precipitation sum (in mm)
#' @param verbose 'logical' scalar that indicates whether or not values of the bioclimatic indices used should be
#'     added to the output.
#' @param filename output filename
#' @param ... additional arguments passed on to \code{\link[raster]{writeRaster}}
#'
#' @details See \code{\link[macroBiome]{cliHoldridgePoints}}.
#'
#' @return Depending on the setting, a RasterStack with one or more layers where the numeric integers encoding the
#'     HLZ type are stored at the last layer, while the additional layers contain the values of bioclimatic indices
#'     used. The meaning of integers is given in the data frame \code{\link[macroBiome]{vegClsNumCodes}}. If
#'     \code{verbose = FALSE}, the return object is a single-layer RasterStack with numeric integers encoding the HLZ
#'     type.
#'
#' @note The objects \code{'rs.temp'} and \code{'rs.prec'} must be 12-layer Raster* objects. These Raster* objects
#'     must have the same bounding box, projection, and resolution.
#'
#' @references
#'
#' \cite{Holdridge LR (1947) Determination of World Plant Formations From Simple Climatic Data. Science
#'     105(2727):367–368. \doi{10.1126/science.105.2727.367}}
#'
#' \cite{Holdridge LR (1967) Life zone ecology. Tropical Science Center, San Jose, Costa Rica}
#'
#' @examples
#' # Loading mandatory data for the Example 'Climate Normal Grid'
#' data(inp_exClnrGrid)
#'
#' # Designate the HLZ types (using the related bioclimatic indices)
#' # for Csongrad-Csanad County (for the normal period 1981-2010)
#' with(inp_exClnrGrid, {
#' rs.HLZ <- cliHoldridgeGrid(temp, prec, verbose = TRUE)
#' rs.HLZ
#' })
#'
#' @import raster
#'
#' @export
#'
cliHoldridgeGrid <- function(rs.temp, rs.prec, verbose = FALSE, filename = "", ...) {

  # ~~~~ FUNCTION WARNINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  cv.arg <- c("rs.temp", "rs.prec")
  for (i in 1 : length(cv.arg)) {
    if (is.null(get(cv.arg[i]))) { stop("Invalid argument: '", cv.arg[i], "' is missing, with no default.") }
  }

  errorCheckingGrid(rs.temp = rs.temp, rs.prec = rs.prec)


  # ~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

  n_lyr <- ifelse(verbose, 4, 1)

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
      df.rslt <- cliHoldridgePoints(temp, prec, verbose = verbose)
      numCode <- hlzDefSubset$Numeric.code[match(df.rslt[["vegCls"]], rownames(hlzDefSubset))]
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
      df.rslt <- cliHoldridgePoints(temp, prec, verbose = verbose)
      numCode <- hlzDefSubset$Numeric.code[match(df.rslt[["vegCls"]], rownames(hlzDefSubset))]
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
