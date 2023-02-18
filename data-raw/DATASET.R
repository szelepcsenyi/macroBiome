## Load raw data from .csv file

# A corrected version of the function 'palinsol::ber90'
pi.ber90 <- palinsol::ber90
body(pi.ber90)[[2]][[3]][[3]] <- substitute(data("BER90", package = "palinsol", envir = .BER90))

# Auxiliary dataframe for checking and correcting objects containing input variables
varFeatures <- read.csv2("data-raw/functions/varFeatures.csv", header = TRUE, row.names = 1,
                         sep = ";", quote = "\"", dec = ".", fill = TRUE, comment.char = "")

# Auxiliary dataframe for defining constants used in the SPLASH algorithm
constants <- read.csv2("data-raw/functions/constants.csv", header = TRUE, sep = ";", quote = "\"",
                       dec = ".", fill = TRUE, comment.char = "")
c <- as.list(setNames(constants$Value, constants$Variable))

# Auxiliary dataframe for selecting outputs of the SPLASH algorithm
opVarChoices <- read.csv2("data-raw/functions/opVarChoices.csv", header = TRUE, sep = ";", quote = "\"",
                          dec = ".", fill = TRUE, comment.char = "")

# An auxiliary dataframe for checking the availability of inputs to calculate bioclimatic indices
ipVarRequirements <- read.csv2("data-raw/functions/ipVarRequirements.csv", header = TRUE, row.names = 1,
                               sep = ";", quote = "\"", dec = ".", fill = TRUE, comment.char = "")

# Auxiliary dataframe for performing vegetation classification that is based on the Holdridge life zone system
hlzDefinitions <- read.csv2("data-raw/functions/vegChlzDefinitions.csv", header = TRUE, row.names = 1,
                            sep = ";", quote = "\"", dec = ".", fill = TRUE, comment.char = "")

# Auxiliary dataframe for the simplified version of the Holdridge life zone system
# Subselect the vegetation classes used
hlzDefSubset <- subset(hlzDefinitions[hlzDefinitions$Type == "core", ], select = -Forest.steppe)
hlzDefSubset <- data.frame(ID = row.names(hlzDefSubset), hlzDefSubset)

# Define some synthetic vegetation classes
synAbtLim <- unique(hlzDefSubset$abt)
intDf <- data.frame()
for (i_lm in 1 : length(synAbtLim)) {
  for (i_sd in 1 : 2) {
    intDf <-
      rbind(intDf,
            cbind(ID = "BaSl", Vegetation.class = "Bare soil and no vegetation", Numeric.code = 39,
                  Type = "syn", abt = synAbtLim[i_lm],
                  tap = get(c("min", "max")[i_sd])(hlzDefSubset$tap[hlzDefSubset$abt == synAbtLim[i_lm]]) *
                    c(0.5, 2.)[i_sd],
                  per = get(c("max", "min")[i_sd])(hlzDefSubset$per[hlzDefSubset$abt == synAbtLim[i_lm]]) *
                    c(2., 0.5)[i_sd]))
  }
}
hlzDefSubset <- rbind(hlzDefSubset,
                      transform(intDf, Numeric.code = as.integer(Numeric.code),
                                abt = as.numeric(abt), tap = as.numeric(tap), per = as.numeric(per)))

# Auxiliary dataframe for the modified version of the Holdridge life zone system
fspChlzDefSubset <- thrSft <- hlzDefinitions[hlzDefinitions$Forest.steppe == TRUE, ]
difAltCrcRad <- 0.5 - (1. / 3. * sqrt(3.) * (0.5 / (0.5 * sqrt(3.))))
thrSft[, c("abt", "tap", "per")] <- NA
thrSft[thrSft$Type == "core", c("abt", "tap", "per")] <- rep(0.5, 3)
thrSft[thrSft$Type == "tran_p", c("abt", "tap", "per")] <- c(1 - difAltCrcRad, difAltCrcRad, difAltCrcRad)
thrSft[thrSft$Type == "tran_m", c("abt", "tap", "per")] <- c(difAltCrcRad, 1. - difAltCrcRad, 1. - difAltCrcRad)

# Auxiliary dataframes for performing vegetation classification by using the BIOME model
bioPFTDefinitions <- read.csv2("data-raw/functions/vegCbioPFTDefinitions.csv", header = TRUE, row.names = 1,
                               sep = ";", quote = "\"", dec = ".", fill = TRUE, comment.char = "")

bioBiomeDefinitions <- read.csv2("data-raw/functions/vegCbioBiomeDefinitions.csv", header = TRUE, row.names = 1,
                                 sep = ";", quote = "\"", dec = ".", fill = TRUE, comment.char = "")

# Auxiliary dataframes for performing vegetation classification by using the Siberian Vegetation Model
svmDefinitions <- read.csv2("data-raw/functions/vegCsvmDefinitions.csv", header = TRUE, row.names = 1,
                            sep = ";", quote = "\"", dec = ".", fill = TRUE, comment.char = "")

# An auxiliary dataframe for checking the availability of vegetation classification schemes
bciRequirements <- read.csv2("data-raw/functions/bciRequirements.csv", header = TRUE, row.names = 1,
                             sep = ";", quote = "\"", dec = ".", fill = TRUE, comment.char = "")

# An auxiliary dataframe decoding outputs of the Köppen-Geiger classification scheme
vegCkgcTypes <- read.csv2("data-raw/functions/vegCkgcTypes.csv", header = TRUE,
                          sep = ";", quote = "\"", dec = ".", fill = TRUE, comment.char = "")

# Auxiliary dataframe for decoding outputs of vegetation classifiers
vegClsNumCodes <- data.frame(matrix(nrow = 39, ncol = 6,
                                    dimnames = list(NULL, c("Name.HLZ", "Code.HLZ", "Name.KGC", "Code.KGC",
                                                            "Name.BIOME", "Code.BIOME"))))
vegClsNumCodes[1, "Name.HLZ"] <- "Polar desert"
vegClsNumCodes[hlzDefSubset$Numeric.code[1:37], "Name.HLZ"] <- hlzDefSubset$Vegetation.class[1:37]
vegClsNumCodes[39, "Name.HLZ"] <- "Bare soil and no vegetation"
vegClsNumCodes[1, "Code.HLZ"] <- "PD"
vegClsNumCodes[hlzDefSubset$Numeric.code[1:37], "Code.HLZ"] <- rownames(hlzDefSubset)[1:37]
vegClsNumCodes[39, "Code.HLZ"] <- "BaSl"
vegClsNumCodes[1:30, "Name.KGC"] <- vegCkgcTypes$NAME.KGC
vegClsNumCodes[1:30, "Code.KGC"] <- vegCkgcTypes$CODE.KGC
vegClsNumCodes[bioBiomeDefinitions$Numeric.code, "Name.BIOME"] <- bioBiomeDefinitions$Vegetation.class
vegClsNumCodes[bioBiomeDefinitions$Numeric.code, "Code.BIOME"] <- rownames(bioBiomeDefinitions)


# Two auxiliary SpatialPolygonsDataFrames to perform the regional classification required
# to select the regional functions used in the BSDF estimation scheme proposed by Yin (1999)
# countriesSP
# countriesSP <- rworldmap::getMap(resolution = 'high')

# islandsSP
mnrIslandsSP <- rnaturalearth::ne_download(scale = 10, type = "minor_islands", category = "physical")

scale <- 10
type <- "geography_regions_polys"
category <- "physical"
file_name <-  paste0('ne_', scale, 'm_', type)
address <- paste0("https://www.naturalearthdata.com/http//",
                  "www.naturalearthdata.com/download/", scale, "m/", category,"/", file_name, ".zip")
utils::download.file(file.path(address), zip_file <- tempfile())
utils::unzip(zip_file, exdir =  tempdir())
landsSP <- rgdal::readOGR(tempdir(), file_name, encoding = 'UTF-8', stringsAsFactors = FALSE, use_iconv = TRUE)
mjrIslandsSP <- landsSP[landsSP$FEATURECLA == c("Island"), ]
names(mjrIslandsSP)[1] <- "featurecla"

islandsSP <- rbind(mjrIslandsSP[,"featurecla"], mnrIslandsSP[,"featurecla"], makeUniqueIDs = TRUE)


# Mandatory data for the Example 'Points'
cv.var <- c("temp", "prec", "bsdf")
inp_exPoints <- setNames(lapply(vector("list", length(cv.var)), function(x) list()), cv.var)
for (i in 1 : length(cv.var)) {
  tmp <- read.csv2(paste0("data-raw/examples/points/", cv.var[i], "_carpatclim_sgly_198101_201012_46.3N_20.2E.csv"),
                   header = TRUE, row.names = 1, sep = ";", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
  inp_exPoints[[i]] <- tmp
  rm(tmp)
}
inp_exPoints <- c(inp_exPoints, lat = 46.3, lon = 20.2, elv = 76.)

# Mandatory data for the Example 'Single-Year Grid'
inp_exSglyGrid <- setNames(lapply(vector("list", (length(cv.var) + 1)), function(x) list()), c(cv.var, "elv"))
for (i in 1 : length(cv.var)) {
  tmp <- read.csv2(paste0("data-raw/examples/grid/", cv.var[i], "_carpatclim_sgly_201001_201012_HUCS.csv"),
                   header = TRUE, sep = ";", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
  inp_exSglyGrid[[i]] <- tmp
  rm(tmp)
}
inp_exSglyGrid$elv <- read.csv2("data-raw/examples/grid/elv_carpatclim_HUCS.csv",
                                header = TRUE, sep = ";", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
inp_exSglyGrid <- lapply(inp_exSglyGrid, raster::rasterFromXYZ, crs = "+proj=longlat +datum=WGS84 +no_defs")

# Mandatory data for the Example 'Climate-Normal Grid'
inp_exClnrGrid <- setNames(lapply(vector("list", (length(cv.var) + 1)), function(x) list()), c(cv.var, "elv"))
for (i in 1 : length(cv.var)) {
  tmp <- read.csv2(paste0("data-raw/examples/grid/", cv.var[i], "_carpatclim_clnr_198101_201012_HUCS.csv"),
                   header = TRUE, sep = ";", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
  inp_exClnrGrid[[i]] <- tmp
  rm(tmp)
}
inp_exClnrGrid$elv <- read.csv2("data-raw/examples/grid/elv_carpatclim_HUCS.csv",
                                header = TRUE, sep = ";", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
inp_exClnrGrid <- lapply(inp_exClnrGrid, raster::rasterFromXYZ, crs = "+proj=longlat +datum=WGS84 +no_defs")


## Save the data in the required R package location

usethis::use_data(pi.ber90, varFeatures, c, opVarChoices, ipVarRequirements, hlzDefSubset, bioPFTDefinitions,
                  bioBiomeDefinitions, fspChlzDefSubset, thrSft, svmDefinitions, bciRequirements, vegCkgcTypes,
                  islandsSP, internal = TRUE, overwrite = TRUE)

usethis::use_data(inp_exPoints, overwrite = T)

usethis::use_data(inp_exSglyGrid, overwrite = T)

usethis::use_data(inp_exClnrGrid, overwrite = T)

usethis::use_data(vegClsNumCodes, overwrite = T)
