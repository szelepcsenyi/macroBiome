% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{inp_exPoints}
\alias{inp_exPoints}
\title{Mandatory Data for the Example 'Points'}
\format{
A list with three data frames and three numeric scalars.

    Data frames contain one-year monthly time series for the period 1981-2010, for the following climate variables:
    \itemize{
      \item{\code{temp}: mean air temperature (in °C)}
      \item{\code{prec}: precipitation sum (in mm)}
      \item{\code{bsdf}: mean relative sunshine duration (dimensionless)}
    }

    The following geographical parameters can be extracted from the numeric scalars:
    \itemize{
      \item{\code{lat}: latitude coordinate (in decimal degrees)}
      \item{\code{lon}: longitude coordinate (in decimal degrees)}
      \item{\code{elv}: elevation (in meters above sea level)}
    }
}
\description{
Data needed to demonstrate the working of functions with a point mode. Measured monthly time series
    of three basic climate variables retrieved from the CarpatClim database (Spinoni et al. 2015), at a grid cell
    near Szeged, Hungary (46.3° N, 20.2° E), for the period 1981-2010; supplemented with the associated
    geographical data. Monthly mean relative sunshine duration has been obtained as a ratio of monthly total of
    sunshine duration and maximum potential number of sunshine hours under clear-sky conditions. Daylength
    (accumulated hours of daylight) was calculated via Eq 1.6.11 in Duffie and Beckman (1991), according to the
    SPLASH radiation algorithm (see \code{\link{cliAvgDlySolIrrPoints}}).
}
\examples{
data(inp_exPoints)
str(inp_exPoints)

}
\references{
\emph{Duffie JA, Beckman WA (1991) Solar Engineering of Thermal Processes. Second Edition. Wiley-Interscience,
    New York, NY}

\emph{Spinoni J, Szalai S, Szentimrey T, Lakatos M, Bihari Z, Nagy A, Németh Á, Kovács T, Mihic D,
    Dacic M, Petrovic P, Kržič A, Hiebl J, Auer I, Milkovic J, Štepánek P, Zahradnícek P, Kilar P, Limanowka D,
    Pyrc R, Cheval S, Birsan M-V, Dumitrescu A, Deak G, Matei M, Antolovic I, Nejedlík P, Štastný P, Kajaba P,
    Bochnícek O, Galo D,  Mikulová K, Nabyvanets Y, Skrynyk O, Krakovska S, Gnatiuk N, Tolasz R, Antofie T, Vogt J
    (2015) Climate  of  the  Carpathian  Region  in  the  period  1961–2010:  climatologies  and  trends  of  10
    variables. Int J Climatol 35(7):1322-1341. \doi{10.1002/joc.4059}}
}
