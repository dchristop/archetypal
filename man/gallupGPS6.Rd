\name{gallupGPS6}
\alias{gallupGPS6}
\docType{data}
\title{
Gallup Global Preferences Study processed data set of six variables
}
\description{
A 76132 x 6 data frame derived from Gallup Global Preferences Study, see [1] and [2] for details.
It can be used as a big data set example. 
}
\usage{data("gallupGPS6")}
\format{
  A data frame with 76132 complete observations on the following 6 variables.
  \describe{
    \item{\code{patience}}{a numeric vector}
    \item{\code{risktaking}}{a numeric vector}
    \item{\code{posrecip}}{a numeric vector}
    \item{\code{negrecip}}{a numeric vector}
    \item{\code{altruism}}{a numeric vector}
    \item{\code{trust}}{a numeric vector}
  }
}
\details{
Data processing:
\enumerate{
  \item The non complete rows have been removed
  \item The duplicated rows have also been removed
  }
}
\source{
Individual data set was downloaded from \url{https://www.briq-institute.org/global-preferences/downloads}, last accessed 2019-12-17.
}
\references{
[1] Falk, A., Becker, A., Dohmen, T., Enke, B., Huffman, D., & Sunde, U. (2018). Global evidence on economic preferences. Quarterly Journal of Economics, 133 (4), 1645-1692.\cr

[2] Falk, A., Becker, A., Dohmen, T. J., Huffman, D., & Sunde, U. (2016). The preference survey module: A validated instrument for measuring risk, time, and social preferences. IZA Discussion Paper No. 9674. \cr
}
\examples{
data(gallupGPS6)
summary(gallupGPS6)
}
\keyword{datasets}
