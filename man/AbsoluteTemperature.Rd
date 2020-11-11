\name{AbsoluteTemperature}
\alias{AbsoluteTemperature}
\docType{data}
\title{
Global Absolute Temperature data set for Northern Hemisphere 1969-2013
}
\description{
It is a subset from the data set which was used for publication [1], i.e. the Global Absolute Temperature for Northern Hemisphere (1800-2013) with only complete yearly observations included. Here we have kept the years 1969-2013.
}
\usage{data("AbsoluteTemperature")}
\format{
  A data frame with 155862 observations on the following 18 variables.
  \describe{
    \item{\code{Year}}{an integer vector of observation years from 1969 to 2013}
    \item{\code{Jan}}{numeric vector of monthly average temperature for January}
    \item{\code{Feb}}{numeric vector of monthly average temperature for February}
    \item{\code{Mar}}{numeric vector of monthly average temperature for March}
    \item{\code{Apr}}{numeric vector of monthly average temperature for April}
    \item{\code{May}}{numeric vector of monthly average temperature for May}
    \item{\code{Jun}}{numeric vector of monthly average temperature for June}
    \item{\code{Jul}}{numeric vector of monthly average temperature for July}
    \item{\code{Aug}}{numeric vector of monthly average temperature for August}
    \item{\code{Sep}}{numeric vector of monthly average temperature for September}
    \item{\code{Oct}}{numeric vector of monthly average temperature for October}
    \item{\code{Nov}}{numeric vector of monthly average temperature for November}
    \item{\code{Dec}}{numeric vector of monthly average temperature for December}
    \item{\code{long}}{a numeric vector for the geographical longitude: positive values for eastings}
    \item{\code{lat}}{a numeric vector for the geographical latitude: positive values for northings}
    \item{\code{h}}{a numeric vector for the altitude in metrs}
    \item{\code{stid}}{an integer vector with the station identity number}
    \item{\code{z}}{an integer vector with the relevant climate zone:
      \itemize{
      \item 1, Tropical Zone
      \item 2, Subtropics
      \item 3, Temperate zone
      \item 4, Cold Zone
      }
    }
  }
}
\details{
That data set was the output of the procedure described in [1]. Initial data set was downloaded from [2]
at 2014-12-17.
}
\examples{
#
######################################
## Load absolute temperature data set:
######################################
#
data("AbsoluteTemperature")
df=AbsoluteTemperature
## Find proportions for climate zones
pcs=table(df$z)/dim(df)[1]
## Choose an approximate size of the new sample and compute resample sizes
N=1000
resamplesizes=as.integer(round(N*pcs))
sum(resamplesizes)
## Create the grouping matrix
groupmat=data.frame("Group_ID"=1:4,"Resample_Size"=resamplesizes)
groupmat
## Simple resampling:
resample_simple <- grouped_resample(in_data = df,grp_vector = "z",
grp_matrix = groupmat,replace = FALSE, option = "Simple", rseed = 20191119)
cat(dim(resample_simple),"\n")
## Dirichlet resampling:
resample_dirichlet <- grouped_resample(in_data = df,grp_vector = "z", 
grp_matrix = groupmat, replace = FALSE, option = "Dirichlet", rseed = 20191119)
cat(dim(resample_dirichlet),"\n")
#
#########################################
## Reproduce the results of 2015 article
#########################################
##
data("AbsoluteTemperature")
dh=AbsoluteTemperature
## Create yearly averages for every station
dh$avg = rowMeans(df[,month.abb[1:12]])
head(dh)
## Compute mean average of every year for all Northern Hemisphere
dagg=data.frame(aggregate(avg~Year,dh,function(x){c(mean(x),sd(x))}))
## Find used stations per year
daggn=aggregate(stid ~ Year,dh,length)
head(daggn)
tail(daggn)
## Combine all in a data frame
dagyears=data.frame(dagg$Year,daggn$stid,dagg$avg[,1],dagg$avg[,2])
colnames(dagyears)=c("Year","Nv","mu","Smu")
head(dagyears)
tail(dagyears)
#
## Compare with Table 7 (Columns: Year, Nv, mu_bar, Smu_bar), page 77 of article
## Extraction of the global absolute temperature for Northern Hemisphere 
## using a set of 6190 meteorological stations from 1800 to 2013
## https://doi.org/10.1016/j.jastp.2015.03.009
## and specifically the years 1969--2013
}
\references{

[1] Demetris T. Christopoulos. Extraction of the global absolute temperature for
Northern Hemisphere using a set of 6190 meteorological stations from 1800 to
2013. Journal of Atmospheric and Solar-Terrestrial Physics, 128:70 - 83, 3 2015.
doi:10.1016/j.jastp.2015.03.009

[2] Met Office Hadley Centre observations datasets, station data sets,\cr
http:///www.metoffice.gov.uk/hadobs/crutem4/data/station_files/CRUTEM.4.2.0.0.station_files.zip \cr
(last visited 17.12.14)

}
\keyword{datasets}
