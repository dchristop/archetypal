\name{align_archetypes_from_list}
\alias{align_archetypes_from_list}
\title{Align archetypes from a list either by the most frequent found or by using a given archetype}
\description{
Align archetypes from a list either by the most frequent or by using a given archetype
}
\usage{
align_archetypes_from_list(archs_list, given_arch = NULL,
  varnames = NULL, ndigits = 0, parallel = FALSE, nworkers = NULL)
}
\arguments{
\item{archs_list}{The list of archetypes that must be aligned}

\item{given_arch}{If it is not NULL, then the given archetype will by used
as guide for aligning other archetypes of list. Otherwise, a heuristic
for finding the most frequent archetype will be used.}

\item{varnames}{The character vector of variable names that must be used.
If it is NULL, then the column names of first archetype will be used.}

\item{ndigits}{The character of digits that will be used for truncation.}

\item{parallel}{If it set to TRUE, then parallel processing will be applied.}

\item{nworkers}{The number of logical processors that will be used for
parallel computing (usually it is the double of available physical cores)}
}
\value{
A list with members:
\enumerate{
  \item arch_guide, the archetype used as guide for aligning others
  \item phrases_most, a table with all rounded phrases from archetypes.
    Frequencies are in decreasing order,so first row indicates the
    most frequent sequence, if exists. Otherwise we take randomly
    a case and proceed.
  \item archs_aa_output, a data frame with rows all given archetypes
  \item archs_aligned, the final list of aligned archetypes    
}
}
\examples{
data("wd2") #2D demo 
df=wd2
# Define 4 archetypes found for it
dalist=list(c(2.172991,3.200754,5.384013,2.579770,4.860343,3.085111),
            c(5.430821,3.128493,2.043495,3.146342,4.781851,2.710885),
            c(5.430752,2.043403,3.128520,3.146252,2.710979,4.781880),
            c(2.043854,5.430890,3.127183,2.710522,3.146432,4.780432))
archslist=lapply(dalist, function(x){matrix(x,ncol=2)}) #not aligned
# Run aligner
yy=align_archetypes_from_list(archs_list = archslist,given_arch = archslist[[1]])
yy$arch_guide
aligned_archs=yy$archs_aligned
aligned_archs #observe thatn they are comparabnle now
}
\references{
This function is a modifcation of "align_arc" function from package "ParetoTI",
see \url{https://github.com/vitkl/ParetoTI} and \url{https://github.com/vitkl/ParetoTI/blob/master/R/align_arc.R}
}
