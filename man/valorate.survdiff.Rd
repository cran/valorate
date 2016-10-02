% File man/valorate.survdiff.Rd
\name{valorate.survdiff}
\alias{valorate.survdiff}
\title{ESTIMATES THE P-VALUE AND STATISTICS OF THE LOG-RANK TEST}
\description{
  Estimates the p-value using the VALORATE calculation. 
}
\usage{
valorate.survdiff(vro, clusters, p.func)
}
\arguments{
  \item{vro}{the valorate object.}
  \item{clusters}{a numerical or logical vector representing the two survival groups encoded in 1/TRUE for those 'mutated' (in the group of interest) or 1/FALSE for those who not. Basically this value is the 'x' vector in the VALORATE re-formulation. See references.}
  \item{p.func}{the function that provides the estimation. The default is \link{valorate.p.value.sampling}. See \link{valorate.p.value}.}
}
\details{
The main function to estimate the p-value of the difference of two survival curves under the VALORATE algorithm. Because the definition of the survival group as '1' or '0' is arbitrary, the actual calculation is performed over the less frequent group. From clusters and s=sum(clusters), n1 is determined as the min(s, length(clusters)-s). Then a call for prepare.n1 is performed and finally the p-value estimations. The statistics are added as attributes.
}
\value{the estimated p-value.
}
\references{
	Trevino et al. 2016 \url{http://bioinformatica.mty.itesm.mx/valorate}
	}
\author{Victor Trevino \email{vtrevino@itesm.mx}}
\seealso{
  \code{\link{new.valorate}}.
  \code{\link{valorate.p.value}}.
  \code{\link{valorate.plot.empirical}}.
}
\examples{
## Create a random population of 100 subjects 
## having 20 events
subjects <- numeric(100)
subjects[sample(100,20)] <- 1
vo <- new.valorate(rank=subjects, sampling.size=100000, verbose=TRUE)

groups <- numeric(100)
groups[sample(100,4)] <- 1  # only 4 subjects are within the 'mutated' group
pvr <- valorate.survdiff(vo, groups) 
print(pvr)

}
