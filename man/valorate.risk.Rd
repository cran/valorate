% File man/valorate.risk.Rd
\name{valorate.risk}
\alias{valorate.risk}
\title{ESTIMATES RISK}
\description{
  Estimates the risk (hazard ratio), and confidence interval of a 'mutated' group.
}
\usage{
valorate.risk(vro, clusters)
}
\arguments{
  \item{vro}{the valorate object.}
  \item{clusters}{a numerical or logical vector representing the two survival groups encoded in 1/TRUE for those 'mutated' (in the group of interest) or 1/FALSE for those who not. Basically this value is the 'x' vector in the VALORATE re-formulation. See references.}
}
\details{
A coxph model depending on clusters is run to establish the risk/hazard ratio.
}
\value{
A number representing the relative risk. The confidence interval, p-value, and the coxph model are included as attributes.
}
\references{
	Trevino et al. 2016 \url{http://bioinformatica.mty.itesm.mx/valorate}
	}
\author{Victor Trevino \email{vtrevino@itesm.mx}}
\seealso{
  \code{\link{new.valorate}}.
  \code{\link{valorate.survdiff}}.
  \code{\link{coxph}} (survival package).
}
\examples{
## Create a random population of 100 subjects 
## having 20 events
subjects <- numeric(100)
subjects[sample(100,20)] <- 1
vo <- new.valorate(rank=subjects, sampling.size=100000, verbose=TRUE)

groups <- numeric(100)
groups[sample(100,20)] <- 1  # 20 to likely see some difference
pvr <- valorate.survdiff(vo, groups) 
print(pvr)

valorate.risk(vo, groups)

}
