% File man/prepare.n1.Rd
\name{prepare.n1}
\alias{prepare.n1}
\title{ESTIMATES THE LOG-RANK DISTRIBUTION AND STORE IT WITHIN A VALORATE OBJECT}
\description{
  This method  estimates the log-rank distribution for a minor size group equal to n1. 
}
\usage{
prepare.n1(vro, n1)
}
\arguments{
  \item{vro}{the valorate object.}
  \item{n1}{a numerical value (perhaps integer) of the size of the survival group.}
}
\details{This method actually estimates the log-rank distribution for a minor size group equal to n1. This method is internally called by all functions to first check and/or compute the log-rank distribution of n1. It is not intended to be used by final users unless it is intended to prepare the VALORATE object before p-value calculations (perhaps in separated threds or jobs and saving/restoring the object). 
}
\value{The updated valorate object. This is a S4 method operating an object, so the valorate object specified in argument will be updated with the estimations of the distribution of the log-rank for n1 subjects.
}
\references{
	Trevino et al. 2016 \url{http://bioinformatica.mty.itesm.mx/valorateR}
	}
\author{Victor Trevino \email{vtrevino@itesm.mx}}
\seealso{
  \code{\link{new.valorate}}.
  \code{\link{valorate.survdiff}}.
}
\examples{
## Create a random population of 100 subjects 
## having 20 events
subjects <- numeric(100)
subjects[sample(100,20)] <- 1
vo <- new.valorate(rank=subjects, sampling.size=100000)
# and with verbose
vo <- new.valorate(rank=subjects, sampling.size=100000, verbose=TRUE)

## print the initial subpop
ls(vo@subpop) # should be character(0)

## calculate
prepare.n1(vo, 10) # should show messages of partial calculations P(L|k)

## print the current subpop
ls(vo@subpop) #should show now: [1] "subpop10"

names(vo@subpop[["subpop10"]]) #should show the internal names of the estimated subpop

}
