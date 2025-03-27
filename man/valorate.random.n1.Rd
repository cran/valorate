% File man/valorate.random.n1.Rd
\name{valorate.random.n1}
\alias{valorate.random.n1}
\title{ESTIMATES THE LOG-RANK DISTRIBUTION OF RANDOM VECTORS}
\description{
  This method estimates the log-rank distribution for random vectors of size n1. 
}
\usage{
valorate.random.n1(vro, n1, simulations=1000)
}
\arguments{
  \item{vro}{the valorate object.}
  \item{n1}{a numerical value (perhaps integer) of the size of the survival group.}
  \item{simulations}{a numerical value of the number of simulations to perform.}
}
\details{This method  estimates the log-rank distribution for a minor size group equal to n1 where the group is randomly and uniformly generated. It is itended to estimate random log-ranks values. 
}
\value{The vector of log-rank values of size equal to simulations.
}
\references{
	Trevino et al. 2016 \url{https://bioinformatics.mx/index.php/bioinfo-tools/}
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

## generate random log-ranks
valorate.random.n1(vo, 10, 100)

}
