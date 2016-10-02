% File man/valorate.p.value.Rd
\name{valorate.p.value}
\alias{valorate.p.value}
\alias{valorate.p.value.sampling}
\alias{valorate.p.value.normal}
\alias{valorate.p.value.chisq}
\alias{valorate.p.value.gaussian}
\alias{valorate.p.value.weibull}
\alias{valorate.p.value.beta}
\alias{valorate.p.value.all}
\title{ESTIMATES THE P-VALUE OF THE LOG-RANK TEST}
\description{
  Estimates the p-value using specific approximations to the log-rank. 
}
\usage{
valorate.p.value.sampling(vro, vrsubo, lrv, z)
valorate.p.value.chisq(vro, vrsubo, lrv, z)
valorate.p.value.normal(vro, vrsubo, lrv, z)
valorate.p.value.gaussian(vro, vrsubo, lrv, z)
valorate.p.value.weibull(vro, vrsubo, lrv, z)
valorate.p.value.beta(vro, vrsubo, lrv, z)
valorate.p.value.all(vro, vrsubo, lrv, z)
}
\arguments{
  \item{vro}{the valorate object.}
  \item{vrsubo}{the subpop list object (see \link{prepare.n1}) or a numeric value representing n1 used to obtain the subpop.}
  \item{lrv}{if provided, the log-rank value. It is needed for .sampling, .gaussian, .weibull, .beta, .normal, and .all .}
  \item{z}{if provided, the log-rank value in z-score (divided by the approximated standard deviation). It is needed for .normal, .chisq, optionally to .all if normal and chisq are required.}
}
\details{
This family of functions estimates the p-value of the log-rank test using specific approximations.
The intended 'user' function in VALORATE is valorate.p.value.sampling, which is the function that is described in the publications. The rest of the functions are complementary for comparison with the classical approximations (chisq and normal) and for experimental purposes fitting each conditional log-rank distribution sampled (conditioned on k co-occurrences) with the specified distribution (gaussian, weibull, and beta). The function valorate.p.value.all is just a proxy to all calculations in the same function.

}
\value{the estimated p-value (times tails).
}
\references{
	Trevino et al. 2016 \url{http://bioinformatica.mty.itesm.mx/valorate}
	}
\author{Victor Trevino \email{vtrevino@itesm.mx}}
\seealso{
  \code{\link{new.valorate}}.
  \code{\link{valorate.survdiff}}.
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

# the same than the value of pvr
valorate.p.value.sampling(vo, vo@subpop[["subpop4"]], attributes(pvr)[[1]]["LR"])

# the same than the value of pvr
valorate.p.value.sampling(vo, 4, attributes(pvr)[[1]]["LR"])

#classical approximations:
valorate.p.value.normal(vo, 4, attributes(pvr)[[1]]["LR"], attributes(pvr)[[1]]["Z"])
valorate.p.value.chisq(vo, 4, attributes(pvr)[[1]]["LR"], attributes(pvr)[[1]]["Z"])

# approximations of the conditional log-rank sampled density
valorate.p.value.gaussian(vo, 4, attributes(pvr)[[1]]["LR"])
valorate.p.value.beta(vo, 4, attributes(pvr)[[1]]["LR"])
valorate.p.value.weibull(vo, 4, attributes(pvr)[[1]]["LR"])

# all above can be get by:
valorate.p.value.all(vo, 4, attributes(pvr)[[1]]["LR"], attributes(pvr)[[1]]["Z"])

# Estimate a p-value a given log-rank
prepare.n1(vo, 50)
valorate.p.value.all(vo, 50, 0, 0) # 0 log-rank, 0 z-score

}
