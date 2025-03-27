% File man/valorate.plot.subpop.empirical.Rd
\name{valorate.plot.subpop}
\alias{valorate.plot.subpop.}
\alias{valorate.plot.subpop.empirical}
\alias{valorate.plot.subpop.empirical.to.0}
\alias{valorate.plot.subpop.empirical.scaled}
\title{PLOT ALL ESTIMATED LOG-RANK DISTRIBUTIONS}
\description{
  Plots all log-rank distributions estimated with the same object (different values of n1). This family of plots is commonly used to compare the estimated distributions.
}
\usage{
valorate.plot.subpop.empirical(vro, which, 
  type, log, xlim, smooth, legends, 
  density, ylim, ...)

valorate.plot.subpop.empirical.to.0(vro, which, 
  type, log, xlim, smooth, legends, density, ylim, ...)

valorate.plot.subpop.empirical.scaled(vro, which, 
  type, log, xlim, 
  smooth, legends, density, ylim, 
  scale.point, ...)
}
\arguments{
  \item{vro}{the valorate object.}
  \item{which}{The values of n1 that will be shown. NULL to plot them all.}
  \item{type}{typical plot parameter: "p"=points, "l"=lines (default), "o"=overlap.}
  \item{log}{typical plot parameter : specify which axis are shown in logarithm base 10.}
  \item{xlim}{typical plot parameter.}
  \item{ylim}{typical plot parameter.}
  \item{smooth}{the strength of density smoothing for display purposes. The default is 10.}
  \item{legends}{the number of columns in legends. 0 to omit legends.}
  \item{density}{indicates whether all curves should represent density (default to TRUE). FALSE to scale to maximum.}
  \item{...}{arguments passed to plot.}
  \item{scale.point}{a double between 0 and 0.5 (exclusive) that determines the two points in quantiles in which all densities will be 'equalized'. The quantiles are scale.point and 1-scale.point. 0.5 should be avoided.}
}
\details{
\code{valorate.plot.subpop.empirical} plots all log-rank distributions estimated with the same object (different values of n1) in raw densities and scales.
\code{valorate.plot.subpop.empirical.to.0} is similar to \code{valorate.plot.subpop.empirical} but shift distributions to 0 and scale horizontal axis to similar limits. 
\code{valorate.plot.subpop.empirical.scaled} is similar to \code{valorate.plot.subpop.empirical} but scales the distributions to have the same scale.point(s) (in x) for all distributions. It also shifts all distribution to zero. This helps to compare the tendencies of the overall distributions. 
}
\value{
Nothing.
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
vo <- new.valorate(rank=subjects, sampling.size=100000, verbose=TRUE)

for (i in c(5,10,20,50)) {
  groups <- numeric(100)
  groups[sample(100,i)] <- 1 
  valorate.survdiff(vo, groups) 
}

\dontrun{valorate.plot.subpop.empirical(vo)}
\dontrun{valorate.plot.subpop.empirical.to.0(vo)}

}
