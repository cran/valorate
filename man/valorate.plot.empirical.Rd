% File man/valorate.plot.empirical.Rd
\name{valorate.plot.empirical}
\alias{valorate.plot.empirical}
\title{PLOT THE SAMPLED (EMPIRICAL) LOG-RANK DISTRIBUTION}
\description{
  Plots the estimated density of the log-rank distribution.
}
\usage{
valorate.plot.empirical(vro, n1, vstat, type, log, add, include, xlab, ylab, 
	main, samp, smooth, legends, shades, transparency, lwd, xlim, 
	minL=NA, minR=NA, ...)
}
\arguments{
  \item{vro}{the valorate object.}
  \item{n1}{the size of the 'mutated' or interested survival group. It can be also the numerical/logical 'x' vector as in \link{valorate.survdiff} in which case n1 and vstat are estimated.}
  \item{vstat}{log-rank statistic. If provided, a vertical mark on this value is shown and both sides of the density are filled to highlight areas.}
  \item{type}{typical plot parameter: "p"=points, "l"=lines (default), "o"=overlap.}
  \item{log}{typical plot parameter : specify which axis are shown in logarithm base 10.}
  \item{add}{typical plot parameter : specify whether the plot is new or added to an existing one.}
  \item{include}{specify which other estimations are included. This is experimental. Default "none", possibles: "none","gaussian","beta","weibull","all".}
  \item{xlab}{typical plot parameter: how the x axis is labelled, the default is "valorate LR".}
  \item{ylab}{typical plot parameter: how the y axis is labelled, the default is "density".}
  \item{main}{typical plot parameter. The default is an expression depending on the parameters. }
  \item{samp}{a numeric value specifying the length of sampling when using other density functions. It is related to the include argument.}
  \item{smooth}{the strength of density smoothing for display purposes. The default is 10.}
  \item{legends}{logical value that defines whether the legends for each curve should be displayed. The default is FALSE.}
  \item{shades}{define de colors used to shade the empirical distribution when the either 'vstat' argument is used or n1 represent the 'x' vector. The default is c(6,8). }
  \item{transparency}{defines the transparency in shades. The default is 0.25.}
  \item{lwd}{typical plot parameter: width of the lines. Default is 2.}
  \item{xlim}{typical plot parameter.}
  \item{minL,minR}{limits to the estimated empirical density.}
  \item{...}{arguments passed to plot.}
}
\details{
Plots the estimated density of the log-rank distribution.
}
\value{
An invisible data frame of the density estimation.
}
\references{
	Trevino et al. 2016 \url{http://bioinformatica.mty.itesm.mx/valorateR}
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

# Crude density
\dontrun{valorate.plot.empirical(vo, 4)}

# Similar but marking the statistic 
# returned by groups and shading
\dontrun{valorate.plot.empirical(vo, groups)}

# Plot density and check symmetry
\dontrun{valorate.plot.empirical(vo, 4, 0)}

# Now should be almost symmetric
\dontrun{valorate.plot.empirical(vo, 50, 0)}

# Crude density plus gaussian, weibull, and beta estimations
\dontrun{valorate.plot.empirical(vo, 4, include="all", legends=TRUE)}

}
