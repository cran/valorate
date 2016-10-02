% File man/valorate.plot.sampling.densities.Rd
\name{valorate.plot.sampling.densities}
\alias{valorate.plot.sampling.densities}
\alias{valorate.plot.sampling.densities.figure}
\title{PLOT CO-OCCURRENCE DENSITIES FORMING A LOG-RANK DISTRIBUTION}
\description{
  Plots the densities of each co-occurrence that shapes the final log-rank distribution for a n1 group size.
}
\usage{
valorate.plot.sampling.densities(vro, n1, 
  type, log, xlim, ylim, ncol, main, 
  rug, add, w.sum, sampling, 
  weighted,legends.cex, weights.cex,
  weights.pos, w.sum.lwd=3,
  ...)

valorate.plot.sampling.densities.figure(vro, n1, type, log, 
  xlim, ylim, main, rug, rug.size, 
  sub, w.sum, sampling, ncol, ...)

}
\arguments{
  \item{vro}{the valorate object.}
  \item{n1}{the size of the 'mutated' or interested survival group. It can be also the numerical/logical 'x' vector as in \link{valorate.survdiff} in which case n1 and vstat are estimated.}
  \item{type}{typical plot parameter: "p"=points, "l"=lines (default), "o"=overlap.}
  \item{log}{typical plot parameter : specify which axis are shown in logarithm base 10.}
  \item{xlim}{typical plot parameter.}
  \item{ylim}{typical plot parameter.}
  \item{ncol}{number of columns for legends.}
  \item{main}{typical plot parameter. The default is an expression depending on the parameters. }
  \item{rug}{if FALSE removes the drawing of rugs. }
  \item{add}{if FALSE assumes plots are added to existing one. Not valid for all functions. }
  \item{w.sum}{if FALSE removes the drawing of the weighted sum distribution (the final log-rank distibution). }
  \item{sampling}{if TRUE includes the drawing of a crude-histogram version of the overall distribution (the final log-rank distibution). }
  \item{weighted}{if TRUE the densities of each co-occurrence is weighted to its overall contribution (proportion of combinations). Nice! to explain the overall distribution.}
  \item{legends.cex}{the character expansion for legends.}
  \item{weights.cex}{the character expansion for legends of the weights}
  \item{weights.pos}{the position of the weights "middle"=above each curve, "left", or "right".}
  \item{w.sum.lwd}{the line width of the weighted sum line.}
  \item{...}{parameters forwarded to plot.}
  \item{rug.size}{determine the size of rugs made in valorate.plot.sampling.densities.figure}
  \item{sub}{specifies the character to include in each plot.}
}
\details{
Plots the densities of each co-occurrence that shapes the final log-rank distribution for a n1 group size.
\code{valorate.plot.sampling.densities} plots all co-occurrences in the same figure whereas \code{valorate.plot.sampling.densities.figure} plots all in separated figures.
}
\value{
Nothing.
}
\references{
	Trevino et al. 2016 \url{http://bioinformatica.mty.itesm.mx/valorate}
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

\dontrun{valorate.plot.sampling.densities(vo, 5)}
\dontrun{valorate.plot.sampling.densities(vo, 10)}
\dontrun{valorate.plot.sampling.densities(vo, 20, weighted=TRUE)}

\dontrun{valorate.plot.sampling.densities.figure(vo, 5)}

}
