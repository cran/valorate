% File man/valorate.plot.kaplan.Rd
\name{valorate.plot.kaplan}
\alias{valorate.plot.kaplan}
\title{PLOT KAPLAN-MEIER CURVES}
\description{
  Plots the Kaplan-Meier Curves from two groups.
}
\usage{
valorate.plot.kaplan(vro, clusters, p=valorate.survdiff(vro, clusters), 
  main, short.names=TRUE, draw.all=FALSE, mark="|", mark.cex=0.75, 
  margins=TRUE, col=2:3, col.all="skyblue")
}
\arguments{
  \item{vro}{the valorate object.}
  \item{clusters}{a numerical or logical vector representing the two survival groups encoded in 1/TRUE for those 'mutated' (in the group of interest) or 1/FALSE for those who not. Basically this value is the 'x' vector in the VALORATE re-formulation. See references.}
  \item{p}{the estimated p-value of the log-rank test. The default is valorate.survdiff(vro, clusters).}
  \item{main}{typical plot parameter. The default is an expression depending on the parameters. }
  \item{short.names}{if TRUE (default) use 'LR' instead or 'Log-Rank' and 'HR' instead or 'Hazard-Ratio' in legends.}
  \item{draw.all}{if TRUE, the plot includes also the survival curve of all subjects before stratification.}
  \item{mark}{character to mark censoring. The default is "|".}
  \item{mark.cex}{the character expansion. The default is 0.75}
  \item{margins}{if TRUE (default) set the margins properly.}
  \item{col}{specifies the colors for survival curves. The default is 2:3 (red for cluster=0, green for cluster=1).}
  \item{col.all}{specifies the color when draw.all is TRUE. The default is "skyblue".}
}
\details{
Plots the estimated Kaplan-Meier survival curves from data.
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

groups <- numeric(100)
groups[sample(100,20)] <- 1  # 20 to likely see some difference
pvr <- valorate.survdiff(vo, groups) 
print(pvr)

\dontrun{valorate.plot.kaplan(vo, groups, main="Two Curves")}

\dontrun{valorate.plot.kaplan(vo, groups, draw.all=TRUE, 
    main="Three Curves (Including All Data)")}

}
