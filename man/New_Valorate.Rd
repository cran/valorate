% File man/New_Valorate.Rd
\name{new.valorate}
\alias{valorate}
\alias{new.valorate}
\title{CREATE A VALORATE OBJECT}
\description{
  Creates a new valorate object from the survival information
  and basic parameters.
  \code{new.valorate}.
}
\usage{
new.valorate(time, status, censored, rank, sampling.size=max(10000, 2e+05/events), 
	min.sampling.size=1000, tails=2, sampling.ties=30, 
	weights.method=c("logrank", "Wilcoxon", "Tarone-Ware", "Peto", 
		"Flemington-Harrington", "Trevino", "user")[1], 
	weights.parameters=list(p=1, q=1, t=3), weights, 
    verbose=FALSE, save.sampling=TRUE, method="C", 
    estimate.distribution.parameters=c("empirical","gaussian","beta","weibull")[1])
}
\arguments{
  \item{time}{character or numeric vector representing the survival time. If character it could be "76+" representing 76 units of time and censored. In this case, the censored parameter should not be provided.}
  \item{status}{numeric or logical vector representing the status (1 for event and 0 for censoring). This should be specified in the same order than the time argument and only if time was specified and does not include the censoring "+" indicator within.}
  \item{censored}{numeric or logical vector representing the censoring status (1 for censored and 0 for event). This should be specified in the same order than the time argument and only if time was specified and does not include the censoring "+" indicator within.}
  \item{rank}{a numeric or logical vector representing time-ordered subjects and whether they are events (1 or TRUE) or censoring observations (0 or FALSE). If rank is provided, time and status, censored should not. This argument is basically the 'c' vector of the log-rank formulation in the VALORATE publication.}
  \item{sampling.size}{a numeric value representing the length of random samples of the survival group vector (basically the 'x' vector, see the publication and references) that will be used to estimate the log-rank distribution.  See the details section. The default is max(10000,200000/events).}
  \item{min.sampling.size}{a numeric value representing the minimum number of random samples of the survival group vector. See the details section. The default is 1000.}
  \item{tails}{a numeric value indicating whether the p-values generated will represent a 1 tail or two-tails (the default is two-tails).}
  \item{sampling.ties}{a numeric value indicating the number of permutations of tie positions used for the estimation of the log-rank distribution. The default is 30.}
  \item{weights.method}{a character specifying the type of log-rank test (see books in references). It can be "logrank" (default), "Wilcoxon", "Tarone-Ware", "Peto", "Flemington-Harrington", "Trevino",and "user". In case of "user", the 'weights' parameter should also be specified }
  \item{weights.parameters}{a list of values for 'weights.method'. "Flemington-Harrington" uses a p and q parameters. "Trevino" uses a t parameter. The default is list(p=1,q=1,t=3).}
  \item{weights}{a numeric vector having order according to the 'time' or 'rank' parameter only for "user" weights.method}
  \item{verbose}{a logical value indicating whether estimation should show messages of partial calculations. The default is FALSE.}
  \item{save.sampling}{a logical value indicating whether all sampling will be saved within the VALORATE object. The default is TRUE. This can be used to avoid saving all sampling and save memory. See details about memory usage.}
  \item{method}{a character value either "R" or "C" that specify the implemented method of calculation. Both should generate same values but "C" is by far faster (default="C"). This can be used in cases where C calculations does not work for any reason or to compare methods and algorithms.}
  \item{estimate.distribution.parameters}{a character vector containing subsets of "empirical","gaussian","beta",and "weibull". The default is "empirical". This has not been extensively explored but attempts to fit the observed log-rank distribution using sums of other distributions whose parameters are estimated after sampling. The "empirical" means nothing basically whereas "gaussian" for example means the estimation of mean and standard deviation of the observed conditional log-ranks. The results of these estimations can be viewd by \link{valorate.plot.empirical} or within the corresponding variables of the @subpop environment. This is experimental and is not intended for most users and applictions.}
}
\details{The values in time and censored arguments do not need to be sorted by time but it is assumed that further specification of survival groups (e. g. calls to valorate.survdiff) will be provided in the corresponding order than that of this parameter. This function generates a VALORATE object prepared to run further analyses on the estimation of log-rank distributions for the specified population. See the 'value' section for details. 

It is critical for computation time and memory the handling of the sampling. This can be managed by the sampling.size and min.sampling.size.

To save memory, save.sampling can be set to FALSE. In this case, valorate will estimate a ~1000 breaks histogram to store each conditional log-rank distribution. However, this strategy will lead in loosing resolution and therefore precision in the estimation of p-values. Thus, save.sampling=FALSE is not recommended for most applications.
}
\value{A valorate object.
  \item{s}{numeric vector representing the subjects ordered in time. This is basically the same than 'rank' argument if specified. So, this is the 'c' vector of the log-rank formulation in the VALORATE publication.}
  \item{n}{the total number of subjects. It should be equal to the length of s.}
  \item{events}{the total number of events observed. It should be equal to the sum of s.}
  \item{parameters}{this is a list of the parameters specified.}
  \item{sampling.size}{the 'total' number of sampling used to estimate the log-rank distribution. The computation time and memory depends largely in this value. Many numeric vectors will be created whose sum of their length will be approximately this value. This may be a copy of the original argument.}
  \item{min.sampling.size}{the minimum number of sampling that will be used to estimate a conditional log-rank distribution (conditional to 'k' co-occurrences, see publication and references below).  This may be a copy of the original argument.}
  \item{wcensored}{a numeric vector denoting the positions of the 's' vector having censored subjects.}
  \item{wevents}{a numeric vector denoting the positions of the 's' vector having events.}
  \item{order}{the index positions of the time/censoring values needed to sort the subjects by time. This will be used in \link{valorate.survdiff} to re-accomodate the data specified.}
  \item{subpop}{an environment of currently estimated log-rank distributions. The names for each item is given by 'subpop#' where '#' is the number of subjects in the survival group of interest (basically the n1 value that is equal to the sum of 1's within the 'x' vector). Each 'subpop' contains a list of values needed for the estimations and many are further indexed by the value of co-occurrences 'k', including 'sampling' which stores all log-rank conditional samplings,  'emp.hist' a tiny 'histogram' version of the estimated conditional distribution, 'combinations' the number of combinations of each co-occurrence, and 'k.density' its corresponding density or weights. See the tutorial within references for details.}
  \item{ties}{a list of vectos having the positions of ties.}
  \item{sampling.ties}{this is a copy of the original argument.}
  \item{tiesame}{equal to ties.}
  \item{tiesame.pos}{all positions having ties.}
  \item{tiesame.sampling}{the actual value of samplings done to ties. 1 means no additional samplings.}
  \item{verbose}{this is a copy of the original argument.}
  \item{save.sampling}{this is a copy of the original argument.}
  \item{time}{this is a copy of the original argument.}
  \item{tails}{this is a copy of the original argument.}
  \item{weights.method}{this is a copy of the original argument.}
  \item{weights}{this is a copy of the original argument.}
  \item{method}{this is a copy of the original argument.}
  \item{samplings}{this is depracated and has been moved to each subpop.}
}
\references{
	Trevino et al. 2016 \url{http://bioinformatica.mty.itesm.mx/valorateR}

	David G. Kleinbaum and Mitchel Klein (2005). \emph{Survival Analysis: A Self-Learning Text}. Second Edition. New York: Springer.

	David Collett (2004). \emph{Modelling survival data in medical research Collett David}. Second Edition. Chapman & Hall-CRC.
	}
\author{Victor Trevino \email{vtrevino@itesm.mx}}
\seealso{
  \code{\link{valorate.survdiff}}.
  \code{\link{valorate.plot.empirical}}.
}
\examples{
## Create a random population of 100 subjects 
## having 20 events
subjects <- numeric(100)
subjects[sample(100,20)] <- 1
vo <- new.valorate(rank=subjects, sampling.size=100000)

## print the structure of properties
str(vo)

## print slots
slotNames(vo)

}
