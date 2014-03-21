\name{SCR.area}
\alias{SCR.area}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Effective Sampling Area for SCR analysis
%%  ~~function to do ... ~~
}
\description{
	Calculates probability of capture, effective sampling area (ESA), expected number of captures
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SCR.area(obj, SO, Mb = 0, Mbvalue = NULL, Msex = 0, Msexsigma = 0, Xsex = NULL, Mss = 0, Xd = NULL, Meff = 0, Xeff = NULL, iter = NULL, scalein = 1000, useSnowfall = FALSE, nprocs = 1, con.type = "SOCK", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  	Object returned by SCRx.fn
  	
%%     ~~Describe \code{obj} here~~
}
  \item{SO}{
  	The number of sampling opportunities used to construct the data for the model fit in obj
%%     ~~Describe \code{SO} here~~
}
  \item{Mb}{
  	Mb is an indicator variable that=1 if behavioral effects on capture histories are  to be included in the estimate of capture probability, effective sampling area and expected number of captures.  For example trap-happiness or trap shyness. 
%%     ~~Describe \code{Mb} here~~
}
  \item{Mbvalue}{
  	Mbvalue gives the value associated with the behavioral effect of capture that should be incorporated into the estimates calculated by SCR.area. For example, Mbvalue=0 would mean no behavioral effect, while Mbvalue=1 would mean the full behavioral effect of capture. 
%%     ~~Describe \code{Mbvalue} here~~
}
  \item{Msex}{
  	Msex is an indicator variable that = 1 if detection probability was estimated seperately by sex
%%     ~~Describe \code{Msex} here~~
}
  \item{Msexsigma}{
  	Msex sigma is an indicatore variable that = 1 if sigma is to be estimated seperately by sex
%%     ~~Describe \code{Msexsigma} here~~
}
  \item{Xsex}{
  	The value of sex to be used for ESA estimation. Xsex = 0 means female exposure probabilities are calculated, Xsex = 1 gives male exposure probabilities. Values between 0 and 1 will use a weighted mean of exposure probability. The default (NULL), will use the weighted mean assuming the estimated sex ratio, psi.sex.
%%     ~~Describe \code{Xsex} here~~
}
  \item{Mss}{
  	Model indicator of statespace covariate. i.e., a covarate that affects
density. Not yet fully integrated or tested.
%%     ~~Describe \code{Mss} here~~
}
  \item{Xd}{
  	Vector of the same length as the number of points in the statespace that supplies the spatial covariate indicated by Mss
%%     ~~Describe \code{Xd} here~~
}
  \item{Meff}{
  	Model indicator of trap-level "effort" covariate.
%%     ~~Describe \code{Meff} here~~
}
  \item{Xeff}{
  	A numeric value, giving the log-transformed effort to be used to estimate exposure probability, ESA, and expected number of captures. Should be of length=1, but could also give effort unique to each trap. The latter is not yet tested.
%%     ~~Describe \code{Xeff} here~~
}
  \item{iter}{
  	character or numeric vector specifying iterations for which to calculate exposure probabilities and other estimates. If a character vector should be function that can be passed to apply (e.g. mean, median, quantile) or "all" (meaning perform the calculations for all iterations of the chain). If numeric, the particular iteration(s) for which parameter estimates are used to compute exposure probability and other estimates.
%%     ~~Describe \code{iter} here~~
}
  \item{scalein}{
  	Scaling factor of statespace and trap coordinates to ESA estimate, e.g. if in meters, 1000 will transform area to kilometers^2
%%     ~~Describe \code{scalein} here~~
}

  \item{useSnowfall}{
	logical indicating whether to use the snowfall package to perform calculations. Can speed up calculations when examining many iterations
%%     ~~Describe \code{useSnowfall} here~~
}
  \item{nprocs}{
  	Number of processors to use for multicore calculations with snowfall, passed to sfInit
%%     ~~Describe \code{nprocs} here~~
}
  \item{con.type}{
  	Connection type to be used with snowfall. Defaults to "SOCK"
%%     ~~Describe \code{con.type} here~~
}
  \item{\dots}{
  	Optional arguments to be passed to apply when iter is a character vector giving a function.
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	\item{mcmchist}{
		a matrix giving the parameter values used to estimate exposure probability with the derived estimates for ESA and expected number of captures appended in the final two columns
	}
	\item{prob.cap}{
		a matrix of the exposure probabilities where each row represents a node of the statespace and each column represents an iteration at which these probabilities were captures
	}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
	Joshua Goldberg
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Load the mountain lion data
data("lions")

## Reformat the "captures" data into the standard encounter data file (EDF) format
newcaptures.lions<-cbind(
   session = rep(1,nrow(captures.lions)),
   individual=captures.lions[,2],
   occasion=captures.lions[,3],
   trapid=captures.lions[,1])

## No animals were removed (dead) so "alive" is a matrix of 1's
alive=matrix(1,nrow=length(unique(newcaptures.lions[,"individual"])),ncol=1)
## State-space covariate: just constant in this case
Xd<- rep(1,nrow(statespace.lions))

### old format
### scrobj<-list(traps=traps.lions,captures=captures.lions,statespace=statespace.lions)
scrobj<-list(
    traps=traps.lions,
    captures=newcaptures.lions,
    statespace=statespace.lions,
    alive=alive,
    Xd=Xd)

test<-SCRh.fn(scrobj,ni=220, burn=20, skip=2,nz=200,
      Msigma=1, Mb=0, Msex=0, Msexsigma=0, thinstatespace=4)

ESA <- SCR.area(obj=test, SO=1, Mb = 0, Mbvalue = NULL, Msex = 0, Msexsigma = 0, Xsex = NULL, Mss = 0, Xd = NULL, Meff = 0, Xeff = NULL, iter = "median", scalein = 1000, useSnowfall = FALSE, nprocs = 1, con.type = "SOCK", ...)

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
