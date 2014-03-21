\name{SCR.bf}
\alias{SCR.bf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Compute posterior marginal model probabilities and Bayes factors
%%  ~~function to do ... ~~
}
\description{
	Uses the Gelfand and Dey (1994) approximation to compute the marginal likelihood of a model within a candidate set. These probabilities are used to calculate Bayes factors and assess the (qualitative) relative support of models using the criteria given by Jeffreys (1961) and Kass and Raftery (1995). This function also computes summary statistics for the number of observations that were estimated to have a likelihood less than computational precision.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SCR.bf(modelsin, refmodel = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelsin}{
  	The output from a set of SCRx.fn models passed to the function in a named list. The names are used within the output for display purposes.
%%     ~~Describe \code{modelsin} here~~
}
  \item{refmodel}{
  	An optional character vector of length=1 giving the name of model or a numeric giving the list index for the model to be used as a reference model in Bayes factor calculation. If refmodel is omitted, only the marginal posterior probabilities are returned.
%%     ~~Describe \code{refmodel} here~~
}
}
\details{
	Some assumptions made in these calculations: (1) Posterior probability density of the model parameters are assumed to be independent; (2) empirical densities are associated with left side of the intervals used for approximation; (3) samples for approximation are drawn from a multivariate t distribution with df = number of MCMC samples, centered at the mean of the posterior samples and a correlation matrix computed from the covariance of the parameter estimates.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	\item{ModelSelection}{
		A data.frame with the marginal model likelihoods for each model, along with Bayes factor results if refmodel is supplied. JeffreysEvidence refers to the qualitative criteria defined by Jeffreys (1961) and KRevidence refers to the support criteria given by Kass and Raftery (1995).
	}
	\item{LowLikelihoodObservations}{
		A data.frame giving the total (sum), mean, median, minimum and maximum number of observations where the likelihood could not be returned within computational precision during MCMC sampling.
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
	Joshua Goldberg, joshua.goldberg@umontana.edu
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
data(wolverine)

## make the encounter data file (EDF)
EDF<-wolverine$wcaps
dimnames(EDF)<-list(NULL,c("session","individual","occasion","trapid"))
## grab the trap deployment file (TDF)
TDF<-wolverine$wtraps

## we will fit a model with sex-specific parameters, so grab the sex variable
wsex<-wolverine$wsex

## bundle these into an "scrdata" object along with the state-space grid
grid<-cbind(wolverine$grid2,rep(1,nrow(wolverine$grid2)))
wolv.scrdata<-scrData(TDF,EDF,grid,Xsex=wsex)


## now fit models using the current development version of
##   the fitting function
# Null SCR model - no covariates
test1<-SCRh.fn( wolv.scrdata,ni= 1200, burn=200, skip=2,nz=100,theta=1,
Msigma=1, Mb=0, Msex=0, Msexsigma=0, 
coord.scale=5000, area.per.pixel=4, thinstatespace=4)

# With sex covariate
test2<-SCRh.fn( wolv.scrdata,ni= 1200, burn=200, skip=2,nz=100,theta=1,
Msigma=1, Mb=0, Msex=1, Msexsigma=0, 
coord.scale=5000, area.per.pixel=4, thinstatespace=4)

# With sexsigma
test3<-SCRh.fn( wolv.scrdata,ni= 1200, burn=200, skip=2,nz=100,theta=1,
Msigma=1, Mb=0, Msex=0, Msexsigma=0, 
coord.scale=5000, area.per.pixel=4, thinstatespace=4)

# With sex and sexsigma
test4<-SCRh.fn( wolv.scrdata,ni= 1200, burn=200, skip=2,nz=100,theta=1,
Msigma=1, Mb=0, Msex=0, Msexsigma=0, 
coord.scale=5000, area.per.pixel=4, thinstatespace=4)

modelsin <- list("Distance"=test1, "Sex"=test2, "SexSigma"=test3,"Sex+SexSigma"=test4)

wolv.bf <- SCR.bf(modelsin, refmodel="Distance")
wolv.bf


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
