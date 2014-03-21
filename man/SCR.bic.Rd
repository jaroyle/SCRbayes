\name{SCR.bic}
\alias{SCR.bic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Compute Bayesian Information Criterion (BIC)
%%  ~~function to do ... ~~
}
\description{
	This function produces a Schwarz (1978) Bayesian information criterion (BIC) table for models with SCRx.fn by the SCRbayes package. The BIC table can be sorted by BIC and includes max(ln(L)), delta BIC and BIC weights for each model in the set. This function also computes summary statistics for the number of observations that were estimated to have a likelihood less than computational precision.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SCR.bic(modelsin, ncaps, sort = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelsin}{
  	The output from a set of SCRx.fn models passed to the function in a named list. The names are used within the output for display purposes.
%%     ~~Describe \code{modelsin} here~~
}
  \item{ncaps}{
  	The number of captures in the dataset(s) used to fit the models. Can be a single number if all models used the same dataset, or a vector with the number of captures for each model.
%%     ~~Describe \code{ncaps} here~~
}
  \item{sort}{
  	logical. Should the BIC table be sorted by BIC values. If false, table is ordered by the order of models in modelsin.
%%     ~~Describe \code{sort} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	\item{ModelSelection}{
		A data.frame with the max(ln(L)), BIC, delta BIC, BIC weight and number of points with likelihood smaller than system precision in the iteration with max(ln(L)).
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

wolv.bic <- SCR.bic(modelsin, nrow(EDF))
wolv.bic
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
