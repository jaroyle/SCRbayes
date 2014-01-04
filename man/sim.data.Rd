\name{sim.data}
\alias{sim.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate some SCR data %%  ~~function to do ... ~~
}
\description{
This function simulates SCR data including auxiliary telemetry data if
desired. 
}
\usage{
sim.data(N = 200, sigma = 0.3, loglam0 = log(0.35), K = 12,
statespace,
traplocs,
Xss = NULL, alpha1 = 0, coord.scale = 5000, Ntel = 2, nfixes = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{loglam0}{
%%     ~~Describe \code{loglam0} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{statespace}{
%%     ~~Describe \code{statespace} here~~
}
  \item{traplocs}{
%%     ~~Describe \code{traplocs} here~~
}
  \item{Xss}{
%%     ~~Describe \code{Xss} here~~
}
  \item{alpha1}{
%%     ~~Describe \code{alpha1} here~~
}
  \item{coord.scale}{
%%     ~~Describe \code{coord.scale} here~~
}
  \item{Ntel}{
%%     ~~Describe \code{Ntel} here~~
}
  \item{nfixes}{
%%     ~~Describe \code{nfixes} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
