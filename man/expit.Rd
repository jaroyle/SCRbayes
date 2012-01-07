\name{expit}
\alias{expit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
inverse-logit transformation 
}
\description{
computes exp(x)/(1+exp(x))
}
\usage{
expit(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 a real number 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  returns a number between 0 and 1 --
  the inverse-logit transform of x. %%  ~Describe the value returned
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

## The function is currently defined as
function (x) 
{
    exp(x)/(1 + exp(x))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
