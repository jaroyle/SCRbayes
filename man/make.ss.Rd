\name{make.ss}
\alias{make.ss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
c%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make.ss(traps, buffer, grid.density)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{traps}{
%%     ~~Describe \code{traps} here~~
}
  \item{buffer}{
%%     ~~Describe \code{buffer} here~~
}
  \item{grid.density}{
%%     ~~Describe \code{grid.density} here~~
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

## The function is currently defined as
function (traps, buffer, grid.density) 
{
    y <- seq(min(traps$x) - buffer, max(traps$x) + buffer, by = grid.density)
    x <- seq(min(traps$y) - buffer, max(traps$y) + buffer, by = grid.density)
    ss.y <- rep(y, length(x))
    ss.x <- rep(x[1], length(y))
    for (i in 2:length(x)) {
        ss.x <- c(ss.x, rep(x[i], length(y)))
    }
    ss.animal <- data.frame(X = ss.x, Y_coord = ss.y)
    list(ss.animal = ss.animal)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
