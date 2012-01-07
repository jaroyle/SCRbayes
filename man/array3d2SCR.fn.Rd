\name{array3d2SCR.fn}
\alias{array3d2SCR.fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
c%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
array3d2SCR.fn(y = tiger.data3d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
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
function (y = tiger.data3d) 
{
    dd <- dim(y)
    nind <- dd[1]
    nrep <- dd[2]
    ntrap <- dd[3]
    indid <- repid <- trapid <- NULL
    for (i in 1:nind) {
        yi <- y[i, , ]
        indid <- c(indid, rep(i, sum(yi > 0, na.rm = TRUE)))
        repid <- c(repid, row(yi)[!is.na(yi) & yi > 0])
        trapid <- c(trapid, col(yi)[!is.na(yi) & yi > 0])
    }
    o <- cbind(trapid, indid, repid)
    dimnames(o) <- list(NULL, c("trapid", "individual", "sample"))
    o <- o[order(o[, 1]), ]
    o
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
