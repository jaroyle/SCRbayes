\name{SCR.gof}
\alias{SCR.gof}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
x%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SCR.gof(out, nx = 20, ny = 20, Xl = NULL, Xu = NULL, Yl = NULL, Yu = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{out}{
%%     ~~Describe \code{out} here~~
}
  \item{nx}{
%%     ~~Describe \code{nx} here~~
}
  \item{ny}{
%%     ~~Describe \code{ny} here~~
}
  \item{Xl}{
%%     ~~Describe \code{Xl} here~~
}
  \item{Xu}{
%%     ~~Describe \code{Xu} here~~
}
  \item{Yl}{
%%     ~~Describe \code{Yl} here~~
}
  \item{Yu}{
%%     ~~Describe \code{Yu} here~~
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
function (out, nx = 20, ny = 20, Xl = NULL, Xu = NULL, Yl = NULL, 
    Yu = NULL) 
{
    S <- out$Sout
    G <- out$G
    Sxout <- Syout <- matrix(NA, nrow = nrow(S), ncol = ncol(S))
    for (i in 1:nrow(S)) {
        Sxout[i, ] <- G[, 1][S[i, ]]
        Syout[i, ] <- G[, 2][S[i, ]]
    }
    z <- out$zout
    niter <- nrow(z)
    if (is.null(Xl)) {
        Xl <- min(Sxout) * 0.999
        Xu <- max(Sxout) * 1.001
        Yl <- min(Syout) * 0.999
        Yu <- max(Syout) * 1.001
    }
    xg <- seq(Xl, Xu, , nx)
    yg <- seq(Yl, Yu, , ny)
    Sxout2 <- cut(Sxout[z == 1], breaks = xg)
    Syout2 <- cut(Syout[z == 1], breaks = yg)
    Dn <- table(Sxout2, Syout2)/niter
    image(xg, yg, Dn, col = terrain.colors(10))
    image.scale(Dn, col = terrain.colors(10))
    stat <- statsim <- rep(NA, niter)
    for (i in 1:niter) {
        Dn <- table(cut(Sxout[i, ][z[i, ] == 1], breaks = xg), 
            cut(Syout[i, ][z[i, ] == 1], breaks = yg))
        Dnv <- Dn[1:length(Dn)]
        stat[i] <- (length(Dnv) - 1) * (var(Dnv)/mean(Dnv))
        Sxsim <- sample(G[, 1], sum(z[i, ]), replace = TRUE)
        Sysim <- sample(G[, 2], sum(z[i, ]), replace = TRUE)
        Dnsim <- table(cut(Sxsim, breaks = xg), cut(Sysim, breaks = yg))
        Dnsimv <- Dnsim[1:length(Dnsim)]
        statsim[i] <- (length(Dnsimv) - 1) * (var(Dnsimv)/mean(Dnsimv))
    }
    out <- cbind(data = stat, newdata = statsim)
    cat("P-value: ", mean(out[, 1] > out[, 2]), fill = TRUE)
    invisible(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
