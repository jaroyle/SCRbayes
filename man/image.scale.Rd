\name{image.scale}
\alias{image.scale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
c%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
image.scale(z, col, x, y = NULL, size = NULL, digits = 2, labels = c("breaks", "ranges"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
%%     ~~Describe \code{z} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{size}{
%%     ~~Describe \code{size} here~~
}
  \item{digits}{
%%     ~~Describe \code{digits} here~~
}
  \item{labels}{
%%     ~~Describe \code{labels} here~~
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
function (z, col, x, y = NULL, size = NULL, digits = 2, labels = c("breaks", 
    "ranges")) 
{
    n <- length(col)
    usr <- par("usr")
    mx <- mean(usr[1:2])
    my <- mean(usr[3:4])
    dx <- diff(usr[1:2])
    dy <- diff(usr[3:4])
    if (missing(x)) 
        x <- mx + 1.05 * dx/2
    else if (is.list(x)) {
        if (length(x$x) == 2) 
            size <- c(diff(x$x), -diff(x$y)/n)
        y <- x$y[1]
        x <- x$x[1]
    }
    else x <- x[1]
    if (is.null(size)) 
        if (is.null(y)) {
            size <- 0.618 * dy/n
            y <- my + 0.618 * dy/2
        }
        else size <- (y - my) * 2/n
    if (length(size) == 1) 
        size <- rep(size, 2)
    if (is.null(y)) 
        y <- my + n * size[2]/2
    i <- seq(along = col)
    rect(x, y - i * size[2], x + size[1], y - (i - 1) * size[2], 
        col = rev(col), xpd = TRUE)
    rng <- range(z, na.rm = TRUE)
    bks <- seq(from = rng[2], to = rng[1], length = n + 1)
    bks <- formatC(bks, format = "f", digits = digits)
    labels <- match.arg(labels)
    if (labels == "breaks") 
        ypts <- y - c(0, i) * size[2]
    else {
        bks <- paste(bks[-1], bks[-(n + 1)], sep = " - ")
        ypts <- y - (i - 0.5) * size[2]
    }
    text(x = x + 1.2 * size[1], y = ypts, labels = bks, adj = ifelse(size[1] > 
        0, 0, 1), xpd = TRUE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
