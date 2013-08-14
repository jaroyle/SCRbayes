\name{print.scrdata}
\alias{print.scrdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
summary method for objects of class "scrdata" %%  ~~function to do ... ~~
}
\description{
prints a summary %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.scrdata(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
must be of class "scrdata" %%     ~~Describe \code{x} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
  objects of class "scrdata" are produced using the scrData() function
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
Andy Royle 
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

library("SCRbayes")
data(lions)

dimnames(captures.lions)<-list(1:nrow(captures.lions),c("trapid","individual","occasion"))

lion.scrdata<- scrData(traps.lions,captures.lions,statespace.lions)

### print(lion.scrdata)
#> lion.scrdata
#Hello, world!
#
#Total encounters of each individual: 
#
# 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 
# 1  1  1  1  2  1  1  1  1  2  1  1  1  1  1  1  2  1  2  2  3  1  1  1  3  2 
#27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 
# 1  4  2  2  1  1  1  1  3  2  1  1  2  1  3  1  3  1  1  1  1  1  1  1 
#
#Summary encounter frequencies: 
#nx
# 1  2  3  4 
#34 10  5  1 
#
#Num. unique traps of capture: 
#
# 1  2  3  4 
#36 10  3  1 


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
