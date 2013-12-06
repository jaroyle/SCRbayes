\name{scrData}
\alias{scrData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  makes an object of class "scrdata" which can be analyzed by the 
  functions of the SCRbayes package
    %%  ~~function to do ... ~~
}
\description{
packages up capture data, trap information and state-space information
to make a data object of class "scrobj".

Right now it doesn't do anything special except package the input as a list
and slightly reformat some things. In the future it will check for
consistency of the different data files, run some error checks, and 
produce some additional output.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scrData(traps, captures, statespace, alive = NULL, Xd = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{traps}{
ntraps x (3+noccasions) matrix: 
column 1 = trap ID, 
column 2 = x coordinate of trap, 
column 3 = y coordinate of trap, 
columns 4:(3+noccasions) = trap operation mask, 1 if operated, 0 if not
%%     ~~Describe \code{traps} here~~
}
  \item{captures}{
Matrix having encounter information about each individual. Some flexibility 
is allowed in data format (more flexibility in the future).
Minimally, the matrix should have columns labeled 
"individual" "occasion" and "trapid".

"individual" = integer from 1:nind

"occasion" = integer from 1:(noccasions)

"trapid" = integer from 1:ntraps mapping to rows of "traps"

Optionally may have a column "session" which is the year or primary
period of the sample. Currently no multi-session functionality is available
but the main functions of SCRbayes accepts this data structure. 
%%     ~~Describe \code{captures} here~~
}
  \item{statespace}{
This is an "ngrid x 3" matrix describing the state-space of the point process
model describing the distribution of individuals.  The column fields are:
column 1 = x-coordinate of potential activity center,
column 2 = y-coordinate of potential activity center,
column 3 = binary indicator of suitable habitat (1) or not (0)
%%     ~~Describe \code{statespace} here~~
}
  \item{alive}{
A matrix describing whether each individual was available for 
sampling (1) or not (0).  Models in SCRbayes assume you know 
whether individual were removed or not, without error. 
The matrix should be "nind x noccasions". 
%%     ~~Describe \code{alive} here~~
}
  \item{Xd}{
A "density covariate" -- should be a vector of "ngrid x 1" where ngrid = nrow(statespace). 
%%     ~~Describe \code{Xd} here~~
}
\item{Ytel}{
 Telemetry data on a sample of individuals. 
 This must be an "ngrix x Ntel" matrix of count frequencies where
 Ytel[j,i] is the frequency of telemetry fixes in pixel j for individual
 i. The number of ROWS of Ytel must be the same dimension as the
 state-space grid and therefore the telemetry data MUST BE
 pre-binned. Obviously covariates for the RSF model have to correspond
 to the state-space grid also. Right now these formatting requirements
  are MANDATORY but only for expediency on my end
 of things. 
  }
}
\details{
more to come .... %%  ~~ If necessary, more details than the description above ~~
}
\value{
an object of class "scrdata"
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
For the use of telemetry data see the paper:


Royle, J.A., R.B. Chandler, C.C. Sun and A.K. Fuller. 
(2013). Integrating resource selection information with spatial
capture-recapture. Methods in Ecology and Evolution, 3:545-554.
}
\author{
Andy Royle %%  ~~who you are~~
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

library("SCRbayes")
data(lions)

dimnames(captures.lions)<-
   list(1:nrow(captures.lions),c("trapid","individual","occasion"))

lion.scrdata <- scrData(traps.lions,captures.lions,statespace.lions)

test <- SCRh.fn(lion.scrdata,ni=220, burn=20,
 skip=2,nz=200,Msigma=1, Mb=0, Msex=0, Msexsigma=0, 
 Xeff=NULL, Xsex=NULL, thinstatespace=2)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
