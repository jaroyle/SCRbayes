\name{wolverine}
\alias{wolverine}
\docType{data}
\title{
wolverine data from Audrey Magoun
}
\description{
  An R list with 5 components (see Section 5.9 of the SCR book):

  wcaps which is the encounter data file (EDF), 

  wtraps is the trap deployment file (TDF) which has x,y coordinates 
of traps and a MASK of which trap was operating (1) or not (0) during
each day of the study.

 3 discrete state-space grids are provided: grid8 grid4 and grid2 which are 8x8
4x4 and 2x2 state-space grids, respectively.


%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(wolverine)}
\format{
  The format is:
> str(wolverine)
List of 5
 $ wcaps : num [1:115, 1:4] 1 1 1 1 1 1 1 1 1 1 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:4] "year" "individual" "day" "trap"
 $ wtraps:'data.frame': 37 obs. of  167 variables:
  ..$ Easting : int [1:37] 632538 634822 638455 634649 637738 625278 631690 632631 631374 634068 ...
  ..$ Northing: int [1:37] 6316012 6316568 6309781 6320016 6313994 6318386 6325157 6316609 6331273 6328575 ...
  ..$ 1       : num [1:37] 0 1 0 0 0 0 0 0 0 0 ...
  ..$ 2       : num [1:37] 0 1 0 0 0 0 0 0 0 0 ...
  ..$ 3       : num [1:37] 0 1 0 0 0 0 0 0 0 0 ...
  ..$ 4       : num [1:37] 0 1 0 0 0 0 0 0 0 0 ...
  ..$ 5       : num [1:37] 0 1 0 0 0 0 0 0 0 0 ...
  ..$ 6       : num [1:37] 0 1 0 0 0 0 0 0 0 0 ...
  ..$ 7       : num [1:37] 0 1 0 0 0 0 0 0 0 0 ...
  ..$ 8       : num [1:37] 0 1 0 0 0 0 0 0 0 0 ...
  ..$ 9       : num [1:37] 0 1 0 0 0 0 0 0 0 0 ...
.... [deleted]......
  ..$ 96      : num [1:37] 1 1 0 0 1 1 1 0 1 1 ...
  ..$ 97      : num [1:37] 1 1 0 0 1 1 1 0 1 1 ...
  .. [list output truncated]
 $ grid8 :'data.frame': 157 obs. of  2 variables:
  ..$ UTMx: num [1:157] 601677 609674 617670 589840 597836 ...
  ..$ UTMy: num [1:157] 6425005 6425162 6425320 6416773 6416930 ...
 $ grid4 :'data.frame': 619 obs. of  2 variables:
  ..$ UTMsX: num [1:619] 607596 611594 615592 601677 605676 ...
  ..$ UTMsY: num [1:619] 6429121 6429200 6429279 6425005 6425084 ...
 $ grid2 :'data.frame': 2466 obs. of  2 variables:
  ..$ Xutm: num [1:2466] 612554 614553 607596 609595 611594 ...
  ..$ Yutm: num [1:2466] 6431219 6431258 6429121 6429161 6429200 ...
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
  Data are described in detail here:
  
  Magoun, A. J., Long, C. D., Schwartz, M. K., Pilgrim, K. L., Lowell,
  R. E., and Valkenburg, P. (2011). Integrating motion-detection cameras
  and hair snags for wolverine identification. The Journal of Wildlife
  Management, 75(3), 731-739.
  
  Royle, J. A., Magoun, A. J., Gardner, B., Valkenburg, P., and Lowell,
  R. E. (2011). Density estimation in a wolverine population using
  spatial capture-recapture models. The Journal of Wildlife Management,
  75(3), 604-611.
  %%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
  SCR book, Section 5.9

  %%  ~~ possibly secondary sources and usages ~~
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

## now fit a model using the current development version of
##   the fitting function. With the 2 km grid the area per pixel is 4 km^2

test<-SCRh.fn( wolv.scrdata,ni= 1200, burn=200, skip=2,nz=100,theta=1,
               Msigma=1, Mb=0, Msex=1, Msexsigma=1, 
               coord.scale=5000, area.per.pixel=4, thinstatespace=4)

## print a summary 
print(test, digits=3)

## maybe str(wolverine) ; plot(wolverine) ...
}
\keyword{datasets}
