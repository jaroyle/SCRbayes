make.statespace <-
function (ll = NA, buffer=.1, minx = NA, maxx = NA, miny = NA, maxy = NA,
    nx = 20, ny = NULL)
{
    ## need to edit this function to accept "average home range size"
    ## argument and that should be converted to a number of
    ## state-space points

    if (!is.na(ll[1,1])) {
        minx <- min(ll[, 1])
        maxx <- max(ll[, 1])
        miny <- min(ll[, 2])
        maxy <- max(ll[, 2])
        bx <- (maxx - minx) * buffer
        by <- (maxy - miny) * buffer
        minx <- minx - bx
        maxx <- maxx + bx
        miny <- miny - by
        maxy <- maxy + by
        r.x<- maxx-minx
        r.y<- maxy-miny
    }

    if (is.null(ny)){
      ny<- round(nx*(r.y/r.x), 0)        # approximately equal area cells
    }

    x <- sort(rep(seq(minx, maxx, , nx), ny))
    y <- rep(seq(maxy, miny, , ny), nx)
    ss<- cbind(x, y)

attr(ss,"area")<- (r.x*r.y)/nrow(ss)
attr(ss,"traps")<- ll

    return(ss)
}
