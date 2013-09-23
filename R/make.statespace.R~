make.statespace <-
function (ll = NA, minx = NA, maxx = NA, miny = NA, maxy = NA,
    nx = 40, ny = NULL, buffer = 0)
{
    if (is.null(ny))
        ny <- nx
    if (!is.na(ll)) {
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
    }
    x <- sort(rep(seq(minx, maxx, , nx), ny))
    y <- rep(seq(maxy, miny, , ny), nx)
    cbind(x, y)
}
