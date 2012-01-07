array3d2SCR.fn <-
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
