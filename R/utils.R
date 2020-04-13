## (c) Vlado Uzunangelov 2016 uzunangelov@soe.ucsc.edu

## expand a list to its component objects - the list entities now become individual entries in the (current) R environment if the names of the output
## variables are not supplied, use names of expanded list/vector variables values is a list or a vector of variables we want to expand
expand <- function(vals, nms = NULL) {
    if (is.null(nms)) {
        nms = names(vals)
    }
    mapply(assign, nms, vals, MoreArgs = list(envir = parent.frame()))
    ## just so it does not return the mapply output
    invisible()
}

## taken from Hadley Wickham's scales package rescale continuous vectors to a specified interval
rescale <- function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {
    if (zero_range(from) || zero_range(to)) 
        return(rep(mean(to), length(x)))
    (x - from[1])/diff(from) * diff(to) + to[1]
}

