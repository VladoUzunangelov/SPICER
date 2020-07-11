## (c) Vlado Uzunangelov 2016 uzunangelov@soe.ucsc.edu

##This is for Roxygen to compile the Rcpp code correctly, needs to be in an .R file
##https://stackoverflow.com/questions/41051003/rcpp-c-functions-dont-work-in-r-package

#' @useDynLib SPICER, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' expand a list to its component objects - the list entities now become individual entries in the R environment from which the function is called.
#'@title expand
#' @param vals a named list or vector of values/objects we want to expand in individual objects in the current environment
#' @param nms optional list of names for the new object in the current environment. If the names of the output variables are not supplied, names of vals argument are used instead
#' @export
expand <- function(vals, nms = NULL) {
    if (is.null(nms)) {
        nms = names(vals)
    }
    mapply(assign, nms, vals, MoreArgs = list(envir = parent.frame()))
    ## just so it does not return the mapply output
    invisible()
}


