## (c) Vlado Uzunangelov 2016 uzunangelov@soe.ucsc.edu

## Evaluates the kernel norm of a given vector u for each kernel Also returns dot of each kernel with vector u (cached because frequently used in other
## computations) K is an NxNxM array where N - no. samples, M- no.kernels u - NxM matrix of M Nx1 vectors, whose kernel norm is evaluated for the kernel
## of matching index activeset defaults to all kernel indices is actiset is a vector of length zero, returns results with one dimension being zero
kernel_norm <- function(K, u, activeset = NULL) {
    if (is.null(activeset)) {
        activeset <- 1:dim(u)[2]
    }


    N <- dim(u)[1]
    M <- length(activeset)
    wdot <- matrix(0, N, M)
    normj <- rep(0, M)

    if (M > 0) {
        for (i in 1:length(activeset)) {
            ind <- activeset[i]
            wdot[, i] <- K[, , ind] %*% u[, ind]
            ## the absolute value here is for small numerical errors (-1* 10e-15 and such)
            normj[i] <- sqrt(max(wdot[, i] %*% u[, ind], 0))
        }
    }
    return(list(wdot = wdot, normj = normj))
}
