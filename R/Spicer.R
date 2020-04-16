#' Port of Tomioka and Suzuki's SpicyMKL to R, expanded for multiclass and probability outputs.


#' @param K : N x N x M array. the (i,j,m)-element contains the (i,j)-element of the m-th kernel gram matrix.
#' @param yapp :  vector of length N with sample labels.  It should be a factor for binary/multiclass classification!
#' @param C : regularization parameter . Large values of C induce strong regularization. For L1 regularization C is a scalar: or elasticnet, C is a vector of length 2: C(1)|x| + C(2)x^2/2
#' @param opt : list of options which control the behavior of Spicer:
#' loss: type of loss function:  'logit' (logistic regression, log(1+exp(- f(x)y))) for classification,
#'  'square' (square loss, 0.5*(y - f(x))^2) for regression
#'  regname: type of regularization: 'l1' (default), 'elasticnet'
#'  outerMaxIter: maximum number of iteration of outer loop. (default 300)
#'  innerMaxIter: maximum number of iteration of inner loop. (default 500) -
#'  stopdualitygap: TRUE/FALSE. If TRUE, Spicer employs duality gap for stopping criterion of outer loop. Default TRUE.  -
#'   stopIneqViolation: TRUE/FALSE. If TRUE, Spicer employs violation of inequality for stopping criterion of outer loop. Default FALSE
#'   tolOuter:tollerance of stopping criteria of outer loop. (default 0.001)
#'   tolInner: tollerance of stopping criteria of inner loop. (default tolOuter/1000) -
#'   calpha: increment factor of gamma: gamma^(t+1)=calpha*gamma^(t).  (default 10)
#'   display: 1:display no progress messages, 2(default):display outer loop progress messages, 3:display inner loop progress messaages.
#'   @return A spicer model with the following components:
#'   comb_alpha : N x 1 coefficient vector.
#'   kern_weight : 1 x M kernel weight vector, scaled to sum to 1
#'   bias : bias term
#'   activeset : indices of kernels that are active ({m : kern_weight[m] is not zero}).
#'   sorted_kern_weight: vector of non-zero kernel weights sorted by magnitude, scaled to sum to 1.
#'   opt : list of Spicer options used in run
#'   history : contains history of primal objective, dual objective, number of active kernels, and duality gap.
#'   Citation: Suzuki,Tomioka.SpicyMKL: a fast algorithm for Multiple Kernel Learning with thousands of kernels. Mach Learn (2011) 85:77–108 (c) Vlado Uzunangelov
#' @export
spicer <- function(K, yapp, C, opt = list()) {



  if (is.null(opt$regname))
    opt$regname <- "l1"


  if (is.null(opt$optname))
    opt$optname <- "Newton"  ## choices are Newton and BFGS

  ## if C has two components(elastic net), but the second one is zero, revert to l1-regularization
  if (length(C) == 2 && C[2] == 0) {
    opt$regname = "l1"
    C = C[1]
  }

  if (is.null(opt$loss)) opt$loss <- if (length(unique(yapp)) == 2) "logit" else "square"
  ## convert labels to -1,1
  if (opt$loss %in% c("logit")) {
    ## change a factor to a numeric vector
    if (is.factor(yapp)) {
      lyapp <- levels(yapp)
      yapp <- as.numeric(yapp)
    }
    ## labels were numbers represented as character strings
    if (is.character(yapp))
      yapp <- as.numeric(yapp)

    uniqYapp <- sort(unique(yapp))
    yapp[yapp == uniqYapp[1]] <- -1
    yapp[yapp == uniqYapp[2]] <- 1
    opt$classes <- if (exists("lyapp"))
      lyapp else uniqYapp
  }


  if (is.null(opt$tolOuter))
    opt$tolOuter <- 0.001

  if (is.null(opt$tolInner))
    opt$tolInner <- opt$tolOuter/1000

  if (is.null(opt$outerMaxIter))
    opt$outerMaxIter <- 300

  if (is.null(opt$innerMaxIter))
    opt$innerMaxIter <- 500

  if (is.null(opt$calpha))
    opt$calpha <- 10

  if (is.null(opt$stopIneqViolation))
    opt$stopIneqViolation <- FALSE

  if (is.null(opt$stopDualityGap))
    opt$stopDualityGap <- TRUE

  if (is.null(opt$minIter))
    opt$minIter <- 30

  if (is.null(opt$tolMinIter))
    opt$tolMinIter <- 0.001

  if (is.null(opt$display))
    opt$display <- 2

  if (is.null(opt$includeSubW))
    opt$includeSubW <- FALSE

  opt$C <- C
  opt$nkern <- dim(K)[3]

    if (is.factor(yapp) && length(levels(yapp)) > 2) {
        res <- spicer.multiclass(K, yapp, C, opt)
    } else {
        res <- spicer.default(K, yapp, C, opt)
    }

    return(res)
}
####################################################
spicer.multiclass <- function(K, yapp, C, opt) {

    opt$loss <- "logit"
    lvls <- levels(yapp)

    combos <- combn(lvls, 2)

    res <- foreach(i = 1:ncol(combos)) %do% {
        idx <- yapp %in% combos[, i]
        spicer.default(K[idx, idx, , drop = FALSE],
                       factor(yapp[idx], levels = combos[, i]), C, opt)

    }


    cw <- rep(0, dim(K)[3])
    names(cw) <- dimnames(K)[[3]]

    for (i in 1:length(res)) {
        cw[names(res[[i]]$sorted_kern_weight)] <- cw[names(res[[i]]$sorted_kern_weight)] + res[[i]]$sorted_kern_weight
    }
    cw <- cw[cw > 0]
    cw <- sort(cw, decreasing = TRUE)/sum(cw)
    ## set attribute of res to the overall kernel weights
    res <- structure(res, sorted_kern_weight = cw)

    class(res) <- c("spicer", class(res))

    return(res)
}
#######################################################

spicer.default <- function(K, yapp, C, opt) {

    ## selection of regularization functions
    expand(get.reg.funcs(opt$regname))

    ## set up trajectory tracking and some defaults
    history <- data.frame(primalobj = double(), dualobj = double(), numActiv = integer(), dualityGap = double(), elapsedTime = double(), stringsAsFactors = FALSE)
    nowTime <- Sys.time()
    oldTime <- nowTime
    elpsdTime <- 0

    N <- length(yapp)
    M <- dim(K)[3]
    oneN <- rep(1, N)  #matrix(1,nrow=N,ncol=1)
    oneM <- rep(1, M)  #matrix(1,nrow=1,ncol=M)

    ## lagrangian for equality constraint on new variable z (rho in SpicyMKL Tomioka Suzuki JMLR 2011)
    rho <- -yapp/2

    ## kernel weights divided by gamma (alpha/gamma in SpicyMKL Tomioka Suzuki JMLR 2011)
    krnlWMod <- matrix(0, nrow = N, ncol = M)
    ## bias term intial estimate
    cb <- 0

    ## elastic net appears to be a lot closer to smoothness under the heuristic below, so no need to start with a very small gamma value
    if (opt$regname == "elasticnet" && C[2]/C[1] > 5/95) {
        cgamma <- rep(1000, M)  #matrix(1000,nrow=M,ncol=1)
        cgammab <- 1000
    } else {
        cgamma <- rep(10, M)  #matrix(10,nrow=M,ncol=1)
        cgammab <- 1
    }

    if (opt$optname == "BFGS") {
        hess <- diag(1, N)

    }
    ## measure of KKT constraint violation
    ck <- Inf
    cbeta <- 1/2

    ## input for proximal update for the kernel weights alpha/gamma + rho in SpicyMKL Tomioka Suzuki JMLR 2011 remember, N=length(yapp)=length(rho)
    krnlWPrIn <- matrix(rho, N, M) + krnlWMod
    expand(kernel_norm_cpp(K, krnlWPrIn), c("wdot", "normj"))
    ## this is the proximal of normj (not proximal of krnlWPrIn, which goes in the krnlWMod update equaiton)
    pr <- prox(normj * cgamma, C, cgamma)
    activeset <- which(pr > 0)

    ## outer loop
    for (l in 1:opt$outerMaxIter) {
        ## inner loop
        for (step in 1:opt$innerMaxIter) {
            sumrho <- sum(rho)
            yrho <- yapp * rho

            fval <- funceval(opt$loss, normj, yapp, rho, yrho, sumrho, cgamma, cgammab, cb, pr, C, reg.func)

            grad <- gradient(opt$loss, yapp, yrho, rho, sumrho, cgamma, cgammab, cb, C, wdot, normj, pr, activeset)

            switch(opt$optname, Newton = {
                hess <- hessian(opt$loss, yapp, yrho, cgamma, cgammab, C, K, normj, wdot, pr, prox.deriv, activeset)

                ## find descent direction


                ## desc <- tryCatch( { -solve(hess,grad) }, error = function(err) { ##if the matrix is ill-conditioned ## add some regularization
                ## as.vector(-qr.solve(hess+diag(1e-06,N,N),grad))

                ## })
                desc <- as.vector(-qr.solve(hess + diag(1e-07, N, N), grad))
                names(desc) <- names(grad)
                ## check descent direction is valid
                graddotd <- grad %*% desc
                if (graddotd > 0) {
                  ## revert to steepest descent
                  desc <- -grad
                }
            }, BFGS = {
                desc <- as.vector(-hess %*% grad)
                graddotd <- grad %*% desc
                if (step > 1) {
                  if (graddotd >= 0) {
                    desc <- -grad
                    hess <- diag(1, N)

                  }

                  deltarho <- rho - prevrho
                  deltagrad <- grad - prevgrad
                  deltarhograd <- (deltagrad %*% deltarho)[, 1]

                  hessmult <- diag(1, N) - tcrossprod(deltagrad, deltarho)/deltarhograd
                  hess <- hessmult %*% hess %*% t(hessmult) + tcrossprod(deltagrad, deltagrad)/deltarhograd


                }
                prevrho <- rho
                prevgrad <- grad

            })

            ## store pre-line search values for later use
            old <- list(fval = fval, grad = grad, rho = rho, normj = normj, wdot = wdot, krnlWPrIn = krnlWPrIn, activeset = activeset, yrho = yrho)

            ## update rho with a full step size
            ss <- 1  #step size
            rho <- old$rho + ss * desc

            ## special case for logit - need to confirm yrho is in the bounds specified for the logit convex conjugate (see SpicyMKL Tomioka Suzuki JMLR 2011 eqn 16)

            if (opt$loss == "logit") {
                yrho <- rho * yapp
                if (any(yrho <= -1 | yrho >= 0)) {
                  yd <- yapp * desc
                  ss_mod <- 0.99 * min(c((-1 - old$yrho[yd < 0])/yd[yd < 0], -old$yrho[yd > 0]/yd[yd > 0]))
                  ss <- min(ss, ss_mod)
                  ## re-compute rho update with new step size
                  rho <- old$rho + ss * desc
                }
            }

            ## update rest of variables this update is needed to get the initial base obj function value
            krnlWPrIn <- matrix(rho, length(rho), M) + krnlWMod
            expand(kernel_norm_cpp(K, krnlWPrIn), c("wdot", "normj"))
            pr <- prox(normj * cgamma, C, cgamma)
            activeset <- which(pr > 0)

            sumrho <- sum(rho)
            yrho <- yapp * rho
            fval <- funceval(opt$loss, normj, yapp, rho, yrho, sumrho, cgamma, cgammab, cb, pr, C, reg.func)

            ## compute step length in descent direction via Armijo's rule can be moved to a separate function, but calculation depends on too many arguments, some of
            ## which large - performance hit for slight increase in modularity
            dir <- rho - old$rho
            tmpActiveset <- union(activeset, old$activeset)
            actdif <- setdiff(tmpActiveset, activeset)

            if (length(actdif) > 0) {
                old$krnlWPrIn[, actdif] <- matrix(old$rho, length(old$rho), length(actdif)) + krnlWMod[, actdif]
                ## need to update portions of old$wdot and old$normj there might be a more efficient method to do it not even sure if these updates are necessary!! - to
                ## be confirmed
                tmp1 <- kernel_norm_cpp(K, old$krnlWPrIn, actdif)
                old$wdot[, actdif] <- tmp1$wdot
                old$normj[actdif] <- tmp1$normj
                rm(tmp1)
            }

            ## calculate normj and wdot including the directional component intermediate calcs for computing normj(alpha+rho+dir*steplen) from normj(alpha+rho) (the
            ## latter is what we had originally)
            dirNorm <- 0 * normj
            dirDotWDot <- 0 * normj
            dirNorm[tmpActiveset] <- dir %*% (wdot[, tmpActiveset] - old$wdot[, tmpActiveset])
            dirDotWDot[tmpActiveset] <- dir %*% old$wdot[, tmpActiveset]
            ## Wolfe conditions
            stepL <- 1


            while (fval > old$fval + stepL * 0.01 * (dir %*% old$grad) && (grad %*% dir < 0.9 * old$grad %*% dir || stepL == 1) && stepL > 0) {
                ## while ( fval> old$fval + stepL*0.01*(dir%*%old$grad) && stepL > 0) {

                stepL <- stepL/2
                rho <- old$rho + stepL * dir
                ## all variables having rho as a component need to be updated prior to the funceval call
                normj <- 0 * normj
                ## max is to avoid sqrt of (slightly) negative numbers
                normj[tmpActiveset] <- sqrt(pmax(0, old$normj[tmpActiveset]^2 + 2 * stepL * dirDotWDot[tmpActiveset] + (stepL^2) * dirNorm[tmpActiveset]))
                pr <- prox(normj * cgamma, C, cgamma)

                fval <- funceval(opt$loss, normj, yapp, rho, yapp * rho, sum(rho), cgamma, cgammab, cb, pr, C, reg.func)
                sumrho <- sum(rho)
                yrho <- yapp * rho
                activeset <- which(pr > 0)
                grad <- gradient(opt$loss, yapp, yrho, rho, sumrho, cgamma, cgammab, cb, C, wdot, normj, pr, activeset)
            }
            ## end of Wolfe conditions line search


            ## in case armijo rule cut step length update rho-dependent variables not used in Armijo step size calculation
            if (stepL != 1) {
                activeset <- which(pr > 0)
                krnlWPrIn <- matrix(rho, length(rho), M) + krnlWMod
                wdot[, activeset] <- kernel_norm_cpp(K, krnlWPrIn, activeset)$wdot
                sumrho <- sum(rho)
                yrho <- yapp * rho
            }

            ## if fval is unexpected, print max singular value of wdot
            if (is.complex(fval)) {
                stop(paste0("Evaluation  of the dual resulted in imaginary numbers. The largest singular value of wdot is ", round(norm(wdot, "2"), digits = 3)))
            }

            if (opt$display >= 3) {
                write(paste0("inner iter: ", step, " fval: ", fval, " step length: ", stepL), stderr())
            }

            if (sqrt((old$rho - rho) %*% (old$rho - rho))/sqrt(old$rho %*% old$rho) <= opt$tolInner) {
                break
            }

            if (step == opt$innerMaxIter) {
                warning(paste0("Inner optimization did not converge to an optimal solution for rho. Increase iterations to more than ", opt$innerMaxIter),
                  immediate. = TRUE)
                return(list(comb_alpha = numeric(), kern_weight = numeric(), sorted_kern_weight = numeric(), bias = numeric(), activeset = numeric(), history = history,
                  opt = opt))
            }
        }
        ## end of inner loop (finding optimal value of rho)



        ## primal vaiable update
        krnlWMod <- krnlWMod * 0
        if (length(activeset) > 0) {
            ## this follows from the discussion in section 4.2 of Suzuki (2011) SpicyMKL see the first pop-up note with the relevant derivations
            krnlWMod[, activeset] <- krnlWPrIn[, activeset] * matrix(((pr[activeset]/normj[activeset])/cgamma[activeset]), N, length(activeset), byrow = TRUE)
        }
        cb <- cb + cgammab * sumrho

        ## Duality Gap calculations ####################### these are the actual vectors of kernel weights now (alpha in in SpicyMKL Tomioka Suzuki JMLR 2011)
        krnlW <- -krnlWMod * matrix(cgamma, N, M, byrow = TRUE)
        krnlWNorm <- kernel_norm_cpp(K, krnlW, activeset)$normj

        ## see second paragraph of section 6.1 in SpicyMKL Tomioka Suzuki JMLR 2011
        rhoMod <- rho - oneN * (sum(rho)/N)
        rhoNorm <- kernel_norm_cpp(K, matrix(rhoMod, N, M), activeset)$normj

        if (length(activeset) > 0 && opt$regname == "l1") {
            ## see second paragraph of section 6.1 in SpicyMKL Tomioka Suzuki JMLR 2011
            rhoMod <- rhoMod * min(1, C[1]/max(rhoNorm))
            rhoNorm <- rhoNorm * min(1, C[1]/max(rhoNorm))
        }

        primalArg <- rep(-cb, N)

        for (i in activeset) {
            primalArg <- primalArg + K[, , i] %*% krnlW[, i]
        }

        primalVal <- primal.obj(opt$loss, yapp, primalArg, length(activeset), reg.func, krnlWNorm, C)

        dualVal <- dual.obj(opt$loss, yapp, rhoMod, rhoNorm, reg.dual, C)

        dualGap <- if (is.infinite(primalVal))
            Inf else abs(primalVal - dualVal)/abs(primalVal)

        nowTime <- Sys.time()
        elpsdTime <- elpsdTime + as.numeric(nowTime - oldTime)
        oldTime <- nowTime

        ## update history data object
        history <- rbind(history, setNames(list(primalVal, dualVal, length(activeset), dualGap, elpsdTime), names(history)))

        if (opt$display > 1) {
            write(paste0("outer iter: ", l, " primal: ", primalVal, " dual: ", dualVal, " duality_gap: ", dualGap), stderr())
        }

        ## Duality Gap stopping criterion#######################
        if (dualGap < opt$tolOuter && opt$stopDualityGap) {
            break
        }
        ## end Duality Gap stopping ###########################
        if (opt$stopDualityGap && nrow(history) > opt$minIter && (history[step - opt$minIter, "dualityGap"] - dualGap)/history[step - opt$minIter, "dualityGap"] <
            opt$tolMinIter) {
            warning("The duality gap has been closing very slowly indicating slow convergence.You should examine your kernels for multicollinearity and or change regularization parameters.Alternatively you can increase minIter or decrease tolMinIter.",
                immediate. = TRUE)
            return(list(comb_alpha = numeric(), kern_weight = numeric(), sorted_kern_weight = numeric(), bias = numeric(), activeset = numeric(), history = history,
                opt = opt))
        }

        ## End of Duality Gap calculations###################

        ## KKT stopping criterion####################### Also Update Augmented Lagrangian coefficients (cgamma and cgammab) ###


        krnlWPrInKKT <- krnlWPrIn
        if (length(activeset) > 0) {
            krnlWPrInKKT[, activeset] <- krnlWPrIn[, activeset] * matrix(pmin(1, C[1]/normj[activeset]), N, length(activeset), byrow = TRUE)
        }

        ## from MATLAB code hresid <- (normj*cgamma - pr).^2 i.e. hresid measures which proximal to the norm of a given kernel's weights show the biggest
        ## difference from the actual norm
        hresid <- sqrt(colSums((krnlWPrInKKT - matrix(rho, N, M))^2))
        maxgap <- 0

        ## tmp is how much cgammab is changing from one iteration to the next, it is used to determine how quickly to ramp up cgammab
        tmp <- abs(sumrho)
        ## determine how quickly to adjust cgammab (the augmented Lagrangian coefficient of cb)
        big_cb <- (tmp > cbeta * ck)

        maxgap <- max(tmp, maxgap)

        ## compute which residuals are large enough ro require faster adjustment of the augmented Lagrangian coefficient cgamma for the respective kernel norm

        largeRes <- which(hresid > cbeta * ck)
        maxgap <- max(max(abs(hresid)), maxgap)

        if (opt$display >= 2 && opt$stopIneqViolation) {
            write(paste0("outer iter: ", l, " maxgap: ", maxgap, " ck: ", ck), stderr())
        }

        if (maxgap <= ck) {
            ck <- maxgap
        }

        if (ck < opt$tolOuter && opt$stopIneqViolation) {
            break
        }

        ## end KKT stopping ###########################

        ## update cgammab according to how close we are to an optimal cgammab (determined by the size of the sumrho update relative to the size of the gaps
        ## between the kernel norm proximals and actual kernel norms - large update means we are still making big steps (i.e. we can move faster towards removing
        ## the AL term); smaller updates mean we are still moving in the vicinity of the current solution, which might be a result of a nearby non-smooth region
        ## -we should therefore maintain higher resolution in our AL elimination)
        cgammab <- ifelse(big_cb, opt$calpha * cgammab, opt$calpha + cgammab)

        # update cgammas
        cgamma[largeRes] = opt$calpha * cgamma[largeRes]
        krnlWMod[, largeRes] = krnlWMod[, largeRes]/opt$calpha

        ## smallRes kernel norms are the complement of the largeRes ones
        smallRes <- setdiff(1:M, largeRes)

        if (length(smallRes) > 0) {
            cgamma[smallRes] = opt$calpha + cgamma[smallRes]
            krnlWMod[, smallRes] = krnlWMod[, smallRes] * matrix((cgamma[smallRes] - opt$calpha)/cgamma[smallRes], N, length(smallRes), byrow = TRUE)
        }
        ## End update Augmented Lagrangian coefficients (cgamma and cgammab) ###


        ## update all main variables
        krnlWPrIn <- matrix(rho, N, M) + krnlWMod
        expand(kernel_norm_cpp(K, krnlWPrIn), c("wdot", "normj"))
        pr <- prox(normj * cgamma, C, cgamma)
        activeset <- which(pr > 0)

        ## make sure user knows if algorithm exceeded outerMaxIter
        if (step == opt$outerMaxIter) {
            warning(paste0("Outer optimization did not converge to an optimal solution for alpha and b. Increase iterations to more than ", opt$outerMaxIter),
                immediate. = TRUE)
            return(list(comb_alpha = numeric(), kern_weight = numeric(), sorted_kern_weight = numeric(), bias = numeric(), activeset = numeric(), history = history,
                opt = opt))
        }


    }
    ## end of outer loop

    tmp <- rep(0, M)
    tmp[activeset] <- pr[activeset]/normj[activeset]
    kern_weight <- tmp/(1 - (tmp/cgamma))
    sum_k_weight <- sum(kern_weight)

    if (sum_k_weight != 0) {
        kern_weight <- kern_weight/sum_k_weight
    }

    ## you rescale rho so that the scaling of the kernel weights is negated and the model can be used directly without consideration about the kernel weight
    ## scaling (which is for visualization/ease of interpretation purposes really)
    comb_alpha <- -rho * sum_k_weight
    bias <- -cb

    krnlW <- krnlWMod[, activeset] * matrix(cgamma[activeset], N, length(activeset), byrow = TRUE)
    rownames(krnlW) <- names(rho)
    colnames(krnlW) <- dimnames(K)[[3]][activeset]
    names(comb_alpha) <- names(rho)
    names(kern_weight) <- dimnames(K)[[3]]
    names(activeset) <- dimnames(K)[[3]][activeset]

    sorted_kern_weight <- sort(kern_weight[activeset], decreasing = T)
    res <- list(comb_alpha = comb_alpha, kern_weight = kern_weight, sorted_kern_weight = sorted_kern_weight, bias = bias, activeset = activeset, idx.train = dimnames(K)[[1]],
        history = history, opt = opt)
    class(res) <- c("spicer", class(res))

    if (opt$includeSubW)
        res$kern_alpha <- krnlW

    return(res)
}
## end of Spicer function


## Spicer Helpers#############################

## Evaluates the augmented dual proximal (in fact the negative of it since we ar eminimizing)
funceval <- function(loss, normj, yapp, rho, yrho, sumrho, cgamma, cgammab, cb, pr, C, reg.func) {
    val <- switch(loss,
                  logit = logit.loss(yrho),
                  square = square.loss(yapp, rho))
    ## they use Proposition 1 (eqn 23) in SpicyMKL Tomioka Suzuki JMLR 2011 to convert Moreau envelope of the convex conjugate into moreau envelope of
    ## norm(alpha+gamma*rho)^2/2 - moreau envelope of regularization function the latter can then be evaluated at prox(norm(alpha+gamma*rho)) as in Boyd
    ## Promixal algorithms monograph, 3.1 , right before eqn 3.2

    val <- val - (sum(reg.func(pr, C)) + sum(pr^2/(2 * cgamma)) - normj %*% pr)
    val <- val + cgammab * sumrho^2/2 + cb * sumrho

    return(val)
}

## Evaluates gradient of the augmented dual proximal objective function

gradient <- function(loss, yapp, yrho, rho, sumrho, cgamma, cgammab, cb, C, wdot, normj, pr, activeset) {
    val <- switch(loss,
                  logit = logit.grad(yrho, yapp),
                  square = square.grad(rho, yapp))

    for (i in activeset) {
        val <- val + wdot[, i] * (pr[i]/normj[i])
    }

    val <- val + (cgammab * sumrho + cb)

    return(val)
}

## Evaluates Hessian on augmented dual proximal objective function

hessian <- function(loss, yapp,  yrho, cgamma, cgammab, C, K, normj, wdot, pr, prox.deriv, activeset) {
    val <- switch(loss,
                  logit = logit.hess(yrho),
                  square = square.hess(length(yapp)))

    dpr <- prox.deriv(normj * cgamma, C, cgamma)
    w1 <- pr/normj
    w2 <- (cgamma * dpr * normj - pr)/(normj^3)

    for (i in activeset) {
        val <- val + w1[i] * K[, , i]
        val <- val + w2[i] * (wdot[, i] %*% t(wdot[, i]))
    }

    val <- val + cgammab

    return(val)

}


dual.obj <- function(loss, yapp, rho, rhoNorm, reg.dual, C) {

    dual <- switch(loss, logit = {
        dual <- -logit.loss(yapp * rho) - sum(reg.dual(rhoNorm, C))
    }, square = {
        dual <- -square.loss(yapp, rho) - sum(reg.dual(rhoNorm, C))
    })

    return(dual)
}


## Primal Obj Funcitons########################### z here has the opposite sign to the MATLAB code (but the right one theoretically) - the opposite sign
## is compensated for by negating z in the argument pre-processing
primal.obj <- function(loss, yapp, z, nactive, reg.func, krnlWNorm, C) {
    primal <- switch(loss,
                     logit = sum(log(1 + exp(-yapp * z))),
                     square = sum((yapp - z)^2) * 0.5)

    if (nactive > 0) {
        primal <- primal + sum(reg.func(krnlWNorm, C))
    }

    return(primal)
}

## for the logit implementation, yrho needs to stay which (0,1), so
##-yrho needs to stay within (-1,0) (all open intervals!!)
check.yrho <- function(yrho) {
    ## this is to make sure yrho falls within the required constraint
    ##-1<-yrho<0,
    ## which is the negated actual constraint 0>yrho>1
    pmin(pmax(yrho, -1 + 1e-07), -1e-07)
}


## Loss functions################################################### Dual Loss Functions############################## all functions below are actually
## related to the dual of the loss function (aka convex conjugate), i.e. *.loss is the evaluation of the dual of the respective loss, *.grad evaluate the
## gradient of the dual of the loss *.hess evaluates the Hessian of the dual of the loss the dual variable is rho, which is the lagrangian of the
## K*alpha+beta (pseudo) constraint - it is Nx1 dimensional (i.e. rho is per-sample, computed across all kernels)

## Logit Loss#######################################

## aka dual/convex conjugate of logit loss see see table 1 in Tomioka 2009 - Dual Augmented Lagrangian also Suzuki 2011 Spicy MKL eqn(16) difference is
## that negative sign of rho is accounted for in function calculation

logit.loss <- function(yrho) {
  yrho<-check.yrho(yrho)
    loss <- sum((1 + yrho) * log(1 + yrho) - yrho * log(-yrho))
    return(loss)
}

## different from table 1 in Tomioka 2009 - Dual Augmented Lagrangian again because of the substitution of rho for negative rho also this is technically
## -logit.grad, even with the variable substitution!!
logit.grad <- function(yrho, yapp) {
  yrho<-check.yrho(yrho)
    grad <- yapp * log((1 + yrho)/(-yrho))
    return(grad)
}

## remember, the logit uses binary {-1,1} y labels, so the y^2 term that should have been in the numerator disappears!!
logit.hess <- function(yrho) {
    hess <- diag(1/(-yrho * (1 + yrho)))
    return(hess)
}

## SVM Loss################################### Suzuki 2011 Spicy MKL eqn(17)
svm.loss <- function(yrho) {
    loss <- sum(yrho)
    return(loss)
}


## Square Loss################
square.loss <- function(yapp, rho) {
    loss <- 0.5 * rho %*% rho + rho %*% yapp
    return(loss)
}

square.grad <- function(rho, yapp) {
    grad <- rho + yapp
    return(grad)
}

square.hess <- function(N) {
    hess <- diag(nrow = N)
    return(hess)
}

## regularization terms############################################ the proximal operators are vectorized

## regularization function factory######
get.reg.funcs <- function(regname = c("l1", "elasticnet")) {
    regname = match.arg(regname)

    switch(regname, l1 = list(reg.func = l1.reg, reg.dual = l1.dual, prox = l1.prox, prox.deriv = l1.prox.deriv), elasticnet = list(reg.func = elas.reg,
        reg.dual = elas.dual, prox = elas.prox, prox.deriv = elas.prox.deriv))

}

## L1-norm#######################################
l1.reg <- function(x, C) {
    reg <- C * abs(x)
    return(reg)
}

## aka convex conjugate see Suzuki (2011) SpicyMKL -eqn (21)
l1.dual <- function(x, C) {
    regDual <- rep(0, length(x))
    regDual[x > C] <- Inf
    return(regDual)
}

## see table 2 in Tomioka 2009 - Dual Augmented Lagrangian can also be derived with Moreau decomposition and projection on the Linf-norm ball - see Boyd
## Proximal Algorithms monograph, section 6.5.1
l1.prox <- function(x, C, eta) {
    prox <- pmax(x - C * eta, 0)
    return(prox)
}


## derivative of l1 proximal function
l1.prox.deriv <- function(x, C, eta) {
    proxDeriv <- as.numeric(x > C * eta)
    return(proxDeriv)
}

## Elastic Net Regularization##################################
elas.reg <- function(x, C) {
    reg <- C[1] * abs(x) + (C[2]/2) * (x^2)
    return(reg)
}

## aka convex conjugate see Suzuki(2011) -SpicyMKL eqn(46)
elas.dual <- function(x, C) {
    regDual <- (0.5/C[2]) * (x^2 - 2 * C[1] * abs(x) + C[1]^2) * (x > C[1])
    return(regDual)
}

## see see Boyd Proximal Algorithms monograph, section 6.5.3
elas.prox <- function(x, C, eta) {
    prox <- pmax(x - C[1] * eta, 0)/(1 + C[2] * eta)
    return(prox)
}

elas.prox.deriv <- function(x, C, eta) {
    proxDeriv <- (1/(1 + C[2] * eta)) * (x > C[1] * eta)
    return(proxDeriv)
}


#####################################################################



