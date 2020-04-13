#' Computes predictions from Spicer model.

#' @param kTest is of dimensions nTrain, nTest, length(model$sorted_kern_weight) (can pass extra kernels so long as ones referenced in model$sorted_kern_weight are included). In the case of multiclass prediction, the third dimension of kTest should include all possible kernels for each pairwise classification task.
#' type - applicable for classification (binary & multiclass) only - "response" returns the predicted class labels, while "probability" returns the probability of the positive class (the second class labels in model$opt$classes)
#' @return Output is a prediction vector of length nTest computed by
#' f(x) = \sum{model$sorted_kern_weigth[i]*KTest[,,i]}%*%model$comb_alpha + model$beta
#' if prediction task is 'regression', output is continuous values
#' if prediction task is 'classification' and type is 'response' - output are predicted labels
#' if prediction task is 'classification' and type is 'probability' - output is probabilities of positive class (model$opt$classes[2])
#' @export
predict.spicer <- function(model, kTest, type = "probability") {

    if (!is.null(model$opt)) {
        res <- predict.spicer.default(model, kTest[, , names(model$sorted_kern_weight), drop = FALSE], type)
    } else {
        ## multiclass option
        res <- predict.spicer.multiclass(model, kTest, type)
    }

    return(res)
}


predict.spicer.default <- function(model, kTest, type = "probability") {

    if (length(model$sorted_kern_weight) != dim(kTest)[3]) {
        stop("The number of weights in your model is different from the number of test kernels!")
    }

    if (length(model$sorted_kern_weight) > 0) {
        combK <- Reduce("+", mapply("*", plyr::alply(kTest, 3), model$sorted_kern_weight, SIMPLIFY = FALSE))
        res <- model$comb_alpha %*% combK + model$bias
        ## convert to vector
        res <- res[1, ]
        names(res) <- dimnames(kTest)[[2]]
        if (model$opt$loss == "logit") {
            switch(type, response = {
                .res <- rep(model$opt$classes[2], length(res))
                .res[res < 0] <- model$opt$classes[1]
                res <- .res
                names(res) <- dimnames(kTest)[[2]]
            }, probability = {
                res <- 1/(1 + exp(-res))
                res <- cbind(1 - res, res)
                rownames(res) <- dimnames(kTest)[[2]]
                colnames(res) <- model$opt$classes
            }, error("This type of output is not yet implemented!"))

        }
    } else {
        res <- rep(NA, dim(kTest)[2])
    }


    return(res)
}


predict.spicer.multiclass <- function(model, kTest, type = "probability") {

    classes <- unique(unlist(lapply(model, function(x) x$opt$classes)))
    comb <- matrix(0, dim(kTest)[2], length(classes))
    rownames(comb) <- dimnames(kTest)[[2]]
    colnames(comb) <- classes
    for (i in 1:length(model)) {
        out <- predict.spicer.default(model[[i]], kTest[model[[i]]$idx.train, , names(model[[i]]$sorted_kern_weight), drop = FALSE], type = "probability")
        comb[, colnames(out)] <- comb[, colnames(out)] + out
    }
    comb <- comb/rowSums(comb)

    switch(type, probability = {
        res <- comb
    }, response = {
        ## scores <- foreach(i=iter(model),.combine=cbind)%docomb%{ predict.spicer.default(i, kTest[i$idx.train,, names(i$sorted_kern_weight),drop=FALSE],
        ## type=type) }

        ## res <- apply(scores,1,function(x) names(which.max(table(x))))
        res <- apply(comb, 1, function(x) names(which.max(x)))
    })

    return(res)
}

