## (c) Vlado Uzunangelov 2016
## uzunangelov@soe.ucsc.edu

##Computes predictions from Spicer model.Minimum required components:
##object of class spicer, with
##sorted_kern_weight - active weights (equivalent to Spicer model$sorted_kern_weight
##bias - offset term (equivalent to Spicer model$bias)
##comb_alpha - weight vector for training samples (equivalent to Spicer model$comb_alpha)

##kTest is of dimensions nTrain, nTest, length(model$sorted_kern_weight) (can pass extra kernels so long as ones referenced in model$sorted_kern_weight are included)
##In the case of multiclass prediction, the third dimensionof kTest should include all possible kernels for each pairwise classification task.
##
##type - applicable for classification (binary & multiclass) only - response returns the predicted class labels, while probability returns the probability of the positive class (the second class labels in model$opt$classes
##output is a prediction vector of length nTest
##if prediction task is "regression", output is continuous prediciton values
##if prediction task is classification and type is response - predictions are the labels of the two classes
##if prediction task is classification and type is probability - predictions are probabilities of positive class (model$opt$classes[2])

predict.spicer <- function(model,kTest,type="probability"){

    if(!is.null(model$opt)){
        res <- predict.spicer.generic(model,kTest[,,
                                                  names(model$sorted_kern_weight),
                                                  drop=FALSE],type)
    } else {
        ##multiclass option
        res <- predict.spicer.multiclass(model,kTest,type)
    }

    return(res)
}


predict.spicer.generic <- function(model,kTest,type="probability") {
    library(plyr)
    
    if(length(model$sorted_kern_weight)!=dim(kTest)[3]){
        stop("The number of weights in your model is different from the number of test kernels!")
    }
    
    if(length(model$sorted_kern_weight)>0){
        combK <- Reduce('+',mapply('*',alply(kTest,3),model$sorted_kern_weight,SIMPLIFY=FALSE))
        res <- model$comb_alpha%*%combK + model$bias
        ##convert to vector
        res <- res[1,]
        names(res) <- dimnames(kTest)[[2]]
        if(model$opt$loss=="logit"){
            switch(type,
                   response= {
                       .res <- rep(model$opt$classes[2],length(res))
                       .res[res<0] <- model$opt$classes[1]
                       res <- .res
                       names(res) <- dimnames(kTest)[[2]]
                   },
                   probability = {
                       res <- 1/(1+exp(-res))
                       res <- cbind(1-res,res)
                       rownames(res) <- dimnames(kTest)[[2]]
                       colnames(res) <- model$opt$classes
                   },
                   error("This type of output is not yet implemented!"))
            
        }
    } else {
        res <- rep(NA,dim(kTest)[2])
    }


    return(res)
}


predict.spicer.multiclass <- function(model,kTest,type="probability"){

    classes <- unique(unlist(lapply(model,function(x) x$opt$classes)))
    comb <- matrix(0,dim(kTest)[2],length(classes))
    rownames(comb) <- dimnames(kTest)[[2]]
    colnames(comb) <- classes
    for(i in 1:length(model)){
        out <- predict.spicer.generic(model[[i]],
                       kTest[model[[i]]$idx.train,,
                             names(model[[i]]$sorted_kern_weight),drop=FALSE],
                       type="probability")
        comb[,colnames(out)] <- comb[,colnames(out)]+out
    }
    comb <- comb/rowSums(comb)
    
    switch(type,
           probability={
               res <- comb
           },
           response={
               ## scores <- foreach(i=iter(model),.combine=cbind)%docomb%{
               ##     predict.spicer.generic(i,
               ##             kTest[i$idx.train,,
               ##                   names(i$sorted_kern_weight),drop=FALSE],
               ##             type=type)
               ## }
    
               ## res <- apply(scores,1,function(x)  names(which.max(table(x))))
               res <- apply(comb,1,function(x)  names(which.max(x)))
           })

    return(res)
    }

