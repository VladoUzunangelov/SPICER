This is an R port of Tomioka and Suzuki's SpicyMKL.  
SpicyMKL solves  both binary classification (hinge & logitstic losses) and regression (squared error loss) problems. SPICER solves binary classification (logistic only), regression (squared error loss) and multiclass problems (logistic loss). The latter is done by expanding the binary SpicyMKL solver via one-versus-one
class pairs training approach. In addition, SPICER outputs class or probability predictions for binary/multiclass classification tasks.  

**Installation**\
Requires devtools package. From an R interpreter do: \
devtools::install_github("VladoUzunangelov/SPICER")

**Training**\
*Inputs*\
*K* : N x N x M array. the (i,j,m)-element contains the (i,j)-element of the m-th kernel gram matrix.\
*yapp* :  vector of length N with sample labels.  It should be a factor for binary/multiclass classification. \
*C* : regularization parameter . Large values of C induce strong regularization. For L1 regularization C is a scalar: for elasticnet, C is a vector of length 2: C(1)|x| + C(2)x^2/2\
*opt* : list of options which control the behavior of SPICER:
* loss: type of loss function:  'logit' (logistic regression, log(1+exp(- f(x)y))) for classification,'square' (square loss, 0.5*(y - f(x))^2) for regression
* regname: type of regularization: 'l1' (default), 'elasticnet'
* optname: optimization solver for dual variable (rho) inner minimization: 'Newton' (default), 'BFGS'
* stop_duality_gap: TRUE/FALSE. If TRUE, SPICER employs duality gap for stopping criterion of outer loop. Default TRUE.  -
* stop_ineq_violation: TRUE/FALSE. If TRUE, SPICER employs violation of inequality for stopping criterion of outer loop. Default FALSE
* outer_maxiter: maximum number of iterations of outer loop. (default 300)
* inner_maxiter: maximum number of iterations of inner loop. (default 500)
* tol_outer:tolerance of stopping criteria of outer loop. (default 0.001)
* tol_inner: tolerance of stopping criteria of inner loop. (default tolOuter/1000)
* miniter: minimum number of iterations of outer loop. (default 30)
* tol_miniter: minimum relative improvement required for each outer iteration after miniter. (default 0.001)
* calpha: increment factor of gamma: gamma^(t+1)=calpha*gamma^(t).  (default 10)
* incl_subw: TRUE/FALSE. If TRUE, Spicer includes the full NxM matrix of alpha coefficients. Default FALSE.
* display: 1:display no progress messages, 2(default):display outer loop progress messages, 3:display inner loop progress messaages.

*Output*\
A SPICER model with the following components:
* *comb_alpha* : N x 1 coefficient vector.
* *kern_weight* : 1 x M kernel weight vector, scaled to sum to 1
* *bias* : bias term
* *activeset* : indices of kernels that are active ({m : kern_weight[m] is not zero}).
* *sorted_kern_weight*: vector of non-zero kernel weights sorted by magnitude, scaled to sum to 1.
* *opt* : list of SPICER options used in run
* *history* : contains history of primal objective, dual objective, number of active kernels, and duality gap.
* *kern_alpha* : If incl_subw is TRUE, the NxM matrix of alpha coefficients.

**Prediction**\
*Inputs*
* *model* SPICER model
* *kTest* a list of test kernels of dimensions nTrain, nTest, length(model\$sorted_kern_weight) (you can pass extra kernels so long as ones referenced in model\$sorted_kern_weight are included). In the case of multiclass prediction, the third dimension of kTest should include all possible kernels for each pairwise classification task.\
* *type* - applicable for classification (binary & multiclass) only - "response" returns the predicted class labels, while "probability" returns the class probability (for classification, positive class is the second class label in model\$opt\$classes)\

*Output* \
A prediction vector of length nTest computed by\
f(x) = \sum{model$sorted_kern_weigth[i]*KTest[,,i]}*model$comb_alpha + model$beta\
If prediction task is 'regression', output is continuous valued predictions\
If prediction task is 'classification' and type is 'response' - output is predicted labels.\
If prediction task is 'classification' and type is 'probability' - output is probabilities of classes.

**Citations** 
1.  V. Uzunangelov, C. K. Wong, and J. Stuart. *Highly Accurate Cancer Phenotype Prediction with AKLIMATE, a Stacked Kernel Learner Integrating Multimodal Genomic Data and Pathway Knowledge.* bioRxiv, July 2020.
2.  T. Suzuki and R. Tomioka.  *SpicyMKL: a fast algorithm for Multiple Kernel Learningwith thousands of kernels.* Machine Learning, 85(1-2):77–108, Oct. 2011. ISSN 0885-6125,1573-0565. 



(c) Vladislav Uzunangelov 2016  
uzunangelov@gmail.com
