# Spicer

This is a port to R of Tomioka and Suzuki's SpicyMKL.  
SpicyMKL solves  both binary classification (hinge & logitstic losses) and regression (squared error loss) problems.  
Spicer solves binary classification (logistic only), regression (squared error loss)
and multiclass problems (logistic loss).  
The latter is done by expanding the binary SpicyMKL solver via one-versus-one
class pairs training approach.  
In addition, Spicer outputs class or probability predictions for binary/multiclass classification tasks.  

Input:   
K        : N x N x M array. the (i,j,m)-element contains the (i,j)-element of the m-th kernel gram matrix.   
yapp     : signal vector of length N (NOT a Nx1 or 1xN df). Should be a factor for binary/multiclass classification!  
C        : regularization parameter. (large: strong, small: weak)  
for l1 regularization C is a scalar: C|x|  
for elasticnet, C is a two dimensional vector: C(1)|x| + C(2)x^2/2  
opt  : list of options which control the behavior of Spicer  
- loss: type of loss function: 'logit', 'square'.  
(default: 'logit' for classification, 'square' for regression)  
logit: logistic regression, log(1+exp(- f(x)y))  
square: square loss, 0.5*(y - f(x))^2  
- regname: type of regularization: 'l1', 'elasticnet' (default:'l1')  
- outerMaxIter: maximum number of iteration of outer loop. (default 300)  
- innerMaxIter: maximum number of iteration of inner loop. (default 500)  
- stopdualitygap: TRUE/FALSE. If TRUE, Spicer employs duality gap for stopping criterion of outer loop. Default TRUE.  
- stopIneqViolation: TRUE/FALSE. If TRUE, Spicer employs violation of inequality for stopping criterion of outer loop. Default FALSE.  
- tolOuter: tollerance of stopping criteria of outer loop. (default 0.001)  
- tolInner: tollerance of stopping criteria of inner loop. (default tolOuter/1000)  
- calpha: increment factor of gamma: gamma^(t+1)=calpha*gamma^(t). (default 10)  
- display: 1:no-display, 2:display outer loop, 3:display inner loop.(default 2)  

Output:  
f(x) = \sum_{m \in activeset}(d(m)*km(x,:)*alpha) + b  
comb_alpha     :  N x 1 coefficient vector.  
kern_weight         :  1 x M kernel weight vector.  
bias         :  bias.  
activeset :  indices of kernels that are active ({m : kern_weight[m] is not zero}).  
sorted_kern_weight: vector of non-zero kernel weights sorted by magnitude.  
opt       : list of Spicer options used in run  
history     :  contains history of primal objective, dual objective, number of active kernels, and duality gap.  
- primalobj: primal objective  
- len_active: number of active kernels  
- dualitygap: duality gap  

Citation:  
1.Suzuki,Tomioka.SpicyMKL: a fast algorithm for Multiple Kernel Learning with thousands of kernels. Mach Learn (2011) 85:77–108  
2. https://github.com/VladoUzunangelov/SPICER


(c) Vlado Uzunangelov 2016  
uzunangelov@soe.ucsc.edu
