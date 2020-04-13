// (c) Vlado Uzunangelov 2016
// uzunangelov@soe.ucsc.edu

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List kernel_norm_cpp ( NumericVector K ,arma::mat u, Nullable<IntegerVector> idx = R_NilValue) {

  // need a separate variable to store seq_len output in case of NULL -
  // can't use idx argument! (convoluted, but ok)
  IntegerVector idxx;
  if (idx.isNull()) {
    // indices are still 1-based (like in the case when they are not null)
    idxx = seq_len(u.n_cols);
  } else {
    idxx = idx;
  }

  // in keeping with convention that when there are matching data types,
  // the arma one is better
  arma::uvec idx_ = as<arma::uvec> (idxx);

  //convert the NumericVector to arma::cube
  //all this should be unnecessary with the latest RcppArmadillo,
  // but that requires g++ >= 4.6
  IntegerVector Kdim = K.attr("dim");
  arma::cube kernels(K.begin(),Kdim(0),Kdim(1),Kdim(2),false);

  int N = u.n_rows;
  int M = idx_.n_elem;

  arma::mat wdot = arma::zeros(N,M);
  arma::vec normj = arma::zeros(M);

  // convert to zero-based indexing
  idx_ = idx_ - 1;
  if (M > 0) {
    for (int i = 0; i < M ; i++) {
      arma::uword id = idx_[i];
      wdot.col(i) = kernels.slice(id) * u.col(id);
      double  dot = arma::as_scalar(wdot.col(i).t()* u.col(id));
      normj(i) = sqrt( max( NumericVector::create ( dot , 0.0 ) ) ) ;
    }
  }

  // need to change normj to a std::vector, otherwise it gets returned as a Nx1 matrix
  return List::create(Named("wdot") = wdot,
                      Named("normj") = as<std::vector<double> >(wrap(normj)));

}


