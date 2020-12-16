#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//[[Rcpp::export]]
List grahm_schimdtCpp(arma::mat A) {
  int n = A.n_cols;
  int m = A.n_rows;
  arma::mat Q(m, n);
  Q.fill(0);
  arma::mat R(n, n);
  R.fill(0);  
  for (int j = 0; j < n; j++) {
    arma::vec v = A.col(j);
    if (j > 0) {
      for(int i = 0; i < j; i++) {
        R(i, j) = arma::as_scalar(Q.col(i).t() *  A.col(j));
        v = v - R(i, j) * Q.col(i);
      }
    }
    R(j, j) = arma::norm(v, 2);
    Q.col(j) = v / R(j, j);
  }
  return List::create(_["Q"] = Q,
                      _["R"] = R
  );
}