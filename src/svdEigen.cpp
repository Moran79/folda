#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Rcpp::List svdEigen(const Eigen::MatrixXd &A, bool uFlag = true, bool vFlag = true) {
  int computeOptions = 0;

  if (uFlag) {
    computeOptions |= Eigen::ComputeThinU;
  }
  if (vFlag) {
    computeOptions |= Eigen::ComputeThinV;
  }

  Eigen::BDCSVD<Eigen::MatrixXd> svd(A, computeOptions);

  Rcpp::List result;

  result["d"] = svd.singularValues();

  if (uFlag) {
    result["u"] = svd.matrixU();
  }
  if (vFlag) {
    result["v"] = svd.matrixV();
  }

  return result;
}
