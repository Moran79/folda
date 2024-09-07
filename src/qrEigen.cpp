#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Rcpp::List qrEigen(const Eigen::MatrixXd &A) {
  Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
  Eigen::MatrixXd Q = qr.householderQ() * Eigen::MatrixXd::Identity(A.rows(), A.cols());
  Eigen::MatrixXd R = qr.matrixQR().topRows(A.cols()).template triangularView<Eigen::Upper>();

  return Rcpp::List::create(Rcpp::Named("Q") = Q, Rcpp::Named("R") = R);
}
