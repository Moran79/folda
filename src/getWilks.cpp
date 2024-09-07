#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

double getWilks(const Eigen::MatrixXd& Sw, const Eigen::MatrixXd& St, double tolerance = 1e-5) {
  try {
    // Perform LLT decomposition of St
    Eigen::LLT<Eigen::MatrixXd> lltSt(St);

    // Check if St decomposition was successful
    if (lltSt.info() != Eigen::Success) {
      return 2.0;  // St is not positive definite
    }

    // Check diagonal elements of St to ensure strict positive definiteness
    const Eigen::MatrixXd& Lst = lltSt.matrixL();
    if (Lst.diagonal().minCoeff() <= tolerance) {
      return 2.0;  // St is not strictly positive definite
    }

    // Perform LLT decomposition of Sw
    Eigen::LLT<Eigen::MatrixXd> lltSw(Sw);

    // Check if Sw decomposition was successful
    if (lltSw.info() != Eigen::Success) {
      return 0.0;  // Sw is not positive definite, return 0
    }

    // Check diagonal elements of Sw to ensure strict positive definiteness
    const Eigen::MatrixXd& Lsw = lltSw.matrixL();
    if (Lsw.diagonal().minCoeff() <= tolerance) {
      return 0.0;  // Sw is not strictly positive definite, return 0
    }

    // Calculate the determinant ratio det(Sw) / det(St)
    double logDetSw = Lsw.diagonal().array().log().sum() * 2;
    double logDetSt = Lst.diagonal().array().log().sum() * 2;

    return std::exp(logDetSw - logDetSt);  // Return the determinant ratio
  } catch (...) {
    return 2.0;  // Return 2 if any error occurs
  }
}
