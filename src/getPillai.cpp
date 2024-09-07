#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

double getPillai(const Eigen::MatrixXd& Sw, const Eigen::MatrixXd& St, double tolerance = 1e-5) {
  try {
    // Perform LLT decomposition of St
    Eigen::LLT<Eigen::MatrixXd> llt(St);

    // Check if decomposition was successful
    if (llt.info() != Eigen::Success) {
      return -1.0;  // St is not positive definite
    }

    // Check the diagonal elements of the lower triangular matrix to ensure strict positive definiteness
    Eigen::MatrixXd L = llt.matrixL(); // Get the lower triangular matrix
    if (L.diagonal().minCoeff() <= tolerance) {
      return -1.0;  // St is not strictly positive definite
    }

    // Compute Sb = St - Sw
    Eigen::MatrixXd Sb = St - Sw;

    // Solve St * X = Sb
    Eigen::MatrixXd result = llt.solve(Sb);

    // Calculate the trace of the result matrix
    double traceResult = result.trace();

    return traceResult;
  } catch (...) {
    return -1.0;  // Return -1 if any error occurs
  }
}
