// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// getPillai
double getPillai(const Eigen::MatrixXd& Sw, const Eigen::MatrixXd& St, double tolerance);
RcppExport SEXP _folda_getPillai(SEXP SwSEXP, SEXP StSEXP, SEXP toleranceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Sw(SwSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type St(StSEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    rcpp_result_gen = Rcpp::wrap(getPillai(Sw, St, tolerance));
    return rcpp_result_gen;
END_RCPP
}
// getWilks
double getWilks(const Eigen::MatrixXd& Sw, const Eigen::MatrixXd& St, double tolerance);
RcppExport SEXP _folda_getWilks(SEXP SwSEXP, SEXP StSEXP, SEXP toleranceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Sw(SwSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type St(StSEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    rcpp_result_gen = Rcpp::wrap(getWilks(Sw, St, tolerance));
    return rcpp_result_gen;
END_RCPP
}
// qrEigen
Rcpp::List qrEigen(const Eigen::MatrixXd& A);
RcppExport SEXP _folda_qrEigen(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(qrEigen(A));
    return rcpp_result_gen;
END_RCPP
}
// svdEigen
Rcpp::List svdEigen(const Eigen::MatrixXd& A, bool uFlag, bool vFlag);
RcppExport SEXP _folda_svdEigen(SEXP ASEXP, SEXP uFlagSEXP, SEXP vFlagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type A(ASEXP);
    Rcpp::traits::input_parameter< bool >::type uFlag(uFlagSEXP);
    Rcpp::traits::input_parameter< bool >::type vFlag(vFlagSEXP);
    rcpp_result_gen = Rcpp::wrap(svdEigen(A, uFlag, vFlag));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_folda_getPillai", (DL_FUNC) &_folda_getPillai, 3},
    {"_folda_getWilks", (DL_FUNC) &_folda_getWilks, 3},
    {"_folda_qrEigen", (DL_FUNC) &_folda_qrEigen, 1},
    {"_folda_svdEigen", (DL_FUNC) &_folda_svdEigen, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_folda(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
