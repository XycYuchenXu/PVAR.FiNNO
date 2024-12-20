// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// lasso_regression
NumericVector lasso_regression(NumericMatrix gram_matrix, NumericVector xy, NumericVector beta, double lambda, NumericVector penalty_factors, int max_iter, double tolerance);
RcppExport SEXP _PVAR_FiNNO_lasso_regression(SEXP gram_matrixSEXP, SEXP xySEXP, SEXP betaSEXP, SEXP lambdaSEXP, SEXP penalty_factorsSEXP, SEXP max_iterSEXP, SEXP toleranceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type gram_matrix(gram_matrixSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xy(xySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type penalty_factors(penalty_factorsSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    rcpp_result_gen = Rcpp::wrap(lasso_regression(gram_matrix, xy, beta, lambda, penalty_factors, max_iter, tolerance));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_PVAR_FiNNO_lasso_regression", (DL_FUNC) &_PVAR_FiNNO_lasso_regression, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_PVAR_FiNNO(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
