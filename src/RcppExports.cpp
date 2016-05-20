// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// logisticMH
NumericMatrix logisticMH(NumericMatrix dataR, NumericVector nsamplesR, NumericVector ini_pars, NumericMatrix priors, int niter, int ninitial, double scale, int nadapt, int nprintsum, double maxscale, double niterdim, int nrand, List randindexesL, IntegerMatrix data_randR, IntegerVector nblock, List blockR, int printini, List noncentreintR, NumericVector noncentreintRE);
RcppExport SEXP BayesLog_logisticMH(SEXP dataRSEXP, SEXP nsamplesRSEXP, SEXP ini_parsSEXP, SEXP priorsSEXP, SEXP niterSEXP, SEXP ninitialSEXP, SEXP scaleSEXP, SEXP nadaptSEXP, SEXP nprintsumSEXP, SEXP maxscaleSEXP, SEXP niterdimSEXP, SEXP nrandSEXP, SEXP randindexesLSEXP, SEXP data_randRSEXP, SEXP nblockSEXP, SEXP blockRSEXP, SEXP printiniSEXP, SEXP noncentreintRSEXP, SEXP noncentreintRESEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type dataR(dataRSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nsamplesR(nsamplesRSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ini_pars(ini_parsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type priors(priorsSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< int >::type ninitial(ninitialSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< int >::type nadapt(nadaptSEXP);
    Rcpp::traits::input_parameter< int >::type nprintsum(nprintsumSEXP);
    Rcpp::traits::input_parameter< double >::type maxscale(maxscaleSEXP);
    Rcpp::traits::input_parameter< double >::type niterdim(niterdimSEXP);
    Rcpp::traits::input_parameter< int >::type nrand(nrandSEXP);
    Rcpp::traits::input_parameter< List >::type randindexesL(randindexesLSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type data_randR(data_randRSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nblock(nblockSEXP);
    Rcpp::traits::input_parameter< List >::type blockR(blockRSEXP);
    Rcpp::traits::input_parameter< int >::type printini(printiniSEXP);
    Rcpp::traits::input_parameter< List >::type noncentreintR(noncentreintRSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type noncentreintRE(noncentreintRESEXP);
    __result = Rcpp::wrap(logisticMH(dataR, nsamplesR, ini_pars, priors, niter, ninitial, scale, nadapt, nprintsum, maxscale, niterdim, nrand, randindexesL, data_randR, nblock, blockR, printini, noncentreintR, noncentreintRE));
    return __result;
END_RCPP
}
// classification
List classification(NumericMatrix pred, IntegerVector obs, NumericVector thresh);
RcppExport SEXP BayesLog_classification(SEXP predSEXP, SEXP obsSEXP, SEXP threshSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type pred(predSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thresh(threshSEXP);
    __result = Rcpp::wrap(classification(pred, obs, thresh));
    return __result;
END_RCPP
}
