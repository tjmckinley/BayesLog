// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// logisticMH
List logisticMH(arma::mat data, arma::vec nsamples, int nrandint, arma::ivec randint, arma::ivec cumrandindex, IntegerVector factindex, IntegerVector cumfactindex, NumericVector ini_pars, NumericVector ini_sigma, double ini_sigmarand, int gen_inits, NumericMatrix priors, int niter, int nitertraining, double scale, int orignpars, int varselect, int ninitial, int random, int nprintsum);
RcppExport SEXP BayesLog_logisticMH(SEXP dataSEXP, SEXP nsamplesSEXP, SEXP nrandintSEXP, SEXP randintSEXP, SEXP cumrandindexSEXP, SEXP factindexSEXP, SEXP cumfactindexSEXP, SEXP ini_parsSEXP, SEXP ini_sigmaSEXP, SEXP ini_sigmarandSEXP, SEXP gen_initsSEXP, SEXP priorsSEXP, SEXP niterSEXP, SEXP nitertrainingSEXP, SEXP scaleSEXP, SEXP orignparsSEXP, SEXP varselectSEXP, SEXP ninitialSEXP, SEXP randomSEXP, SEXP nprintsumSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP );
        Rcpp::traits::input_parameter< arma::vec >::type nsamples(nsamplesSEXP );
        Rcpp::traits::input_parameter< int >::type nrandint(nrandintSEXP );
        Rcpp::traits::input_parameter< arma::ivec >::type randint(randintSEXP );
        Rcpp::traits::input_parameter< arma::ivec >::type cumrandindex(cumrandindexSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type factindex(factindexSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type cumfactindex(cumfactindexSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type ini_pars(ini_parsSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type ini_sigma(ini_sigmaSEXP );
        Rcpp::traits::input_parameter< double >::type ini_sigmarand(ini_sigmarandSEXP );
        Rcpp::traits::input_parameter< int >::type gen_inits(gen_initsSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type priors(priorsSEXP );
        Rcpp::traits::input_parameter< int >::type niter(niterSEXP );
        Rcpp::traits::input_parameter< int >::type nitertraining(nitertrainingSEXP );
        Rcpp::traits::input_parameter< double >::type scale(scaleSEXP );
        Rcpp::traits::input_parameter< int >::type orignpars(orignparsSEXP );
        Rcpp::traits::input_parameter< int >::type varselect(varselectSEXP );
        Rcpp::traits::input_parameter< int >::type ninitial(ninitialSEXP );
        Rcpp::traits::input_parameter< int >::type random(randomSEXP );
        Rcpp::traits::input_parameter< int >::type nprintsum(nprintsumSEXP );
        List __result = logisticMH(data, nsamples, nrandint, randint, cumrandindex, factindex, cumfactindex, ini_pars, ini_sigma, ini_sigmarand, gen_inits, priors, niter, nitertraining, scale, orignpars, varselect, ninitial, random, nprintsum);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
