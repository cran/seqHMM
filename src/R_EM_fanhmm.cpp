#include "config.h"
#include "EM_nhmm.h"
#include "fanhmm.h"
#include "create_Q.h"
#include "eta_to_gamma.h"
#include "list_to_field.h"

// [[Rcpp::export]]
Rcpp::List Rcpp_EM_LBFGS_fanhmm(
    const arma::field<arma::umat>& obs,
    const arma::uvec& Ti,
    const arma::uvec& M,
    const arma::mat& X_pi,
    const arma::field<arma::mat>& X_A,
    const Rcpp::List& X_B,
    const bool icpt_only_pi,
    const bool icpt_only_A,
    const arma::uvec& icpt_only_B,
    const bool iv_A,
    const arma::uvec& iv_B,
    const bool tv_A,
    const arma::uvec& tv_B,
    const arma::mat& eta_pi,
    const arma::cube& eta_A,
    const arma::field<arma::cube>& eta_B,
    const arma::vec& prior_y,
    const Rcpp::List& W_X_B, 
    const double lambda,
    const arma::uword maxeval, 
    const double ftol_abs, 
    const double ftol_rel, 
    const double xtol_abs, 
    const double xtol_rel, 
    const arma::uword print_level,
    const arma::uword maxeval_m, 
    const double ftol_abs_m, 
    const double ftol_rel_m, 
    const double xtol_abs_m, 
    const double xtol_rel_m, 
    const arma::uword print_level_m,
    const double bound) {
  
  arma::uword S = eta_A.n_slices;
  arma::uword C = obs(0).n_rows;
  arma::mat Qs = create_Q(S);
  arma::field<arma::mat> Qm(C);
  for (arma::uword c = 0; c < C; ++c) {
    Qm(c) = create_Q(M(c));
  }
  
  fanhmm model(
      obs, Ti, M, X_pi, X_A, matlist_to_2d_field(X_B), 
      icpt_only_pi, icpt_only_A, icpt_only_B, 
      iv_A, iv_B, tv_A, tv_B, eta_to_gamma(eta_pi, Qs), 
      eta_to_gamma(eta_A, Qs), eta_to_gamma(eta_B, Qm), prior_y, W_X_B
  );

  EM_nhmm EM(
      model, Qs, Qm, lambda,
      maxeval, ftol_abs, ftol_rel, xtol_abs, xtol_rel, print_level,
      maxeval_m, ftol_abs_m, ftol_rel_m, xtol_abs_m, xtol_rel_m, 
      print_level_m, bound
  );
  return EM.run();
}
