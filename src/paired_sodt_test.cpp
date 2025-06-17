#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List sodt_permutation_loop(Rcpp::List G0_list, Rcpp::List G1_list, Rcpp::List HY_list) {
  int nperm = G0_list.size();

  arma::mat H0 = Rcpp::as<arma::mat>(HY_list[0]);
  int n = H0.n_rows;
  int d = H0.n_cols;

  Rcpp::NumericVector delta_SSB(nperm);
  Rcpp::NumericVector delta_SSW(nperm);
  Rcpp::NumericVector delta_F(nperm);
  Rcpp::NumericMatrix delta_SSW_g(nperm, d);

  for (int b = 0; b < nperm; ++b) {
    arma::mat G0 = Rcpp::as<arma::mat>(G0_list[b]);
    arma::mat G1 = Rcpp::as<arma::mat>(G1_list[b]);
    arma::mat H  = Rcpp::as<arma::mat>(HY_list[b]);

    arma::mat I = arma::eye(n, n);
    arma::mat R = I - H;

    double trace_G0 = arma::trace(G0);
    double SSB0 = arma::accu(H % G0);
    double SSW0 = trace_G0 - SSB0;
    arma::mat GR0 = R * G0;
    arma::vec SSW0_g(d);
    for (int j = 0; j < d; ++j) {
      arma::uvec idx = arma::find(H.col(j) > 0.5);
      arma::mat Rg = R.rows(idx);
      arma::mat GRg = GR0.rows(idx);
      SSW0_g(j) = arma::trace(GRg * Rg.t());
    }
    double F0 = (SSB0 / (d - 1)) / (SSW0 / (n - d));

    double trace_G1 = arma::trace(G1);
    double SSB1 = arma::accu(H % G1);
    double SSW1 = trace_G1 - SSB1;
    arma::mat GR1 = R * G1;
    arma::vec SSW1_g(d);
    for (int j = 0; j < d; ++j) {
      arma::uvec idx = arma::find(H.col(j) > 0.5);
      arma::mat Rg = R.rows(idx);
      arma::mat GRg = GR1.rows(idx);
      SSW1_g(j) = arma::trace(GRg * Rg.t());
    }
    double F1 = (SSB1 / (d - 1)) / (SSW1 / (n - d));

    delta_SSB[b] = SSB1 - SSB0;
    delta_SSW[b] = SSW1 - SSW0;
    delta_F[b] = F1 - F0;

    for (int j = 0; j < d; ++j) {
      delta_SSW_g(b, j) = SSW1_g(j) - SSW0_g(j);
    }
  }

  return Rcpp::List::create(
    Rcpp::_["delta_SSB"] = delta_SSB,
    Rcpp::_["delta_SSW"] = delta_SSW,
    Rcpp::_["delta_F"] = delta_F,
    Rcpp::_["delta_SSW_g"] = delta_SSW_g
  );
}


