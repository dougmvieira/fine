#include <complex>
#include <tuple>
#include <vector>

#include <math.h>

#include "fine.hpp"


extern "C" {
  void c_formula(double *price, int n, double *maturities,
                 double underlying_price, double volatility, double *strikes,
                 double interest_rate, double kappa, double theta, double nu,
                 double rho);

  void c_deltafunc(double *delta, int n, double *maturities,
                   double underlying_price, double volatility, double *strikes,
                   double interest_rate, double kappa, double theta, double nu,
                   double rho);

  void c_vegafunc(double *vega, int n, double *maturities,
                  double underlying_price, double volatility, double *strikes,
                  double interest_rate, double kappa, double theta, double nu,
                  double rho);

  void c_calibrate(double *params, double *prices, double *t, double *k,
                   int d, double underlying_price, double interest_rate);

  void c_calibrate_ts(double *params, double *prices, double *t, double *k,
                      int d, double *underlying_price, int n,
                      double interest_rate);

  double c_calibrate_vol(double *prices, double *t, double *k, int d,
                         double underlying_price, double interest_rate,
                         double kappa, double theta, double nu, double rho);}

Eigen::VectorXd fine::heston_formula(std::vector<double> strikes,
                                     std::vector<double> maturities,
                                     double underlying_price, double volatility,
                                     double interest_rate, double kappa,
                                     double theta, double nu, double rho) {
  int d = strikes.size();
  Eigen::VectorXd prices(d);

  c_formula(&prices(0), strikes.size(), &maturities[0], underlying_price,
            volatility, &strikes[0], interest_rate, kappa, theta, nu, rho);

  return prices;}

Eigen::VectorXd fine::heston_delta(std::vector<double> strikes,
                                   std::vector<double> maturities,
                                   double underlying_price, double volatility,
                                   double interest_rate, double kappa,
                                   double theta, double nu, double rho) {
  int d = strikes.size();
  Eigen::VectorXd deltas(d);

  c_deltafunc(&deltas(0), strikes.size(), &maturities[0], underlying_price,
              volatility, &strikes[0], interest_rate, kappa, theta, nu, rho);

  return deltas;}

Eigen::VectorXd fine::heston_vega(std::vector<double> strikes,
                                  std::vector<double> maturities,
                                  double underlying_price, double volatility,
                                  double interest_rate, double kappa,
                                  double theta, double nu, double rho) {
   int d = strikes.size();
   Eigen::VectorXd vegas(d);

   c_vegafunc(&vegas(0), strikes.size(), &maturities[0], underlying_price,
              volatility, &strikes[0], interest_rate, kappa, theta, nu, rho);

   return vegas;}

quintet fine::heston_calibration(std::vector<double> prices,
                                 std::vector<double> strikes,
                                 std::vector<double> maturities,
                                 double underlying_price,
                                 double interest_rate) {
  double params[5];
  c_calibrate(&params[0], &prices[0], &maturities[0], &strikes[0],
              prices.size(), underlying_price, interest_rate);
  return std::make_tuple(params[0], params[1], params[2], params[3],
                         params[4]);}

vquintet fine::heston_calibration_ts(std::vector<double> prices,
                                     std::vector<double> strikes,
                                     std::vector<double> maturities,
                                     std::vector<double> underlying_prices,
                                     double interest_rate) {
  int n = underlying_prices.size();
  std::vector<double> params(n + 4);
  double kappa, theta, nu, rho;
  c_calibrate_ts(&params[0], &prices[0], &maturities[0], &strikes[0],
                 strikes.size(), &underlying_prices[0],
                 underlying_prices.size(), interest_rate);
  kappa = params[n    ];
  theta = params[n + 1];
  nu    = params[n + 2];
  rho   = params[n + 3];
  params.resize(n);

  return std::make_tuple(params, kappa, theta, nu, rho);}

double fine::heston_calibrate_vol(std::vector<double> prices,
                                  std::vector<double> strikes,
                                  std::vector<double> maturities,
                                  double underlying_price, double interest_rate,
                                  double kappa, double theta, double nu,
                                  double rho) {
  return c_calibrate_vol(&prices[0], &maturities[0], &strikes[0],
                         strikes.size(), underlying_price, interest_rate, kappa,
                         theta, nu, rho);}
