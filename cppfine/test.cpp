#include <cassert>
#include <string>

#include <Eigen/Dense>
#include "fine.hpp"

using namespace std;


int is_approx_equal(double a, double b, double eps) {
  return std::abs(b-a) < eps;}

double interest_rate = 0.02;

std::vector<double> strikes = {
   93.71,  86.03,  81.12,  77.60,  74.70,  72.16,  66.99,  61.37,
   99.56,  98.68,  97.28,  95.88,  94.64,  93.58,  91.75,  90.25,
  104.27, 104.63, 104.99, 105.30, 105.62, 105.93, 106.63, 107.66,
  122.87, 123.99, 124.85, 126.59, 126.46, 127.15, 128.59, 130.46,
  139.39, 141.02, 142.91, 144.56, 146.03, 147.36, 150.05, 153.28};

std::vector<double> maturities = {
  0.1190, 0.2381, 0.3571, 0.4762, 0.5952, 0.7143, 1.0714, 1.4286,
  0.1190, 0.2381, 0.3571, 0.4762, 0.5952, 0.7143, 1.0714, 1.4286,
  0.1190, 0.2381, 0.3571, 0.4762, 0.5952, 0.7143, 1.0714, 1.4286,
  0.1190, 0.2381, 0.3571, 0.4762, 0.5952, 0.7143, 1.0714, 1.4286,
  0.1190, 0.2381, 0.3571, 0.4762, 0.5952, 0.7143, 1.0714, 1.4286};

int test() {
  double underlying_price = 100;

  double volatility =  0.08;
  double kappa      =  3.0;
  double theta      =  0.1;
  double nu         =  0.25;
  double rho        = -0.8;

  Eigen::VectorXd prices_e(40);
  prices_e = fine::heston_formula(strikes, maturities, underlying_price,
                                  volatility, interest_rate, kappa, theta, nu,
                                  rho);
  std::vector<double> prices(prices_e.data(), prices_e.data() + prices_e.size());

  double volatility_hat, kappa_hat, theta_hat, nu_hat, rho_hat;
  std::tie(volatility_hat, kappa_hat, theta_hat, nu_hat,
           rho_hat) = fine::heston_calibration(prices, strikes, maturities,
                                               underlying_price, interest_rate);

  assert(is_approx_equal(volatility, volatility_hat, 0.01));
  assert(is_approx_equal(kappa,      kappa_hat,      0.01));
  assert(is_approx_equal(theta,      theta_hat,      0.01));
  assert(is_approx_equal(nu,         nu_hat,         0.01));
  assert(is_approx_equal(rho,        rho_hat,        0.01));

  return 0;}

int test_ts() {
  std::vector<double> underlying_price = {80, 100, 120};

  double volatility[] = {0.08, 0.12, 0.10};
  double kappa        =  3.0;
  double theta        =  0.1;
  double nu           =  0.25;
  double rho          = -0.8;

  std::vector<double> prices(120);

  Eigen::VectorXd::Map(&prices[ 0], 40) =
    fine::heston_formula(strikes, maturities, underlying_price[0],
                         volatility[0], interest_rate, kappa, theta, nu, rho);
  Eigen::VectorXd::Map(&prices[40], 40) =
    fine::heston_formula(strikes, maturities, underlying_price[1],
                         volatility[1], interest_rate, kappa, theta, nu, rho);
  Eigen::VectorXd::Map(&prices[80], 40) =
    fine::heston_formula(strikes, maturities, underlying_price[2],
                         volatility[2], interest_rate, kappa, theta, nu, rho);

  double kappa_hat, theta_hat, nu_hat, rho_hat;
  std::vector<double> volatility_hat;
  std::tie(volatility_hat, kappa_hat, theta_hat, nu_hat, rho_hat)
    = fine::heston_calibration_ts(prices, strikes, maturities, underlying_price,
                                  interest_rate);

  assert(is_approx_equal(volatility[0], volatility_hat[0], 0.01));
  assert(is_approx_equal(volatility[1], volatility_hat[1], 0.01));
  assert(is_approx_equal(volatility[2], volatility_hat[2], 0.01));
  assert(is_approx_equal(kappa,         kappa_hat,         0.01));
  assert(is_approx_equal(theta,         theta_hat,         0.01));
  assert(is_approx_equal(nu,            nu_hat,            0.01));
  assert(is_approx_equal(rho,           rho_hat,           0.01));

  return 0;}

int test_vol() {
  double underlying_price = 100;

  double volatility =  0.08;
  double kappa      =  3.0;
  double theta      =  0.1;
  double nu         =  0.25;
  double rho        = -0.8;

  Eigen::VectorXd prices_e(40);
  prices_e = fine::heston_formula(strikes, maturities, underlying_price,
                                  volatility, interest_rate, kappa, theta, nu,
                                  rho);
  std::vector<double> prices(prices_e.data(),
                             prices_e.data() + prices_e.size());

  // Calibrate using analytical gradient
  double volatility_hat =
    fine::heston_calibrate_vol(prices, strikes, maturities, underlying_price,
                               interest_rate, kappa, theta, nu, rho);

  assert(is_approx_equal(volatility, volatility_hat, 0.01));

  return 0;}

int test_delta() {
  double underlying_price = 100;

  double volatility =  0.08;
  double kappa      =  3.0;
  double theta      =  0.1;
  double nu         =  0.25;
  double rho        = -0.8;

  Eigen::VectorXd analytic(40);
  analytic = fine::heston_delta(strikes, maturities, underlying_price,
                                volatility, interest_rate, kappa, theta, nu,
                                rho);

  Eigen::VectorXd finite_diff(40);
  double h = 1;
  finite_diff = (
    fine::heston_formula(strikes, maturities, underlying_price + h/2,
                         volatility, interest_rate, kappa, theta, nu, rho)
    - fine::heston_formula(strikes, maturities, underlying_price - h/2,
                           volatility, interest_rate, kappa, theta, nu, rho))/h;

  assert(analytic.isApprox(finite_diff, 0.0001));

  return 0;}

int test_vega() {
  double underlying_price = 100;

  double volatility =  0.08;
  double kappa      =  3.0;
  double theta      =  0.1;
  double nu         =  0.25;
  double rho        = -0.8;

  Eigen::VectorXd analytic(40);
  analytic = fine::heston_vega(strikes, maturities, underlying_price,
                               volatility, interest_rate, kappa, theta, nu,
                               rho);

  double h = 0.0001;
  Eigen::VectorXd finite_diff(40);
  finite_diff = (
    fine::heston_formula(strikes, maturities, underlying_price,
                         volatility + h/2, interest_rate, kappa, theta, nu, rho)
    - fine::heston_formula(strikes, maturities, underlying_price,
                           volatility - h/2, interest_rate, kappa, theta, nu,
                           rho))/h;

  assert(analytic.isApprox(finite_diff, 0.0000001));

  return 0;}

int main() {
  return test() || test_ts() || test_vol() || test_delta() || test_vega();}
