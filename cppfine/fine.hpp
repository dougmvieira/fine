#pragma once
#include <tuple>
#include <vector>
#include <Eigen/Dense>


typedef std::tuple<double, double, double, double, double> quintet;
typedef std::tuple<std::vector<double>, double, double, double, double>
        vquintet;


// user interface

namespace fine {
Eigen::VectorXd heston_formula(std::vector<double> strikes,
                               std::vector<double> maturities,
                               double underlying_price, double volatility,
                               double interest_rate, double kappa, double theta,
                               double nu, double rho);

Eigen::VectorXd heston_delta(std::vector<double> strikes,
                             std::vector<double> maturities,
                             double underlying_price, double volatility,
                             double interest_rate, double kappa, double theta,
                             double nu, double rho);

Eigen::VectorXd heston_vega(std::vector<double> strikes,
                            std::vector<double> maturities,
                            double underlying_price, double volatility,
                            double interest_rate, double kappa, double theta,
                            double nu, double rho);

quintet heston_calibration(std::vector<double> prices,
                           std::vector<double> strikes,
                           std::vector<double> maturities,
                           double underlying_price,
                           double interest_rate);

vquintet heston_calibration_ts(std::vector<double> prices,
                               std::vector<double> strikes,
                               std::vector<double> maturities,
                               std::vector<double> underlying_prices,
                               double interest_rate);

double heston_calibrate_vol(std::vector<double> prices,
                            std::vector<double> strikes,
                            std::vector<double> maturities,
                            double underlying_price, double interest_rate,
                            double kappa, double theta, double nu, double rho);
}
