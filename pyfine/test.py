import numpy as np
from fine import (heston_formula, heston_delta, heston_vega,
                  heston_calibration, heston_calibration_ts,
                  heston_calibrate_vol)


interest_rate = 0.02

strikes = np.array([
     93.71,  86.03,  81.12,  77.60,  74.70,  72.16,  66.99,  61.37,
     99.56,  98.68,  97.28,  95.88,  94.64,  93.58,  91.75,  90.25,
    104.27, 104.63, 104.99, 105.30, 105.62, 105.93, 106.63, 107.66,
    122.87, 123.99, 124.85, 126.59, 126.46, 127.15, 128.59, 130.46,
    139.39, 141.02, 142.91, 144.56, 146.03, 147.36, 150.05, 153.28])

maturities = np.array([
    0.1190, 0.2381, 0.3571, 0.4762, 0.5952, 0.7143, 1.0714, 1.4286,
    0.1190, 0.2381, 0.3571, 0.4762, 0.5952, 0.7143, 1.0714, 1.4286,
    0.1190, 0.2381, 0.3571, 0.4762, 0.5952, 0.7143, 1.0714, 1.4286,
    0.1190, 0.2381, 0.3571, 0.4762, 0.5952, 0.7143, 1.0714, 1.4286,
    0.1190, 0.2381, 0.3571, 0.4762, 0.5952, 0.7143, 1.0714, 1.4286])

def test():
    underlying_price = 100
    volatility =  0.08
    kappa      =  3.0
    theta      =  0.1
    nu         =  0.25
    rho        = -0.8

    prices = heston_formula(strikes, maturities, underlying_price, volatility,
                            interest_rate, kappa, theta, nu, rho)

    (volatility_hat, kappa_hat, theta_hat, nu_hat,
     rho_hat) = heston_calibration(prices, strikes, maturities,
                                   underlying_price, interest_rate)

    assert np.abs(volatility - volatility_hat) < 0.01
    assert np.abs(kappa      - kappa_hat)      < 0.01
    assert np.abs(theta      - theta_hat)      < 0.01
    assert np.abs(nu         - nu_hat)         < 0.01
    assert np.abs(rho        - rho_hat)        < 0.01

def test_ts():
    underlying_prices = [80, 100, 120]
    volatilities =  [0.08, 0.12, 0.10]
    kappa        =  3.0
    theta        =  0.1
    nu           =  0.25
    rho          = -0.8

    prices = np.concatenate([heston_formula(strikes, maturities,
                                            underlying_price, volatility,
                                            interest_rate, kappa, theta, nu,
                                            rho)
                             for underlying_price, volatility
                             in zip(underlying_prices, volatilities)])
    (volatilities_hat, kappa_hat, theta_hat, nu_hat,
     rho_hat) = heston_calibration_ts(prices, strikes, maturities,
                                      underlying_prices, interest_rate)

    assert np.abs(volatilities[0] - volatilities_hat[0]) < 0.01
    assert np.abs(volatilities[1] - volatilities_hat[1]) < 0.01
    assert np.abs(volatilities[2] - volatilities_hat[2]) < 0.01
    assert np.abs(kappa           - kappa_hat)           < 0.01
    assert np.abs(theta           - theta_hat)           < 0.01
    assert np.abs(nu              - nu_hat)              < 0.01
    assert np.abs(rho             - rho_hat)             < 0.01

def test_vol():
    underlying_price = 100
    volatility =  0.08
    kappa      =  3.0
    theta      =  0.1
    nu         =  0.25
    rho        = -0.8

    prices = heston_formula(strikes, maturities, underlying_price, volatility,
                            interest_rate, kappa, theta, nu, rho)
    volatility_hat = heston_calibrate_vol(prices, strikes, maturities,
                                          underlying_price, interest_rate,
                                          kappa, theta, nu, rho)

    assert np.abs(volatility - volatility_hat) < 0.01

def test_delta():
    underlying_price = 100
    volatility =  0.08
    kappa      =  3.0
    theta      =  0.1
    nu         =  0.25
    rho        = -0.8

    analytic = heston_delta(strikes, maturities, underlying_price, volatility,
                            interest_rate, kappa, theta, nu, rho)

    h = 0.1
    finite_diff = (heston_formula(strikes, maturities, underlying_price + h/2,
                                  volatility, interest_rate, kappa, theta, nu,
                                  rho)
                   - heston_formula(strikes, maturities, underlying_price - h/2,
                                    volatility, interest_rate, kappa, theta, nu,
                                    rho))/h

    assert np.all(np.abs(finite_diff - analytic) < 0.00001)

def test_vega():
    underlying_price = 100
    volatility =  0.08
    kappa      =  3.0
    theta      =  0.1
    nu         =  0.25
    rho        = -0.8

    analytic = heston_vega(strikes, maturities, underlying_price, volatility,
                           interest_rate, kappa, theta, nu, rho)

    h = 0.0001
    finite_diff = (heston_formula(strikes, maturities, underlying_price,
                                  volatility + h/2, interest_rate, kappa, theta,
                                  nu, rho)
                   - heston_formula(strikes, maturities, underlying_price,
                                    volatility - h/2, interest_rate, kappa, theta,
                                    nu, rho))/h

    assert np.all(np.abs(finite_diff - analytic) < 0.001)

test()
test_ts()
test_vol()
test_delta()
test_vega()
