module fine

  use iso_c_binding, only: c_double, c_double_complex, c_int
  use heston, only: deltafunc, vegafunc, formula, calibrate, calibrate_ts, &
                    calibrate_vol

  implicit none

contains

  subroutine c_formula(price, n, t, S, v, K, r, kappa, theta, nu, rho) bind(c)
    integer(c_int), value :: n
    real(c_double), intent(in) :: t(n), K(n)
    real(c_double), value :: S, v, r, kappa, theta, nu, rho
    real(c_double), intent(out) :: price(n)

    price = formula(t, S, v, K, r, kappa, theta, nu, rho)
  end subroutine c_formula

  subroutine c_deltafunc(delta, n, t, S, v, K, r, kappa, theta, nu, rho) bind(c)
    integer(c_int), value :: n
    real(c_double), intent(in) :: t(n), K(n)
    real(c_double), value :: S, v, r, kappa, theta, nu, rho
    real(c_double), intent(out) :: delta(n)

    delta = deltafunc(t, S, v, K, r, kappa, theta, nu, rho)
  end subroutine c_deltafunc

  subroutine c_vegafunc(vega, n, t, S, v, K, r, kappa, theta, nu, rho) bind(c)
    integer(c_int), value :: n
    real(c_double), intent(in) :: t(n), K(n)
    real(c_double), value :: S, v, r, kappa, theta, nu, rho
    real(c_double), intent(out) :: vega(n)

    vega = vegafunc(t, S, v, K, r, kappa, theta, nu, rho)
  end subroutine c_vegafunc

  subroutine c_calibrate(params, prices, t, K, d, S, r) bind(c)
    integer(c_int), value :: d
    real(c_double), value :: S, r
    real(c_double), intent(in) :: prices(d), t(d), k(d)
    real(c_double), intent(out) :: params(5)

    params = calibrate(prices, t, K, S, r)
  end subroutine c_calibrate

  subroutine c_calibrate_ts(params, prices, t, K, d, S, n, r) bind(c)
    integer(c_int), value :: d, n
    real(c_double), value :: r
    real(c_double), intent(in) :: prices(d, n), t(d), K(d), S(n)
    real(c_double), intent(out) :: params(n + 4)

    call calibrate_ts(params, prices, t, K, S, r)
  end subroutine c_calibrate_ts

  function c_calibrate_vol(prices, t, K, d, S, r, kappa, theta, nu, &
                           rho) result(v) bind(c)
    integer(c_int), value :: d
    real(c_double), value :: S, r, kappa, theta, nu, rho
    real(c_double), intent(in) :: prices(d), t(d), K(d)
    real(c_double) :: v

    v = calibrate_vol(prices, t, K, S, r, kappa, theta, nu, rho)
  end function c_calibrate_vol

end module fine
