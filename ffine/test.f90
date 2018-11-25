program test
  use numerics, only: dp
  use fine, only: formula, deltafunc, vegafunc, calibrate, calibrate_ts, &
                  calibrate_vol

  implicit none

  real(dp) :: r = 0.02

  real(dp) :: K(40) = [ &
     93.71,  86.03,  81.12,  77.60,  74.70,  72.16,  66.99,  61.37, &
     99.56,  98.68,  97.28,  95.88,  94.64,  93.58,  91.75,  90.25, &
    104.27, 104.63, 104.99, 105.30, 105.62, 105.93, 106.63, 107.66, &
    122.87, 123.99, 124.85, 126.59, 126.46, 127.15, 128.59, 130.46, &
    139.39, 141.02, 142.91, 144.56, 146.03, 147.36, 150.05, 153.28]

  real(dp) :: t(40) = [ &
    0.1190, 0.2381, 0.3571, 0.4762, 0.5952, 0.7143, 1.0714, 1.4286, &
    0.1190, 0.2381, 0.3571, 0.4762, 0.5952, 0.7143, 1.0714, 1.4286, &
    0.1190, 0.2381, 0.3571, 0.4762, 0.5952, 0.7143, 1.0714, 1.4286, &
    0.1190, 0.2381, 0.3571, 0.4762, 0.5952, 0.7143, 1.0714, 1.4286, &
    0.1190, 0.2381, 0.3571, 0.4762, 0.5952, 0.7143, 1.0714, 1.4286]

  call test_calibrate()
  call test_calibrate_ts()
  call test_calibrate_vol()
  call test_deltafunc()
  call test_vegafunc()

contains

  subroutine assert(condition)
    logical :: condition

    if(.not. condition) stop 1
  end subroutine assert

  subroutine test_calibrate()
    real(dp) :: s = 100
    real(dp) :: v     =  0.08
    real(dp) :: kappa =  3.0
    real(dp) :: theta =  0.1
    real(dp) :: nu    =  0.25
    real(dp) :: rho   = -0.8
    real(dp) :: prices(40), params_hat(5)

    prices = formula(t, s, v, k, r, kappa, theta, nu, rho)
    params_hat = calibrate(prices, t, k, s, r)

    call assert(abs(v     - params_hat(1)) < 0.01)
    call assert(abs(kappa - params_hat(2)) < 0.01)
    call assert(abs(theta - params_hat(3)) < 0.01)
    call assert(abs(nu    - params_hat(4)) < 0.01)
    call assert(abs(rho   - params_hat(5)) < 0.01)
  end subroutine test_calibrate

  subroutine test_calibrate_ts()
    real(dp) :: S(3) = [80, 100, 120]
    real(dp) :: v(3) = [0.08, 0.12, 0.10]
    real(dp) :: kappa =  3.0
    real(dp) :: theta =  0.1
    real(dp) :: nu    =  0.25
    real(dp) :: rho   = -0.8
    real(dp) :: prices(40, 3), params_hat(7)
    integer :: i

    do i = 1, 3
       prices(:, i) = formula(t, S(i), v(i), K, r, kappa, theta, nu, rho)
    end do

    call calibrate_ts(params_hat, prices, t, K, S, r)

    call assert(abs(v(1)  - params_hat(1)) < 0.01)
    call assert(abs(v(2)  - params_hat(2)) < 0.01)
    call assert(abs(v(3)  - params_hat(3)) < 0.01)
    call assert(abs(kappa - params_hat(4)) < 0.01)
    call assert(abs(theta - params_hat(5)) < 0.01)
    call assert(abs(nu    - params_hat(6)) < 0.01)
    call assert(abs(rho   - params_hat(7)) < 0.01)
  end subroutine test_calibrate_ts

  subroutine test_calibrate_vol()
    real(dp) :: S = 100
    real(dp) :: v     =  0.08
    real(dp) :: kappa =  3.0
    real(dp) :: theta =  0.1
    real(dp) :: nu    =  0.25
    real(dp) :: rho   = -0.8
    real(dp) :: prices(40), v_hat

    prices = formula(t, S, v, K, r, kappa, theta, nu, rho)
    v_hat = calibrate_vol(prices, t, K, S, r, kappa, theta, nu, rho)

    call assert(abs(v - v_hat) < 0.01)
  end subroutine test_calibrate_vol

  subroutine test_deltafunc()
    real(dp) :: S = 100
    real(dp) :: v     =  0.08
    real(dp) :: kappa =  3.0
    real(dp) :: theta =  0.1
    real(dp) :: nu    =  0.25
    real(dp) :: rho   = -0.8
    real(dp) :: analytic(40), finite_diff(40)
    real(dp) :: h = 0.1
    integer :: i

    analytic = deltafunc(t, S, v, K, r, kappa, theta, nu, rho)
    finite_diff = (formula(t, S + h/2, v, K, r, kappa, theta, nu, rho) &
                  -formula(t, S - h/2, v, K, r, kappa, theta, nu, rho))/h

    do i = 1, 40
      call assert(abs(finite_diff(i) - analytic(i)) < 0.00001)
    end do
  end subroutine test_deltafunc

  subroutine test_vegafunc()
    real(dp) :: S = 100
    real(dp) :: v     =  0.08
    real(dp) :: kappa =  3.0
    real(dp) :: theta =  0.1
    real(dp) :: nu    =  0.25
    real(dp) :: rho   = -0.8
    real(dp) :: analytic(40), finite_diff(40)
    real(dp) :: h = 0.0001
    integer :: i

    analytic = vegafunc(t, S, v, K, r, kappa, theta, nu, rho)
    finite_diff = (formula(t, S, v + h/2, K, r, kappa, theta, nu, rho) &
                  -formula(t, S, v - h/2, K, r, kappa, theta, nu, rho))/h

    do i = 1, 40
       call assert(abs(finite_diff(i) - analytic(i)) < 0.00001)
    end do
  end subroutine test_vegafunc

end program test
