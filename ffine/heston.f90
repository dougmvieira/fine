module heston
  use numerics, only: dp, unused, quad, quadvec, evalunique, matevalunique, &
                      newtonstep
  use ieee_arithmetic, only: ieee_value,  ieee_positive_inf

  implicit none

  private
  public formula, deltafunc, vegafunc, calibrate, calibrate_ts, calibrate_vol

  complex(dp), parameter :: i = cmplx(0, 1, kind=dp)
  complex(dp), parameter :: zero = cmplx(0, 0, kind=dp)
  real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)

contains

  function psigradfunc(u, t, kappa, a, nu, rho) result(psigrad)
    complex(dp), intent(in) :: u
    real(dp), intent(in) :: t, kappa, a, nu, rho
    complex(dp) :: psigrad(5, 2), psigrad_volnorm(5, 2)

    psigrad_volnorm = psigradfunc_volnorm(u, nu*t, kappa/nu, a/nu**2, rho)
    psigrad(1, :) = psigrad_volnorm(1, :)
    psigrad(2, :) = psigrad_volnorm(3, :)/nu
    psigrad(3, :) = psigrad_volnorm(4, :)/nu**2
    psigrad(4, :) = t*psigrad_volnorm(2, :) &
                    - kappa*psigrad_volnorm(3, :)/nu**2 &
                    - 2*a*psigrad_volnorm(4, :)/nu**3
    psigrad(4, 2) = psigrad(4, 2) - psigrad_volnorm(1, 2)/nu
    psigrad(5, :) = psigrad_volnorm(5, :)
    psigrad(:, 2) = psigrad(:, 2)/nu
  end function psigradfunc

  function psigradfunc_volnorm(u, t, kappa, a, rho) result(psigrad)
    complex(dp), intent(in) :: u
    real(dp), intent(in) :: t, kappa, a, rho
    complex(dp) :: psigrad(5, 2)

    include 'heston_psi_grad.f90'

    psigrad(1, 1) = psi_1
    psigrad(1, 2) = psi_2
    psigrad(2, 1) = psi_1_t_chain
    psigrad(2, 2) = psi_2_t_chain
    psigrad(3, 1) = psi_1_kappa_chain
    psigrad(3, 2) = psi_2_kappa_chain
    psigrad(4, 1) = psi_1_a_chain
    psigrad(4, 2) = psi_2_a_chain
    psigrad(5, 1) = psi_1_rho_chain
    psigrad(5, 2) = psi_2_rho_chain
  end function psigradfunc_volnorm

  pure function psifunc_volnorm(u, t, kappa, a, rho) result(psi)
    complex(dp), intent(in) :: u
    real(dp), intent(in) :: t, kappa, a, rho
    complex(dp) :: psi(2)

    include 'heston_psi.f90'

    psi(1) = psi_1
    psi(2) = psi_2
  end function psifunc_volnorm

  pure function psifunc(u, t, kappa, a, nu, rho) result(psi)
    complex(dp), intent(in) :: u
    real(dp), intent(in) :: t, kappa, a, nu, rho
    complex(dp) :: psi(2)

    psi = psifunc_volnorm(u, nu*t, kappa/nu, a/nu**2, rho)
    psi(2) = psi(2)/nu
  end function psifunc

  function integrandfunc(u, t, v, k, kappa, a, nu, rho, &
                         alpha) result(integrand)
    real(dp), intent(in) :: u
    real(dp), intent(in) :: t, v, k(:), kappa, a, nu, rho, alpha
    real(dp) :: integrand(size(K))
    complex(dp) :: logvarphi, psi(2), num(size(k)), den

    psi = psifunc(u - (alpha + 1)*i, t, kappa, a, nu, rho)
    logvarphi = psi(1) + v*psi(2)
    num = exp(logvarphi - i*u*k)
    den = cmplx(u*u - alpha*(alpha + 1), -2*u*alpha - u, kind=dp)
    integrand = -(real(num)*real(den) + aimag(num)*aimag(den)) &
                /(real(den)*real(den) + aimag(den)*aimag(den))
  end function integrandfunc

  function pricegradintegrandfunc(u, t, v, k, kappa, a, nu, rho, &
                                  alpha) result(integrand)
    real(dp), intent(in) :: u
    real(dp), intent(in) :: t, v, k(:), kappa, a, nu, rho, alpha
    real(dp) :: integrand(size(k), 6)
    complex(dp) :: psigrad(5, 2), logvarphigrad(5), num(size(k), 6), den

    psigrad = psigradfunc(u - (alpha + 1)*i, t, kappa, a, nu, rho)
    logvarphigrad = psigrad(:, 1) + v*psigrad(:, 2)
    num(:, 1) = exp(logvarphigrad(1) - i*u*k)
    num(:, 2) = psigrad(1, 2)*num(:, 1)
    num(:, 3) = logvarphigrad(2)*num(:, 1)
    num(:, 4) = logvarphigrad(3)*num(:, 1)
    num(:, 5) = logvarphigrad(4)*num(:, 1)
    num(:, 6) = logvarphigrad(5)*num(:, 1)
    den = cmplx(u*u - alpha*(alpha + 1), -2*u*alpha - u, kind=dp)
    integrand = -(real(num)*real(den) + aimag(num)*aimag(den)) &
                /(real(den)*real(den) + aimag(den)*aimag(den))
  end function pricegradintegrandfunc

  elemental function deltaintegrandfunc(u, t, v, k, kappa, a, nu, &
                                        rho) result(integrand)
    real(dp), intent(in) :: u
    real(dp), intent(in) :: t, v, k, kappa, a, nu, rho
    real(dp) :: integrand
    complex(dp) :: psi(2), logvarphi

    psi = psifunc(u - i, t, kappa, a, nu, rho)
    logvarphi = psi(1) + v*psi(2)

    integrand = real(exp(-u*i*k + logvarphi)/(u*i))
  end function deltaintegrandfunc

  elemental function vegaintegrandfunc(u, t, v, k, kappa, a, nu, rho, &
                                       alpha) result(integrand)
    real(dp), intent(in) :: u
    real(dp), intent(in) :: t, v, k, kappa, a, nu, rho, alpha
    complex(dp) :: psi(2), logvarphi, num, den
    real(dp) :: integrand

    psi = psifunc(u - (alpha + 1)*i, t, kappa, a, nu, rho)
    logvarphi = psi(1) + v*psi(2)
    num = psi(2)*exp(logvarphi - i*u*k)
    den = cmplx(u*u - alpha*(alpha + 1), -2*u*alpha - u, kind=dp)
    integrand = -(real(num)*real(den) + aimag(num)*aimag(den)) &
                /(real(den)*real(den) + aimag(den)*aimag(den))
  end function vegaintegrandfunc

  function deltafunc_normalised(t, v, K, kappa, theta, nu, rho) result(delta)
    real(dp), intent(in) :: t, v, K, kappa, theta, nu, rho
    real(dp) :: delta, logK

    logK = log(K)
    delta = 0.5 + quad(closure)/pi
  contains
    function closure(u) result(integrand)
      real(dp), intent(in) :: u
      real(dp) :: integrand

      integrand = deltaintegrandfunc(u, t, v, logK, kappa, kappa*theta, nu, rho)
    end function closure
  end function deltafunc_normalised

  function deltafunc(t, S, v, K, r, kappa, theta, nu, rho) result(delta)
    real(dp), intent(in) :: t(:), S, v, K(:), r, kappa, theta, nu, rho
    real(dp) :: delta(size(K)), Ktilde(size(K))
    integer :: i

    Ktilde = K/(exp(r*t)*S)

    do i = 1, size(K)
       delta(i) = deltafunc_normalised(t(i), v, Ktilde(i), kappa, theta, nu, &
                                       rho)
    end do
  end function deltafunc

  function vegafunc_normalised(t, v, K, kappa, theta, nu, rho) result(vega)
    real(dp), intent(in) :: t, v, K, kappa, theta, nu, rho
    real(dp) :: vega, logK

    logK = log(K)
    vega = quad(closure)/pi
  contains
    function closure(u) result(integrand)
      real(dp), intent(in) :: u
      real(dp) :: integrand

      integrand = vegaintegrandfunc(u, t, v, logK, kappa, kappa*theta, nu, &
                                    rho, 0._dp)
    end function closure
  end function vegafunc_normalised

  function vegafunc(t, S, v, K, r, kappa, theta, nu, rho) result(vega)
    real(dp), intent(in) :: t(:), S, v, K(:), r, kappa, theta, nu, rho
    real(dp) :: vega(size(K)), Ktilde(size(K))
    integer :: i

    Ktilde = K/(exp(r*t)*S)

    do i = 1, size(K)
       vega(i) = S*vegafunc_normalised(t(i), v, Ktilde(i), kappa, theta, nu, &
                                       rho)
    end do
  end function vegafunc

  function smilepricegradfunc(t, v, K, kappa, a, nu, rho) result(pricegrad)
    real(dp), intent(in) :: t, v, K(:), kappa, a, nu, rho
    real(dp) :: pricegrad(size(K), 6), integral(size(K)*6), logK(size(K))

    logK = log(K)
    call quadvec(closure, integral)
    pricegrad = reshape(integral, [size(K), 6])/pi
    pricegrad(:, 1) = 0.5_dp + pricegrad(:, 1)
  contains
    subroutine closure(u, integrand)
      real(dp), intent(in) :: u
      real(dp), intent(out) :: integrand(:)

      integrand = reshape(pricegradintegrandfunc(u, t, v, logK, kappa, a, nu, &
                                                 rho, 0._dp), [6*size(K)])
    end subroutine closure
  end function smilepricegradfunc

  function pricegradfunc(t, v, K, kappa, a, nu, rho) result(pricegrad)
    real(dp), intent(in) :: t(:), v, K(:), kappa, a, nu, rho
    real(dp) :: pricegrad(size(K), 6)

    pricegrad = matevalunique(matevalunique_wrap, K, t, 6)
  contains
    function matevalunique_wrap(K, t, d) result(price)
      integer, intent(in) :: d
      real(dp), intent(in) :: K(:), t
      real(dp) :: price(size(K), d)
      price = smilepricegradfunc(t, v, K, kappa, a, nu, rho)
    end function matevalunique_wrap
  end function pricegradfunc

  function smileformula(t, v, K, kappa, a, nu, rho) result(price)
    real(dp), intent(in) :: t, v, K(:), kappa, a, nu, rho
    real(dp) :: price(size(K)), integral(size(K))

    call quadvec(closure, integral)
    price = 0.5_dp + integral/pi
  contains
    subroutine closure(u, integrand)
      real(dp), intent(in) :: u
      real(dp), intent(out) :: integrand(:)

      integrand = integrandfunc(u, t, v, log(K), kappa, a, nu, rho, 0._dp)
    end subroutine closure
  end function smileformula

  function formula_normalised(t, v, K, kappa, theta, nu, rho) result(price)
    real(dp), intent(in) :: t(:), v, K(:), kappa, theta, nu, rho
    real(dp) :: price(size(K))

    price = evalunique(evalunique_wrap, K, t)
  contains
    function evalunique_wrap(K, t) result(price)
      real(dp), intent(in) :: K(:), t
      real(dp) :: price(size(K))
      price = smileformula(t, v, K, kappa, kappa*theta, nu, rho)
    end function evalunique_wrap
  end function formula_normalised

  function formula(t, S, v, K, r, kappa, theta, nu, rho) result(price)
    real(dp), intent(in) :: t(:), S, v, K(:), r, kappa, theta, nu, rho
    real(dp) :: price(size(K))

    price = formula_normalised(t, v, K/(exp(r*t)*S), kappa, theta, nu, rho)
    price = price*S
  end function formula

  function calibrate(prices, t, K, S, r) result(params)
    integer, parameter :: p = 5
    real(dp), parameter :: tol = 0.001_dp
    real(dp), intent(in) :: prices(:), t(:), K(:), S, r
    integer :: info, ipvt(p)
    real(dp) :: wa((p + 1)*size(prices) + p*5)
    real(dp) :: pricestilde(size(prices)), Ktilde(size(prices))
    real(dp) :: params_cache(p), pricegrad(size(prices), 6)
    real(dp) :: params(p), fvec(size(prices))

    params = [0.2_dp,        &! v0
              1.2_dp,        &! kappa
              1.2_dp*0.2_dp, &! a=kappa*theta
              0.3_dp,        &! nu
             -0.6_dp]         ! rho

    pricestilde = prices/S
    Ktilde = K/(exp(r*t)*S)

    call dnls1e(lmderwrap, 2, size(prices), p, params, fvec, tol, 0, info, &
                ipvt, wa, size(wa))
    call positivenu(params(4), params(5))
    params(3) = params(3)/params(2)

  contains
    subroutine lmderwrap(iflag, m, p, params, fvec, fjac, ldfjac)
      integer, intent(in) :: iflag, m, p, ldfjac
      real(dp), intent(in) :: params(p)
      real(dp), intent(out) :: fvec(m), fjac(m, p)

      call unused(ldfjac)

      if (any(params /= params_cache)) then
         pricegrad = pricegradfunc(t, params(1), Ktilde, params(2), params(3), &
                                   params(4), params(5))
         params_cache = params
      end if

      if (iflag == 1) then
         fvec = pricegrad(:, 1) - pricestilde
      else
         fjac = pricegrad(:, 2:)
      end if
    end subroutine lmderwrap
  end function calibrate

  subroutine calibrate_ts(params, prices, t, K, S, r)
    real(dp), parameter :: tol = 0.001_dp
    real(dp), intent(in) :: prices(:, :), t(:), K(:), S(:), r
    real(dp), intent(out) :: params(:)
    integer :: d, m, n, p, info, ipvt(size(S) + 4)
    real(dp) :: wa((size(S) + 4)*(size(prices) + 5) + size(prices))
    real(dp) :: params_cache(size(S) + 4), pricegrad(size(K), size(S), 6)
    real(dp) :: fvec(size(prices))

    d = size(K)
    n = size(S)
    p = n + 4
    m = size(prices)

    params(:n)     =  0.2_dp           ! v0
    params(n + 1:) = [1.2_dp, &        ! kappa
                      1.2_dp*0.2_dp, & ! a = kappa*theta
                      0.3_dp, &        ! nu
                     -0.6_dp]          ! rho

    call dnls1e(lmderwrap, 2, m, p, params, fvec, tol, 0, info, ipvt, wa, &
                size(wa))
    call positivenu(params(n + 3), params(n + 4))
    params(n + 2) = params(n + 2)/params(n + 1)

  contains
    subroutine lmderwrap(iflag, m, p, params, fvec, fjac, ldfjac)
      integer, intent(in) :: iflag, m, p, ldfjac
      real(dp), intent(in) :: params(p)
      real(dp), intent(out) :: fvec(m), fjac(m, p)
      real(dp) :: fjac_r(d, n, p)
      integer :: i

      call unused(ldfjac)

      if (any(params /= params_cache)) then
         do i = 1, n
            pricegrad(:, i, :) = &
                 pricegradfunc(t, params(i), K/(exp(r*t)*S(i)), params(n + 1), &
                               params(n + 2), params(n + 3), params(n + 4))
            pricegrad(:, i, 1) = pricegrad(:, i, 1) - prices(:, i)/S(i)
        end do
        params_cache = params
      end if

      if (iflag == 1) then
        fvec = reshape(pricegrad(:, :, 1), [m])
      else
        fjac = 0
        fjac_r = reshape(fjac, [d, n, p])
        do i = 1, n
           fjac_r(:, i, n + 1:) = pricegrad(:, i, 3:)
           fjac_r(:, i,      i) = pricegrad(:, i, 2)
        end do
        fjac = reshape(fjac_r, [m, p])
      end if
    end subroutine lmderwrap
  end subroutine calibrate_ts

  function calibrate_vol(prices, t, K, S, r, kappa, theta, nu, rho) result(v)
    integer, parameter :: p = 1
    real(dp), parameter :: tol = 0.001_dp
    real(dp), intent(in) :: prices(:), t(:), K(:), S, r, kappa, theta, nu, rho
    integer :: info, ipvt(p)
    real(dp) :: wa((p + 1)*size(prices) + p*5)
    real(dp) :: pricestilde(size(prices)), Ktilde(size(prices))
    real(dp) :: v, v_vec(p), fvec(size(prices))

    v_vec = 0.2_dp
    pricestilde = prices/S
    Ktilde = K/(exp(r*t)*S)

    call dnls1e(lmderwrap, 2, size(prices), p, v_vec, fvec, tol, 0, info, &
                ipvt, wa, size(wa))
    v = v_vec(1)

  contains
    subroutine lmderwrap(iflag, m, p, v_vec, fvec, fjac, ldfjac)
      integer, intent(in) :: iflag, m, p, ldfjac
      real(dp), intent(in) :: v_vec(p)
      real(dp), intent(out) :: fvec(m), fjac(m, p)
      integer :: i

      call unused(ldfjac)

      if (iflag == 1) then
        fvec = formula_normalised(t, v_vec(1), Ktilde, kappa, theta, nu, rho) &
               - pricestilde
      else
        do i = 1, m
           fjac(i, 1) = vegafunc_normalised(t(i), v_vec(1), Ktilde(i), kappa, &
                                            theta, nu, rho)
        end do
      end if
    end subroutine lmderwrap
  end function calibrate_vol

  subroutine positivenu(nu, rho)
    real(dp), intent(inout) :: nu, rho
    if (nu < 0) then
       nu = -nu
       rho = -rho
    end if
  end subroutine positivenu

end module heston
