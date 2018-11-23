module numerics
  use iso_c_binding, only: c_double, c_int, c_ptr, c_loc, c_funptr, c_funloc, &
                           c_f_pointer

  implicit none

  private
  public dp, unused, quad, quadvec, evalunique, matevalunique, newtonstep

  integer, parameter :: dp = kind(0.d0) ! double precision

  type opaque
  end type opaque

  interface
     function c_memoise_project(x, k, f, vecsize, cache) bind(c) result(y)
       import c_double, c_funptr, c_int, c_ptr
       real(c_double), value :: x
       integer(c_int), value :: k, vecsize
       type(c_funptr), value :: f
       type(c_ptr), value :: cache
       real(c_double) :: y
     end function c_memoise_project

     function c_build_cache() bind(c) result(cache)
       import c_ptr
       type(c_ptr) :: cache
     end function c_build_cache

     subroutine c_destroy_cache(cache) bind(c)
       import c_ptr
       type(c_ptr), value :: cache
     end subroutine c_destroy_cache

     subroutine c_vecfunc(x, y) bind(c)
       import c_double
       real(c_double), value :: x
       real(c_double), intent(out) :: y(:)
     end subroutine c_vecfunc

     function matrealparam(x, p, d) result(y)
       import dp
       integer, intent(in) :: d
       real(dp), intent(in) :: x(:), p
       real(dp) :: y(size(x), d)
     end function matrealparam

     function vecrealparam(x, p) result(y)
       import dp
       real(dp), intent(in) :: x(:), p
       real(dp) :: y(size(x))
     end function vecrealparam

     subroutine vecrealsub(x, y)
       import dp
       real(dp), intent(in) :: x
       real(dp), intent(out) :: y(:)
     end subroutine vecrealsub
  end interface

contains
  subroutine unused(var)
    integer, intent(in) :: var
  end subroutine unused

  subroutine buildcache(cache)
    type(opaque), pointer :: cache
    type(c_ptr) :: c_cache
    c_cache = c_build_cache()
    call c_f_pointer(c_cache, cache)
  end subroutine buildcache

  subroutine destroycache(cache)
    type(opaque), pointer :: cache
    call c_destroy_cache(c_loc(cache))
  end subroutine destroycache

  function memoiseproject(x, k, vecsub, vecsize, cache) result(y)
    real(dp) :: x
    integer :: k, vecsize
    procedure(vecrealsub) :: vecsub
    type(opaque), pointer :: cache
    real(dp) :: y

    y = c_memoise_project(x, k - 1, c_funloc(c_vecfun), vecsize, c_loc(cache))
  contains
    subroutine c_vecfun(x, y)
      real(dp), value :: x
      real(dp), intent(out) :: y(vecsize)
      call vecsub(x, y)
    end subroutine c_vecfun
  end function memoiseproject

  function quad(func) result(integral)
    real(dp), parameter :: ulim = 200
    real(dp), parameter :: epsabs = 0.0001
    real(dp), parameter :: epsrel = 0.0001
    integer, parameter :: key = 2
    integer, parameter :: limit = 1000
    integer, parameter :: lenw = 4*limit

    real(dp), external :: func
    real(dp) :: integral, abserr, last, work(4*limit)
    integer :: neval, ier, iwork(lenw)

    call dqag(func, 0, ulim, epsabs, epsrel, key, integral, abserr, neval, &
              ier, limit, lenw, last, iwork, work)
  end function quad

  function quadvec_core(f, k, n, cache) result(y)
    procedure(vecrealsub) :: f
    integer, intent(in) :: k, n
    type(opaque), pointer :: cache
    real(dp) :: y

    y = quad(integrandclosure)
  contains
    function integrandclosure(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      y = memoiseproject(x, k, f, n, cache)
    end function integrandclosure
  end function quadvec_core

  subroutine quadvec(f, y)
    procedure(vecrealsub) :: f
    real(dp), intent(out) :: y(:)
    integer :: k
    type(opaque), pointer :: cache

    call buildcache(cache)
    do k=1,size(y)
       y(k) = quadvec_core(f, k, size(y), cache)
    end do
    call destroycache(cache)
  end subroutine quadvec

  recursive function evalunique(f, x, p) result(y)
    procedure(vecrealparam) :: f
    real(dp), intent(in) :: x(:), p(:)
    real(dp) :: y(size(x))
    logical :: mask(size(x))

    if (size(x) > 0) then
       mask = abs(p - p(1)) > tiny(1._dp)
       y = unpack(evalunique(f, pack(x, mask), pack(p, mask)), mask, &
                  unpack(f(pack(x, .not. mask), p(1)), .not. mask, 0._dp))
    end if
  end function evalunique

  recursive function matevalunique(f, x, p, d) result(y)
    procedure(matrealparam) :: f
    integer, intent(in) :: d
    real(dp), intent(in) :: x(:), p(:)
    real(dp) :: y(size(x), d)
    logical :: mask(size(x))
    real(dp), allocatable :: y_rec(:, :), y_step(:, :)
    integer :: i

    if (size(x) > 0) then
       mask = abs(p - p(1)) > tiny(1._dp)
       y_rec = matevalunique(f, pack(x, mask), pack(p, mask), d)
       y_step = f(pack(x, .not. mask), p(1), d)
       do i = 1, d
          y(:, i) = unpack(y_rec(:, i), mask, &
               unpack(y_step(:, i), .not. mask, 0._dp))
       end do
    end if
  end function matevalunique
  function newtonstep(x, y, dy) result(xprime)
    real(dp), intent(in) :: x, y, dy
    real(dp) :: xprime

    xprime = x - y/dy
  end function newtonstep
end module numerics
