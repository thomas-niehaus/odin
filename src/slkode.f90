subroutine slkode(dum, i, j, em, iovpar)
  implicit none
  include 'maxima.h'
  integer :: iovpar
  integer :: i, j, l, k, maxmax, minmax
  real(8) :: dum(3), em(LDIM, LDIM), dummy(LDIM, LDIM)
  real(8) :: x(6), x2(6), r2, r2i, ri
  integer :: lmax(MAXTYP)

  external :: iovpar
  ! external routines assumed defined elsewhere
  external :: skss, skpp, sksp, skdd, sksd, skpd, selfs, selfp, selfd

  ! common block (consider replacing with a module)
  common /lmax/ lmax

  r2 = 0.0d0
  do l = 1, 3
    x(l) = dum(l)
    x2(l) = x(l) * x(l)
    r2 = r2 + x2(l)
  end do

  if (r2 >= 1.0d-8) then
    r2i = 1.0d0 / r2
    ri = sqrt(r2i)
    do l = 1, 3
      x(l) = x(l) * ri
      x(l + 3) = x(l)
      x2(l) = x2(l) * r2i
      x2(l + 3) = x2(l)
    end do

    maxmax = max(lmax(i), lmax(j))
    minmax = min(lmax(i), lmax(j))

    ! s interaction
    call skss(x, x2, i, j, r2, iovpar, em, LDIM)
    if (maxmax <= 1) return

    ! p interaction
    if (minmax >= 2) then
      call skpp(x, x2, i, j, r2, iovpar, em(2,2), LDIM)
      call sksp(x, x2, i, j, r2, iovpar, em(1,2), em(2,1), LDIM)
      if (i /= j) then
        call sksp(x, x2, j, i, r2, iovpar, dummy, em(2,1), LDIM)
      end if
    else if (lmax(j) >= 2) then
      call sksp(x, x2, i, j, r2, iovpar, em(1,2), em(2,1), LDIM)
    else
      call sksp(x, x2, j, i, r2, iovpar, dummy, em(2,1), LDIM)
    end if

    if (maxmax <= 2) return

    ! d interaction
    if (minmax == 3) then
      call skdd(x, x2, i, j, r2, iovpar, em(5,5), LDIM)
      call sksd(x, x2, i, j, r2, iovpar, em(1,5), em(5,1), LDIM)
      call skpd(x, x2, i, j, r2, iovpar, em(2,5), em(5,2), LDIM)
      if (i /= j) then
        call sksd(x, x2, j, i, r2, iovpar, dummy, em(5,1), LDIM)
        call skpd(x, x2, j, i, r2, iovpar, dummy, em(5,2), LDIM)
      end if
    else if (lmax(i) == 1) then
      call sksd(x, x2, i, j, r2, iovpar, em(1,5), em(5,1), LDIM)
    else if (lmax(i) == 2) then
      call sksd(x, x2, i, j, r2, iovpar, em(1,5), em(5,1), LDIM)
      call skpd(x, x2, i, j, r2, iovpar, em(2,5), em(5,2), LDIM)
    else if (lmax(j) == 1) then
      call sksd(x, x2, j, i, r2, iovpar, dummy, em(5,1), LDIM)
    else
      call sksd(x, x2, j, i, r2, iovpar, dummy, em(5,1), LDIM)
      call skpd(x, x2, j, i, r2, iovpar, dummy, em(5,2), LDIM)
    end if

  else
    do k = 1, LDIM
      do l = 1, LDIM
        em(k, l) = 0.0d0
      end do
    end do

    call selfs(i, j, r2, iovpar, em, LDIM)
    if (lmax(i) <= 1) return
    call selfp(i, j, r2, iovpar, em(2,2), LDIM)
    if (lmax(i) <= 2) return
    call selfd(i, j, r2, iovpar, em(5,5), LDIM)
  end if

end subroutine slkode

