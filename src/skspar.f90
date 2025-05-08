function skspar(i, j, r2, dd) result(res)
  implicit none

  ! Arguments
  integer, intent(in) :: i, j
  real(8), intent(in) :: r2
  real(8), intent(out) :: dd(13)
  integer :: res

  ! Parameters
  include 'maxima.h'

  ! Local variables
  integer :: maxmax, minmax, inu, ind, in
  real(8) :: r, h, rch, rch3
  real(8) :: dx1, dx2, dx3, dx4
  real(8) :: fac14, fac23
  real(8) :: f1, f2, f3, f4

  ! Shared data structures
  integer :: dim(MAXTYP, MAXTYP)
  real(8) :: skhtab(10, MAXTAB, MAXTYP, MAXTYP)
  real(8) :: skstab(10, MAXTAB, MAXTYP, MAXTYP)
  real(8) :: skself(3, MAXTYP), dr(MAXTYP, MAXTYP), rcdr(MAXTYP, MAXTYP)
  real(8) :: rslako
  integer :: lmax(MAXTYP)

  ! Common blocks
  common /sktab/ skhtab, skstab, skself, dr, rcdr, dim
  common /range/ rslako
  common /lmax/ lmax

  ! Initialize return value
  res = 0

  ! Determine interaction type
  maxmax = max(lmax(i), lmax(j))
  minmax = min(lmax(i), lmax(j))

  if (maxmax <= 1) then
    inu = 10
  else if (maxmax <= 2) then
    if (minmax <= 1) then
      inu = 9
    else
      inu = 6
    end if
  else
    if (minmax <= 1) then
      inu = 8
    else if (minmax <= 2) then
      inu = 4
    else
      inu = 1
    end if
  end if

  ! Compute distance and index
  r = sqrt(r2)
  if (r > rslako) then
    ind = dim(i, j) + 1
  else
    ind = int(r * rcdr(i, j) + 1.0d0)
  end if

  if (r2 < 1.0d-8) then
    do in = 1, 3
      dd(in + 10) = 1.0d0
    end do
  else if (ind + 3 > dim(i, j)) then
    do in = inu, 10
      dd(in) = 0.0d0
    end do
  else
    h = dr(i, j)
    rch = rcdr(i, j)
    rch3 = rch * rch * rch

    dx1 = r - (ind - 1) * h
    dx2 = dx1 - h
    dx3 = dx2 - h
    dx4 = dx3 - h

    fac14 = dx2 * dx3 * rch3 / 6.0d0
    fac23 = dx1 * dx4 * rch3 / 2.0d0

    do in = inu, 10
      f1 = skstab(in, ind,     i, j)
      f2 = skstab(in, ind + 1, i, j)
      f3 = skstab(in, ind + 2, i, j)
      f4 = skstab(in, ind + 3, i, j)

      dd(in) = fac14 * (f4 * dx1 - f1 * dx4) + fac23 * (f2 * dx3 - f3 * dx2)
    end do
  end if

end function skspar
