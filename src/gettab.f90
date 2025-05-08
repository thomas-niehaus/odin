subroutine gettab(ntype, skfn)
  implicit none
  include 'maxima.h'
  integer, intent(in) :: ntype
  real(kind(0.0d0)) :: skhtab(10,MAXTAB,MAXTYP,MAXTYP),skstab(10,MAXTAB,MAXTYP,MAXTYP)
  real(kind(0.0d0)) :: espin(MAXTYP)
  real(kind(0.0d0)) :: skself(3,MAXTYP),dr(MAXTYP,MAXTYP)
  real(kind(0.0d0)) :: coeff(6,MAXINT,MAXTYP,MAXTYP),xr(2,MAXINT,MAXTYP,MAXTYP)
  real(kind(0.0d0)) :: efkt(3,MAXTYP,MAXTYP),cutoff(MAXTYP,MAXTYP)
  real(kind(0.0d0)) :: qzero(3,MAXTYP),uhubb(3,MAXTYP),uhbuu(MAXTYP),uhbud(MAXTYP),dipol(3)
  real(kind(0.0d0)) :: rcdr(MAXTYP,MAXTYP), rslako
  integer :: numint(MAXTYP,MAXTYP)
  integer :: dim(MAXTYP,MAXTYP),ppp
  integer :: i, j, k, l, ios
  character(len=128) :: chdummy
  character(len=128) :: skfn(MAXTYP,MAXTYP)

  common /sktab/ skhtab,skstab,skself,dr,rcdr,dim
  common /mcharge/ qzero, uhubb, dipol, uhbuu, uhbud
  common /range/ rslako
  common /spltab/ coeff,xr,efkt,cutoff,numint

  rslako = 0.0d0

  do i = 1, ntype
     do j = 1, ntype
        open (unit=3, file=trim(skfn(i,j)), status='old', action='read') 
        rewind (unit=3)
        read (3, *) dr(i, j), dim(i, j)
        rslako = max(rslako, dr(i, j) * dim(i, j))
        if (dr(i, j) <= 0.0d0) then
           write (*, *) 'radial increment for SK-data must be > 0'
           stop
        end if
      rcdr(i, j) = 1.0d0 / dr(i, j)
      if (i == j) then
         read (3, *) (skself(l, i), l=1,3), espin(i), &
              (uhubb(l, i), l=1,3), (qzero(4-l, i), l=1,3)
         uhbuu(i) = uhubb(1,i)
         uhbud(i) = uhubb(2,i)
         uhubb(1,i)=uhubb(3,i)
         uhubb(2,i)=uhubb(3,i)
      end if

      do k = 1, dim(i, j)
         read (3, *) (skhtab(l, k, i, j), l=1,10), &
              (skstab(l, k, i, j), l=1,10)
      end do
      
      do
         read (3, '(A)', iostat=ios) chdummy
         if (ios /= 0) exit ! Exit loop on end-of-file or error
         if (trim(chdummy) == 'Spline') exit
      end do
      
      read (3, *) numint(i, j), cutoff(i, j)
      if (numint(i, j) > MAXINT) then
         write (*, *) 'Too many intervalls!'
         stop
      end if
      read (3, *) (efkt(ppp, i, j), ppp=1, 3)
      
      do ppp = 1, numint(i, j)
         if (ppp < numint(i, j)) then
            read (3, *) xr(1, ppp, i, j), xr(2, ppp, i, j), &
                 coeff(1, ppp, i, j), coeff(2, ppp, i, j), &
                 coeff(3, ppp, i, j), coeff(4, ppp, i, j)
         else
            read (3, *) xr(1, ppp, i, j), xr(2, ppp, i, j), &
                       coeff(1, ppp, i, j), coeff(2, ppp, i, j), &
                       coeff(3, ppp, i, j), coeff(4, ppp, i, j), &
                       coeff(5, ppp, i, j), coeff(6, ppp, i, j)
         end if
      end do
      
      if (abs(xr(2, numint(i, j), i, j) - cutoff(i, j)) > 1.0d-10) then ! Using tolerance for comparison
         write (*, *) 'Error in Datafile for pair ', i, j
         stop
      end if
      write (*, *) 'skfile for pair ', i, j, ' :', trim(skfn(i,j))
      close (unit=3)
   end do
end do

end subroutine gettab
