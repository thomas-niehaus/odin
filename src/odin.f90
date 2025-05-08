program odin

  include 'maxima.h'
  !!
  integer :: dim(MAXTYP,MAXTYP), lmax(MAXTYP)
  !! The order of basis functions in the slater-koster files changed 
  integer :: icv(9) = (/ 1, 4, 2, 3, 5, 6, 8, 9, 7 /) 
  integer, allocatable :: izp(:), ind(:)
  integer :: nn, izpk, izpj, k, j, ndim
  integer :: skspar
  character(len=128) :: infile, strForm, prefix, sep, suffix, line
  character(len=128) :: skfn(MAXTYP,MAXTYP)
  character(len=2) :: tmpchr, atomsym(MAXTYP)
  real(kind(0.0d0)), allocatable :: x(:,:), overl(:,:)
  real(kind(0.0d0)) :: skhtab(10,MAXTAB,MAXTYP,MAXTYP), skstab(10,MAXTAB,MAXTYP,MAXTYP)
  real(kind(0.0d0)) :: skself(3,MAXTYP), dr(MAXTYP,MAXTYP), rcdr(MAXTYP,MAXTYP)
  real(kind(0.0d0)) :: xhlp(3), dif(3), au(LDIM,LDIM)
  common /sktab/ skhtab,skstab,skself,dr,rcdr,dim
  common /lmax/ lmax
  external gettab, slkode, skspar

  print *,'** odin (t.niehaus, based on dylcao by seifert,porezag,blaudeck) **'
  print *,'Creates dftb overlap matrix and stores the result in oversqr.dat'
  print *,'** Version 01  02.05.2025 **'
  print *,'enter filename for input structure'
  read *,infile
  print *,'infile :', trim(infile)

  open (1,file=infile,status='old'); rewind 1
  read (1,*) nn,tmpchr
  allocate(x(3,nn))
  allocate(izp(nn))
  allocate(ind(nn+1))

  read(1, '(A)') line
  call get_strings_from_line(line, atomsym, nel)

  do n = 1,nn 
     read(1,*) j,izp(n),(xhlp(j),j=1,3)
     do j=1,3 
        x(j,n)=xhlp(j)/0.529177
     enddo
  enddo
  close (1)

 !! Check if ntype in limits of MAXTYP 
  ntype = 1
  do n = 1,nn 
     if (izp(n) > ntype) ntype = izp(n)
  enddo
  if (ntype > MAXTYP) then
     print *,'odyn: ntype >', MAXTYP
     stop
  endif

  print *,'enter prefix for Slater-Koster path' 
  read *, prefix
  print *,'Prefix is ', trim(prefix)
  print *,'enter separator of Slater-Koster file names'
  read *, sep
  print *,'Separator is ', trim(sep)
  print *,'enter suffix of Slater-Koster file'
  read *, suffix
  print *,'Suffix is ', trim(suffix)

  do i = 1, ntype
     do j = 1, ntype
       skfn(i,j) = trim(prefix) // trim(atomsym(i)) // trim(sep) // trim(atomsym(j)) // trim(suffix)
     enddo
  enddo

  !! Read in maximal angular momentum for each type
  !! 1,2,3  for  s,p,d, respectively
  write (*,*) 'enter ',ntype,' * lmax'
  read *,(lmax(i),i = 1,ntype)
  call gettab(ntype, skfn)

  !! calculation of indices for matrices S 
  ind(1) = 0
  do j = 1,nn 
     izpj = izp(j)
     ind(j+1) = ind(j)+lmax(izpj)**2
  enddo 
  ndim = ind(nn+1) ! actual dimension of matrix
  allocate(overl(ndim,ndim))

  !! setup of overlap
  do j = 1,nn 
     izpj = izp(j)
     do k = 1,j 
        izpk = izp(k)
        do n = 1,ind(k+1)-ind(k) 
           do m = 1,ind(j+1)-ind(j) 
              overl(ind(j)+m,ind(k)+n) = 0.0d0 
           enddo
        enddo   
        do i = 1,3 
           dif(i) = x(i,k) - x(i,j) 
        enddo
        call slkode(dif,izpj,izpk,au,skspar)

        do n = 1,ind(k+1)-ind(k) 
           do m = 1,ind(j+1)-ind(j) 
              overl(ind(j)+icv(m),ind(k)+icv(n)) = overl(ind(j)+icv(m),ind(k)+icv(n)) + au(m,n)
           enddo
        enddo
        do n = 1,ind(k+1)-ind(k) 
           do m = 1,ind(j+1)-ind(j) 
              overl(ind(k)+icv(n),ind(j)+icv(m)) = overl(ind(j)+icv(m),ind(k)+icv(n))
           enddo
        enddo
     enddo
  enddo

  open(12, file='oversqr.dat')
  write(12, "(A1,A10,A10,A10,A10)") "#", "REAL", "NALLORB", "NKPOINT"
  write(12, "(1X,L10,I10,I10,I10)") .true., ndim, 1
  write (strForm, "(A,I0,A)") "(", ndim, "ES24.15)"
  write(12, "(A1,A10,A10)") "#", "IKPOINT"
  write(12, "(1X,I10,I10)") 1
  write(12, "(A1,A)") "#", " MATRIX"
  write(12, strForm) overl
  close(12)
  print *,'Results in oversqr.dat'
  print *,'Enjoy!'

  deallocate(overl)
  deallocate(izp)
  deallocate(ind)
  deallocate(x)
end program odin


subroutine get_strings_from_line(input_string, output_array, string_count)
  implicit none
  include 'maxima.h'
  character(len=*), intent(in) :: input_string
  character(len=*) :: output_array(MAXTYP)
  integer, intent(out) :: string_count
  integer :: start, end_pos, i, array_index
  logical :: in_string
  
  string_count = 0
  array_index = 1
  start = 1
  in_string = .false.
  
  do i = 1, len(input_string)
     if (input_string(i:i) /= ' ') then
        if (.not. in_string) then
           start = i
           in_string = .true.
        end if
     else
        if (in_string) then
           end_pos = i - 1
           if (array_index <= size(output_array)) then
              output_array(array_index) = trim(adjustl(input_string(start:end_pos)))
              array_index = array_index + 1
              string_count = string_count + 1
           else
              exit
           end if
           in_string = .false.
        end if
     end if
  end do
  
  ! Handle the last string
  if (in_string) then
     end_pos = len(input_string)
     if (array_index <= size(output_array)) then
        output_array(array_index) = trim(adjustl(input_string(start:end_pos)))
        string_count = string_count + 1
     end if
  end if

end subroutine get_strings_from_line



