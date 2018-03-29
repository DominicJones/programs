module constants_m
! use prec_m, sp => rsp, dp => rdp, rp => rsp, mpi_rp => mpi_rsp
use prec_m, sp => rsp, dp => rdp, rp => rdp, mpi_rp => mpi_rdp
implicit none


! computational
integer,parameter::sl=256
integer::proc=0
integer::n_proc=1
type wrk_t
  integer,dimension(:),allocatable::i
  real(sp),dimension(:),allocatable::s
  real(dp),dimension(:),allocatable::d
end type
type(wrk_t),dimension(:),allocatable::wrk


! numerical
real(rp),parameter::zero=0
real(rp),parameter::one=1
real(rp),parameter::ten=10

real(rp),parameter::half=one/2
real(rp),parameter::third=one/3
real(rp),parameter::fifth=one/5
real(rp),parameter::sevth=one/7
real(rp),parameter::tenth=one/10

real(rp),parameter::small=ten**(-20)
real(rp),parameter::large=one/small


! mathematical
real(rp),parameter::pi=3.14159265


! physical
real(rp),parameter::c_gas=8.314472
real(rp),parameter::c_grv=-9.81

contains




function rad_to_deg(rad) result(deg)
  real(rp),intent(in)::rad
  real(rp)::deg
  deg = (180/pi)*rad
end function




function deg_to_rad(deg) result(rad)
  real(rp),intent(in)::deg
  real(rp)::rad
  rad = (pi/180)*deg
end function



function fix_prec(x,p) result(y)
  real(rp),intent(in)::x
  integer,intent(in)::p
  real(rp)::y

  character(20)::ystr
  character(10)::pfmt

  pfmt = "(g20."//trim(itoa(p))//")"

  write(ystr,pfmt) x
  read(ystr,pfmt) y
end function




function atoi(a) result(i)
  character(*)::a
  integer::i
  read(a,"(i10)") i
end function




function itoa(i) result(a)
  integer::i
  character(10)::a
  write(a,"(i10)") i
  a=adjustl(trim(a))
end function




function dir (filename)
  character(*),intent(in)::filename
  character(sl)::dir
  integer::n,i,j

  n=len_trim(filename); j=0
  do i=n,1,-1
    if(filename(i:i)=="/")then
      j=i; exit
    end if
  end do
  dir=filename(:j)
end function




function notdir (filename)
  character(*),intent(in)::filename
  character(sl)::notdir
  integer::n,i,j

  n=len_trim(filename); j=0
  do i=n,1,-1
    if(filename(i:i)=="/")then
      j=i; exit
    end if
  end do
  notdir=filename(j+1:n)
end function




function file_handle(rel) result(ios)
  integer,optional::rel
  integer::i,ios

  integer,parameter::n_unit=5
  integer,dimension(n_unit),parameter::unit=(/10,20,30,40,50/)
  logical,dimension(n_unit),save::in_use=.false.

  ios=-1

  if(present(rel))then
    do i=1,n_unit
      select case(rel-unit(i))
      case(0)
        in_use(i)=.false.
        return
      end select
    end do
  else
    do i=1,n_unit
      select case(in_use(i))
      case(.false.)
        ios=unit(i)
        in_use(i)=.true.
        return
      end select
    end do
  end if
end function
end module
