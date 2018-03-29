module const_m
implicit none


! computational
integer,parameter::rp=kind(1.d0)
real(rp),parameter::zero=0,one=1,ten=10
real(rp),parameter::half=one/2,third=one/3
real(rp),parameter::fifth=one/5,sevth=one/7
real(rp),parameter::small=ten**(-20)
real(rp),parameter::large=one/small
! character(256)::input_path=""
! character(256)::output_path=""


! mathematical
real(rp),parameter::pi=3.14159265
real(rp),parameter::deg2rad=pi/180
real(rp),parameter::rad2deg=1/deg2rad


! physical
real(rp),parameter::c_gas=8.314472
real(rp),parameter::c_grv=-9.81


contains

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

end module
