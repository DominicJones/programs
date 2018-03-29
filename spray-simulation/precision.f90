module precision_m
  integer,parameter::wp=selected_real_kind(10,30)
  real(wp),parameter::zero=0._wp
  real(wp),parameter::one=1._wp
  real(wp),parameter::ten=10._wp
  real(wp),parameter::pi=4._wp*atan(one)
  real(wp),parameter::rad=pi/180._wp
  real(wp),parameter::small=ten**(-20)
  real(wp),parameter::large=one/small
end module
