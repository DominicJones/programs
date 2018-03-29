module grd_math_m
  use grd_data_m
  implicit none
  contains

  ! vector magnitude
  function vec_mag(vec)
    real(wp),dimension(:),intent(in)::vec
    real(wp)::vec_mag
    vec_mag=sqrt(sum(vec**2))
  end function

  ! unit vector
  function unit_vec(vec)
    real(wp),dimension(:),intent(in)::vec
    real(wp),dimension(size(vec))::unit_vec
    unit_vec=vec/(vec_mag(vec)+small)
  end function

  ! cross product
  function cross_prod(vec1,vec2)
    real(wp),dimension(3),intent(in)::vec1,vec2
    real(wp),dimension(3)::cross_prod
    integer,dimension(3),parameter::c1=(/2,3,1/),c2=(/3,1,2/)
    cross_prod=vec1(c1)*vec2(c2)-vec2(c1)*vec1(c2)
  end function
end module
