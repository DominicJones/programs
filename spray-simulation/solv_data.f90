module solv_data_m
  use precision_m
  implicit none

  ! list of faces and neigbours for each volume
  type::v_lst_t
    integer,dimension(:),pointer::fce_lst=>null()
    integer,dimension(:),pointer::ng_lst=>null()
  end type
  type(v_lst_t),dimension(:),allocatable::v_lst


  ! neigbour indicies and coefficients (off-diagonal)
  type::ng_t
    integer::ip=0,in=0
    real(wp)::p=0.,n=0.
  end type
  type(ng_t),dimension(:),allocatable::ng


  ! pole coefficients (diagonal)
  real(wp),dimension(:),allocatable::pl


  ! phi type (locally a vector)
  type::phi_t
    real(wp),dimension(:),pointer::vc=>null()
  end type
end module
