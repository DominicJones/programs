module mesh_data_m
use list_m
implicit none

type::msh_t
  character(50)::path,fnm

  integer::key=0,n_dx=0,ty=0,ph_ty=0
  integer::n_ja=0,n_ia=0

  integer::n_vrt=0,n_fce=0,n_elm=0
  integer::n_bvrt=0,n_bedg=0,n_bfce=0,n_ebf=0

  integer,dimension(:),allocatable::part
  integer,dimension(:),allocatable::elm_idx

  integer,dimension(:),allocatable::vrt_ja
  integer,dimension(:),allocatable::vrt_ia

  integer,dimension(:),allocatable::fce_ja
  integer,dimension(:),allocatable::fce_ia

  integer,dimension(:),allocatable::elm_ja
  integer,dimension(:),allocatable::elm_ia

  type(lst_t),dimension(:),allocatable::vrt_tag
  type(lst_t),dimension(:),allocatable::vrt_edg
  type(lst_t),dimension(:),allocatable::vrt_fce
  type(lst_t),dimension(:),allocatable::vrt_elm

  type(lst_t),dimension(:),allocatable::edg_tag
  type(lst_t),dimension(:),allocatable::edg_vrt
  type(lst_t),dimension(:),allocatable::edg_fce
  type(lst_t),dimension(:),allocatable::edg_elm

  type(lst_t),dimension(:),allocatable::fce_tag
  type(lst_t),dimension(:),allocatable::fce_vrt
  type(lst_t),dimension(:),allocatable::fce_edg
  type(lst_t),dimension(:),allocatable::fce_elm

  type(lst_t),dimension(:),allocatable::elm_tag
  type(lst_t),dimension(:),allocatable::elm_vrt
  type(lst_t),dimension(:),allocatable::elm_edg
  type(lst_t),dimension(:),allocatable::elm_fce
end type

end module
