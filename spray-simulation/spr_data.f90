module spr_data_m
  use precision_m
  implicit none
  integer,dimension(:),allocatable::mom_pow
  integer::mp1=1,mpr=3,pdf_mth=5,vdf_mth=1
  integer::inj_st=1,reb_st=2,spl_st=3

  logical::hyd=.true.,imp=.true.
  logical::bre=.true.,col=.true.,drg=.true.
  integer::bre_mod=1
  real(wp)::vel_exp=0.4,n_spl=2.

  real(wp)::r_orif=0.00025,sheet_th=0.
  real(wp)::inj_dur=1.
  real(wp)::vof_inj=0.6,r32_inj=20.e-6,skw_inj=7.
  real(wp)::spd_inj=100.,ang_inj=10.*(pi/180.)

  real(wp)::vof_wall=0.01
  real(wp),dimension(1)::bre_we=(/6./)
  real(wp),dimension(1)::bre_ts=(/20./)
  real(wp),dimension(2)::wall_we=(/1320.,2634./)

  contains
  function mom_id(imom)
    integer,intent(in)::imom
    integer::mom_id,id

    mom_id=0
    do id=1,size(mom_pow)
      if(mom_pow(id)==imom)then
        mom_id=id
        exit
      end if
    end do
  end function
end module
