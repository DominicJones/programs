module grd_solv_m
  use grd_data_m
  use solv_data_m
  implicit none
  contains

  subroutine grd_solv(igrd)
    integer,intent(in)::igrd
    integer::ip,in,ivf,ifce

    allocate(v_lst(grd(igrd)%n_vol))

    allocate(pl(grd(igrd)%n_vol))
    do ip=1,size(pl)
      v_lst(ip)%fce_lst=>grd(igrd)%vol(ip)%fce_lst
      v_lst(ip)%ng_lst=>grd(igrd)%vol(ip)%ng_lst
    end do

    allocate(ng(grd(igrd)%n_fce))
    do ifce=1,size(ng)
      ng(ifce)%ip=grd(igrd)%fce(ifce)%ip
      ng(ifce)%in=grd(igrd)%fce(ifce)%in
    end do
  end subroutine
end module
