module grd_conn_m
  use grd_read_m
  implicit none
  contains

  ! link faces to associated volumes
  subroutine grd_conn(igrd)
    integer,intent(in)::igrd
    integer::i,j,ip,in,ivf,ivrt,ivol,ifce,ifce1,ifce2
    integer::n_fce,n_vrt,n_vrt_fce,ivol_ty,n_mtch

    type::fce_lst2_t
      integer::n_vrt,ip=0,in=0
      integer,dimension(:),allocatable::vrt_lst
      logical::mtchd=.false.,copy=.false.
    end type
    integer::n_fce_lst=0
    type(fce_lst2_t),dimension(:),allocatable::fce_lst

    type(grd_t),pointer::grd_p
    type(vol_t),pointer::vol_p,vol_n
    type(fce_t),pointer::fce_p
    type(vrt_t),pointer::vrt_p
    type(vol_ty_t),pointer::vol_ty_p
    type(fce_lst_t),pointer::fce_lst_p


    print*,'finding neighbour volumes...'
    grd_p=>grd(igrd)

    do ivol=1,grd(igrd)%n_vol
      vol_p=>grd_p%vol(ivol)
      n_fce_lst=n_fce_lst+vol_p%n_fce
    end do
    allocate(fce_lst(n_fce_lst))

!     if(n_dx==2)n_vrt_fce=8
!     if(n_dx==3)n_vrt_fce=16
!     do ivrt=1,grd_p%n_vrt
!       vrt_p=>grd_p%vrt(ivrt)
!       allocate(vrt_p%fce_lst(n_vrt_fce))
!       vrt_p%fce_lst=0
!     end do

    ifce=0
    do ivol=1,grd_p%n_vol
      vol_p=>grd_p%vol(ivol)

      do ivf=1,vol_p%n_fce
        ifce=ifce+1
        fce_lst(ifce)%ip=ivol
        ivol_ty=vol_p%v_ty
        n_vrt=vol_ty(ivol_ty)%fce_lst(ivf)%n_vrt
        fce_lst(ifce)%n_vrt=n_vrt
        allocate(fce_lst(ifce)%vrt_lst(n_vrt))

        do i=1,n_vrt
          ivrt=vol_ty(ivol_ty)%fce_lst(ivf)%vrt_lst(i)
          fce_lst(ifce)%vrt_lst(i)=vol_p%vrt_lst(ivrt)

!           ivrt=fce_lst(ifce)%vrt_lst(i)
!           do j=1,n_vrt_fce
!             vrt_p=>grd_p%vrt(ivrt)
!             if(vrt_p%fce_lst(j)==0)then
!               vrt_p%fce_lst(j)=ifce
!               vrt_p%n_fce=vrt_p%n_fce+1
!               exit
!             end if
!           end do

        end do

      end do
    end do


    ! find matching faces (slow)
    n_mtch=0
    do ifce1=1,n_fce_lst
      if(fce_lst(ifce1)%mtchd)cycle
      do ifce2=1,n_fce_lst
        if(ifce1==ifce2)cycle
        if(fce_lst(ifce2)%mtchd)cycle
        if(match_vrt(fce_lst(ifce1)%vrt_lst,fce_lst(ifce2)%vrt_lst))then
          fce_lst(ifce1)%in=fce_lst(ifce2)%ip
          fce_lst(ifce1)%mtchd=.true.
          fce_lst(ifce2)%mtchd=.true.
          fce_lst(ifce2)%copy=.true.
          n_mtch=n_mtch+1
          exit
        end if
      end do
    end do

    print*,'t. no. faces:',n_fce_lst
    print*,'no.  matches:',n_mtch


    ! set pole and neighbour indices
    ifce2=0
    do ifce=1,n_fce_lst
      if(.not.fce_lst(ifce)%copy)then
        ifce2=ifce2+1
        fce_p=>grd_p%fce(ifce2)
        fce_p%n_vrt=fce_lst(ifce)%n_vrt
        allocate(fce_p%vrt_lst(fce_p%n_vrt))
        fce_p%vrt_lst=fce_lst(ifce)%vrt_lst
        fce_p%ip=fce_lst(ifce)%ip
        if(fce_lst(ifce)%in>0)then
          fce_p%in=fce_lst(ifce)%in
        end if
      end if
    end do


    ! deallocate arrays
!     do ivrt=1,grd_p%n_vrt
!       vrt_p=>grd_p%vrt(ivrt)
!       deallocate(vrt_p%fce_lst)
!     end do

    do ifce=1,n_fce_lst
      deallocate(fce_lst(ifce)%vrt_lst)
    end do
    deallocate(fce_lst)


    ! create volume face lists
    do ifce=1,grd_p%n_fce
      fce_p=>grd_p%fce(ifce)
      ip=fce_p%ip; in=fce_p%in
      vol_p=>grd_p%vol(ip)
      n_fce=vol_p%n_fce
      do ivf=1,n_fce
        ifce1=vol_p%fce_lst(ivf)
        if(ifce1==0)then
          vol_p%fce_lst(ivf)=ifce
          exit
        end if
      end do
      if(in>0)then
        vol_n=>grd_p%vol(in)
        n_fce=vol_n%n_fce
        do ivf=1,n_fce
          ifce1=vol_n%fce_lst(ivf)
          if(ifce1==0)then
            vol_n%fce_lst(ivf)=ifce
            exit
          end if
        end do
      end if
    end do

    call c_vfce_lst(igrd)
    call c_neig_lst(igrd)
  end subroutine



  ! compare two lists of verticies
  function match_vrt(vrt_lst1,vrt_lst2) result(match)
    integer,dimension(:),intent(in)::vrt_lst1,vrt_lst2
    integer::idx1,idx2,n_vrt1,n_vrt2,sum_vrt1,sum_vrt2,n_vrt
    logical::match
    match=.false.
    n_vrt1=size(vrt_lst1)
    n_vrt2=size(vrt_lst2)
    sum_vrt1=sum(vrt_lst1)
    sum_vrt2=sum(vrt_lst2)
    if(n_vrt1==n_vrt2.and.sum_vrt1==sum_vrt2)then
      n_vrt=0
      do idx1=1,n_vrt1
        do idx2=1,n_vrt1
          if(vrt_lst1(idx1)==vrt_lst2(idx2))then
            n_vrt=n_vrt+1; exit
          end if
        end do
      end do
      if(n_vrt==n_vrt1)match=.true.
    end if
  end function



  ! construct volume face lists
  subroutine c_vfce_lst(igrd)
    integer,intent(in)::igrd
    integer::ivol,ip,in,ivf,ifce
    type(grd_t),pointer::grd_p
    type(vol_t),pointer::vol_p,vol_n
    type(fce_t),pointer::fce_p

    grd_p=>grd(igrd)
    do ifce=1,grd_p%n_fce
      fce_p=>grd_p%fce(ifce)
      ip=fce_p%ip
      vol_p=>grd_p%vol(ip)
      do ivf=1,vol_p%n_fce
        if(vol_p%fce_lst(ivf)==0)then
          vol_p%fce_lst(ivf)=ifce
        end if
      end do
      if(fce_p%in>0)then
        in=fce_p%in
        vol_n=>grd_p%vol(in)
        do ivf=1,vol_n%n_fce
          if(vol_n%fce_lst(ivf)==0)then
            vol_n%fce_lst(ivf)=ifce
          end if
        end do
      end if
    end do
  end subroutine



  ! construct volume neighbour lists
  subroutine c_neig_lst(igrd)
    integer,intent(in)::igrd
    integer::ivol,ip,in,ivf,ifce
    type(grd_t),pointer::grd_p
    type(vol_t),pointer::vol_p,vol_n
    type(fce_t),pointer::fce_p

    grd_p=>grd(igrd)
    do ivol=1,grd_p%n_vol
      vol_p=>grd_p%vol(ivol)
      do ivf=1,vol_p%n_fce
        ip=grd_p%fce(vol_p%fce_lst(ivf))%ip
        in=grd_p%fce(vol_p%fce_lst(ivf))%in
        if(ip>0.and.ip/=ivol)then
          vol_p%ng_lst(ivf)=ip
        else if(in>0.and.in/=ivol)then
          vol_p%ng_lst(ivf)=in
        end if
      end do
    end do
  end subroutine
end module
