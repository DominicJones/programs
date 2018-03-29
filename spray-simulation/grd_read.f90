module grd_read_m
  use grd_constr_m
  implicit none
  contains



  ! read in Gambit neutral grid and store data
  subroutine grd_read(igrd,file_nm)
    integer,intent(in)::igrd
    character(len=*),intent(in)::file_nm
    character(len=32)::str,bnd_nm
    integer::ivrt,ivol,ifce,ivf,idx,ibs=0,i,n,ty
    integer::iunit=100
    type(grd_t),pointer::grd_p
    type(vol_t),pointer::vol_p
    type(fce_t),pointer::fce_p
    type(vrt_t),pointer::vrt_p

    grd_p=>grd(igrd)
    print*,'opening: ',trim(file_nm)
    open(iunit,file=trim(file_nm),status='old')
    do
      read(iunit,*)str


      if(trim(str)=='NUMNP')then
        read(iunit,'(6(1x,i9))') &
          grd_p%n_vrt,grd_p%n_vol, &
          grd_p%n_grp,grd_p%n_bs, &
          n_dx,grd_p%n_vc

        allocate(grd_p%vrt(grd_p%n_vrt))
        allocate(grd_p%vol(grd_p%n_vol))
      end if


      if(trim(str)=='NODAL')then
        ivrt=0
        do
          ivrt=ivrt+1
          read(iunit,'(i10,1x,3e20.12)')idx
          vrt_p=>grd(igrd)%vrt(idx)

          allocate(vrt_p%x(n_dx))

          backspace(iunit)
          read(iunit,'(i10,1x,3e20.12)')idx,vrt_p%x
          vrt_p%x=vrt_p%x*l_sc
          if(ivrt==grd_p%n_vrt)exit
        end do
      end if


      if(trim(str)=='ELEMENTS')then
        ivol=0
        do
          ivol=ivol+1
          read(iunit,'(i8,1x,i2,1x,i2,1x,7i8:/(15x,7i8:))')idx, &
            grd(igrd)%vol(idx)%v_ty,grd(igrd)%vol(idx)%n_vrt

          vol_p=>grd_p%vol(idx)
          allocate(vol_p%vrt_lst(vol_p%n_vrt)); vol_p%vrt_lst=0

          backspace(iunit)
          read(iunit,'(i8,1x,i2,1x,i2,1x,7i8:/(15x,7i8:))')idx, &
            vol_p%v_ty,vol_p%n_vrt,vol_p%vrt_lst

          vol_p%n_fce=vol_ty(vol_p%v_ty)%n_fce
          allocate(vol_p%fce_lst(vol_p%n_fce)); vol_p%fce_lst=0
          allocate(vol_p%bnd_lst(vol_p%n_fce)); vol_p%bnd_lst=0
          allocate(vol_p%ng_lst(vol_p%n_fce)); vol_p%ng_lst=0

          grd_p%n_fce=grd_p%n_fce+vol_p%n_fce

          if(ivol==grd_p%n_vol)exit
        end do
      end if


      if(trim(str)=='BOUNDARY')then
        ibs=ibs+1
        read(iunit,'(a32,8i10)')bnd_nm,ty,n
        do i=1,n
          read(iunit,'(i10,3i5)')ivol,ty,ivf
          vol_p=>grd_p%vol(ivol)
          vol_p%bnd_v=.true.
          vol_p%bnd_lst(ivf)=bnd_ty(bnd_nm)
          vol_p%n_bnd=vol_p%n_bnd+1
        end do
      end if


      if(ibs==grd(igrd)%n_bs)then
        do ivol=1,grd_p%n_vol
          vol_p=>grd_p%vol(ivol)
          grd_p%n_fce=grd_p%n_fce+vol_p%n_bnd
        end do
        grd_p%n_fce=grd_p%n_fce/2

        allocate(grd_p%fce(grd_p%n_fce))
        exit
      end if
    end do
    close(iunit)

    print*,'no. volumes: ',grd_p%n_vol
    print*,'no. faces:   ',grd_p%n_fce
    print*,'no. vertices:',grd_p%n_vrt
  end subroutine



  ! convert boundary name to boundary index
  function bnd_ty(bnd_nm) result(bt)
    character(len=*),intent(in)::bnd_nm
    character(len=len(bnd_nm))::nm
    integer::i,bt

    nm=trim(adjustl(bnd_nm))
    i=index(nm,'_')+1
    read(nm(i:i),'(i1)')bt
    if(nm(1:5)=='inlet')bt=bt+10
    if(nm(1:6)=='outlet')bt=bt+20
    if(nm(1:4)=='wall')bt=bt+30
    if(nm(1:8)=='symmetry')bt=bt+40
    if(nm(1:8)=='pressure')bt=bt+50
  end function
end module
