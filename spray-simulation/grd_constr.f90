module grd_constr_m
  use grd_math_m
  implicit none
  contains



  ! define element types
  subroutine grd_prepr(file_nm)
    character(len=*),intent(in)::file_nm
    character(len=32)::str
    integer::iunit=100
    integer::ivol_ty,ifce,i1
    type(vol_ty_t),pointer::vol_ty_p
    type(fce_lst_t),pointer::fce_lst_p

    print*,'opening: ',trim(file_nm)
    open(iunit,file=trim(file_nm),status='old')
    do
      read(iunit,*)str

      if(trim(str)=='num_elements')then
        read(iunit,*)n_vol_ty,i1
        ivol_ty=i1-1

        allocate(vol_ty(i1:i1+n_vol_ty))
      end if

      if(trim(str)=='new_element')then
        ivol_ty=ivol_ty+1
        vol_ty_p=>vol_ty(ivol_ty)

        read(iunit,*) &
          vol_ty_p%n_dx,vol_ty_p%n_vrt, &
          vol_ty_p%n_fce,vol_ty_p%v_ty

        allocate(vol_ty_p%vrt_lst(vol_ty_p%n_vrt))
        allocate(vol_ty_p%fce_lst(vol_ty_p%n_fce))

        do ifce=1,vol_ty_p%n_fce
          fce_lst_p=>vol_ty_p%fce_lst(ifce)
          read(iunit,*)fce_lst_p%n_vrt

          allocate(fce_lst_p%vrt_lst(fce_lst_p%n_vrt))

          backspace(iunit)
          read(iunit,*) &
            fce_lst_p%n_vrt,fce_lst_p%vrt_lst
        end do
        if(ivol_ty==n_vol_ty)exit
      end if
    end do
    close(iunit)
  end subroutine



  ! vertex ordering for post-processing program
  subroutine grd_postpr(file_nm,vrt_lst1,vrt_lst2)
    character(len=*),intent(in),optional::file_nm
    integer,dimension(:),intent(in),optional::vrt_lst1
    integer,dimension(:),intent(inout),optional::vrt_lst2
    integer::iunit=100
    integer,dimension(2,4),save::arr_2d=0
    integer,dimension(4,8),save::arr_3d=0
    integer::n_vrt

    if(present(file_nm))then
      print*,'opening: ',trim(file_nm)
      open(iunit,file=trim(file_nm),status='old')
      read(iunit,*)arr_2d(1,:)  ! square
      read(iunit,*)arr_2d(2,:)  ! triangle
      read(iunit,*)arr_3d(1,:)  ! brick
      read(iunit,*)arr_3d(2,:)  ! wedge
      read(iunit,*)arr_3d(3,:)  ! tetrahedron
      read(iunit,*)arr_3d(4,:)  ! pyramid
      close(iunit)
    else
      n_vrt=size(vrt_lst1)
      if(n_dx==2.and.n_vrt==4)vrt_lst2=vrt_lst1(arr_2d(1,:))
      if(n_dx==2.and.n_vrt==3)vrt_lst2=vrt_lst1(arr_2d(2,:))
      if(n_dx==3.and.n_vrt==8)vrt_lst2=vrt_lst1(arr_3d(1,:))
      if(n_dx==3.and.n_vrt==6)vrt_lst2=vrt_lst1(arr_3d(2,:))
      if(n_dx==3.and.n_vrt==4)vrt_lst2=vrt_lst1(arr_3d(3,:))
      if(n_dx==3.and.n_vrt==5)vrt_lst2=vrt_lst1(arr_3d(4,:))
    end if
  end subroutine
end module
