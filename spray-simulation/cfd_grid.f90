module cfd_grid_m
  use cfd_math_m
  use cfd_alloc_m
  implicit none
  contains
!
! convert boundary name to boundary index
!
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
!
! read in Gambit neutral grid and store data
!
  subroutine read_grid(igrd,file_nm)
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
            grd(igrd)%vol(idx)%vol_ty,grd(igrd)%vol(idx)%n_vrt
          vol_p=>grd_p%vol(idx)
          allocate(vol_p%vrt_lst(vol_p%n_vrt)); vol_p%vrt_lst=0
          backspace(iunit)
          read(iunit,'(i8,1x,i2,1x,i2,1x,7i8:/(15x,7i8:))')idx, &
            vol_p%vol_ty,vol_p%n_vrt,vol_p%vrt_lst
          vol_p%n_fce=vol_ty(vol_p%vol_ty)%n_fce
          allocate(vol_p%fce_lst(vol_p%n_fce)); vol_p%fce_lst=0
          allocate(vol_p%bnd_lst(vol_p%n_fce)); vol_p%bnd_lst=0
          allocate(vol_p%neig_lst(vol_p%n_fce)); vol_p%neig_lst=0
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
          vol_p%bnd_vol=.true.
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


    call find_neig(igrd)
    call c_neig_lst(igrd)
    call grid_geom(igrd)

    allocate(pole(grd(igrd)%n_vol))
    allocate(neig(grd(igrd)%n_fce))
    allocate(phi(grd(igrd)%n_vbf))
  end subroutine
!
! compare two lists of verticies
!
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
!
! link faces to associated volumes
!
  subroutine find_neig(igrd)
    integer,intent(in)::igrd
    integer::i,j,ip,in,ivf,ivrt,ivol,ifce,ifce1,ifce2
    integer::n_fce,n_vrt,n_vrt_fce,ivol_ty,n_mtch
    real::t1,t2
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
    if(n_dx==2)n_vrt_fce=8
    if(n_dx==3)n_vrt_fce=16
    do ivrt=1,grd_p%n_vrt
      vrt_p=>grd_p%vrt(ivrt)
      allocate(vrt_p%fce_lst(n_vrt_fce))
      vrt_p%fce_lst=0
    end do

    ifce=0
    do ivol=1,grd_p%n_vol
      vol_p=>grd_p%vol(ivol)
      do ivf=1,vol_p%n_fce
        ifce=ifce+1
        fce_lst(ifce)%ip=ivol
        ivol_ty=vol_p%vol_ty
        n_vrt=vol_ty(ivol_ty)%fce_lst(ivf)%n_vrt
        fce_lst(ifce)%n_vrt=n_vrt
        allocate(fce_lst(ifce)%vrt_lst(n_vrt))

        do i=1,n_vrt
          ivrt=vol_ty(ivol_ty)%fce_lst(ivf)%vrt_lst(i)
          fce_lst(ifce)%vrt_lst(i)=vol_p%vrt_lst(ivrt)
          ivrt=fce_lst(ifce)%vrt_lst(i)
          do j=1,n_vrt_fce
            vrt_p=>grd_p%vrt(ivrt)
            if(vrt_p%fce_lst(j)==0)then
              vrt_p%fce_lst(j)=ifce
              vrt_p%n_fce=vrt_p%n_fce+1
              exit
            end if
          end do
        end do
      end do
    end do

    ! find matching faces
    n_mtch=0
    do ifce=1,n_fce_lst
      if(.not.fce_lst(ifce)%mtchd)then
        do ifce1=1,n_vrt_fce
          vrt_p=>grd_p%vrt(fce_lst(ifce)%vrt_lst(1))
          ifce2=vrt_p%fce_lst(ifce1)
          if(ifce2==0.or.ifce2==ifce)cycle
          if(.not.fce_lst(ifce2)%mtchd.and.&
            match_vrt(fce_lst(ifce)%vrt_lst, &
            fce_lst(ifce2)%vrt_lst))then
            fce_lst(ifce)%in=fce_lst(ifce2)%ip
            fce_lst(ifce)%mtchd=.true.
            fce_lst(ifce2)%copy=.true.
            n_mtch=n_mtch+1
            exit
          end if
        end do
      end if
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
    do ivrt=1,grd_p%n_vrt
      vrt_p=>grd_p%vrt(ivrt)
      deallocate(vrt_p%fce_lst)
    end do
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
  end subroutine
!
! construct volume face lists
!
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
!
! construct volume neighbour lists
!
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
          vol_p%neig_lst(ivf)=ip
        else if(in>0.and.in/=ivol)then
          vol_p%neig_lst(ivf)=in
        end if
      end do
    end do
  end subroutine
!
! calculate grid geometric quantities
!
  subroutine grid_geom(igrd)
    integer,intent(in)::igrd
    real(kind=wp)::r_p,v_p,t_vol=0.,t_s_f=0.
    real(kind=wp),dimension(:,:),allocatable::vec
    real(kind=wp),dimension(:),allocatable::sc
    real(wp),dimension(n_dx)::n_sum
    integer::i,ivol,ifce,ifce1,ivrt,ivrt1,ivrt2,ibnd,ip,in
    integer::ivf,ifv,ivv,n_fv,ivf2
    integer,dimension(:),allocatable::vrt_lst
    type(grd_t),pointer::grd_p
    type(vol_t),pointer::vol_p,vol_n
    type(fce_t),pointer::fce_p
    type(vrt_t),pointer::vrt_p

    print*,'calculating geometry...'
    grd_p=>grd(igrd)
    allocate(vec(1,n_dx))
    do ivol=1,grd_p%n_vol
      vol_p=>grd_p%vol(ivol)

      ! volume centre
      vec(1,:)=0.
      do ivrt=1,vol_p%n_vrt
        ivrt2=vol_p%vrt_lst(ivrt)
        vec(1,:)=vec(1,:)+grd_p%vrt(ivrt2)%x
      end do
      allocate(vol_p%x(n_dx))
      vol_p%x=vec(1,:)/vol_p%n_vrt

      ! set boundary types
      if(vol_p%bnd_vol)then
      do ivf=1,vol_p%n_fce
        if(vol_p%bnd_lst(ivf)>0)then
          n_fv=vol_ty(vol_p%vol_ty)%fce_lst(ivf)%n_vrt
          allocate(vrt_lst(n_fv))
          do ifv=1,n_fv
            ivv=vol_ty(vol_p%vol_ty)%fce_lst(ivf)%vrt_lst(ifv)
            vrt_lst(ifv)=vol_p%vrt_lst(ivv)
          end do

          do ivf2=1,vol_p%n_fce
            ifce=vol_p%fce_lst(ivf2)
            fce_p=>grd_p%fce(ifce)
            if(match_vrt(fce_p%vrt_lst,vrt_lst))then
              fce_p%bnd_ty=vol_p%bnd_lst(ivf)
              grd_p%n_bf=grd_p%n_bf+1
              fce_p%ib=grd_p%n_vol+grd_p%n_bf
            end if
          end do
          deallocate(vrt_lst)
        end if
      end do
      end if
    end do

    grd_p%n_vbf=grd_p%n_vol+grd_p%n_bf
    print*,'no. b. faces:',grd_p%n_bf
    print*,'no. i. faces:',grd_p%n_fce-grd_p%n_bf
    print*,'no. vol+b.f.:',grd_p%n_vbf


    do ifce=1,grd_p%n_fce
      fce_p=>grd_p%fce(ifce)
      ip=fce_p%ip
      vol_p=>grd_p%vol(ip)

      allocate(fce_p%x(n_dx),fce_p%n(n_dx))
      allocate(fce_p%x_p_aux(n_dx),fce_p%x_n_aux(n_dx))
      allocate(fce_p%dx_p(n_dx),fce_p%dx_n(n_dx))


      ! 2d: normal, face centre & area
      r_p=0.
      if(n_dx==2)then
        if(allocated(vec))deallocate(vec)
        allocate(vec(2,3))
        if(allocated(sc))deallocate(sc)
        allocate(sc(0:1))

        ivrt1=fce_p%vrt_lst(1)
        ivrt2=fce_p%vrt_lst(2)
        vec(1,1:2)=grd_p%vrt(ivrt2)%x-grd_p%vrt(ivrt1)%x
        vec(1,3)=0.
        sc(0)=vec_mag(vec(1,:))

        vec(2,:)=(/0.,0.,1./)
        vec(1,:)=cross_prod(vec(1,:),vec(2,:))
        fce_p%n=vec(1,1:2)/vec_mag(vec(1,1:2))
        fce_p%s=sc(0)

        if(cyl)then
          r_p=0.5*(grd_p%vrt(ivrt1)%x(2) &
            +grd_p%vrt(ivrt2)%x(2))
          fce_p%s=fce_p%s*r_p
        end if

        vec(1,1:2)=grd_p%vrt(ivrt1)%x+grd_p%vrt(ivrt2)%x
        fce_p%x=0.5*vec(1,1:2)
      end if


      ! 3d: normal, face centre & area
      if(n_dx==3)then
        if(allocated(vec))deallocate(vec)
        allocate(vec(0:fce_p%n_vrt,3)); vec=0.
        if(allocated(sc))deallocate(sc)
        allocate(sc(0:1)); sc=0.

        ivrt1=fce_p%vrt_lst(1)
        do ivrt=2,fce_p%n_vrt
          ivrt2=fce_p%vrt_lst(ivrt)
          vec(ivrt,:)=grd_p%vrt(ivrt2)%x-grd_p%vrt(ivrt1)%x

          if(ivrt>2)then
            vec(1,:)=0.5*cross_prod(vec(ivrt-1,:),vec(ivrt,:))
            sc(1)=vec_mag(vec(1,:))
            vec(0,:)=vec(0,:)+vec(1,:)/sc(1)
            sc(0)=sc(0)+sc(1)
          end if
        end do

        vec(0,:)=vec(0,:)/(fce_p%n_vrt-2)
        fce_p%n=vec(0,:)
        fce_p%s=sc(0)

        vec(0,:)=0.
        do ivrt=1,fce_p%n_vrt
          ivrt1=fce_p%vrt_lst(ivrt)
          vec(0,:)=vec(0,:)+grd_p%vrt(ivrt1)%x
        end do
        fce_p%x=vec(0,:)/fce_p%n_vrt
      end if


      ! volume
      v_p=fce_p%x(1)*fce_p%s*fce_p%n(1)
      vol_p%v=vol_p%v+v_p

      ! front face area
      if(n_dx==2.and.cyl)then
        vol_p%s_f=vol_p%s_f+v_p/(r_p+small)
      end if

      ! aux. pole coordinate
      fce_p%x_p_aux=fce_p%x-fce_p%n &
        *dot_product(fce_p%x-vol_p%x,fce_p%n)

      ! pole to aux. pole displacement
      fce_p%dx_p=fce_p%x_p_aux-vol_p%x

!       print*,'f',ifce,fce_p%x,fce_p%s,fce_p%n
!       print*,'f',ifce,fce_p%x_p_aux,fce_p%dx_p

      ! internal faces
      if(fce_p%in>0)then
        in=fce_p%in
        vol_n=>grd_p%vol(in)


        ! volume
        v_p=fce_p%x(1)*fce_p%s*fce_p%n(1)
        vol_n%v=vol_n%v-v_p

        ! front face area
        if(n_dx==2.and.cyl)then
          vol_n%s_f=vol_n%s_f-v_p/(r_p+small)
        end if

        ! aux. neig coordinate
        fce_p%x_n_aux=fce_p%x-fce_p%n &
          *dot_product(fce_p%x-vol_n%x,fce_p%n)

        ! neig. to aux. neig. displacement
        fce_p%dx_n=fce_p%x_n_aux-vol_n%x


        ! pole to neig. length
        fce_p%l_pn=vec_mag(vol_n%x-vol_p%x)
        fce_p%l_pn_aux=vec_mag(fce_p%x_n_aux-fce_p%x_p_aux)

        ! interpolation factor
        fce_p%w_fp=one-vec_mag(fce_p%x-vol_p%x) &
          /fce_p%l_pn
        fce_p%w_fp_aux=one-vec_mag(fce_p%x-fce_p%x_p_aux) &
          /fce_p%l_pn_aux

!         print*,'f',fce_p%l_pn,fce_p%l_pn_aux,fce_p%w_fp,fce_p%w_fp_aux
      ! boundary faces
      else

        ! pole to face length
        fce_p%l_pn=vec_mag(fce_p%x-vol_p%x)
        fce_p%l_pn_aux=vec_mag(fce_p%x-fce_p%x_p_aux)

        ! link face to volume
        fce_p%bnd_vol=ip
      end if
    end do


    ! volume reciprocal
    do ivol=1,grd_p%n_vol
      vol_p=>grd_p%vol(ivol)
      vol_p%v_r=1./vol_p%v
      t_vol=t_vol+vol_p%v
      t_s_f=t_s_f+vol_p%s_f

!       print*,'v',ivol,vol_p%x,vol_p%v,vol_p%s_f
    end do

    print*,'volume:',t_vol
    if(n_dx==2.and.cyl)then
      print*,'frontal area:',t_s_f
    end if
  end subroutine
end module
