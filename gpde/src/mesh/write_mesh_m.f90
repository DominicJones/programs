module write_mesh_m
implicit none
contains




subroutine write_vertices(ios,fnm,n_vrt,key,x)
  use const_m
  use mesh_data_m

  integer, intent(in)::ios
  character(*), intent(in)::fnm
  integer, intent(in)::n_vrt
  integer,dimension(:), intent(in)::key
  real(rp),dimension(:,:), intent(in)::x

  real(rp),dimension(3)::x3
  integer::i,j,n_dim


  print "(1x,a,a)","creating: ",trim(fnm)
  open(ios,file=trim(fnm))

!   write(ios,"(a)")"$MeshFormat"
!   write(ios,"(3i7)")2,0,8
!   write(ios,"(a)")"$EndMeshFormat"

!   write(ios,"(a)")"$Nodes"

!   write(ios,"(i7)") n_vrt
  write(ios,*) n_vrt

  n_dim = size(x,1)

  do i=1,n_vrt
    j = key(i)
    x3 = zero
    x3(:n_dim) = x(:,j)
    where (abs(x3) < small) x3 = zero
!     write(ios,"(i7,3es14.4)") j, x3
    write(ios,*) j, x3
  end do

!   write(ios,"(a)")"$EndNodes"

  close(ios)
  print "(1x,a,a)","closed:   ",trim(fnm)
end subroutine




subroutine write_mesh(ios,fnm,msh,x)
  use const_m
  use mesh_data_m

  integer,intent(in)::ios
  character(*),intent(in)::fnm
  type(msh_t),intent(inout)::msh
  real(rp),dimension(:,:)::x

  real(rp),dimension(3)::x3
  integer::i,key,n_dx,ty,ph_ty


  print "(1x,a,a)","creating: ",trim(fnm)
  open(ios,file=trim(fnm))


  write(ios,"(a)")"$MeshFormat"
  write(ios,"(3i7)")2,0,8
  write(ios,"(a)")"$EndMeshFormat"


  write(ios,"(a)")"$Nodes"
  write(ios,"(i7)")msh%n_vrt

  do i=1,msh%n_vrt
    key=msh%vrt_tag(i)%lst(1)
    x3=0; x3(:msh%n_dx)=x(:,i)
    where(abs(x3)<small)x3=0
    write(ios,"(i7,3es14.4)")key,x3
  end do

  write(ios,"(a)")"$EndNodes"


  write(ios,"(a)")"$Elements"
  write(ios,"(i7)")msh%n_ebf

  do i=1,msh%n_bfce
    key=msh%fce_tag(i)%lst(1)
    ty=msh%fce_tag(i)%lst(3)
    ph_ty=msh%fce_tag(i)%lst(4)
    write(ios,"(20i7)")key,ty,1,ph_ty,msh%fce_vrt(i)%lst
  end do

  do i=1,msh%n_elm
    key=msh%elm_tag(i)%lst(1)
    ty=msh%elm_tag(i)%lst(3)
    ph_ty=msh%elm_tag(i)%lst(4)
    write(ios,"(20i7)")key,ty,1,ph_ty,msh%elm_vrt(i)%lst
  end do

  write(ios,"(a)")"$EndElements"


  close(ios)
  print "(1x,a,a)","closed:   ",trim(fnm)
end subroutine



! NOT CORRECT IN 3D
subroutine write_background_mesh(ios,fnm,msh,x,phi_vrt)
  use const_m
  use mesh_data_m

  character(*),intent(in)::fnm
  integer,intent(in)::ios
  type(msh_t),intent(in)::msh
  real(rp),dimension(:,:)::x
  real(rp),dimension(:),intent(in)::phi_vrt

  integer::ni,nj,i,j,k
  real(rp),dimension(3)::x3

  print "(1x,a,a)","creating: ",trim(fnm)
  open(ios,file=trim(fnm))
  write(ios,"(a)")'View "background mesh" {'

  ni=msh%n_elm
  do i=1,ni
    write(ios,"(a)",advance="no")"ST("

    nj=size(msh%elm_vrt(i)%lst)
    do j=1,nj
      k=msh%elm_vrt(i)%lst(j)
      x3=0; x3(:msh%n_dx)=x(:,k)
      if(j==nj)then
        write(ios,"(es12.4,a,es12.4,a,es12.4)",advance="no") &
          x3(1),",",x3(2),",",x3(3)
      else
        write(ios,"(es12.4,a,es12.4,a,es12.4,a)",advance="no") &
          x3(1),",",x3(2),",",x3(3),","
      end if
    end do

    write(ios,"(a)",advance="no")")  {"

    do j=1,nj
      k=msh%elm_vrt(i)%lst(j)
      if(j==nj)then
        write(ios,"(es12.4)",advance="no")phi_vrt(k)
      else
        write(ios,"(es12.4,a)",advance="no")phi_vrt(k),","
      end if
    end do

    write(ios,"(a)")"};"
  end do
  write(ios,"(a)")"};"

  close(ios)
  print "(1x,a,a)","closed:  ",trim(fnm)
end subroutine
end module
