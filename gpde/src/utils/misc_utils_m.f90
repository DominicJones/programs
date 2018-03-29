module misc_utils_m
implicit none


contains



! split "path/file" into "path/" and "file"
subroutine split_fnm(filename,path,file)
  character(*),intent(in)::filename
  character(*),intent(inout)::path,file
  integer::n,i,j

  n=len_trim(filename); j=0
  do i=n,1,-1
    if(filename(i:i)=="/")then
      j=i; exit
    end if
  end do
  path=filename(:j)
  file=filename(j+1:n)
end subroutine



! print to screen the licence file
subroutine header(fnm,ios)
  character(*),intent(in)::fnm
  integer,intent(in)::ios
  character(80)::str
  integer::istat

  print *,"opening: ",trim(fnm)
  open(ios,file=trim(fnm),status="old")
  do
    read(ios,"(a80)",iostat=istat)str
    if(trim(str)=="EOH")exit
    if(istat/=0)exit
    write(*,"(a)")str(:len_trim(str))
  end do
  close(ios)
  print *,"closed:  ",trim(fnm)
end subroutine



subroutine last_fixed_point_iter(iter,n_iter,norm,cutoff,last_iter)
  use const_m

  integer,intent(in)::iter,n_iter
  real(rp),intent(in)::norm,cutoff
  logical,intent(out)::last_iter

  last_iter = .false.

  if(iter == n_iter)then
    last_iter = .true.
  else if(norm < cutoff)then
    last_iter = .true.
  end if

  if(last_iter)print *, "Last iteration"
end subroutine


! fixed-point iteration control:
! ctrl = F for primal or tangent (Forward) iterations
! ctrl = R for adjoint (Reverse) iterations
! id indicates which system of equations is being converged; mesh deformation or flow eqn
! iter is the iteration count going from 1 to n_iter in Forward mode and from n_iter to 1 in Reverse mode
! k_iter is how frequently to print the iteration details
! fp_exit is a flag to tell the iteration loop when to break
! cutoff is the minimum tolerance of the global norm
! norm is the residual 1-norm from dx,dy,dz, and u,v,w,p
! norm_gl is the global maximum normalized norm
subroutine fp_iter_ctrl(ctrl,id,iter,n_iter,k_iter,fp_exit,cutoff,norm,norm_gl)
  use const_m
  use pde_data_m, only: last_fp_iter

  integer,intent(in)::id
  integer,intent(in)::iter,n_iter,k_iter
  real(rp),dimension(:),intent(inout)::norm
  real(rp),intent(in)::cutoff
  integer,intent(inout)::fp_exit
  real(rp),intent(inout)::norm_gl
  character,intent(in)::ctrl

  real(rp),save::norm_mx
  integer::i,ios


  last_fp_iter=.false.

  if((ctrl=='F'.and.iter==1) .or. (ctrl=='R'.and.iter==n_iter))then
    fp_exit=0
    norm_mx=0.d0
  end if


  norm_gl=maxval(norm)
  norm=0.d0
  norm_mx=max(norm_mx,norm_gl)
  norm_gl=norm_gl/(norm_mx+1.d-20)

  if(norm_gl<cutoff .or. (ctrl=='F'.and.iter>n_iter-2))then
    fp_exit=fp_exit+1
    if(fp_exit==1)last_fp_iter=.true.
  end if


  if (ctrl=='F') ios=100+id
  if (ctrl=='R') ios=200+id

  if(k_iter>=0 .and. mod(iter,k_iter)==0)then
    print "(1x,a,i1,a,i6,a,es12.4)",ctrl,id," iteration:",iter,",  Norm:",norm_gl
    ! if (id == 2) write(ios,"(i6,1x,es12.4)")iter,norm_gl
  end if
end subroutine

end module
