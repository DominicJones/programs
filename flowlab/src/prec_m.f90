module prec_m
#ifdef USE_MPI
use mpi
#endif

implicit none

integer, parameter :: rsp = kind (1.e0)
integer, parameter :: rdp = kind (1.d0)

#ifdef USE_MPI
integer, parameter :: mpi_rsp = mpi_real
integer, parameter :: mpi_rdp = mpi_double_precision
#else
integer, parameter :: mpi_rsp = rsp
integer, parameter :: mpi_rdp = rdp
#endif
end module
