#include "machine.h"

struct MPI_Comm;
struct MPI_Status;
struct MPI_Request;

#ifdef USE_MPI
#include "mpi.h"
#include "parmetis.h"
#endif

// mpi_comm_world = MPI_COMM_WORLD;
// mpi_stat = (MPI_Status*) malloc (mpi_size * sizeof(MPI_Status));
// mpi_rqst = (MPI_Request*) malloc (mpi_size * sizeof(MPI_Request));

void machine_t::initialize(int argc, char** argv) {
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif
}

void machine_t::finalize() {
#ifdef USE_MPI
  MPI_Finalize();
#endif
}

int machine_t::rank() {
  int value = 0;
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &value);
#endif
  return value;
}

int machine_t::size() {
  int value = 1;
#ifdef USE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &value);
#endif
  return value;
}

void machine_t::all_gather(const int * sendbuf, int sendcount,
                           int * recvbuf, int recvcount) {
#ifdef USE_MPI
  MPI_Allgather(sendbuf, sendcount, MPI_INT,
                recvbuf, recvcount, MPI_INT, MPI_COMM_WORLD);
#else
  *recvbuf = *sendbuf;
#endif
}
