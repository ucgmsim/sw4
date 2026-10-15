// PMPI wrappers that charge the cost of each MPI call to the running
// sw4timer region (see sw4timer.h). Each wrapper times the call and records
// the payload this rank sends: count * type size for point-to-point sends and
// for the contribution to a collective. Exchanges with MPI_PROC_NULL (domain
// edges) cost time but count no message or bytes. Receives and waits record
// time only, so bytes are never counted twice.
//
// Only the functions SW4 uses in its time loop are wrapped; anything else goes
// straight to the MPI library. Compiled out by SW4_NO_TIMERS or
// SW4_NO_MPI_TIMERS.
#include "sw4timer.h"

#if !defined(SW4_NO_TIMERS) && !defined(SW4_NO_MPI_TIMERS)

#include <mpi.h>

namespace
{
inline double payload( MPI_Datatype type, int count )
{
   int size = 0;
   PMPI_Type_size( type, &size );
   return static_cast<double>( size ) * count;
}

inline void charge( int fn, double t0, double bytes, int msgs )
{
   sw4timer::add_mpi( fn, PMPI_Wtime() - t0, bytes, msgs );
}
} // namespace

using namespace sw4timer;

extern "C" {

int MPI_Sendrecv( const void* sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag,
                  void* recvbuf, int recvcount, MPI_Datatype recvtype, int source, int recvtag,
                  MPI_Comm comm, MPI_Status* status )
{
   const double t0 = PMPI_Wtime();
   const int rc = PMPI_Sendrecv( sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount,
                                 recvtype, source, recvtag, comm, status );
   const bool real = dest != MPI_PROC_NULL;
   charge( MPIFN_SENDRECV, t0, real ? payload( sendtype, sendcount ) : 0.0, real ? 1 : 0 );
   return rc;
}

int MPI_Send( const void* buf, int count, MPI_Datatype type, int dest, int tag, MPI_Comm comm )
{
   const double t0 = PMPI_Wtime();
   const int rc = PMPI_Send( buf, count, type, dest, tag, comm );
   const bool real = dest != MPI_PROC_NULL;
   charge( MPIFN_SEND, t0, real ? payload( type, count ) : 0.0, real ? 1 : 0 );
   return rc;
}

int MPI_Recv( void* buf, int count, MPI_Datatype type, int source, int tag, MPI_Comm comm,
              MPI_Status* status )
{
   const double t0 = PMPI_Wtime();
   const int rc = PMPI_Recv( buf, count, type, source, tag, comm, status );
   charge( MPIFN_RECV, t0, 0.0, 0 );
   return rc;
}

int MPI_Isend( const void* buf, int count, MPI_Datatype type, int dest, int tag, MPI_Comm comm,
               MPI_Request* request )
{
   const double t0 = PMPI_Wtime();
   const int rc = PMPI_Isend( buf, count, type, dest, tag, comm, request );
   const bool real = dest != MPI_PROC_NULL;
   charge( MPIFN_ISEND, t0, real ? payload( type, count ) : 0.0, real ? 1 : 0 );
   return rc;
}

int MPI_Irecv( void* buf, int count, MPI_Datatype type, int source, int tag, MPI_Comm comm,
               MPI_Request* request )
{
   const double t0 = PMPI_Wtime();
   const int rc = PMPI_Irecv( buf, count, type, source, tag, comm, request );
   charge( MPIFN_IRECV, t0, 0.0, 0 );
   return rc;
}

int MPI_Wait( MPI_Request* request, MPI_Status* status )
{
   const double t0 = PMPI_Wtime();
   const int rc = PMPI_Wait( request, status );
   charge( MPIFN_WAIT, t0, 0.0, 0 );
   return rc;
}

int MPI_Waitall( int count, MPI_Request requests[], MPI_Status statuses[] )
{
   const double t0 = PMPI_Wtime();
   const int rc = PMPI_Waitall( count, requests, statuses );
   charge( MPIFN_WAITALL, t0, 0.0, 0 );
   return rc;
}

int MPI_Barrier( MPI_Comm comm )
{
   const double t0 = PMPI_Wtime();
   const int rc = PMPI_Barrier( comm );
   charge( MPIFN_BARRIER, t0, 0.0, 1 );
   return rc;
}

int MPI_Bcast( void* buf, int count, MPI_Datatype type, int root, MPI_Comm comm )
{
   const double t0 = PMPI_Wtime();
   const int rc = PMPI_Bcast( buf, count, type, root, comm );
   charge( MPIFN_BCAST, t0, payload( type, count ), 1 );
   return rc;
}

int MPI_Reduce( const void* sendbuf, void* recvbuf, int count, MPI_Datatype type, MPI_Op op, int root,
                MPI_Comm comm )
{
   const double t0 = PMPI_Wtime();
   const int rc = PMPI_Reduce( sendbuf, recvbuf, count, type, op, root, comm );
   charge( MPIFN_REDUCE, t0, payload( type, count ), 1 );
   return rc;
}

int MPI_Allreduce( const void* sendbuf, void* recvbuf, int count, MPI_Datatype type, MPI_Op op,
                   MPI_Comm comm )
{
   const double t0 = PMPI_Wtime();
   const int rc = PMPI_Allreduce( sendbuf, recvbuf, count, type, op, comm );
   charge( MPIFN_ALLREDUCE, t0, payload( type, count ), 1 );
   return rc;
}

int MPI_Gather( const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount,
                MPI_Datatype recvtype, int root, MPI_Comm comm )
{
   const double t0 = PMPI_Wtime();
   const int rc = PMPI_Gather( sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm );
   charge( MPIFN_GATHER, t0, payload( sendtype, sendcount ), 1 );
   return rc;
}

int MPI_Gatherv( const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
                 const int recvcounts[], const int displs[], MPI_Datatype recvtype, int root,
                 MPI_Comm comm )
{
   const double t0 = PMPI_Wtime();
   const int rc = PMPI_Gatherv( sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype,
                                root, comm );
   charge( MPIFN_GATHERV, t0, payload( sendtype, sendcount ), 1 );
   return rc;
}

int MPI_Allgather( const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
                   int recvcount, MPI_Datatype recvtype, MPI_Comm comm )
{
   const double t0 = PMPI_Wtime();
   const int rc = PMPI_Allgather( sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm );
   charge( MPIFN_ALLGATHER, t0, payload( sendtype, sendcount ), 1 );
   return rc;
}

int MPI_Allgatherv( const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
                    const int recvcounts[], const int displs[], MPI_Datatype recvtype, MPI_Comm comm )
{
   const double t0 = PMPI_Wtime();
   const int rc = PMPI_Allgatherv( sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs,
                                   recvtype, comm );
   charge( MPIFN_ALLGATHERV, t0, payload( sendtype, sendcount ), 1 );
   return rc;
}

} // extern "C"

#endif
