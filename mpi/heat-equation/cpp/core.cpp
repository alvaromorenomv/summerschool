// Main solver routines for heat equation solver

#include <mpi.h>

#include "heat.hpp"

// Exchange the boundary values
void exchange(Field& field, const ParallelData parallel)
{

    double* sbuf;
    double* rbuf;
    // TODO start: implement halo exchange

    // You can utilize the data() method of the Matrix class to obtain pointer
    // to element, e.g. field.temperature.data(i, j)

    // Send to up, receive from down
    sbuf = field.temperature.data(1,0);
    rbuf = field.temperature.data(field.nx+1,0);
    MPI_Sendrecv(sbuf,field.ny+2,MPI_DOUBLE,parallel.nup,11,rbuf,field.ny+2,MPI_DOUBLE,
    parallel.ndown, 11, parallel.comunnicator, MPI_STATUS_IGNORE);

    // Send to down, receive from up
    sbuf = field.temperature.data(field.nx,0);
    rbuf = field.temperature.data(0,0);
    MPI_Sendrecv(sbuf,field.ny+2,MPI_DOUBLE,parallel.ndown,12,rbuf,field.ny+2,MPI_DOUBLE,
    parallel.nup, 12, parallel.comunnicator, MPI_STATUS_IGNORE);

    // TODO end
}

// Exchange the boundary values
void startExchangeNB(Field& field, const ParallelData parallel, std::vector<MPI_Request>& request)
{

    double* sbuf;
    MPI_Request sendToUp, sendToDown;

    // You can utilize the data() method of the Matrix class to obtain pointer
    // to element, e.g. field.temperature.data(i, j)

    // Send to up, receive from down
    sbuf = field.temperature.data(1,0);
    MPI_Isend(sbuf,field.ny+2, MPI_DOUBLE, parallel.nup, 11,
    parallel.comunnicator, &sendToUp);

    // Send to down, receive from up
    sbuf = field.temperature.data(field.nx,0);
    MPI_Isend(sbuf,field.ny+2, MPI_DOUBLE, parallel.ndown, 12,
    parallel.comunnicator, &sendToDown);
    request.push_back(sendToUp);
    request.push_back(sendToDown);
}

// Update the temperature values using five-point stencil */
void computeInnerValues(Field& curr, const Field& prev, const double a, const double dt)
{

  // Compilers do not necessarily optimize division to multiplication, so make it explicit
  auto inv_dx2 = 1.0 / (prev.dx * prev.dx);
  auto inv_dy2 = 1.0 / (prev.dy * prev.dy);

  // Determine the temperature field at next time step
  // As we have fixed boundary conditions, the outermost gridpoints
  // are not updated.
  for (int i = 2; i < curr.nx; i++) {
    for (int j = 1; j < curr.ny+1; j++) {
            curr(i, j) = prev(i, j) + a * dt * (
	       ( prev(i + 1, j) - 2.0 * prev(i, j) + prev(i - 1, j) ) * inv_dx2 +
	       ( prev(i, j + 1) - 2.0 * prev(i, j) + prev(i, j - 1) ) * inv_dy2
               );
    }
  }

}

void endExchangeNB(Field& field, const ParallelData parallel, std::vector<MPI_Request>& request)
{

    double* rbuf;
    MPI_Request recvFromUp, recvFromDown;

    // You can utilize the data() method of the Matrix class to obtain pointer
    // to element, e.g. field.temperature.data(i, j)

    // Send to up, receive from down
    rbuf = field.temperature.data(field.nx+1,0);
    MPI_Irecv(rbuf,field.ny+2, MPI_DOUBLE, parallel.ndown, 11,
    parallel.comunnicator, &recvFromDown);

    // Send to down, receive from up
    rbuf = field.temperature.data(0,0);
    MPI_Irecv(rbuf,field.ny+2, MPI_DOUBLE, parallel.nup, 12,
    parallel.comunnicator, &recvFromUp);
    request.push_back(recvFromUp);
    request.push_back(recvFromDown);
    MPI_Waitall(4,&request[0],MPI_STATUS_IGNORE);
}

// Update the temperature values using five-point stencil */
void computeBorderValues(Field& curr, const Field& prev, const double a, const double dt)
{

  // Compilers do not necessarily optimize division to multiplication, so make it explicit
  auto inv_dx2 = 1.0 / (prev.dx * prev.dx);
  auto inv_dy2 = 1.0 / (prev.dy * prev.dy);

  // Determine the temperature field at next time step
  // As we have fixed boundary conditions, the outermost gridpoints
  // are not updated.

  for (int j = 1; j < curr.ny + 1; j++) {
      curr(1, j) = prev(1, j) + a * dt * (
	       ( prev(2, j) - 2.0 * prev(1, j) + prev(1 - 1, j) ) * inv_dx2 +
	       ( prev(1, j + 1) - 2.0 * prev(1, j) + prev(1, j - 1) ) * inv_dy2
      );

      curr(curr.nx, j) = prev(curr.nx, j) + a * dt * (
	       ( prev(curr.nx + 1, j) - 2.0 * prev(curr.nx, j) + prev(curr.nx - 1, j) ) * inv_dx2 +
	       ( prev(curr.nx, j + 1) - 2.0 * prev(curr.nx, j) + prev(curr.nx, j - 1) ) * inv_dy2
      );

  }

}


// Update the temperature values using five-point stencil */
void evolve(Field& curr, const Field& prev, const double a, const double dt)
{

  // Compilers do not necessarily optimize division to multiplication, so make it explicit
  auto inv_dx2 = 1.0 / (prev.dx * prev.dx);
  auto inv_dy2 = 1.0 / (prev.dy * prev.dy);

  // Determine the temperature field at next time step
  // As we have fixed boundary conditions, the outermost gridpoints
  // are not updated.
  for (int i = 1; i < curr.nx + 1; i++) {
    for (int j = 1; j < curr.ny + 1; j++) {
            curr(i, j) = prev(i, j) + a * dt * (
	       ( prev(i + 1, j) - 2.0 * prev(i, j) + prev(i - 1, j) ) * inv_dx2 +
	       ( prev(i, j + 1) - 2.0 * prev(i, j) + prev(i, j - 1) ) * inv_dy2
               );
    }
  }

}
