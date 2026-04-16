/*****************************************************
    AUTHOR  : Sébastien Valat
    MAIL    : sebastien.valat@univ-grenoble-alpes.fr
    LICENSE : BSD
    YEAR    : 2021
    COURSE  : Parallel Algorithms and Programming
*****************************************************/

//////////////////////////////////////////////////////
//
// Goal: Implement 2D grid communication scheme with
//       8 neighbors using manual copy for non
//       contiguous side and blocking communications
//
// SUMMARY:
//     - 2D splitting along X and Y
//     - 8 neighbors communications
//     - Blocking communications
//     - Manual copy for non continguous cells
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"
#include <math.h>
#include <stdlib.h>

/****************************************************/
void lbm_comm_init_ex4(lbm_comm_t * comm, int total_width, int total_height)
{
    int rank;
    int comm_size;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &comm_size );
 
    if (comm_size == 1)
        fatal("Invalid communicator size, should not be 1 !");
 
    int s = sqrt(comm_size);
    if ((s * s) != comm_size)
        fatal("Invalid communicator size, should be an integer sqrt!");
 
    // Number of tasks along each axis
    comm->nb_x = s;
    comm->nb_y = s;
 
    if (total_width % comm->nb_x != 0)
        fatal("nb_x does not divide total width!");
    if (total_height % comm->nb_y != 0)
        fatal("nb_y does not divide total height!");
 
    // Create a 2D Cartesian communicator (non-periodic)
    int dims[2]    = { comm->nb_x, comm->nb_y };
    int periods[2] = { 0, 0 };
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm->communicator);
 
    // Get position of current task in the Cartesian grid
    int coords[2];
    MPI_Cart_coords(comm->communicator, rank, 2, coords);
    comm->rank_x = coords[0];
    comm->rank_y = coords[1];
 
    // Local subdomain size (with ghost cells)
    comm->width  = (total_width  / comm->nb_x) + 2;
    comm->height = (total_height / comm->nb_y) + 2;
 
    // Absolute position in global mesh (ignoring ghost cells)
    comm->x = comm->rank_x * (total_width  / comm->nb_x);
    comm->y = comm->rank_y * (total_height / comm->nb_y);
 
	// Buffers for non-contiguous rows (Full width: width * DIRECTIONS)
    size_t row_size = comm->width * DIRECTIONS * sizeof(double);
    comm->buffer_send_up   = (double*)malloc(row_size);
    comm->buffer_recv_up   = (double*)malloc(row_size);
    comm->buffer_send_down = (double*)malloc(row_size);
    comm->buffer_recv_down = (double*)malloc(row_size);
 
    if (!comm->buffer_send_up || !comm->buffer_recv_up ||
        !comm->buffer_send_down || !comm->buffer_recv_down)
        fatal("Failed to allocate communication buffers!");
 
    // lbm_comm_print(comm);
}

/****************************************************/
void lbm_comm_release_ex4(lbm_comm_t * comm)
{
	//free allocated ressources
	free(comm->buffer_send_up);
    free(comm->buffer_recv_up);
    free(comm->buffer_send_down);
    free(comm->buffer_recv_down);
    MPI_Comm_free(&comm->communicator);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex4(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
    int rank_left, rank_right, rank_up, rank_down;
    // Dimension 0 = X, Dimension 1 = Y
    MPI_Cart_shift(comm->communicator, 0, 1, &rank_left,  &rank_right);
    MPI_Cart_shift(comm->communicator, 1, 1, &rank_up,    &rank_down);

    // ----------------------------------------------------------------
    // STEP 1: LEFT / RIGHT (Contiguous columns)
    // ----------------------------------------------------------------
    // Send right border (width-2), receive into left ghost (0)
    MPI_Sendrecv(
        lbm_mesh_get_cell(mesh, comm->width - 2, 0), comm->height * DIRECTIONS, MPI_DOUBLE, rank_right, 0,
        lbm_mesh_get_cell(mesh, 0,                0), comm->height * DIRECTIONS, MPI_DOUBLE, rank_left,  0,
        comm->communicator, MPI_STATUS_IGNORE
    );
 
    // Send left border (1), receive into right ghost (width-1)
    MPI_Sendrecv(
        lbm_mesh_get_cell(mesh, 1,               0), comm->height * DIRECTIONS, MPI_DOUBLE, rank_left,  1,
        lbm_mesh_get_cell(mesh, comm->width - 1, 0), comm->height * DIRECTIONS, MPI_DOUBLE, rank_right, 1,
        comm->communicator, MPI_STATUS_IGNORE
    );

    // ----------------------------------------------------------------
    // STEP 2: TOP / BOTTOM (Non-contiguous rows)
    // ----------------------------------------------------------------
    
    // PACK: Only bother packing if we actually have somewhere to send it
    for (int x = 0; x < comm->width; x++) {
        double *src_up   = lbm_mesh_get_cell(mesh, x, 1);
        double *src_down = lbm_mesh_get_cell(mesh, x, comm->height - 2);
        for (int d = 0; d < DIRECTIONS; d++) {
            comm->buffer_send_up[x * DIRECTIONS + d]   = src_up[d];
            comm->buffer_send_down[x * DIRECTIONS + d] = src_down[d];
        }
    }
 
    // EXCHANGE
    // Send to UP, receive from DOWN
    MPI_Sendrecv(
        comm->buffer_send_up,   comm->width * DIRECTIONS, MPI_DOUBLE, rank_up,   2,
        comm->buffer_recv_down, comm->width * DIRECTIONS, MPI_DOUBLE, rank_down, 2,
        comm->communicator, MPI_STATUS_IGNORE
    );
 
    // Send to DOWN, receive from UP
    MPI_Sendrecv(
        comm->buffer_send_down, comm->width * DIRECTIONS, MPI_DOUBLE, rank_down, 3,
        comm->buffer_recv_up,   comm->width * DIRECTIONS, MPI_DOUBLE, rank_up,   3,
        comm->communicator, MPI_STATUS_IGNORE
    );
 
    // UNPACK: ONLY if the neighbor exists!
    // If rank_down is NULL, we are at the global bottom; don't touch the ghost cells.
    if (rank_down != MPI_PROC_NULL) {
        for (int x = 0; x < comm->width; x++) {
            double *src = comm->buffer_recv_down + x * DIRECTIONS;
            double *dst = lbm_mesh_get_cell(mesh, x, comm->height - 1);
            for (int d = 0; d < DIRECTIONS; d++) dst[d] = src[d];
        }
    }
 
    // If rank_up is NULL, we are at the global top; don't touch the ghost cells.
    if (rank_up != MPI_PROC_NULL) {
        for (int x = 0; x < comm->width; x++) {
            double *src = comm->buffer_recv_up + x * DIRECTIONS;
            double *dst = lbm_mesh_get_cell(mesh, x, 0);
            for (int d = 0; d < DIRECTIONS; d++) dst[d] = src[d];
        }
    }
}