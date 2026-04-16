/*****************************************************
    AUTHOR  : Sébastien Valat
    MAIL    : sebastien.valat@univ-grenoble-alpes.fr
    LICENSE : BSD
    YEAR    : 2021
    COURSE  : Parallel Algorithms and Programming
*****************************************************/

//////////////////////////////////////////////////////
//
// Goal: Implement 2D grid communication with non-blocking
//       messages.
//
// SUMMARY:
//     - 2D splitting along X and Y
//     - 8 neighbors communications
//     - MPI type for non contiguous cells
// NEW:
//     - Non-blocking communications
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"

/****************************************************/
void lbm_comm_init_ex6(lbm_comm_t * comm, int total_width, int total_height)
{
	//we use the same implementation than ex5
	lbm_comm_init_ex5(comm, total_width, total_height);
}

/****************************************************/
void lbm_comm_release_ex6(lbm_comm_t * comm)
{
	//we use the same implementation than ext 5
	lbm_comm_release_ex5(comm);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex6(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	int rank_left, rank_right, rank_up, rank_down;
	MPI_Request requests[4];
	
	// Get neighbor ranks
	MPI_Cart_shift(comm->communicator, 0, 1, &rank_left,  &rank_right);
	MPI_Cart_shift(comm->communicator, 1, 1, &rank_up,    &rank_down);

	// ----------------------------------------------------------------
	// BATCH 1: LEFT / RIGHT (Contiguous columns)
	// ----------------------------------------------------------------
	
	// Prepare receptions into ghost columns
	MPI_Irecv(lbm_mesh_get_cell(mesh, 0, 0), comm->height * DIRECTIONS, MPI_DOUBLE, 
	          rank_left, 0, comm->communicator, &requests[0]);
	MPI_Irecv(lbm_mesh_get_cell(mesh, comm->width - 1, 0), comm->height * DIRECTIONS, MPI_DOUBLE, 
	          rank_right, 1, comm->communicator, &requests[1]);

	// Prepare sends from border columns
	MPI_Isend(lbm_mesh_get_cell(mesh, 1, 0), comm->height * DIRECTIONS, MPI_DOUBLE, 
	          rank_left, 1, comm->communicator, &requests[2]);
	MPI_Isend(lbm_mesh_get_cell(mesh, comm->width - 2, 0), comm->height * DIRECTIONS, MPI_DOUBLE, 
	          rank_right, 0, comm->communicator, &requests[3]);

	// Wait for Left/Right exchange to finish so corners are ready for Batch 2
	MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);

	// ----------------------------------------------------------------
	// BATCH 2: TOP / BOTTOM (Non-contiguous rows using MPI Datatype)
	// ----------------------------------------------------------------

	// Prepare receptions into ghost rows
	MPI_Irecv(lbm_mesh_get_cell(mesh, 0, 0), 1, comm->type, 
	          rank_up, 3, comm->communicator, &requests[0]);
	MPI_Irecv(lbm_mesh_get_cell(mesh, 0, comm->height - 1), 1, comm->type, 
	          rank_down, 2, comm->communicator, &requests[1]);

	// Prepare sends from border rows
	MPI_Isend(lbm_mesh_get_cell(mesh, 0, 1), 1, comm->type, 
	          rank_up, 2, comm->communicator, &requests[2]);
	MPI_Isend(lbm_mesh_get_cell(mesh, 0, comm->height - 2), 1, comm->type, 
	          rank_down, 3, comm->communicator, &requests[3]);

	// Final wait to ensure Top/Bottom exchange is complete before physics starts
	MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);
}