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
//      8 neighbors using MPI types for non contiguous
//      side.
//
// SUMMARY:
//     - 2D splitting along X and Y
//     - 8 neighbors communications
//     - Blocking communications
// NEW:
//     - >>> MPI type for non contiguous cells <<<
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"

/****************************************************/
void lbm_comm_init_ex5(lbm_comm_t * comm, int total_width, int total_height)
{
	// 1. Reuse the Cartesian topology and buffer allocation from Ex 4
	lbm_comm_init_ex4(comm, total_width, total_height);

	// 2. Create MPI type for non-contiguous rows
	// A row consists of 'width' blocks.
	// Each block contains 'DIRECTIONS' (9) doubles.
	// The stride between the start of one cell and the next cell in a row 
	// is (height * DIRECTIONS) because of column-major layout.
	MPI_Type_vector(
		comm->width,               // count: number of cells in a row
		DIRECTIONS,                // blocklength: doubles per cell
		comm->height * DIRECTIONS, // stride: distance between start of cells
		MPI_DOUBLE,                // oldtype
		&comm->type                // newtype
	);

	// 3. Commit the type so it can be used in communication
	MPI_Type_commit(&comm->type);
}

/****************************************************/
void lbm_comm_release_ex5(lbm_comm_t * comm)
{
	// 1. Release the MPI Datatype
	MPI_Type_free(&comm->type);

	// 2. Reuse the release function from Ex 4 to free buffers and communicator
	lbm_comm_release_ex4(comm);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex5(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	int rank_left, rank_right, rank_up, rank_down;
	
	// Get neighbors from the Cartesian communicator
	MPI_Cart_shift(comm->communicator, 0, 1, &rank_left,  &rank_right);
	MPI_Cart_shift(comm->communicator, 1, 1, &rank_up,    &rank_down);

	// ----------------------------------------------------------------
	// STEP 1: LEFT / RIGHT (Contiguous columns)
	// ----------------------------------------------------------------
	// Send right border (x=width-2), receive into left ghost (x=0)
	MPI_Sendrecv(
		lbm_mesh_get_cell(mesh, comm->width - 2, 0), comm->height * DIRECTIONS, MPI_DOUBLE, rank_right, 0,
		lbm_mesh_get_cell(mesh, 0,                0), comm->height * DIRECTIONS, MPI_DOUBLE, rank_left,  0,
		comm->communicator, MPI_STATUS_IGNORE
	);

	// Send left border (x=1), receive into right ghost (x=width-1)
	MPI_Sendrecv(
		lbm_mesh_get_cell(mesh, 1,               0), comm->height * DIRECTIONS, MPI_DOUBLE, rank_left,  1,
		lbm_mesh_get_cell(mesh, comm->width - 1, 0), comm->height * DIRECTIONS, MPI_DOUBLE, rank_right, 1,
		comm->communicator, MPI_STATUS_IGNORE
	);

	// ----------------------------------------------------------------
	// STEP 2: TOP / BOTTOM (Non-contiguous rows using MPI Datatype)
	// By exchanging full width (0 to width-1), corners are handled 
	// automatically because Step 1 is already finished.
	// ----------------------------------------------------------------

	// Exchange Top: Send real row (y=1) UP, receive into BOTTOM ghost (y=height-1)
	MPI_Sendrecv(
		lbm_mesh_get_cell(mesh, 0, 1),                1, comm->type, rank_up,   2,
		lbm_mesh_get_cell(mesh, 0, comm->height - 1), 1, comm->type, rank_down, 2,
		comm->communicator, MPI_STATUS_IGNORE
	);

	// Exchange Bottom: Send real row (y=height-2) DOWN, receive into TOP ghost (y=0)
	MPI_Sendrecv(
		lbm_mesh_get_cell(mesh, 0, comm->height - 2), 1, comm->type, rank_down, 3,
		lbm_mesh_get_cell(mesh, 0, 0),                1, comm->type, rank_up,   3,
		comm->communicator, MPI_STATUS_IGNORE
	);
}
