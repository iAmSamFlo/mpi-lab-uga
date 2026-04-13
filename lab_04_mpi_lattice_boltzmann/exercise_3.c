/*****************************************************
    AUTHOR  : Sébastien Valat
    MAIL    : sebastien.valat@univ-grenoble-alpes.fr
    LICENSE : BSD
    YEAR    : 2021
    COURSE  : Parallel Algorithms and Programming
*****************************************************/

//////////////////////////////////////////////////////
//
// Goal: Implement non-blocking 1D communication scheme
//       along X axis.
//
// SUMMARY:
//     - 1D splitting along X
// NEW:
//     - >>> Non-blocking communications <<<
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"

/****************************************************/
void lbm_comm_init_ex3(lbm_comm_t * comm, int total_width, int total_height)
{
	//we use the same implementation then ex1
	lbm_comm_init_ex1(comm, total_width, total_height);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex3(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	//
	// TODO: Implement the 1D communication with non-blocking MPI functions.
	//
	// To be used:
	//    - DIRECTIONS: the number of doubles composing a cell
	//    - double[9] lbm_mesh_get_cell(mesh, x, y): function to get the address of a particular cell.
	//    - comm->width : The with of the local sub-domain (containing the ghost cells)
	//    - comm->height : The height of the local sub-domain (containing the ghost cells)
	
	//example to access cell
	//double * cell = lbm_mesh_get_cell(mesh, local_x, local_y);
	//double * cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);

	int rank;
	int comm_size;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	MPI_Comm_size( MPI_COMM_WORLD, &comm_size );
	MPI_Request send, req;

	for (int i = 0; i < comm->height; i++){
        if (rank < comm_size - 1) {
            MPI_Isend(lbm_mesh_get_cell(mesh, comm->width - 2, i), DIRECTIONS, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &send);
            MPI_Wait(&send, MPI_STATUS_IGNORE);
        }
        if (rank > 0){
            MPI_Irecv(lbm_mesh_get_cell(mesh, 0, i), DIRECTIONS, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &req);
            MPI_Wait(&req, MPI_STATUS_IGNORE);
        }
    }

    for (int i = 0; i < comm->height; i++){
        if (rank > 0){
            MPI_Isend(lbm_mesh_get_cell(mesh, 1, i), DIRECTIONS, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &send);
            MPI_Wait(&send, MPI_STATUS_IGNORE);
        }
        if (rank < comm_size - 1){
            MPI_Irecv(lbm_mesh_get_cell(mesh, comm->width - 1, i), DIRECTIONS, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &req);
            MPI_Wait(&req, MPI_STATUS_IGNORE);
        }
    }
}
