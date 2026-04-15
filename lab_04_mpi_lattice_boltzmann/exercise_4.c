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

/****************************************************/
void lbm_comm_init_ex4(lbm_comm_t * comm, int total_width, int total_height)
{
	//
	// TODO: calculate the splitting parameters for the current task.
	//

	//get infos
	int rank;
	int comm_size;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	MPI_Comm_size( MPI_COMM_WORLD, &comm_size );

	//check
	if (comm_size == 1)
		fatal("Invalid communicator size, should not be 1 !");

	int s = sqrt(comm_size);
	if ((s * s) != comm_size) 
		fatal("Invalid communicator size, should be an integer sqrt!");

	// TODO: calculate the number of tasks along X axis and Y axis.
	comm->nb_x = s;
	comm->nb_y = s;

	if(total_width % comm->nb_x != 0){
		fatal("nb_x does not divide total width!");
	}
	if(total_height % comm->nb_y != 0){
		fatal("nb_y does not divide total width!");
	}

	// TODO: calculate the current task position in the splitting
	comm->rank_x = rank % comm->nb_x;
    comm->rank_y = rank % comm->nb_y;

	// TODO : calculate the local sub-domain size (do not forget the 
	//        ghost cells). Use total_width & total_height as starting 
	//        point.
	comm->width = (total_width / comm->nb_x) + 2;
	comm->height = (total_height / comm->nb_y) + 2;

	// TODO : calculate the absolute position  (in cell number) in the global mesh.
	//        without accounting the ghost cells
	//        (used to setup the obstable & initial conditions).
	comm->x = comm->rank_x * (total_width / comm->nb_x);
	comm->y = comm->rank_y * (total_height / comm->nb_y);

	//OPTIONAL : if you want to avoid allocating temporary copy buffer
	//           for every step :
	//comm->buffer_recv_down, comm->buffer_recv_up, comm->buffer_send_down, comm->buffer_send_up

	//if debug print comm
	//lbm_comm_print(comm);
}

/****************************************************/
void lbm_comm_release_ex4(lbm_comm_t * comm)
{
	//free allocated ressources
}

/****************************************************/
void lbm_comm_ghost_exchange_ex4(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	//
	// TODO: Implement the 2D communication with :
	//         - blocking MPI functions
	//         - manual copy in temp buffer for non contiguous side 
	//
	// To be used:
	//    - DIRECTIONS: the number of doubles composing a cell
	//    - double[9] lbm_mesh_get_cell(mesh, x, y): function to get the address of a particular cell.
	//    - comm->width : The with of the local sub-domain (containing the ghost cells)
	//    - comm->height : The height of the local sub-domain (containing the ghost cells)
	//
	// TIP: create a function to get the target rank from x,y task coordinate. 
	// TIP: You can use MPI_PROC_NULL on borders.
	// TIP: send the corner values 2 times, with the up/down/left/write communication
	//      and with the diagonal communication in a second time, this avoid
	//      special cases for border tasks.

	//example to access cell
	//double * cell = lbm_mesh_get_cell(mesh, local_x, local_y);
	//double * cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);

	//TODO:
	//   - implement left/write communications
	//   - implement top/bottom communication (non contiguous)
	//   - implement diagonal communications

	// int rank;
	// int comm_size;
	// MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	// MPI_Comm_size( MPI_COMM_WORLD, &comm_size );
	// MPI_Comm cart_comm;

	// int dim[2] = {comm->nb_x, comm->nb_y};
	// int periods[2] = {0, 0};

	// MPI_Cart_create(MPI_COMM_WORLD, 2, dim, periods, 0, &cart_comm);




}
