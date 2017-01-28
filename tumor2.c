#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <math.h>  
#include <time.h>  
#include "alloc.h"

void swap (int *a, int *b)
{
	// Swap arguments in array.
    int temp = *a;
    *a = *b;
    *b = temp;
}

void shuffle(int arr[], int n) {
	// Shuffles n elements of array arr.
	
	// Start from the last element and swap one by one. We don't
    // need to run for the first element that's why i > 0
	int i;
    for (i = n-1; i > 0; i--)
    {
        // Pick a random index from 0 to i
        int j = rand() % (i+1);
 
        // Swap arr[i] with the element at random index
        swap(&arr[i], &arr[j]);
	}
}

double rand2()
{
	// Returns random number between 0 and 1 
    return (double)rand() / (double)(RAND_MAX);
}

int main(int argc, char **argv) {

	int rank, size, provided, idle = 0;
	int lx, ly, lz;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);	/* starts MPI */
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);  /* get current process id */
	MPI_Comm_size(MPI_COMM_WORLD, &size);  /* get number of processes */
	
	/* Setup variables for life function */
	int ii, jj, kk, n, nsum, noccupied, nOB, liveCells = 0, index = 0;
	int maxCells;
	int ***matrix_a, ***matrix;
	int *matrix_addresses_x, *matrix_addresses_y, *matrix_addresses_z, *index_list;
	int *available_x, *available_y, *available_z, *available_indices;
	double o2;
	int next_x, prev_x;
	double rm, rp, rd;
	int migration_target, prolif_target, check, freeCount, random_order, space, done;
	int xx, yy, zz, xm, xp, ym, yp, zm, zp, x_migrate, y_migrate, z_migrate, x_prolif, y_prolif, z_prolif;
	int m[26];
	int rule_order[3];
	rule_order[0] = 0;
	rule_order[1] = 1;
	rule_order[2] = 2;
	lx = 100;
	ly = 100;
	lz = 100;
	
	// Create dummy seed data.
	int i;
	int * id_array;
	id_array = (int*)malloc(lx*ly*lz * sizeof(int));  //memory allocated using malloc
	
	for (i = 0; i < lx*ly*lz; i++) {
			id_array[i] = 0;
	}
	
	id_array[50+50*lz+50*ly*lz] = 2;
	
	// setup grid dimensions
	int my_xmin, my_xmax;
	my_xmin = rank * round(lx / size);
	my_xmax = (rank + 1) * round(lx / size);
	if (rank == size - 1) {
		if (my_xmax != lx) {
			my_xmax = lx;
		}
	}
	/* Set up Neighbors */
	int prev, next, x_expansion, add_to_x;
	int my_lx = my_xmax - my_xmin;
	double dx = 4.0 / ((double)(lx - 1));
	double dt = dx;
	if (rank == 0 && rank != size - 1) 
	{ 	// initialization specific to first core
		prev = -1;
		next = rank + 1;
		add_to_x = 0;
		x_expansion = 1;
	}
	else if (rank != 0 && rank == size - 1) 
	{ 	// initialization specific to last core
		prev = rank - 1;
		next = -1;
		add_to_x = 1;
		x_expansion = 1;
	}
	else if (rank != 0 && rank != size - 1) 
	{ 	// initialization specific to middle cores
		prev = rank - 1;
		next = rank + 1;
		add_to_x = 1;
		x_expansion = 2;
	}
	else 
	{ 	// serial case
		prev = -1;
		next = -1;
		add_to_x = 0;
		x_expansion = 0;
	}
	
	//Allocate game of life matrices
	matrix_a = malloc((my_lx+2)*sizeof(int**));
	matrix = malloc((my_lx+2)*sizeof(int**));
	
	for(ii=0; ii<my_lx+2; ii++){
		matrix_a[ii] = malloc((ly)*sizeof(int*));
		matrix[ii] = malloc((ly)*sizeof(int*));
		for (kk = 0; kk < ly; kk++) {
			matrix_a[ii][kk] = malloc((lz)*sizeof(int));
			matrix[ii][kk] = malloc((lz)*sizeof(int));
		}
	}	
	
	// Initialize matrices to 0
	for(ii = 0; ii <my_lx+2; ii++){
		for(jj = 0; jj < ly; jj++){
			for (kk = 0; kk < lz; kk++) {
				matrix_a[ii][jj][kk] = 0;
				matrix[ii][jj][kk] = 0;
			}
		}
	}
	
	// Initialize locations of live cells
	for(ii=1;ii<my_lx+1;ii++)
	{	// Leave one slice on each side in x-direction for ghost cells
		for(jj=1;jj<ly-1;jj++)
		{	// Similarly to above
			for(kk=1;kk<lz-1;kk++)
			{	// Similarly to above
				if (id_array[(kk) + (jj)*lz + ((ii-1) + my_xmin)*ly*lz] == 2)
				{
					matrix_a[ii][jj][kk] = 1;
					matrix[ii][jj][kk] = 1;
					liveCells++;
				}
			}
		}
	}
	
	// Initialize matrices to hold neighbor information and communication buffers
	available_x = malloc(26*sizeof(int));
	available_y = malloc(26*sizeof(int));
	available_z = malloc(26*sizeof(int));
	available_indices = malloc(26*sizeof(int));
	
	int **fsend_buf = alloc2DContigInt(ly, lz);
	int **bsend_buf = alloc2DContigInt(ly, lz);
	int **frecv_buf = alloc2DContigInt(ly, lz);
	int **brecv_buf = alloc2DContigInt(ly, lz);
	
	int **fsend_upd = alloc2DContigInt(ly, lz);
	int **bsend_upd = alloc2DContigInt(ly, lz);
	int **frecv_upd = alloc2DContigInt(ly, lz);
	int **brecv_upd = alloc2DContigInt(ly, lz);
	
	/* Set cell parameters for simulation */
	double r_migr = 0.5;
	double r_prolif = 0.95;
	double r_death = 0.13;
	srand(rank*time(NULL));
	
	// setup time
	int n_step;
	int n_step_max = 100;
	int n_step_out = 100;
	
	for (n_step = 1; n_step <= n_step_max; n_step++) 
	{
		// If desired, can make cell cycle periods longer than a single time step. Otherwise, set modulo value to 1.
		if (n_step % 2 == 0)
		{
			for(ii = 0;ii < ly;ii++)
			{
				for(jj=0; jj<lz; jj++)
				{
					// Write information on the edges of the x-direction to communication buffers
					bsend_buf[ii][jj] = matrix_a[0][ii][jj];
					fsend_buf[ii][jj] = matrix_a[my_lx-1][ii][jj];
				}
			}
			
			int count = 0;
			MPI_Request reqs[8];
			MPI_Status status;
			
			if (size>1)
			{	// Send all messages to appropriate neighbors. The end of the domain sends to the beginning
				// of the domain, and vice versa.
				if (prev != -1)
				{
					MPI_Isend(&bsend_buf[0][0], ly*lz, MPI_INT, prev, 0, MPI_COMM_WORLD, reqs + count);
					count = count + 1;
				}
				
				if (next != -1)
				{
					MPI_Isend(&fsend_buf[0][0], ly*lz, MPI_INT, next, 0, MPI_COMM_WORLD, reqs + count);
					count = count + 1;
				}
				
				if (prev != -1)
				{
					MPI_Irecv(&brecv_buf[0][0], ly*lz, MPI_INT, prev, 0, MPI_COMM_WORLD, reqs + count);
					count = count + 1;
				}
				
				if (next != -1)
				{
					MPI_Irecv(&frecv_buf[0][0], ly*lz, MPI_INT, next, 0, MPI_COMM_WORLD, reqs + count);
					count = count + 1;
				}
				
				if (prev == -1)
				{
					MPI_Isend(&bsend_buf[0][0], ly*lz, MPI_INT, size-1, 0, MPI_COMM_WORLD, reqs + count);
					count = count + 1;
				}
				
				if (next == -1)
				{
					MPI_Isend(&fsend_buf[0][0], ly*lz, MPI_INT, 0, 0, MPI_COMM_WORLD, reqs + count);
					count = count + 1;
				}
				
				if (prev == -1)
				{
					MPI_Irecv(&brecv_buf[0][0], ly*lz, MPI_INT, size-1, 0, MPI_COMM_WORLD, reqs + count);
					count = count + 1;
				}
				
				if (next == -1)
				{
					MPI_Irecv(&frecv_buf[0][0], ly*lz, MPI_INT, 0, 0, MPI_COMM_WORLD, reqs + count);
					count = count + 1;
				}
				
				if (count > 0) 
				{
					MPI_Waitall(count, reqs, MPI_STATUSES_IGNORE);
				}
				
				for(ii=0; ii<ly; ii++){
					for (jj = 0; jj < lz; jj++) {
						// Fill out ghost slices with buffer information.
						matrix_a[0][ii][jj] = brecv_buf[ii][jj];
						matrix_a[my_lx+1][ii][jj] = frecv_buf[ii][jj];
					}
				}
			}
			else
			{
				for(ii=0; ii<ly; ii++){
					for (jj = 0; jj < lz; jj++) {
						// If running serially, use periodic boundary condition in x.
						matrix_a[0][ii][jj] = matrix_a[my_lx][ii][jj];
						matrix_a[my_lx+1][ii][jj] = matrix_a[1][ii][jj];
					}
				}
			}
		
			/* Z boundary conditions */
			for(ii=0; ii<my_lx+2; ii++){
				for (jj = 0; jj < ly; jj++) {
					// Periodic boundary condition in z.
					matrix_a[ii][jj][0] = matrix_a[ii][jj][lz-2];
					matrix_a[ii][jj][lz-1] = matrix_a[ii][jj][1];
				}
			}
			
			/* Y boundary conditions */
			for(ii = 0; ii < my_lx+2; ii++) {
				for (jj = 0; jj < lz; jj++) {
					// Periodic boundary condition in y.
					matrix_a[ii][0][jj] = matrix_a[ii][ly-2][jj];
					matrix_a[ii][ly-1][jj] = matrix_a[ii][1][jj];
				}
			}
			
			/* Random Ordering */
			// Count live cells.
			liveCells = 0;
			for (ii = 1; ii < my_lx+1; ii++) {
				for (jj = 1; jj < ly-1; jj++) {
					for (kk = 1; kk < lz-1; kk++) {
						if (matrix_a[ii][jj][kk] == 1) {
							liveCells++;
						}
					}
				}
			}
			
			// This was a tricky one. At first I was looping through liveCells on every 
			// process, but I realized that not all processes would have the same
			// number of live cells. Therefore, some processes were not entering the part of the
			// for loop below that handles communication and the program was getting stuck.
			// So I used MPI Reduce to get the max number of cells on any process and broadcasted 
			// that number to every rank. Then I loop through using max cells as the limit for the
			// for loop and for each process, only do stuff if n, the loop index, is less than the number
			// of live cells for that process.
			MPI_Reduce(&liveCells,&maxCells,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
			MPI_Bcast(&maxCells,1,MPI_INT,0,MPI_COMM_WORLD);
			
			// Allocate arrays according to the number of live cells. 
			matrix_addresses_x = malloc(liveCells * sizeof(int*));
			matrix_addresses_y = malloc(liveCells * sizeof(int*));
			matrix_addresses_z = malloc(liveCells * sizeof(int*));
			index_list = malloc(liveCells * sizeof(int*));
			printf("There are %d live cells.\n",liveCells);
			
			// Store the addresses of all live cells and their indices so they can be shuffled.
			index = 0;
			for (ii = 1; ii < my_lx+1; ii++) {
				for (jj = 1; jj < ly-1; jj++) {
					for (kk = 1; kk < lz-1; kk++) {
						if (matrix_a[ii][jj][kk] == 1) {
							// Create lists of grid indices
							matrix_addresses_x[index] = ii;
							matrix_addresses_y[index] = jj;
							matrix_addresses_z[index] = kk;
							index_list[index] = index;
							index++;
						}
					}
				}
			}
			
			// Shuffle cell indices
			shuffle(index_list,index-1);
			
			for(n=0;n<maxCells;n++)
			{
				if (n<liveCells)
				{
					done = 0;
					// Get n-th index of shuffled grid. This ensures we check all live cells in random order.
					check = index_list[n];
					xx = matrix_addresses_x[check];
					yy = matrix_addresses_y[check];
					zz = matrix_addresses_z[check];
					if(matrix_a[xx][yy][zz] == 1)
					{
						// Define neighbors
						xm = xx - 1;
						xp = xx + 1;
						ym = yy - 1;
						yp = yy + 1;
						zm = zz - 1;
						zp = zz + 1;
						
						// Get number of live cells in immediate vicinity.
						nsum = matrix_a[xm][ym][zm] + matrix_a[xm][ym][zz] + matrix_a[xm][ym][zp] + matrix_a[xm][yy][zm] + matrix_a[xm][yy][zz] 
						+ matrix_a[xm][yy][zp] + matrix_a[xm][yp][zm] + matrix_a[xm][yp][zz] + matrix_a[xm][yp][zp] + matrix_a[xx][ym][zm]
						+ matrix_a[xx][ym][zz] + matrix_a[xx][ym][zp] + matrix_a[xx][yy][zm] + matrix_a[xx][yy][zp] + matrix_a[xx][yp][zm]
						+ matrix_a[xx][yp][zz] + matrix_a[xx][yp][zp] + matrix_a[xp][ym][zm] + matrix_a[xp][ym][zz] + matrix_a[xp][ym][zp]
						+ matrix_a[xp][yy][zm] + matrix_a[xp][yy][zz] + matrix_a[xp][yy][zp] + matrix_a[xp][yp][zm] + matrix_a[xp][yp][zz]
						+ matrix_a[xp][yp][zp];
						
						// Get number of grid points occupied by LBM in immediate vicinity.
						noccupied = 0;
						for(ii=xm+my_xmin;ii<=xp+my_xmin;ii++)
						{
							for(jj=ym;jj<=yp;jj++)
							{
								for(kk=zm;kk<=zp;kk++)
								{
									if (xx>1 && xx < my_lx+1)
									{
										if((id_array[kk + jj*lz + (ii-1)*ly*lz] == 1) && (ii != xx) && (jj !=yy) && (kk != zz) && (jj != 0) )
										{
											noccupied++;
										}
									}
								}
							}
						}
						
						
						// Decide if there is space for division / migration
						if ((nsum + noccupied) >= 25)
						{
							// printf("There is no space. nsum = %d, noccupied = %d\n",nsum,noccupied,nOB);
							space = 0;
						}
						else
						{
							// printf("There is space. nsum = %d, noccupied = %d\n",nsum,noccupied,nOB);
							space = 1;
						}
					}
					
					// Randomize order of checking
					// rule_order is an array containing 0, 1, and 2. Shuffle the matrix each time
					// so we check for proliferation, migration, and death in random orders for each cell.
					shuffle(rule_order,3);
					// Roll RNG for rates.
					rd = rand2();
					rm = rand2();
					rp = rand2();
				}
				
				for(ii=0;ii<ly;ii++)
				{
					for(jj=0;jj<lz;jj++)
					{
						// Initialize matrix update communication buffers to sentinel value.
						bsend_upd[ii][jj] = 1738;
						fsend_upd[ii][jj] = 1738;
					}
				}
				
				if (n<liveCells)
				{
					// There is a series of 3 of these switch statements. On each, either 0, 1 or 2 will be
					// chosen. The order of these has been randomized above. Once an action has been taken by the cell
					// no further action can be taken that time step for that cell.
					switch(rule_order[0]) // Trigger the first check
					{
						case 2: // Check for migration
						if((space == 1) && (done == 0))
						{	// If there's space, and no action has been taken yet...
							if (rm < r_migr)
							{	// If the RNG value is less than the rate of migration...
								done = 1; // An action has now been taken
								freeCount = 0;
								for(ii=xm;ii<=xp;ii++)
								{
									for(jj=ym;jj<=yp;jj++)
									{
										for(kk=zm;kk<=zp;kk++)
										{
											if((id_array[kk + jj*lz + ((ii-1)+my_xmin)*ly*lz] == 0))
											{
												if (matrix_a[ii][jj][kk] == 0)
												{	
													// Get indices of available neighboring spaces.
													available_x[freeCount] = ii;
													available_y[freeCount] = jj;
													available_z[freeCount] = kk;
													freeCount++;
												}
											}
										}
									}
								}
								
								// Randomly choose one to migrate to.
								migration_target = rand() % freeCount;
								x_migrate = available_x[migration_target];
								y_migrate = available_y[migration_target];
								z_migrate = available_z[migration_target];
								
								if (size>1)
								{
									// If running in parallel...
									if(x_migrate < 0 )
									{
										// Write to backward sending buffer if cell migrates out of bounds toward 0.
										bsend_upd[y_migrate][z_migrate] = 1;
									}
									
									else if(x_migrate > my_lx)
									{
										// Write to forward sending buffer if cell migrates out of bounds toward lx.
										fsend_upd[y_migrate][z_migrate] = 1;
									}
									else
									{
										matrix[x_migrate][y_migrate][z_migrate] = 1;
									}
								}
								else
								{
									// If running serially...
									if (x_migrate < 0)
									{
										// Pass back to the end of the domain if cell migrates out of bounds toward 0.
										matrix[my_lx][y_migrate][z_migrate] = 1;
									}
									else if(x_migrate > 0)
									{
										// Pass back to the beginning of the domain if cell migrates out of bounds toward lx.
										matrix[1][y_migrate][z_migrate] = 1;
									}
									else
									{
										matrix[x_migrate][y_migrate][z_migrate] = 1;
									}
								}
								// Cell has left current spot. Set it to zero.
								matrix[xx][yy][zz] = 0;
							}
						}
					
					
						case 1: // Check for death
						if(done==0)
						{	// Check if action has been taken.
							if(rd < r_death)
							{
								// Rather self explanatory.
								done = 1; // An action has been taken.
								matrix[xx][yy][zz] = 0;
							}
						}

					
						case 0: // Check for proliferation
						if((space == 1) && (done == 0))
						{
							if (rm < r_prolif)
							{	// Similarly to above, RNG check.
								done = 1;
								freeCount = 0;
								for(ii=xm;ii<=xp;ii++)
								{
									for(jj=ym;jj<=yp;jj++)
									{
										for(kk=zm;kk<=zp;kk++)
										{
											if((id_array[kk + jj*lz + ((ii-1)+my_xmin)*ly*lz] == 0))
											{
												if (matrix_a[ii][jj][kk] == 0)
												{
													available_x[freeCount] = ii;
													available_y[freeCount] = jj;
													available_z[freeCount] = kk;
													freeCount++;
												}
											}
										}
									}
								}
								
								// Check available spaces and randomly choose one.
							
								prolif_target = rand() % freeCount;
								x_prolif = available_x[prolif_target];
								y_prolif = available_y[prolif_target];
								z_prolif = available_z[prolif_target];
								
								if (size>1)
								{
									// Write to communication if out of bounds in parallel.
									if(x_prolif < 0 )
									{
										bsend_upd[y_prolif][z_prolif] = 1;
									}
									
									else if(x_migrate > my_lx)
									{
										fsend_upd[y_prolif][z_prolif] = 1;
									}
									else
									{
										matrix[x_prolif][y_prolif][z_prolif] = 1;
									}
								}
								else
								{
									// If running serially, periodic boundary conditions.
									if (x_migrate < 0)
									{
										matrix[my_lx][y_prolif][z_prolif] = 1;
									}
									else if(x_migrate > 0)
									{
										matrix[1][y_prolif][z_prolif] = 1;
									}
									else
									{
										matrix[x_prolif][y_prolif][z_prolif] = 1;
									}
								}
							}
						}
					}
				}
				
				if(size>1)
				{
					// Handle all communications.
					count = 0;
					MPI_Request reqs[8];
					MPI_Status status;
					
					if (prev != -1)
					{
						MPI_Isend(&bsend_upd[0][0], ly*lz, MPI_INT, prev, 0, MPI_COMM_WORLD, reqs + count);
						count = count + 1;
					}
					
					if (next != -1)
					{
						MPI_Isend(&fsend_upd[0][0], ly*lz, MPI_INT, next, 0, MPI_COMM_WORLD, reqs + count);
						count = count + 1;
					}
					
					if (prev != -1)
					{
						MPI_Irecv(&brecv_upd[0][0], ly*lz, MPI_INT, prev, 0, MPI_COMM_WORLD, reqs + count);
						count = count + 1;
					}
					
					if (next != -1)
					{
						MPI_Irecv(&frecv_upd[0][0], ly*lz, MPI_INT, next, 0, MPI_COMM_WORLD, reqs + count);
						count = count + 1;
					}
					
					if (prev == -1)
					{
						MPI_Isend(&bsend_upd[0][0], ly*lz, MPI_INT, size-1, 0, MPI_COMM_WORLD, reqs + count);
						count = count + 1;
					}
					
					if (next == -1)
					{
						MPI_Isend(&fsend_upd[0][0], ly*lz, MPI_INT, 0, 0, MPI_COMM_WORLD, reqs + count);
						count = count + 1;
					}
					
					if (prev == -1)
					{
						MPI_Irecv(&brecv_upd[0][0], ly*lz, MPI_INT, size-1, 0, MPI_COMM_WORLD, reqs + count);
						count = count + 1;
					}
					
					if (next == -1)
					{
						MPI_Irecv(&frecv_upd[0][0], ly*lz, MPI_INT, 0, 0, MPI_COMM_WORLD, reqs + count);
						count = count + 1;
					}
					
					if (count > 0) 
					{
						MPI_Waitall(count, reqs, MPI_STATUSES_IGNORE);
					}
					
					// Update "update matrix" with information from communications.
					for(ii=1;ii<ly-1;ii++)
					{
						for(jj=1;jj<lz-1;jj++)
						{
							if(brecv_upd[ii][jj] != 1738)
							{
								matrix[1][ii][jj] = brecv_upd[ii][jj];
							}
							
							if(frecv_upd[ii][jj] != 1738)
							{
								matrix[my_lx][ii][jj] = frecv_upd[ii][jj];
							}
						}
					}
				}
				
				
				// Update final matrix
				for(ii=1; ii<my_lx+1; ii++)
				{
					for(jj=1; jj<ly-1; jj++)
					{
						for(kk=1; kk<lz-1; kk++)
						{
							matrix_a[ii][jj][kk] = matrix[ii][jj][kk];
						}
					}
				}
				
				// And we repeat.
				if (n<liveCells)
				{
					switch(rule_order[1]) // Trigger the second check
					{
						case 2: // Check for migration
						if((space == 1) && (done == 0))
						{
							if (rm < r_migr)
							{
								done = 1;
								freeCount = 0;
								for(ii=xm;ii<=xp;ii++)
								{
									for(jj=ym;jj<=yp;jj++)
									{
										for(kk=zm;kk<=zp;kk++)
										{
											if((id_array[kk + jj*lz + ((ii-1)+my_xmin)*ly*lz] == 0))
											{
												if (matrix_a[ii][jj][kk] == 0)
												{
													
													available_x[freeCount] = ii;
													available_y[freeCount] = jj;
													available_z[freeCount] = kk;
													freeCount++;
												}
											}
										}
									}
								}
								
							
								migration_target = rand() % freeCount;
								x_migrate = available_x[migration_target];
								y_migrate = available_y[migration_target];
								z_migrate = available_z[migration_target];
								
								if (size>1)
								{
									if(x_migrate < 0 )
									{
										bsend_upd[y_migrate][z_migrate] = 1;
									}
									
									else if(x_migrate > my_lx)
									{
										fsend_upd[y_migrate][z_migrate] = 1;
									}
									else
									{
										matrix[x_migrate][y_migrate][z_migrate] = 1;
									}
								}
								else
								{
									if (x_migrate < 0)
									{
										matrix[my_lx][y_migrate][z_migrate] = 1;
									}
									else if(x_migrate > 0)
									{
										matrix[1][y_migrate][z_migrate] = 1;
									}
									else
									{
										matrix[x_migrate][y_migrate][z_migrate] = 1;
									}
								}
								matrix[xx][yy][zz] = 0;
							}
						}
					
					
						case 1: // Check for death
						if(done==0)
						{
							if(rd < r_death)
							{
								done = 1;
								matrix[xx][yy][zz] = 0;
							}
						}

					
						case 0: // Check for proliferation
						if((space == 1) && (done == 0))
						{
							if (rm < r_prolif)
							{
								done = 1;
								freeCount = 0;
								for(ii=xm;ii<=xp;ii++)
								{
									for(jj=ym;jj<=yp;jj++)
									{
										for(kk=zm;kk<=zp;kk++)
										{
											if((id_array[kk + jj*lz + ((ii-1)+my_xmin)*ly*lz] == 0))
											{
												if (matrix_a[ii][jj][kk] == 0)
												{
													available_x[freeCount] = ii;
													available_y[freeCount] = jj;
													available_z[freeCount] = kk;
													freeCount++;
												}
											}
										}
									}
								}
						
								prolif_target = rand() % freeCount;
								x_prolif = available_x[prolif_target];
								y_prolif = available_y[prolif_target];
								z_prolif = available_z[prolif_target];
								
								if (size>1)
								{
									if(x_prolif < 0 )
									{
										bsend_upd[y_prolif][z_prolif] = 1;
									}
									
									else if(x_migrate > my_lx)
									{
										fsend_upd[y_prolif][z_prolif] = 1;
									}
									else
									{
										matrix[x_prolif][y_prolif][z_prolif] = 1;
									}
								}
								else
								{
									if (x_migrate < 0)
									{
										matrix[my_lx][y_prolif][z_prolif] = 1;
									}
									else if(x_migrate > 0)
									{
										matrix[1][y_prolif][z_prolif] = 1;
									}
									else
									{
										matrix[x_prolif][y_prolif][z_prolif] = 1;
									}
								}
							}
						}
					}
				}
				
				if(size>1)
				{
					count = 0;
					MPI_Request reqs[8];
					MPI_Status status;
					
					if (prev != -1)
					{
						MPI_Isend(&bsend_upd[0][0], ly*lz, MPI_INT, prev, 0, MPI_COMM_WORLD, reqs + count);
						count = count + 1;
					}
					
					if (next != -1)
					{
						MPI_Isend(&fsend_upd[0][0], ly*lz, MPI_INT, next, 0, MPI_COMM_WORLD, reqs + count);
						count = count + 1;
					}
					
					if (prev != -1)
					{
						MPI_Irecv(&brecv_upd[0][0], ly*lz, MPI_INT, prev, 0, MPI_COMM_WORLD, reqs + count);
						count = count + 1;
					}
					
					if (next != -1)
					{
						MPI_Irecv(&frecv_upd[0][0], ly*lz, MPI_INT, next, 0, MPI_COMM_WORLD, reqs + count);
						count = count + 1;
					}
					
					if (prev == -1)
					{
						MPI_Isend(&bsend_upd[0][0], ly*lz, MPI_INT, size-1, 0, MPI_COMM_WORLD, reqs + count);
						count = count + 1;
					}
					
					if (next == -1)
					{
						MPI_Isend(&fsend_upd[0][0], ly*lz, MPI_INT, 0, 0, MPI_COMM_WORLD, reqs + count);
						count = count + 1;
					}
					
					if (prev == -1)
					{
						MPI_Irecv(&brecv_upd[0][0], ly*lz, MPI_INT, size-1, 0, MPI_COMM_WORLD, reqs + count);
						count = count + 1;
					}
					
					if (next == -1)
					{
						MPI_Irecv(&frecv_upd[0][0], ly*lz, MPI_INT, 0, 0, MPI_COMM_WORLD, reqs + count);
						count = count + 1;
					}
					
					if (count > 0) 
					{
						MPI_Waitall(count, reqs, MPI_STATUSES_IGNORE);
					}
					
					for(ii=1;ii<ly-1;ii++)
					{
						for(jj=1;jj<lz-1;jj++)
						{
							if(brecv_upd[ii][jj] != 1738)
							{
								matrix[1][ii][jj] = brecv_upd[ii][jj];
							}
							
							if(frecv_upd[ii][jj] != 1738)
							{
								matrix[my_lx][ii][jj] = frecv_upd[ii][jj];
							}
						}
					}
				}
				
				
				// Update matrix
				for(ii=1; ii<my_lx+1; ii++)
				{
					for(jj=1; jj<ly-1; jj++)
					{
						for(kk=1; kk<lz-1; kk++)
						{
							matrix_a[ii][jj][kk] = matrix[ii][jj][kk];
						}
					}
				}
				// And again for the last time.
				
				if (n<liveCells)
				{
					switch(rule_order[2]) // Trigger the third check
					{
						case 2: // Check for migration
						if((space == 1) && (done == 0))
						{
							if (rm < r_migr)
							{
								done = 1;
								freeCount = 0;
								for(ii=xm;ii<=xp;ii++)
								{
									for(jj=ym;jj<=yp;jj++)
									{
										for(kk=zm;kk<=zp;kk++)
										{
											if((id_array[kk + jj*lz + ((ii-1)+my_xmin)*ly*lz] == 0))
											{
												if (matrix_a[ii][jj][kk] == 0)
												{
													
													available_x[freeCount] = ii;
													available_y[freeCount] = jj;
													available_z[freeCount] = kk;
													freeCount++;
												}
											}
										}
									}
								}
							
								migration_target = rand() % freeCount;
								x_migrate = available_x[migration_target];
								y_migrate = available_y[migration_target];
								z_migrate = available_z[migration_target];
								
								if (size>1)
								{
									if(x_migrate < 0 )
									{
										bsend_upd[y_migrate][z_migrate] = 1;
									}
									
									else if(x_migrate > my_lx)
									{
										fsend_upd[y_migrate][z_migrate] = 1;
									}
									else
									{
										matrix[x_migrate][y_migrate][z_migrate] = 1;
									}
								}
								else
								{
									if (x_migrate < 0)
									{
										matrix[my_lx][y_migrate][z_migrate] = 1;
									}
									else if(x_migrate > 0)
									{
										matrix[1][y_migrate][z_migrate] = 1;
									}
									else
									{
										matrix[x_migrate][y_migrate][z_migrate] = 1;
									}
								}
								matrix[xx][yy][zz] = 0;
							}
						}
					
					
						case 1: // Check for death
						if(done==0)
						{
							if(rd < r_death)
							{
								done = 1;
								matrix[xx][yy][zz] = 0;
							}
						}

					
						case 0: // Check for proliferation
						if((space == 1) && (done == 0))
						{
							if (rm < r_prolif)
							{
								done = 1;
								freeCount = 0;
								for(ii=xm;ii<=xp;ii++)
								{
									for(jj=ym;jj<=yp;jj++)
									{
										for(kk=zm;kk<=zp;kk++)
										{
											if((id_array[kk + jj*lz + ((ii-1)+my_xmin)*ly*lz] == 0))
											{
												if (matrix_a[ii][jj][kk] == 0)
												{
													available_x[freeCount] = ii;
													available_y[freeCount] = jj;
													available_z[freeCount] = kk;
													freeCount++;
												}
											}
										}
									}
								}
						
								prolif_target = rand() % freeCount;
								x_prolif = available_x[prolif_target];
								y_prolif = available_y[prolif_target];
								z_prolif = available_z[prolif_target];
								
								if (size>1)
								{
									if(x_prolif < 0 )
									{
										bsend_upd[y_prolif][z_prolif] = 1;
									}
									
									else if(x_migrate > my_lx)
									{
										fsend_upd[y_prolif][z_prolif] = 1;
									}
									else
									{
										matrix[x_prolif][y_prolif][z_prolif] = 1;
									}
								}
								else
								{
									if (x_migrate < 0)
									{
										matrix[my_lx][y_prolif][z_prolif] = 1;
									}
									else if(x_migrate > 0)
									{
										matrix[1][y_prolif][z_prolif] = 1;
									}
									else
									{
										matrix[x_prolif][y_prolif][z_prolif] = 1;
									}
								}
							}
						}
					}
				}
				
				if(size>1)
				{
					count = 0;
					MPI_Request reqs[8];
					MPI_Status status;
					
					if (prev != -1)
					{
						MPI_Isend(&bsend_upd[0][0], ly*lz, MPI_INT, prev, 0, MPI_COMM_WORLD, reqs + count);
						count = count + 1;
					}
					
					if (next != -1)
					{
						MPI_Isend(&fsend_upd[0][0], ly*lz, MPI_INT, next, 0, MPI_COMM_WORLD, reqs + count);
						count = count + 1;
					}
					
					if (prev != -1)
					{
						MPI_Irecv(&brecv_upd[0][0], ly*lz, MPI_INT, prev, 0, MPI_COMM_WORLD, reqs + count);
						count = count + 1;
					}
					
					if (next != -1)
					{
						MPI_Irecv(&frecv_upd[0][0], ly*lz, MPI_INT, next, 0, MPI_COMM_WORLD, reqs + count);
						count = count + 1;
					}
					
					if (prev == -1)
					{
						MPI_Isend(&bsend_upd[0][0], ly*lz, MPI_INT, size-1, 0, MPI_COMM_WORLD, reqs + count);
						count = count + 1;
					}
					
					if (next == -1)
					{
						MPI_Isend(&fsend_upd[0][0], ly*lz, MPI_INT, 0, 0, MPI_COMM_WORLD, reqs + count);
						count = count + 1;
					}
					
					if (prev == -1)
					{
						MPI_Irecv(&brecv_upd[0][0], ly*lz, MPI_INT, size-1, 0, MPI_COMM_WORLD, reqs + count);
						count = count + 1;
					}
					
					if (next == -1)
					{
						MPI_Irecv(&frecv_upd[0][0], ly*lz, MPI_INT, 0, 0, MPI_COMM_WORLD, reqs + count);
						count = count + 1;
					}
					
					if (count > 0) 
					{
						MPI_Waitall(count, reqs, MPI_STATUSES_IGNORE);
					}
					
					for(ii=1;ii<ly-1;ii++)
					{
						for(jj=1;jj<lz-1;jj++)
						{
							if(brecv_upd[ii][jj] != 1738)
							{
								matrix[1][ii][jj] = brecv_upd[ii][jj];
							}
							
							if(frecv_upd[ii][jj] != 1738)
							{
								matrix[my_lx][ii][jj] = frecv_upd[ii][jj];
							}
						}
					}
				}
				
				
				// Update matrix
				for(ii=1; ii<my_lx+1; ii++)
				{
					for(jj=1; jj<ly-1; jj++)
					{
						for(kk=1; kk<lz-1; kk++)
						{
							matrix_a[ii][jj][kk] = matrix[ii][jj][kk];
						}
					}
				}
				
				
			}
			// The address matrices must change size every time, because the number of live cells
			// will presumably change at every time step. Therefore, I free at the end of every iteration of
			// the time loop (or cell cycle period, if so inclined), and reallocate at the beginning of
			// the next iteration.
			free(matrix_addresses_x);
			free(matrix_addresses_y);
			free(matrix_addresses_z);
			free(index_list);
		
		}
	}
	
	// Write outputs
	if(size==1)
	{
		FILE *vtk_fp;
		char fname[15];
		sprintf(fname, "tumor_%d.vtk", rank);
		vtk_fp = fopen(fname,"w");
		
		if (vtk_fp == NULL) {
			fprintf(stderr, "Error: cannot open file life.vtk. \n");
			return 1;
		}
		fprintf(vtk_fp, "# vtk DataFile Version 3.0");
		fprintf(vtk_fp, "\nvtk global output\n");
		fprintf(vtk_fp, "ASCII\n");
		fprintf(vtk_fp, "DATASET STRUCTURED_POINTS\n");
		fprintf(vtk_fp, "DIMENSIONS %i %i %i\n", lx, ly, lz);
		fprintf(vtk_fp, "ORIGIN 0. 0. 0.\n");
		fprintf(vtk_fp, "SPACING 1 1 1\n");
		fprintf(vtk_fp, "\nPOINT_DATA %d\n",lx*ly*lz);
		fprintf(vtk_fp, "SCALARS state double \n");
		fprintf(vtk_fp, "LOOKUP_TABLE default\n");
		
		for(kk = 0; kk < lz; kk++) {
			for(jj = 0; jj < ly; jj++) {
				for(ii = 1; ii < my_lx+1; ii++) {

					if(matrix_a[ii][jj][kk]==1)
					{
						
						fprintf(vtk_fp, "%d\n",1);
					}
					else
					{
						fprintf(vtk_fp, "%d\n",0);
					}
				}
			}
		}

		fprintf(vtk_fp, "VECTORS coordinates double \n");
		for(kk = 0; kk < lz; kk++) {
			for(jj = 0; jj < ly; jj++) {
				for(ii = 1; ii < my_lx+1; ii++) {

					fprintf(vtk_fp, "%f %f %f\n",(ii-1)*dx,jj*dx,kk*dx);
				}
			}
		}
		fclose(vtk_fp);
	}
	
	// Finalize.
	MPI_Finalize();
	return 0;
}
