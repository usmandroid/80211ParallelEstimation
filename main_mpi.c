#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "inputs.h"
#include "utils.h"

void WiFi_channel_estimation_LT_LS(long double complex tx_pre[], long double complex rx_pre[],long double complex H_EST[]);
void WiFi_channel_estimation_PS_Linear(long double complex tx_symbols[], long double complex rx_symbols[],long double complex H_EST[]);
void WiFi_channel_estimation_PS_Cubic(long double complex tx_symbols[], long double complex rx_symbols[],long double complex H_EST[]);
void WiFi_channel_estimation_PS_Sinc(long double complex tx_symbols[], long double complex rx_symbols[],long double complex H_EST[]);
void WiFi_channel_estimation_PS_MMSE(long double complex tx_symbols[], long double complex rx_symbols[], long double complex **F, double ow2, long double complex H_EST_LS[], long double complex H_EST[]);

int main(int argc, char *argv[]) {
	int numprocs, rank, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Get_processor_name(processor_name, &namelen);

	clock_t start, stop, start_tot, stop_tot;
	long double complex tx_symb_vec[SAMPUTIL],rx_symb_vec[SAMPUTIL];
	long double complex H_EST_LT_LS[SAMPUTIL],H_EST_PS_Linear[SAMPUTIL],H_EST_PS_Cubic[SAMPUTIL];
	long double complex H_EST_PS_Sinc[SAMPUTIL], H_EST_PS_MMSE[SAMPUTIL];
	int OFDM_block = 0;

	if(rank == 0){
		/* One OFDM Symbol isextracted to perform channel estimation */
		printf("Proc %d Processing Block %d\n", rank, OFDM_block);
		for(int r=0 ; r<SAMPUTIL ; r++){
			tx_symb_vec[r] = tx_symb[SAMPUTIL*OFDM_block + r];
			rx_symb_vec[r] = rx_symb[SAMPUTIL*OFDM_block + r];
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);	

	/* ----------------------------------------------------------------------------------------*/
	/* ----------------------------------- COMMON VARIABLES -----------------------------------*/
	/* ----------------------------------------------------------------------------------------*/

	long double complex conj1, conj2, H_PILOTS[4];
	long double complex f0, f01, f12, f23, f012, f123, f0123;
	long double complex myres;

	double complex sinc1[SAMPUTIL],sinc2[SAMPUTIL],sinc3[SAMPUTIL],sinc4[SAMPUTIL];

	long double res_LS[4];
	long double tx_real[4], tx_imag[4], rx_real[4], rx_imag[4];
	long double H_PILOTS_real[4],H_PILOTS_imag[4], res_Linear[6];
	long double alpha, delta = P1 - P0;
	long double f0_r, f0_i, f01_r, f01_i;
	long double f12_r, f12_i, f23_r, f23_i, f012_r, f012_i;
	long double f123_r, f123_i, f0123_r, f0123_i;
	long double res_part_real, res_part_imag, myres_real, myres_imag;
	
	int chunk[4], tag1 = 1, tag2 = 2, tag3 = 3, tag = 3;
	int k;

	double a, b, c, d;

	/* Communicators creation. 5 groups of 4 processes are created
	*  in order to perform a reduction operation. A communicator is
	*  assigned to each group, allowing them to share the results after
	*  the reduction.
	*/
	MPI_Comm comm1, comm2, comm3, comm4, comm5;
	MPI_Group gr1, gr2, gr3, gr4, gr5, orig_group;
	int *ranks = (int *) malloc(4*sizeof(rank));
	MPI_Comm_group(MPI_COMM_WORLD, &orig_group);
	MPI_Comm_group(MPI_COMM_WORLD, &orig_group);
	ranks[0]=0;ranks[1]=1;ranks[2]=2;ranks[3]=3;
	MPI_Group_incl(orig_group, 4, ranks, &gr1);
	MPI_Comm_create(MPI_COMM_WORLD, gr1, &comm1);
	ranks[0]=4;ranks[1]=5;ranks[2]=6;ranks[3]=7;
	MPI_Group_incl(orig_group, 4, ranks, &gr2);
	MPI_Comm_create(MPI_COMM_WORLD, gr2, &comm2);
	ranks[0]=8;ranks[1]=9;ranks[2]=10;ranks[3]=11;
	MPI_Group_incl(orig_group, 4, ranks, &gr3);
	MPI_Comm_create(MPI_COMM_WORLD, gr3, &comm3);
	ranks[0]=12;ranks[1]=13;ranks[2]=14;ranks[3]=15;
	MPI_Group_incl(orig_group, 4, ranks, &gr4);
	MPI_Comm_create(MPI_COMM_WORLD, gr4, &comm4);
	ranks[0]=16;ranks[1]=17;ranks[2]=18;ranks[3]=19;
	MPI_Group_incl(orig_group, 4, ranks, &gr5);
	MPI_Comm_create(MPI_COMM_WORLD, gr5, &comm5);
	free(ranks);


	if(rank == 0){
		long double complex tx_pilots[4] = {tx_symb_vec[P0],tx_symb_vec[P1],tx_symb_vec[P2],tx_symb_vec[P3]};
		long double complex rx_pilots[4] = {rx_symb_vec[P0],rx_symb_vec[P1],rx_symb_vec[P2],rx_symb_vec[P3]};

		for(int i=0 ; i<4 ; i++){
			H_PILOTS_real[i] = creal(rx_pilots[i] / tx_pilots[i]);
			H_PILOTS_imag[i] = cimag(rx_pilots[i] / tx_pilots[i]);
			H_PILOTS[i] = rx_pilots[i] / tx_pilots[i];
		}
	}

	MPI_Bcast(H_PILOTS_real, 4, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(H_PILOTS_imag, 4, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);

	for (int r = 0; r < 4; r++)
        H_PILOTS[r] = H_PILOTS_real[r] + I*H_PILOTS_imag[r];

	/* ----------------------------------------------------------------------------------------*/
	/* ----------------------------------------------------------------------------------------*/

	/* ----------------------------------------------------------------------------------------*/
	/* -------------------------------- LEAST SQUARE ESTIMATOR --------------------------------*/
	/* ----------------------------------------------------------------------------------------*/
	MPI_Barrier(MPI_COMM_WORLD);

	if(rank == 0){

		printf("Processing LT Least Square...\n"); start = clock();

		/* Processor 0 Distributes 27 tasks among the 20 processors.
		*  Proc 0-5: Perform 2 tasks. Proc 6-19: Perform 1 task.
		*/
		for(int dest=1 ; dest<numprocs; dest++){
			chunk[0] = dest;
			chunk[1] = dest + 27;
			if(dest<6){
				chunk[2] = dest + numprocs;
				chunk[3] = dest + numprocs + 27;
			}
			MPI_Send(chunk, 4, MPI_INT, dest, tag1, MPI_COMM_WORLD);
		}

		/* Gather results from workers */
		for (int src=1; src<numprocs; src++) {
			MPI_Recv(res_LS, 4, MPI_LONG_DOUBLE, src, tag1, MPI_COMM_WORLD, &status);
			H_EST_LT_LS[src] = res_LS[0] + I*res_LS[1];
			H_EST_LT_LS[src + 27] = res_LS[2] + I*res_LS[3];
			if(src < 6) {
				MPI_Recv(res_LS, 4, MPI_LONG_DOUBLE, src, tag2, MPI_COMM_WORLD, &status);
				H_EST_LT_LS[src + numprocs] = res_LS[0] + I*res_LS[1];
				H_EST_LT_LS[src + numprocs + 27] = res_LS[2] + I*res_LS[3];
			}
		}

		/* Process 0 does its processing part */
		conj1 = creal(tx_preamble_fft[0]) - cimag(tx_preamble_fft[0]);
		conj2 = creal(tx_preamble_fft[0+27]) - cimag(tx_preamble_fft[0+27]);
		H_EST_LT_LS[0] = ( conj1*rx_preamble_fft[0] ) / ( conj1*tx_preamble_fft[0] );
		H_EST_LT_LS[0+27] = ( conj1*rx_preamble_fft[0+27] ) / ( conj1*tx_preamble_fft[0+27] );
		H_EST_LT_LS[26] = 0.0;

	} else {

		MPI_Recv(chunk, 4, MPI_INT, 0, tag1, MPI_COMM_WORLD, &status);

		conj1 = creal(tx_preamble_fft[chunk[0]]) - cimag(tx_preamble_fft[chunk[0]]);
		conj2 = creal(tx_preamble_fft[chunk[1]]) - cimag(tx_preamble_fft[chunk[1]]);
		res_LS[0] = creal(( conj1*rx_preamble_fft[chunk[0]] ) / ( conj1*tx_preamble_fft[chunk[0]] ));
		res_LS[1] = cimag(( conj1*rx_preamble_fft[chunk[0]] ) / ( conj1*tx_preamble_fft[chunk[0]] ));
		res_LS[2] = creal(( conj2*rx_preamble_fft[chunk[1]] ) / ( conj2*tx_preamble_fft[chunk[1]] ) );
		res_LS[3] = cimag(( conj2*rx_preamble_fft[chunk[1]] ) / ( conj2*tx_preamble_fft[chunk[1]] ) );
		MPI_Send(res_LS, 4, MPI_LONG_DOUBLE, 0, tag1, MPI_COMM_WORLD);

		if(rank<6){
			conj1 = creal(tx_preamble_fft[chunk[2]]) - cimag(tx_preamble_fft[chunk[2]]);
			conj2 = creal(tx_preamble_fft[chunk[3]]) - cimag(tx_preamble_fft[chunk[3]]);
			res_LS[0] = creal(( conj1*rx_preamble_fft[chunk[2]] ) / ( conj1*tx_preamble_fft[chunk[2]] ));
			res_LS[1] = cimag(( conj1*rx_preamble_fft[chunk[2]] ) / ( conj1*tx_preamble_fft[chunk[2]] ));
			res_LS[2] = creal(( conj2*rx_preamble_fft[chunk[3]] ) / ( conj2*tx_preamble_fft[chunk[3]] ) );
			res_LS[3] = cimag(( conj2*rx_preamble_fft[chunk[3]] ) / ( conj2*tx_preamble_fft[chunk[3]] ) );
			MPI_Send(res_LS, 4, MPI_LONG_DOUBLE, 0, tag2, MPI_COMM_WORLD);
		}
	}

	/* End of LT Estimator */
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){
		stop = clock();
		printf("LT - Elapsed time %f\n",(double) (stop - start));
	}
	
	/* ----------------------------------------------------------------------------------------*/
	/* ----------------------------------------------------------------------------------------*/

	/* ----------------------------------------------------------------------------------------*/
	/* -------------------------------- PS LINEAR INTERPOLATION -------------------------------*/
	/* ----------------------------------------------------------------------------------------*/
	MPI_Barrier(MPI_COMM_WORLD);

	if(rank == 0){

		printf("Processing PS Linear Interpolation...\n"); start = clock(); 

		/* param1 is a matrix that contains the parameter being used in each processor
		*  for all the iteratiors. In this case, in order to cover 52 samples with 20 processors,
		*  we need 3 iterations in total (ceil(SAMPUTIL/numprocs)).
		*/
		int param1[3][19] = {
			{P0,P0,P0,P0,P0,P0,P0,P0,P0,P0,P0,P0,P0,P0,P0,P0,P0,P0,P0},
			{P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P1,P2,P2,P2,P2,P2},
			{P2,P2,P2,P2,P2,P2,P2,P2,P2,P3,P3,P3,P3,P3,P3,NULL,NULL,NULL,NULL} 
		};

		/* Proc 0 distributes work and data being used by all the processes */
		for(int dest=1 ; dest<numprocs ; dest++){
			MPI_Send(param1, numprocs*3, MPI_INT, dest, tag1, MPI_COMM_WORLD);
		}

		/* Waits and stores value */
		for(int src=1 ; src<numprocs ; src++){
			MPI_Recv(res_Linear, 6, MPI_LONG_DOUBLE, src, tag2, MPI_COMM_WORLD, &status);
			H_EST_PS_Linear[src] = res_Linear[0] + I*res_Linear[1];
			H_EST_PS_Linear[src + numprocs] = res_Linear[2] + I*res_Linear[3]; 	
			if(src < 16)
				H_EST_PS_Linear[src + 2*numprocs] = res_Linear[4] + I*res_Linear[5];
		}

		/* Proc0 does also its part of the job */
		alpha = (rank-P1)/delta;
		H_EST_PS_Linear[0] = (H_PILOTS_real[0]+( (H_PILOTS_real[1]-H_PILOTS_real[0] )*alpha )) + \
			I*(H_PILOTS_imag[0]+( (H_PILOTS_imag[1]-H_PILOTS_imag[0] )*alpha ));

	} else {
		long double alpha;
		int param1[3][19];

		/* Each processor knows the part of the work assigned to it */
		MPI_Recv(param1, numprocs*3, MPI_INT, 0, tag1, MPI_COMM_WORLD, &status);

		alpha = (rank-param1[0][rank-1])/delta;
		res_Linear[0] = H_PILOTS_real[0]+( (H_PILOTS_real[1]-H_PILOTS_real[0] )*alpha );
		res_Linear[1] = H_PILOTS_imag[0]+( (H_PILOTS_imag[1]-H_PILOTS_imag[0] )*alpha );
		alpha = ((rank+numprocs)-param1[1][rank-1])/delta;
		res_Linear[2] = H_PILOTS_real[0]+( (H_PILOTS_real[1]-H_PILOTS_real[0] )*alpha );
		res_Linear[3] = H_PILOTS_imag[0]+( (H_PILOTS_imag[1]-H_PILOTS_imag[0] )*alpha );
		if(rank < 16){
			alpha = ((rank+numprocs)-param1[2][rank-1])/delta;
			res_Linear[4] = H_PILOTS_real[0]+( (H_PILOTS_real[1]-H_PILOTS_real[0] )*alpha );
			res_Linear[5] = H_PILOTS_imag[0]+( (H_PILOTS_imag[1]-H_PILOTS_imag[0] )*alpha );
		}

		/* Send results to Processor 0 */
		MPI_Send(res_Linear, 6, MPI_LONG_DOUBLE, 0, tag2, MPI_COMM_WORLD);
	}	

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){
		stop = clock();
		printf("PS Linear Interpolation - Elapsed time %f\n",(double) (stop - start));
	}

	/* ----------------------------------------------------------------------------------------*/
	/* ----------------------------------------------------------------------------------------*/

	/* ----------------------------------------------------------------------------------------*/
	/* -------------------------------- PS CUBIC INTERPOLATION --------------------------------*
	*
	*  This function considers that there is 20 processes available. 5 different groups are created,
	*  each of them containing 4 processors. Each group will be assigned a value k, representing
	*  the index to compute. Each process within a group will compute 1 out of the 4 tasks that 
	*  each index requires to perform. Afterwars, a reduction operation will store the result on 
	*  the process whose rank is the lowest. Finally, this result will be sent to the master node.
	*  This process will repeat 11 times, so that the entire length of the vector is covered.
	*
	/* ----------------------------------------------------------------------------------------*/
	
	
	MPI_Barrier(MPI_COMM_WORLD);

	/* Each process uses some common variables, which need to be
	*  initialized by the master process.
	*/
	if(rank == 0){

		printf("Processing PS Cubic Interpolation...\n"); start = clock(); 

		f0_r = creal(H_PILOTS[0]);
		f0_i = cimag(H_PILOTS[0]);
		f01_r   = creal((H_PILOTS[1]-H_PILOTS[0]) / delta);
		f01_i   = cimag((H_PILOTS[1]-H_PILOTS[0]) / delta);
		f12_r   = creal((H_PILOTS[2]-H_PILOTS[1]) / delta);
		f12_i   = cimag((H_PILOTS[2]-H_PILOTS[1]) / delta);
	    f23_r   = creal((H_PILOTS[3]-H_PILOTS[2]) / delta);
	    f23_i   = cimag((H_PILOTS[3]-H_PILOTS[2]) / delta);
	    f012_r  = creal((f12_r-f01_r) / delta);
	    f012_i  = cimag((f12_i-f01_i) / delta);
	    f123_r  = creal((f23_r-f12_r) / delta);
	    f123_i  = cimag((f23_i-f12_i) / delta);
	    f0123_r = creal((f123_r-f012_r) / delta);
	    f0123_i = cimag((f123_i-f012_i) / delta);

	}

	/* The Common variables are Broadcasted to all the nodes by
	*  the master node 
	*/
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&f0_r, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&f0_i, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&f01_r, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&f01_i, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&f12_r, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&f12_i, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&f23_r, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&f23_i, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&f012_r, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&f012_i, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&f123_r, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&f123_i, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&f0123_r, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&f0123_i, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	f0 = f0_r + f0_i*I;
	f01 = f01_r + f01_i*I;
	f12 = f12_r + f12_i*I;
	f23 = f23_r + f23_i*I;
	f012 = f012_r + f012_i*I;
	f123 = f123_r + f123_i*I;
	f0123 = f0123_r + f0123_i*I;

	for(int i=0; i<11; i++){

		k = floor(rank/4) + 5*i;

		if(k<SAMPUTIL){

			if(rank%4 == 0){
				myres = f0;
			}else if(rank%4 == 1){
				myres = f01*(k-P0);
			}else if(rank%4 == 2){
				myres = f012*(k-P0)*(k-P1);
			}else if(rank%4 == 3){
				myres = f0123*(k-P0)*(k-P1)*(k-P2);
			}

			myres_real = creal(myres); myres_imag = cimag(myres);	
		} else {
			myres_real = 0.0; myres_imag = 0.0;
		}
		
		/* Each process sends the real and imaginary part the following way:
		*  tag = 0 for the Real and tag = 1 for the imag.
		*  Since we are using communicators, the ranks of the processes are local.
		*  The results are always sent to the process with the lowest rank. Afterwards,
		*  this processes will send the result to the master, who will store it.
		*/
		if(rank<4){
			MPI_Reduce(&myres_real,&res_part_real,1,MPI_LONG_DOUBLE,MPI_SUM,0,comm1);
			MPI_Reduce(&myres_imag,&res_part_imag,1,MPI_LONG_DOUBLE,MPI_SUM,0,comm1);
		} else if(rank < 8){
			MPI_Reduce(&myres_real,&res_part_real,1,MPI_LONG_DOUBLE,MPI_SUM,0,comm2);
			MPI_Reduce(&myres_imag,&res_part_imag,1,MPI_LONG_DOUBLE,MPI_SUM,0,comm2);
		} else if(rank < 12){
			MPI_Reduce(&myres_real,&res_part_real,1,MPI_LONG_DOUBLE,MPI_SUM,0,comm3);
			MPI_Reduce(&myres_imag,&res_part_imag,1,MPI_LONG_DOUBLE,MPI_SUM,0,comm3);
		} else if(rank < 16){
			MPI_Reduce(&myres_real,&res_part_real,1,MPI_LONG_DOUBLE,MPI_SUM,0,comm4);
			MPI_Reduce(&myres_imag,&res_part_imag,1,MPI_LONG_DOUBLE,MPI_SUM,0,comm4);
		} else{
			MPI_Reduce(&myres_real,&res_part_real,1,MPI_LONG_DOUBLE,MPI_SUM,0,comm5);
			MPI_Reduce(&myres_imag,&res_part_imag,1,MPI_LONG_DOUBLE,MPI_SUM,0,comm5);
		}

		if(rank == 0){
			// Rank 0 does its part of the task
			H_EST_PS_Cubic[i*5] = res_part_real + I*res_part_imag;

			// Rank 0 gathers the results of the rest of the master nodes
			for(int gr=1; gr<5; gr++){
				MPI_Recv(&res_part_real, 1, MPI_LONG_DOUBLE, gr*4, 0, MPI_COMM_WORLD, &status);
				MPI_Recv(&res_part_imag, 1, MPI_LONG_DOUBLE, gr*4, 1, MPI_COMM_WORLD, &status);
				H_EST_PS_Cubic[i*5 + gr] = res_part_real + I*res_part_imag;
			}
		} else if(rank%4==0 && rank!=0){
			MPI_Send(&res_part_real, 1, MPI_LONG_DOUBLE, 0, 0, MPI_COMM_WORLD);
			MPI_Send(&res_part_imag, 1, MPI_LONG_DOUBLE, 0, 1, MPI_COMM_WORLD);
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){
		stop = clock();
		printf("PS Cubic Interpolation - Elapsed time %f\n",(double) (stop - start));
	}

	/* ----------------------------------------------------------------------------------------*/
	/* ----------------------------------------------------------------------------------------*/
	
	/* ----------------------------------------------------------------------------------------*/
	/* -------------------------------- PS SINC INTERPOLATION ---------------------------------*
	*
	*  This function considers that there is 20 processes available. 5 different groups are created,
	*  each of them containing 4 processors. Each group will be assigned a value k, representing
	*  the index to compute. Each process within a group will compute 1 out of the 4 tasks that 
	*  each index requires to perform. Afterwars, a reduction operation will store the result on 
	*  the process whose rank is the lowest. Finally, this result will be sent to the master node.
	*  This process will repeat 11 times, so that the entire length of the vector is covered.
	*
	/* ----------------------------------------------------------------------------------------*/	

	MPI_Barrier(MPI_COMM_WORLD);

	if(rank == 0){
		printf("Processing PS Sinc Interpolation...\n"); start = clock(); 
	}

	for(int i=0; i<11; i++){

		k = floor(rank/4) + 5*i;

		if(k<SAMPUTIL){

			if(rank%4 == 0){
				a = (k-P0) / delta;
				myres = H_PILOTS[0]*sinc(a);
			}else if(rank%4 == 1){
				b = (k-P1) / delta;
				myres = H_PILOTS[1]*sinc(b);
			}else if(rank%4 == 2){
				c = (k-P2) / delta;
				myres = H_PILOTS[2]*sinc(c);
			}else if(rank%4 == 3){
				d = (k-P3) / delta;
				myres = H_PILOTS[3]*sinc(d);
			}

			myres_real = creal(myres); myres_imag = cimag(myres);
		} else {
			myres_real = 0.0; myres_imag = 0.0;
		}

		/* Each process sends the real and imaginary part the following way:
		*  tag = 0 for the Real and tag = 1 for the imag.
		*  Since we are using communicators, the ranks of the processes are local.
		*  The results are always sent to the process with the lowest rank. Afterwards,
		*  this processes will send the result to the master, who will store it.
		*/
		if(rank<4){
			MPI_Reduce(&myres_real,&res_part_real,1,MPI_LONG_DOUBLE,MPI_SUM,0,comm1);
			MPI_Reduce(&myres_imag,&res_part_imag,1,MPI_LONG_DOUBLE,MPI_SUM,0,comm1);
		} else if(rank < 8){
			MPI_Reduce(&myres_real,&res_part_real,1,MPI_LONG_DOUBLE,MPI_SUM,0,comm2);
			MPI_Reduce(&myres_imag,&res_part_imag,1,MPI_LONG_DOUBLE,MPI_SUM,0,comm2);
		} else if(rank < 12){
			MPI_Reduce(&myres_real,&res_part_real,1,MPI_LONG_DOUBLE,MPI_SUM,0,comm3);
			MPI_Reduce(&myres_imag,&res_part_imag,1,MPI_LONG_DOUBLE,MPI_SUM,0,comm3);
		} else if(rank < 16){
			MPI_Reduce(&myres_real,&res_part_real,1,MPI_LONG_DOUBLE,MPI_SUM,0,comm4);
			MPI_Reduce(&myres_imag,&res_part_imag,1,MPI_LONG_DOUBLE,MPI_SUM,0,comm4);
		} else{
			MPI_Reduce(&myres_real,&res_part_real,1,MPI_LONG_DOUBLE,MPI_SUM,0,comm5);
			MPI_Reduce(&myres_imag,&res_part_imag,1,MPI_LONG_DOUBLE,MPI_SUM,0,comm5);
		}

		if(rank == 0){
			// Rank 0 does its part of the task
			H_EST_PS_Sinc[i*5] = res_part_real + I*res_part_imag;

			// Rank 0 gathers the results of the rest of the master nodes
			for(int gr=1; gr<5; gr++){
				MPI_Recv(&res_part_real, 1, MPI_LONG_DOUBLE, gr*4, 0, MPI_COMM_WORLD, &status);
				MPI_Recv(&res_part_imag, 1, MPI_LONG_DOUBLE, gr*4, 1, MPI_COMM_WORLD, &status);
				H_EST_PS_Sinc[i*5 + gr] = res_part_real + I*res_part_imag;
			}

		} else if(rank%4==0 && rank!=0){
			MPI_Send(&res_part_real, 1, MPI_LONG_DOUBLE, 0, 0, MPI_COMM_WORLD);
			MPI_Send(&res_part_imag, 1, MPI_LONG_DOUBLE, 0, 1, MPI_COMM_WORLD);
		}

		MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){
		stop = clock();
		printf("PS Sinc Interpolation - Elapsed time %f\n",(double) (stop - start));
	}

	/* ----------------------------------------------------------------------------------------*/
	/* ----------------------------------------------------------------------------------------*/

	/* ----------------------------------------------------------------------------------------*/
	/* -------------------------------- PS MMSE INTERPOLATION ---------------------------------*
	* 
	*  Levels of dependency are analyzed among the functions that are integrated in the MMSE Estimator.
	*  Each level contain several functions that are independent. Thus, they can be executed in 
	*  parallel. Each function is executed in only 1 processor. The dependency mas as well as a 
	*  execution flow graph is shown below. Functions on the same row can be executed in parallel.
	*  After their execution, the results are sent to another processor, which will gather the 
	*  results from other processors to perform its own task. This system can be seen as a chain of 
	*  dependencies, where each level (row) cannot be executed without the completion of the previous
	*  
	*  levels |       Processor1	   |	Processor2	|	Processor3	|	Processor4       |
	* -----   |------------------------------------------------------------------------------|
	*  01st   | inv(F)   	(1)		   |    ow2*I (2)   |	 F' (3)	 	|	X4*F (4)         |
	*  02nd   | inv(F)*Hest (1) 	   |			    |	F'*X4 (5)   |  (X4*F)' (6)       |
	*  03rd   | inv(inv(F)*Hest) (1)   |                |               |                    |
	*         | (inv(F)*Hest)*         |                |               |                    |
	*  04th   |    inv(inv(F)*Hest)    |                |               |                    |
	*         | = Rhh (1)              |                |               |                    |
	*  05th   | Rhh*F'*X4 = Rhy	(7)	   |				|				|	Rhh*(X4*F)' (8)  |
	*  06th   | X4*F*Rhh*(X4*F)' (9)   |   F*Rhy (10)   |               |                    |
	*         | ow2*I +                |                |               |                    |
	*  07th   |   X4*F*Rhh*(X4*F)'     |                |               |                    |
	*         |   =Ryy (11)            |                |               |                    |
	*  08th   | inv(Ryy) (11)		   |	            |               |                    |
	*  09th   | F*Rhy*inv(Ryy) (12)    |                |               |                    |
	*  10th   | F*Rhy*inv(Ryy)*Rx_symbol (0)            |               |                    |
	*/
	/* ----------------------------------------------------------------------------------------*/

	long double complex **Fmatrix; malloc2dLongDoubleComplex(&Fmatrix, SAMPUTIL, SAMPUTIL);
	long double complex **FHermitian; malloc2dLongDoubleComplex(&FHermitian, SAMPUTIL, SAMPUTIL);
	long double complex **X4Hermitian; malloc2dLongDoubleComplex(&X4Hermitian, SAMPUTIL, SAMPUTIL);
	long double complex **X4; malloc2dLongDoubleComplex(&X4, SAMPUTIL, SAMPUTIL);
	long double complex **invRyy; malloc2dLongDoubleComplex(&invRyy, SAMPUTIL, SAMPUTIL);
	long double complex **invF; malloc2dLongDoubleComplex(&invF, SAMPUTIL, SAMPUTIL);
	long double complex **temp1; malloc2dLongDoubleComplex(&temp1, SAMPUTIL, SAMPUTIL);
	long double complex **temp2; malloc2dLongDoubleComplex(&temp2, SAMPUTIL, SAMPUTIL);
	long double complex **temp3; malloc2dLongDoubleComplex(&temp3, SAMPUTIL, SAMPUTIL);
	long double complex **Id; malloc2dLongDoubleComplex(&Id, SAMPUTIL, SAMPUTIL);

	long double complex **rx_symbols1; malloc2dLongDoubleComplex(&rx_symbols1, SAMPUTIL, 1);
	long double complex **H_EST1; malloc2dLongDoubleComplex(&H_EST1, SAMPUTIL, 1);

	long double **Fmatrix_Re; malloc2dLongDouble(&Fmatrix_Re, SAMPUTIL, SAMPUTIL);
	long double **Fmatrix_Im; malloc2dLongDouble(&Fmatrix_Im, SAMPUTIL, SAMPUTIL);
	long double **FHermitian_Re; malloc2dLongDouble(&FHermitian_Re, SAMPUTIL, SAMPUTIL);
	long double **FHermitian_Im; malloc2dLongDouble(&FHermitian_Im, SAMPUTIL, SAMPUTIL);
	long double **X4Hermitian_Re; malloc2dLongDouble(&X4Hermitian_Re, SAMPUTIL, SAMPUTIL);
	long double **X4Hermitian_Im; malloc2dLongDouble(&X4Hermitian_Im, SAMPUTIL, SAMPUTIL);
	long double **X4_Re; malloc2dLongDouble(&X4_Re, SAMPUTIL, SAMPUTIL);
	long double **X4_Im; malloc2dLongDouble(&X4_Im, SAMPUTIL, SAMPUTIL);
	long double **invF_Re; malloc2dLongDouble(&invF_Re, SAMPUTIL, SAMPUTIL);
	long double **invF_Im; malloc2dLongDouble(&invF_Im, SAMPUTIL, SAMPUTIL);
	long double **Id_Re; malloc2dLongDouble(&Id_Re, SAMPUTIL, SAMPUTIL);
	long double **Id_Im; malloc2dLongDouble(&Id_Im, SAMPUTIL, SAMPUTIL);
	long double **temp1_Re; malloc2dLongDouble(&temp1_Re, SAMPUTIL, SAMPUTIL);
	long double **temp1_Im; malloc2dLongDouble(&temp1_Im, SAMPUTIL, SAMPUTIL);
	long double **temp2_Re; malloc2dLongDouble(&temp2_Re, SAMPUTIL, SAMPUTIL);
	long double **temp2_Im; malloc2dLongDouble(&temp2_Im, SAMPUTIL, SAMPUTIL);
	long double **temp3_Re; malloc2dLongDouble(&temp3_Re, SAMPUTIL, SAMPUTIL);
	long double **temp3_Im; malloc2dLongDouble(&temp3_Im, SAMPUTIL, SAMPUTIL);

	if(rank==0){

		printf("Processing PS MMSE Interpolation...\n"); 
		start = clock();

		for(int i=0 ; i<SAMPUTIL; i++){
			rx_symbols1[i][0] = rx_symb_vec[i];
		}

		for(int r=0; r<SAMPUTIL; r++){
			for(int c=0; c<SAMPUTIL; c++){
				if((r==P0)&&(c==P0))
					X4[r][c] = tx_symb_vec[P0];
				else if ((r==P1)&&(c==P1))
					X4[r][c] = tx_symb_vec[P1];
				else if ((r==P2)&&(c==P2))
					X4[r][c] = tx_symb_vec[P2];
				else if ((r==P3)&&(c==P3))
					X4[r][c] = tx_symb_vec[P3];
				else
					X4[r][c] = 0.0;

				Fmatrix[c][r] = cexp(-2*I*PI*c*r/SAMPUTIL);
			}
		}

		complexToDouble(SAMPUTIL,X4,X4_Re,X4_Im);
		complexToDouble(SAMPUTIL,Fmatrix,Fmatrix_Re,Fmatrix_Im);
	}

	// Process 0 Broadcast Matrices X4 and Fmatrix
	MPI_Bcast(&(X4_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&(X4_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&(Fmatrix_Re[0][0]),SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&(Fmatrix_Im[0][0]),SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);

	// Each Process creates its own X4 and FMatrix
	doubleToComplex(SAMPUTIL,X4,X4_Re,X4_Im);
	doubleToComplex(SAMPUTIL,Fmatrix,Fmatrix_Re,Fmatrix_Im);

	if(rank==1){
		for(int r=0 ; r<SAMPUTIL ; r++)
			H_EST1[r][0] = H_EST_LT_LS[r];
		printf("Proc 1 working on inverse...\n");
		inverse(Fmatrix,SAMPUTIL,invF);
		printf("Proc 1 done with inverse...\n");
		multiply(invF,SAMPUTIL,SAMPUTIL,H_EST1,SAMPUTIL,SAMPUTIL,temp1);		//temp1 = invF*H_EST
		hermitian(temp1,SAMPUTIL,SAMPUTIL,temp2);								//temp2 = (invF*H_EST)'
		multiplyVxVeqM(temp1,SAMPUTIL,SAMPUTIL,temp2,SAMPUTIL,SAMPUTIL,temp3);	//temp3=Rhh
		complexToDouble(SAMPUTIL,temp3,temp3_Re,temp3_Im);
		MPI_Send(&(temp3_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 7, tag1, MPI_COMM_WORLD);
		MPI_Send(&(temp3_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 7, tag2, MPI_COMM_WORLD);
		MPI_Send(&(temp3_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 8, tag1, MPI_COMM_WORLD);
		MPI_Send(&(temp3_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 8, tag2, MPI_COMM_WORLD);
		printf("Proc %d sent its job to 7\n",rank);
		printf("Proc %d sent its job to 8\n",rank);
	} else if(rank==2){
		identity(Id,SAMPUTIL,OW2);
		complexToDouble(SAMPUTIL,Id,Id_Re,Id_Im);		
		MPI_Send(&(Id_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 11, tag1, MPI_COMM_WORLD);
		MPI_Send(&(Id_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 11, tag2, MPI_COMM_WORLD);
		printf("Proc %d sent its job to 11\n",rank);
	} else if(rank==3){
		hermitian(Fmatrix,SAMPUTIL,SAMPUTIL,FHermitian);
		complexToDouble(SAMPUTIL,FHermitian,FHermitian_Re,FHermitian_Im);  
		MPI_Send(&(FHermitian_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 5, tag1, MPI_COMM_WORLD);
		MPI_Send(&(FHermitian_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 5, tag2, MPI_COMM_WORLD);
		printf("Proc %d sent its job to 5\n",rank);
	} else if(rank==4){
		multiply(X4,SAMPUTIL,SAMPUTIL,Fmatrix,SAMPUTIL,SAMPUTIL,temp1);
		complexToDouble(SAMPUTIL,temp1,temp1_Re,temp1_Im);
		MPI_Send(&(temp1_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 6, tag1, MPI_COMM_WORLD);
		MPI_Send(&(temp1_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 6, tag2, MPI_COMM_WORLD);
		MPI_Send(&(temp1_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 9, tag1, MPI_COMM_WORLD);
		MPI_Send(&(temp1_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 9, tag2, MPI_COMM_WORLD);
		printf("Proc %d sent its job to 6\n",rank);
		printf("Proc %d sent its job to 9\n",rank);
	} else if(rank==5){
		printf("Proc %d waits...\n",rank);
		MPI_Recv(&(temp1_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 3, tag1, MPI_COMM_WORLD, &status);
		MPI_Recv(&(temp1_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 3, tag2, MPI_COMM_WORLD, &status);
		printf("Proc %d recv job from 3\n",rank);
		doubleToComplex(SAMPUTIL,temp1,temp1_Re,temp1_Im);
		multiply(temp1,SAMPUTIL,SAMPUTIL,X4,SAMPUTIL,SAMPUTIL,temp2);			//temp2=F'*X4
		complexToDouble(SAMPUTIL,temp2,temp2_Re,temp2_Im);
		MPI_Send(&(temp2_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 7, tag1, MPI_COMM_WORLD);
		MPI_Send(&(temp2_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 7, tag2, MPI_COMM_WORLD);
		printf("Proc %d sent its job to 7\n",rank);
	} else if(rank==6){
		printf("Proc %d waits...\n",rank);
		MPI_Recv(&(temp1_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 4, tag1, MPI_COMM_WORLD, &status);
		MPI_Recv(&(temp1_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 4, tag2, MPI_COMM_WORLD, &status);
		printf("Proc %d recv job from 4\n",rank);
		doubleToComplex(SAMPUTIL,temp1,temp1_Re,temp1_Re);						//temp1=X4*F
		hermitian(temp1,SAMPUTIL,SAMPUTIL,temp2);								//temp2=(X4*F)'
		complexToDouble(SAMPUTIL,temp2,temp2_Re,temp2_Im);
		MPI_Send(&(temp2_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 8, tag1, MPI_COMM_WORLD);
		MPI_Send(&(temp2_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 8, tag2, MPI_COMM_WORLD);
		MPI_Send(&(temp2_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 9, tag1, MPI_COMM_WORLD);
		MPI_Send(&(temp2_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 9, tag2, MPI_COMM_WORLD);
		printf("Proc %d sent its job to 8\n",rank);
		printf("Proc %d sent its job to 9\n",rank);
	} else if(rank==7){
		printf("Proc %d waits...\n",rank);
		MPI_Recv(&(temp1_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 5, tag1, MPI_COMM_WORLD, &status);
		MPI_Recv(&(temp1_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 5, tag2, MPI_COMM_WORLD, &status);
		printf("Proc %d recv job from 5\n",rank);
		MPI_Recv(&(temp2_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 1, tag1, MPI_COMM_WORLD, &status);
		MPI_Recv(&(temp2_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 1, tag2, MPI_COMM_WORLD, &status);
		printf("Proc %d recv job from 1\n",rank);
		doubleToComplex(SAMPUTIL,temp1,temp1_Re,temp1_Re);						//temp1=F'*X4
		doubleToComplex(SAMPUTIL,temp2,temp2_Re,temp2_Re);						//temp2=Rhh
		multiply(temp2,SAMPUTIL,SAMPUTIL,temp1,SAMPUTIL,SAMPUTIL,temp3);		//temp3=Rhy=Rhh*F'*X4
		complexToDouble(SAMPUTIL,temp3,temp3_Re,temp3_Im);
		MPI_Send(&(temp3_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 10, tag1, MPI_COMM_WORLD);
		MPI_Send(&(temp3_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 10, tag2, MPI_COMM_WORLD);
		printf("Proc %d sent its job to 10\n",rank);
	} else if(rank==8){
		printf("Proc %d waits...\n",rank);
		MPI_Recv(&(temp1_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 6, tag1, MPI_COMM_WORLD, &status);
		MPI_Recv(&(temp1_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 6, tag2, MPI_COMM_WORLD, &status);
		printf("Proc %d recv job from 6\n",rank);
		printf("Proc %d waits...\n",rank);
		MPI_Recv(&(temp2_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 1, tag1, MPI_COMM_WORLD, &status);
		MPI_Recv(&(temp2_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 1, tag2, MPI_COMM_WORLD, &status);
		printf("Proc %d recv job from 1\n",rank);
		doubleToComplex(SAMPUTIL,temp1,temp1_Re,temp1_Re);						//temp1=(X4*F)'
		doubleToComplex(SAMPUTIL,temp2,temp2_Re,temp2_Re);						//temp2=Rhh
		multiply(temp2,SAMPUTIL,SAMPUTIL,temp1,SAMPUTIL,SAMPUTIL,temp3);		//temp3=Rhh*(X4*F)'
		complexToDouble(SAMPUTIL,temp3,temp3_Re,temp3_Im);
		MPI_Send(&(temp3_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 9, tag1, MPI_COMM_WORLD);
		MPI_Send(&(temp3_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 9, tag2, MPI_COMM_WORLD);
		printf("Proc %d sent its job to 9\n",rank);
	} else if(rank==9){
		printf("Proc %d waits...\n",rank);
		MPI_Recv(&(temp1_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 4, tag1, MPI_COMM_WORLD, &status);
		MPI_Recv(&(temp1_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 4, tag2, MPI_COMM_WORLD, &status);
		printf("Proc %d recv job from 4\n",rank);
		printf("Proc %d waits...\n",rank);
		MPI_Recv(&(temp2_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 8, tag1, MPI_COMM_WORLD, &status);
		MPI_Recv(&(temp2_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 8, tag2, MPI_COMM_WORLD, &status);
		printf("Proc %d recv job from 8\n",rank);
		doubleToComplex(SAMPUTIL,temp1,temp1_Re,temp1_Re);						//temp1=X4*F
		doubleToComplex(SAMPUTIL,temp2,temp2_Re,temp2_Re);						//temp2=Rhh*(X4*F)'
		multiply(temp1,SAMPUTIL,SAMPUTIL,temp2,SAMPUTIL,SAMPUTIL,temp3);		//temp3=X4*F*Rhh*(X4*F)'
		complexToDouble(SAMPUTIL,temp3,temp3_Re,temp3_Im);
		MPI_Send(&(temp3_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 11, tag1, MPI_COMM_WORLD);
		MPI_Send(&(temp3_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 11, tag2, MPI_COMM_WORLD);
		printf("Proc %d sent its job to 11\n",rank);
	} else if(rank==10){
		printf("Proc %d waits...\n",rank);
		MPI_Recv(&(temp1_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 7, tag1, MPI_COMM_WORLD, &status);
		MPI_Recv(&(temp1_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 7, tag2, MPI_COMM_WORLD, &status);
		printf("Proc %d recv job from 7\n",rank);
		doubleToComplex(SAMPUTIL,temp1,temp1_Re,temp1_Re);						//temp1=Rhy
		multiply(Fmatrix,SAMPUTIL,SAMPUTIL,temp1,SAMPUTIL,SAMPUTIL,temp2);		//temp2=F*Rhy
		complexToDouble(SAMPUTIL,temp2,temp2_Re,temp2_Im);
		MPI_Send(&(temp2_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 12, tag1, MPI_COMM_WORLD);
		MPI_Send(&(temp2_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 12, tag2, MPI_COMM_WORLD);
		printf("Proc %d sent its job to 12\n",rank);
	} else if(rank==11){
		printf("Proc %d waits...\n",rank);
		MPI_Recv(&(temp1_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 2, tag1, MPI_COMM_WORLD, &status);
		MPI_Recv(&(temp1_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 2, tag2, MPI_COMM_WORLD, &status);
		printf("Proc %d recv job from 2\n",rank);
		printf("Proc %d waits...\n",rank);
		MPI_Recv(&(temp2_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 9, tag1, MPI_COMM_WORLD, &status);
		MPI_Recv(&(temp2_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 9, tag2, MPI_COMM_WORLD, &status);
		printf("Proc %d recv job from 9\n",rank);
		doubleToComplex(SAMPUTIL,temp1,temp1_Re,temp1_Re);						//temp1=ow2*I
		doubleToComplex(SAMPUTIL,temp2,temp2_Re,temp2_Re);						//temp2=X4*F*Rhh*(X4*F)'
		addition(temp1,SAMPUTIL,SAMPUTIL,temp2,SAMPUTIL,SAMPUTIL,temp3);		//temp3=ow2*I+X4*F*Rhh*(X4*F)'=Ryy
		printf("Proc 11 working on inverse...\n");
		inverse(temp3,SAMPUTIL,temp1);											//temp1=inv(Ryy)
		printf("Proc 11 done with inverse...\n");
		complexToDouble(SAMPUTIL,temp1,temp1_Re,temp1_Im);
		MPI_Send(&(temp1_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 12, tag1, MPI_COMM_WORLD);
		MPI_Send(&(temp1_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 12, tag2, MPI_COMM_WORLD);
		printf("Proc %d sent its job to 12\n",rank);
	} else if(rank==12){
		printf("Proc %d waits...\n",rank);
		MPI_Recv(&(temp1_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 10, tag1, MPI_COMM_WORLD, &status);
		MPI_Recv(&(temp1_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 10, tag2, MPI_COMM_WORLD, &status);
		printf("Proc %d recv job from 10\n",rank);
		printf("Proc %d waits...\n",rank);
		MPI_Recv(&(temp2_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 11, tag1, MPI_COMM_WORLD, &status);
		MPI_Recv(&(temp2_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 11, tag2, MPI_COMM_WORLD, &status);
		printf("Proc %d recv job from 11\n",rank);
		doubleToComplex(SAMPUTIL,temp1,temp1_Re,temp1_Re);						//temp1=F*Rhy
		doubleToComplex(SAMPUTIL,temp2,temp2_Re,temp2_Re);						//temp2=inv(Ryy)
		multiply(temp1,SAMPUTIL,SAMPUTIL,temp2,SAMPUTIL,SAMPUTIL,temp3);		//temp3=F*Rhy*inv(Ryy)
		complexToDouble(SAMPUTIL,temp3,temp3_Re,temp3_Im);
		MPI_Send(&(temp3_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 0, tag1, MPI_COMM_WORLD);
		MPI_Send(&(temp3_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 0, tag2, MPI_COMM_WORLD);
		printf("Proc %d sent its job to 0\n",rank);
	} else if(rank == 0){
		printf("Proc %d waits...\n",rank);
		MPI_Recv(&(temp1_Re[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 12, tag1, MPI_COMM_WORLD, &status);
		MPI_Recv(&(temp1_Im[0][0]), SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 12, tag2, MPI_COMM_WORLD, &status);
		printf("Proc %d recv job from 12\n",rank);
		doubleToComplex(SAMPUTIL,temp1,temp1_Re,temp1_Re);						//temp1=F*Rhy*inv(Ryy)
		multiply(temp1,SAMPUTIL,SAMPUTIL,rx_symbols1,SAMPUTIL,SAMPUTIL,temp2);	//temp3=F*Rhy*inv(Ryy)*rx_symbols1

		for(int r=0 ; r<SAMPUTIL ; r++)
			H_EST_PS_MMSE[r] = temp2[r][0];
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){
		stop = clock();
		printf("PS MMSE - Elapsed time %f\n",(double) (stop - start));
	}

	/* ----------------------------------------------------------------------------------------*/
	/* ----------------------------------------------------------------------------------------*/

	MPI_Group_free(&gr1);MPI_Group_free(&gr2);MPI_Group_free(&gr3);MPI_Group_free(&gr4);MPI_Group_free(&gr5);
	MPI_Comm_free(&comm1);MPI_Comm_free(&comm2);MPI_Comm_free(&comm3);MPI_Comm_free(&comm4);MPI_Comm_free(&comm5);

	MPI_Finalize();

	return 0;
}