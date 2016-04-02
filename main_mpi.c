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

	long double complex **Fmatrix = new long double complex*[SAMPUTIL];
	long double **Fmatrix_real = new long double *[SAMPUTIL];
	long double **Fmatrix_imag = new long double *[SAMPUTIL];
	for (int i = 0; i < SAMPUTIL; i++) {
		Fmatrix[i] = new long double complex[SAMPUTIL];
		Fmatrix_real[i] = new long double[SAMPUTIL];
		Fmatrix_imag[i] = new long double[SAMPUTIL];
	}

	if(rank == 0){
		for (int f=0; f<SAMPUTIL; f++){
	        for (int t=0; t<SAMPUTIL; t++){
	            Fmatrix_real[t][f] = creal(cexp(-2*I*PI*t*f/SAMPUTIL));
	            Fmatrix_imag[t][f] = cimag(cexp(-2*I*PI*t*f/SAMPUTIL));
		 	}
    	}

		/* One OFDM Symbol isextracted to perform channel estimation */
		printf("Proc %d Processing Block %d\n", rank, OFDM_block);
		for(int r=0 ; r<SAMPUTIL ; r++){
			tx_symb_vec[r] = tx_symb[SAMPUTIL*OFDM_block + r];
			rx_symb_vec[r] = rx_symb[SAMPUTIL*OFDM_block + r];
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&(Fmatrix_real[0][0]),SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&(Fmatrix_imag[0][0]),SAMPUTIL*SAMPUTIL, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	// Each Process creates the complex Fourier Matrix
	for(int i=0 ; i<SAMPUTIL ; i++){
		for(int j=0 ; j<SAMPUTIL ; j++){
			Fmatrix[i][j] = Fmatrix_real[i][j] + I*Fmatrix_imag[i][j];
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);	

	/* ----------------------------------------------------------------------------------------*/
	/* -------------------------------- LEAST SQUARE ESTIMATOR --------------------------------*/
	/* ----------------------------------------------------------------------------------------*/
	
	long double complex conj1, conj2;
	long double res_LS[4];
	int chunk[4], tag1 = 1, tag2 = 2, tag3 = 3;
	long double tx_real[4], tx_imag[4], rx_real[4], rx_imag[4];
	if(rank == 0){

		printf("Processing LT Least Square...\n");
		start = clock();

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
	/* -------------------------------- PS COMMON VARIABLES -----------------------------------*/
	/* ----------------------------------------------------------------------------------------*/

	long double complex H_PILOTS[4];
	long double H_PILOTS_real[4],H_PILOTS_imag[4], res_Linear[6];
	long double alpha, delta = P1 - P0;
	tag1 = 1, tag2 = 2;
	int tag = 3;
	long double complex f0, f01, f12, f23, f012, f123, f0123;

	if(rank == 0){
		long double complex tx_pilots[4] = {tx_symb_vec[P0],tx_symb_vec[P1],tx_symb_vec[P2],tx_symb_vec[P3]};
		long double complex rx_pilots[4] = {rx_symb_vec[P0],rx_symb_vec[P1],rx_symb_vec[P2],rx_symb_vec[P3]};

		for(int i=0 ; i<4 ; i++){
			H_PILOTS_real[i] = creal(rx_pilots[i] / tx_pilots[i]);
			H_PILOTS_imag[i] = cimag(rx_pilots[i] / tx_pilots[i]);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	/* ----------------------------------------------------------------------------------------*/
	/* ----------------------------------------------------------------------------------------*/

	/* ----------------------------------------------------------------------------------------*/
	/* -------------------------------- PS LINEAR INTERPOLATION -------------------------------*/
	/* ----------------------------------------------------------------------------------------*/

	if(rank == 0){

		printf("Processing PS Linear Interpolation...\n"); 
		start = clock(); 

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
			MPI_Send(H_PILOTS_real, 4, MPI_LONG_DOUBLE, dest, tag2, MPI_COMM_WORLD);
			MPI_Send(H_PILOTS_imag, 4, MPI_LONG_DOUBLE, dest, tag3, MPI_COMM_WORLD);
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
		MPI_Recv(H_PILOTS_real, 4, MPI_LONG_DOUBLE, 0, tag2, MPI_COMM_WORLD, &status);
		MPI_Recv(H_PILOTS_imag, 4, MPI_LONG_DOUBLE, 0, tag3, MPI_COMM_WORLD, &status);

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
		printf("PS Linear - Elapsed time %f\n",(double) (stop - start));
	}

	/* ----------------------------------------------------------------------------------------*/
	/* ----------------------------------------------------------------------------------------*/

	/* ----------------------------------------------------------------------------------------*/
	/* -------------------------------- PS CUBIC INTERPOLATION --------------------------------*/
	/* ----------------------------------------------------------------------------------------*/

	long double complex myres, res_part_real, res_part_imag;
	long double myres_real, myres_imag;
	int k;

	/* Each process uses some common variables, which need to be
	*  initialized by the master process.
	*/
	if(rank == 0){

		printf("Processing PS Cubic Interpolation...\n"); start = clock(); 

		f0 = H_PILOTS[0];
		f01   = (H_PILOTS[1]-H_PILOTS[0]) / delta;
		f12   = (H_PILOTS[2]-H_PILOTS[1]) / delta;
	    f23   = (H_PILOTS[3]-H_PILOTS[2]) / delta;
	    f012  = (f12-f01) / delta;
	    f123  = (f23-f12) / delta;
	    f0123 = (f123-f012) / delta;
	}

	/* The Common variables are Broadcasted to all the nodes by
	*  the master node 
	*/
	MPI_Bcast(&f0,1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&f01,1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&f12,1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&f23,1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&f012,1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&f123,1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&f0123,1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	/* Communicators creation. 5 groups of 4 processes are created
	*  in order to perform a reduction operation. A communicator is
	*  assigned to each group, allowing them to share the results after
	*  the reduction.
	*/
	MPI_Comm MPI_GROUP1, MPI_GROUP2, MPI_GROUP3, MPI_GROUP4, MPI_GROUP5;
	MPI_Comm gr1, gr2, gr3, gr4, gr5;
	MPI_Comm comm1, comm2, comm3, comm4, comm5;
	static int group1[] = {0,1,2,3};
	static int group2[] = {4,5,6,7};
	static int group3[] = {8,9,10,11};
	static int group4[] = {12,13,14,15};
	static int group5[] = {16,17,18,19};
	MPI_Comm_group(MPI_COMM_WORLD,&MPI_GROUP1);
	MPI_Comm_group(MPI_COMM_WORLD,&MPI_GROUP2);
	MPI_Comm_group(MPI_COMM_WORLD,&MPI_GROUP3);
	MPI_Comm_group(MPI_COMM_WORLD,&MPI_GROUP4);
	MPI_Comm_group(MPI_COMM_WORLD,&MPI_GROUP5);
	MPI_Group_incl(MPI_GROUP1,4,group1,&gr1);
	MPI_Group_incl(MPI_GROUP2,4,group2,&gr2);
	MPI_Group_incl(MPI_GROUP3,4,group3,&gr3);
	MPI_Group_incl(MPI_GROUP4,4,group4,&gr4);
	MPI_Group_incl(MPI_GROUP5,4,group5,&gr5);

	MPI_Comm_create(MPI_COMM_WORLD, gr1, &comm1);
	MPI_Comm_create(MPI_COMM_WORLD, gr2, &comm2);
	MPI_Comm_create(MPI_COMM_WORLD, gr3, &comm3);
	MPI_Comm_create(MPI_COMM_WORLD, gr4, &comm4);
	MPI_Comm_create(MPI_COMM_WORLD, gr5, &comm5);

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0) printf("caca\n");

	for(int i=0; i<10; i++){
		/*   0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | .. | 18 | 19 |
		* k= 0   0   0   0   1   1   1   1   2        4    4
		* k= 5   5   5   5   6   6   6   6   7        9    9
		* ...
		* k= 50  50  50  50  51  51  51  51  52 .. 
		* 
		* Processes are divided into groups of four, where each component
		* of the group will perform one of the four operations needed 
		* in each channel subcarrier. At the end, the reduction sums up
		* the value obtained and its assigned to the vector in the correct
		* location.
		*/

		k = floor((rank + 1)/5) + 5*i;
		
		if(rank%4 == 0){
			myres = f0;
		}else if(rank%4 == 1){
			myres = f01*(k-P0);
		}else if(rank%4 == 2){
			myres = f012*(k-P0)*(k-P1);
		}else if(rank%4 == 3){
			myres = f0123*(k-P0)*(k-P1)*(k-P2);
		}

		myres_real = creal(myres); 
		myres_imag = cimag(myres);

		/* Each process send the real and the imaginary part in the following way:
		*  tag = 0 for the Real and tag = 1 for the imag
		*/
		if(rank<4 && rank!=0){
			MPI_Reduce(&myres_real,&res_part_real,1,MPI_LONG_DOUBLE,MPI_SUM,0*4,comm1);
			MPI_Reduce(&myres_imag,&res_part_imag,1,MPI_LONG_DOUBLE,MPI_SUM,0*4,comm1);
		} else if(rank < 8){
			MPI_Reduce(&myres_real,&res_part_real,1,MPI_LONG_DOUBLE,MPI_SUM,1*4,comm2);
			MPI_Reduce(&myres_imag,&res_part_imag,1,MPI_LONG_DOUBLE,MPI_SUM,1*4,comm2);
		} else if(rank < 8){
			MPI_Reduce(&myres_real,&res_part_real,1,MPI_LONG_DOUBLE,MPI_SUM,2*4,comm3);
			MPI_Reduce(&myres_imag,&res_part_imag,1,MPI_LONG_DOUBLE,MPI_SUM,2*4,comm3);
		} else if(rank < 12){
			MPI_Reduce(&myres_real,&res_part_real,1,MPI_LONG_DOUBLE,MPI_SUM,3*4,comm4);
			MPI_Reduce(&myres_imag,&res_part_imag,1,MPI_LONG_DOUBLE,MPI_SUM,3*4,comm4);
		} else if(rank < 16){
			MPI_Reduce(&myres_real,&res_part_real,1,MPI_LONG_DOUBLE,MPI_SUM,4*4,comm5);
			MPI_Reduce(&myres_imag,&res_part_imag,1,MPI_LONG_DOUBLE,MPI_SUM,4*4,comm5);
		} else if(rank == 0){
			/* Rank 0 does its part of the task */
			H_EST_PS_Cubic[i*5] = res_part_real + I*res_part_imag;
			/* Rank 0 gathers the results of the rest of the master nodes */
			for(int gr=1; gr<5; gr++){
				MPI_Recv(&res_part_real, 1, MPI_LONG_DOUBLE, 1*4, 0, MPI_COMM_WORLD, &status);
				MPI_Recv(&res_part_imag, 1, MPI_LONG_DOUBLE, 1*4, 1, MPI_COMM_WORLD, &status);
				H_EST_PS_Cubic[i*5 + gr] = res_part_real + I*res_part_imag;
			}
		}

		if(rank%4==0 && rank!=0){
			MPI_Send(&res_part_real, 1, MPI_LONG_DOUBLE, 0, 0, MPI_COMM_WORLD);
			MPI_Send(&res_part_real, 1, MPI_LONG_DOUBLE, 0, 1, MPI_COMM_WORLD);
		}

		/* Wait until Process 0 has received all the results */
		MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Group_free(&MPI_GROUP1); MPI_Group_free(&MPI_GROUP2);
	MPI_Group_free(&MPI_GROUP3); MPI_Group_free(&MPI_GROUP4); MPI_Group_free(&MPI_GROUP5);
	MPI_Group_free(&gr1);MPI_Group_free(&gr2);MPI_Group_free(&gr3);MPI_Group_free(&gr4);MPI_Group_free(&gr5);
	MPI_Comm_free(&comm1);MPI_Comm_free(&comm2);MPI_Comm_free(&comm3);MPI_Comm_free(&comm4);MPI_Comm_free(&comm5);


	/* ----------------------------------------------------------------------------------------*/
	/* ----------------------------------------------------------------------------------------*/
	
	/* ----------------------------------------------------------------------------------------*/
	/* -------------------------------- PS SINC INTERPOLATION ---------------------------------*/
	/* ----------------------------------------------------------------------------------------*/

	// HERE MY CODE

	/* ----------------------------------------------------------------------------------------*/
	/* ----------------------------------------------------------------------------------------*/

	/* ----------------------------------------------------------------------------------------*/
	/* -------------------------------- PS MMSE INTERPOLATION ---------------------------------*/
	/* ----------------------------------------------------------------------------------------*/

	// HERE MY CODE

	/* ----------------------------------------------------------------------------------------*/
	/* ----------------------------------------------------------------------------------------*/

	MPI_Finalize();

	return 0;
}