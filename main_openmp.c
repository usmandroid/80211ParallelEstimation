#include <omp.h>
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

int main(void) {

	printf( "OMP_STACKSIZE=%s \n", getenv( "OMP_STACKSIZE" ) );

	clock_t start, stop, start_tot, stop_tot;
	long double complex tx_symb_vec[SAMPUTIL],rx_symb_vec[SAMPUTIL];
	long double complex H_EST_LT_LS[SAMPUTIL],H_EST_PS_Linear[SAMPUTIL],H_EST_PS_Cubic[SAMPUTIL];
	long double complex H_EST_PS_Sinc[SAMPUTIL], H_EST_PS_MMSE[SAMPUTIL];
	int OFDM_block = 0;

	long double complex **Fmatrix = new long double complex*[SAMPUTIL];
	for (int i = 0; i < SAMPUTIL; i++) {
		Fmatrix[i] = new long double complex[SAMPUTIL];
	}
	for (int f=0; f<SAMPUTIL; f++){
        for (int t=0; t<SAMPUTIL; t++){
            Fmatrix[t][f] = cexp(-2*I*PI*t*f/SAMPUTIL);
		}
    }

	/* One OFDM Symbol isextracted to perform channel estimation */
	printf("**** Processing Block %d\n",OFDM_block);
	for(int r=0 ; r<SAMPUTIL ; r++){
		tx_symb_vec[r] = tx_symb[SAMPUTIL*OFDM_block + r];
		rx_symb_vec[r] = rx_symb[SAMPUTIL*OFDM_block + r];
	}

	start_tot = clock();
	printf("CLOCKS_PER_SEC = %d\n",CLOCKS_PER_SEC);
	printf("Processing LT Least Squares...\n"); start = clock(); 
	WiFi_channel_estimation_LT_LS(tx_symb_vec,rx_symb_vec,H_EST_LT_LS); stop = clock();
	printf("Elapsed time: %f\n",(double) (stop - start));		
	printf("Processing PS Linear Interpolation...\n"); start = clock(); 
	WiFi_channel_estimation_PS_Linear(tx_symb_vec,rx_symb_vec,H_EST_PS_Linear); stop = clock();
	printf("Elapsed time: %f\n",(double) (stop - start));		
	printf("Processing PS Cubic Interpolation...\n"); start = clock();
	WiFi_channel_estimation_PS_Cubic(tx_symb_vec,rx_symb_vec,H_EST_PS_Cubic); stop = clock();
	printf("Elapsed time: %f\n",(double) (stop - start));		
	printf("Processing PS Sinc Interpolation...\n"); start = clock();
	WiFi_channel_estimation_PS_Sinc(tx_symb_vec,rx_symb_vec,H_EST_PS_Sinc); stop = clock();
	printf("Elapsed time: %f\n",(double) (stop - start));		
	printf("Processing PS MMSE...\n"); start = clock();
	WiFi_channel_estimation_PS_MMSE(tx_symb_vec,rx_symb_vec,Fmatrix,OW2,H_EST_LT_LS,H_EST_PS_MMSE); stop = clock();
	printf("Elapsed time: %f\n",(double) (stop - start));		

	stop_tot = clock();
	printf("**** Total Elapsed time: %f\n",(double) (stop - start_tot));

	// printVect(H_EST_LT_LS,SAMPUTIL,"H_EST_LT_LS");
	// printVect(H_EST_PS_Linear,SAMPUTIL,"H_EST_PS_Linear");
	// printVect(H_EST_PS_Cubic,SAMPUTIL,"H_EST_PS_Cubic");
	// printVect(H_EST_PS_Sinc,SAMPUTIL,"H_EST_PS_Sinc");
	// printVect(H_EST_PS_MMSE,SAMPUTIL,"H_EST_PS_MMSE");

	free(Fmatrix);
	
	return 0;
}

void WiFi_channel_estimation_LT_LS(long double complex tx_pre[], long double complex rx_pre[],long double complex H_EST[]){
	long double complex conj;
	#pragma omp parallel num_threads(SAMPUTIL) private(conj)
	#pragma omp for schedule(static,1)
	for(int i=0 ; i<SAMPUTIL ; i++){
		conj = creal(tx_pre[i]) - cimag(tx_pre[i]);
		H_EST[i] = ( conj*rx_pre[i] ) / ( conj*tx_pre[i] );
	} 
	H_EST[26] = 0.0;
	
}

void WiFi_channel_estimation_PS_Linear(long double complex tx_symbols[], long double complex rx_symbols[],long double complex H_EST[]){
	long double complex H_PILOTS[4];
	long double complex tx_pilots[4] = {tx_symbols[P0],tx_symbols[P1],tx_symbols[P2],tx_symbols[P3]};
	long double complex rx_pilots[4] = {rx_symbols[P0],rx_symbols[P1],rx_symbols[P2],rx_symbols[P3]};
	int tid;
	long double alpha, delta = P1 - P0;
	for(int i=0 ; i<4 ; i++){
		H_PILOTS[i] = rx_pilots[i] / tx_pilots[i];
	}
	
	omp_set_num_threads(SAMPUTIL);
	#pragma omp parallel private(tid,alpha)
	{
		tid = omp_get_thread_num();
		if(tid<P1){
			alpha = (tid-P0)/delta;
			H_EST[tid] = H_PILOTS[0]+( (H_PILOTS[1]-H_PILOTS[0] )*alpha );
		} else if ((tid>=P1) && (tid<P2)){
            alpha = (tid-P1)/delta;
            H_EST[tid] = H_PILOTS[1]+( (H_PILOTS[2]-H_PILOTS[1] )*alpha );
        } else if ((tid>=P2) && (tid<P3)) {
            alpha = (tid-P2)/delta;
            H_EST[tid] = H_PILOTS[2]+( (H_PILOTS[3]-H_PILOTS[2] )*alpha );
        } else if (tid>=P3) {
            alpha = (tid-P2)/delta;
            H_EST[tid] = H_PILOTS[2]+( (H_PILOTS[3]-H_PILOTS[2] )*alpha );
        }
	}
	
}

void WiFi_channel_estimation_PS_Cubic(long double complex tx_symbols[], long double complex rx_symbols[],long double complex H_EST[]){
	long double complex H_PILOTS[4];
	long double complex tx_pilots[4] = {tx_symbols[P0],tx_symbols[P1],tx_symbols[P2],tx_symbols[P3]};
	long double complex rx_pilots[4] = {rx_symbols[P0],rx_symbols[P1],rx_symbols[P2],rx_symbols[P3]};
	long double complex f0, f01, f12, f23, f012, f123, f0123;
	long double alpha, delta = P1 - P0;
	int tid;
	
	for(int i=0 ; i<4 ; i++){
		H_PILOTS[i] = rx_pilots[i] / tx_pilots[i];
	}
	f0 	  = H_PILOTS[0];
	f01   = (H_PILOTS[1]-H_PILOTS[0]) / delta;
	f12   = (H_PILOTS[2]-H_PILOTS[1]) / delta;
    f23   = (H_PILOTS[3]-H_PILOTS[2]) / delta;
    f012  = (f12-f01) / delta;
    f123  = (f23-f12) / delta;
    f0123 = (f123-f012) / delta;
	
	#pragma omp parallel for num_threads(SAMPUTIL)
	for(int k=0 ; k<SAMPUTIL ; k++){
        H_EST[k] = f0 + f01*(k-P0) + f012*(k-P0)*(k-P1) + f0123*(k-P0)*(k-P1)*(k-P2);
    }
	
}

void WiFi_channel_estimation_PS_Sinc(long double complex tx_symbols[], long double complex rx_symbols[],long double complex H_EST[]){
	double complex H_PILOTS[4],sinc1[SAMPUTIL],sinc2[SAMPUTIL],sinc3[SAMPUTIL],sinc4[SAMPUTIL];
	long double complex tx_pilots[4] = {tx_symbols[P0],tx_symbols[P1],tx_symbols[P2],tx_symbols[P3]};
	long double complex rx_pilots[4] = {rx_symbols[P0],rx_symbols[P1],rx_symbols[P2],rx_symbols[P3]};
	long double delta = P1 - P0;
	double a;

	for(int i=0 ; i<4 ; i++){
		H_PILOTS[i] = rx_pilots[i] / tx_pilots[i];
	}

	#pragma omp parallel for num_threads(SAMPUTIL) private(a)
	for(int k=0 ; k<SAMPUTIL ; k++){
		a = (k-P0) / delta;
		sinc1[k] = H_PILOTS[0]*sinc(a);
	}
    #pragma omp parallel for num_threads(SAMPUTIL) private(a)
	for(int k=0 ; k<SAMPUTIL ; k++){
		a = (k-P1) / delta;
    	sinc2[k] = H_PILOTS[1]*sinc(a);
	}
    #pragma omp parallel for num_threads(SAMPUTIL) private(a)
	for(int k=0 ; k<SAMPUTIL ; k++) {
		a = (k-P2) / delta;
		sinc3[k] = H_PILOTS[2]*sinc(a);
	}
	#pragma omp parallel for num_threads(SAMPUTIL) private(a)
	for(int k=0 ; k<SAMPUTIL ; k++){
		a = (k-P3) / delta;
		sinc4[k] = H_PILOTS[3]*sinc(a);
	}
	#pragma omp parallel for num_threads(SAMPUTIL)
	for(int k=0 ; k<SAMPUTIL ; k++){
		H_EST[k] = sinc1[k] + sinc2[k] + sinc3[k] + sinc4[k];
	}

}

void WiFi_channel_estimation_PS_MMSE(long double complex tx_symbols[], long double complex rx_symbols[], long double complex **F, double ow2, long double complex H_EST_LS[], long double complex H_EST_MMSE[]){
	long double complex **FHermitian = new long double complex*[SAMPUTIL];
	long double complex **X4Hermitian = new long double complex*[SAMPUTIL];
	long double complex **X4 = new long double complex*[SAMPUTIL];
	long double complex **Rhh = new long double complex*[SAMPUTIL];
	long double complex **Rhy = new long double complex*[SAMPUTIL];
	long double complex **Ryy = new long double complex*[SAMPUTIL];
	long double complex **invRyy = new long double complex*[SAMPUTIL];
	long double complex **invF = new long double complex*[SAMPUTIL];
	long double complex **temp1 = new long double complex*[SAMPUTIL];
	long double complex **temp2 = new long double complex*[SAMPUTIL];
	long double complex **temp3 = new long double complex*[SAMPUTIL];
	long double complex **Id = new long double complex*[SAMPUTIL];
	long double complex **rx_symbols1 = new long double complex*[SAMPUTIL];
	long double complex **H_EST1 = new long double complex*[SAMPUTIL];

	for (int i = 0; i < SAMPUTIL; i++) {
		/* MATRICES */
		FHermitian[i] = new long double complex[SAMPUTIL];
		X4Hermitian[i] = new long double complex[SAMPUTIL];
		X4[i] = new long double complex[SAMPUTIL];
		Rhh[i] = new long double complex[SAMPUTIL];
		Rhy[i] = new long double complex[SAMPUTIL];
		Ryy[i] = new long double complex[SAMPUTIL];
		invRyy[i] = new long double complex[SAMPUTIL];
		invF[i] = new long double complex[SAMPUTIL];
		temp1[i] = new long double complex[SAMPUTIL];
		temp2[i] = new long double complex[SAMPUTIL];
		temp3[i] = new long double complex[SAMPUTIL];
		Id[i] = new long double complex[SAMPUTIL];
		/* VECTORS */ 
		rx_symbols1[i] = new long double complex[1];
		H_EST1[i] = new long double complex[1];
	}

	for(int i=0 ; i<SAMPUTIL; i++){
		rx_symbols1[i][0] = rx_symbols[i];
	}

	#pragma omp parallel num_threads(SAMPUTIL)
	{
		#pragma omp for schedule(static)
		for(int r=0; r<SAMPUTIL; r++){
			for(int c=0; c<SAMPUTIL; c++){
				if((r==P0)&&(c==P0))
					X4[r][c] = tx_symbols[P0];
				else if ((r==P1)&&(c==P1))
					X4[r][c] = tx_symbols[P1];
				else if ((r==P2)&&(c==P2))
					X4[r][c] = tx_symbols[P2];
				else if ((r==P3)&&(c==P3))
					X4[r][c] = tx_symbols[P3];
				else
					X4[r][c] = 0.0;
			}
		}
	}

/*
	int size = 54;
	long double complex **Matrix = new long double complex*[size];
	long double complex **invMatrix = new long double complex*[size];
	long double complex **resMatrix = new long double complex*[size];
	for(int i=0 ; i<size ; i++){
		Matrix[i] = new long double complex[size];
		invMatrix[i] = new long double complex[size];
		resMatrix[i] = new long double complex[size];
		for(int j=0 ; j<size ; j++){
			Matrix[i][j] = ((int) rand()%20) + I*((int) rand()%20);
		}
	}
	double complex det = (double complex) determinant_impl(Matrix,size);
	//printf("%lf + i%lf\n", creal(det), cimag(det));	

	inverse(Matrix,size,invMatrix);											//invMatrix
	multiply(Matrix,size,size,Matrix,size,size,resMatrix);
*/

	for(int r=0 ; r<SAMPUTIL ; r++)
		H_EST1[r][0] = H_EST_LS[r];

	hermitian_omp(F,SAMPUTIL,SAMPUTIL,FHermitian);								// FHermitian = F'
	hermitian_omp(X4,SAMPUTIL,SAMPUTIL,X4Hermitian);							//X4Hermitian = X4'

	inverse_omp(F,SAMPUTIL,invF);												//invF
	multiply_omp(invF,SAMPUTIL,SAMPUTIL,H_EST1,SAMPUTIL,1,temp1);				//temp1 = invF*H_EST
	hermitian_omp(temp1,SAMPUTIL,SAMPUTIL,temp2);								//temp2 = (invF*H_EST)'
	multiplyVxVeqM_omp(temp1,SAMPUTIL,SAMPUTIL,temp2,SAMPUTIL,SAMPUTIL,Rhh);	//Rhh = (invF*H_EST)*(invF*H_EST)'

	multiply_omp(Rhh,SAMPUTIL,SAMPUTIL,FHermitian,SAMPUTIL,SAMPUTIL,temp1);		//temp1 = Rhh*F'
	multiply_omp(temp1,SAMPUTIL,SAMPUTIL,X4,SAMPUTIL,SAMPUTIL,Rhy);				//Rhy

	multiply_omp(temp1,SAMPUTIL,SAMPUTIL,X4Hermitian,SAMPUTIL,SAMPUTIL,temp2);	//temp2 = Rhh*F'*X4'
	multiply_omp(F,SAMPUTIL,SAMPUTIL,temp1,SAMPUTIL,SAMPUTIL,temp1);  			//temp1 = F*Rhh*F'*X4'
	multiply_omp(X4,SAMPUTIL,SAMPUTIL,temp1,SAMPUTIL,SAMPUTIL,temp2);			//temp2 = X4*F*Rhh*F'*X4'

	identity_omp(Id,SAMPUTIL,ow2);
	addition_omp(Id,SAMPUTIL,SAMPUTIL,temp2,SAMPUTIL,SAMPUTIL,Ryy);				//Ryy

	inverse_omp(Ryy,SAMPUTIL,invRyy);											//invRyy

 	multiply_omp(F,SAMPUTIL,SAMPUTIL,Rhy,SAMPUTIL,SAMPUTIL,temp1);				//temp1 = F*Rhy
	multiply_omp(invRyy,SAMPUTIL,SAMPUTIL,rx_symbols1,SAMPUTIL,1,temp3);		//temp3 = invRyy*rx_symbols
	multiply_omp(temp1,SAMPUTIL,SAMPUTIL,temp3,SAMPUTIL,1,temp2);				//H_EST = F*Rhy*invRyy*rx_symbols

	for(int r=0 ; r<SAMPUTIL ; r++)
		H_EST_MMSE[r] = temp2[r][0];

	free(FHermitian);free(X4Hermitian);free(X4);free(Rhh);free(Rhy);free(Ryy);free(invRyy);
	free(invF);free(temp1);free(temp2);free(Id);free(temp3);free(rx_symbols1);free(H_EST1);
}