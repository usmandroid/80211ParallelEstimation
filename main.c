#include "inputs.h"
#include "utils.h"

void WiFi_channel_estimation_LT_LS(long double complex tx_pre[], long double complex rx_pre[],long double complex H_EST[]);
void WiFi_channel_estimation_PS_Linear(long double complex tx_symbols[], long double complex rx_symbols[],long double complex H_EST[]);
void WiFi_channel_estimation_PS_Cubic(long double complex tx_symbols[], long double complex rx_symbols[],long double complex H_EST[]);
void WiFi_channel_estimation_PS_Sinc(long double complex tx_symbols[], long double complex rx_symbols[],long double complex H_EST[]);
void WiFi_channel_estimation_PS_MMSE(long double complex tx_symbols[], long double complex rx_symbols[], long double complex **F, double ow2, long double complex H_EST_LS[], long double complex H_EST[]);

int main(void) {

	clock_t start, stop, start_tot;
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

	// printf("Processing LT Least Square...\n"); start = clock(); 
	// WiFi_channel_estimation_LT_LS(tx_preamble_fft,rx_preamble_fft,H_EST_LT_LS); stop = clock();
	// printf("Elapsed time: %f\n",(double) (stop - start));
	printf("Processing PS Linear Interpolation...\n"); start = clock(); 
	WiFi_channel_estimation_PS_Linear(tx_symb_vec,rx_symb_vec,H_EST_PS_Linear); stop = clock();
	for(int i=0; i<SAMPUTIL; i++){
    	printf("H_EST[%d] = %f + %fi\n",i,creal(H_EST_PS_Linear[i]),cimag(H_EST_PS_Linear[i]));
    }
	printf("Elapsed time: %f\n",(double) (stop - start));
	// printf("Processing PS Cubic Interpolation...\n"); start = clock();
	// WiFi_channel_estimation_PS_Cubic(tx_symb_vec,rx_symb_vec,H_EST_PS_Cubic); stop = clock();
	// printf("Elapsed time: %f\n",(double) (stop - start));
	// printf("Processing PS Sinc Interpolation...\n"); start = clock();
	// WiFi_channel_estimation_PS_Sinc(tx_symb_vec,rx_symb_vec,H_EST_PS_Sinc); stop = clock();
	// printf("Elapsed time: %f\n",(double) (stop - start));
	// printf("Processing PS MMSE...\n"); start = clock();
	// WiFi_channel_estimation_PS_MMSE(tx_symb_vec,rx_symb_vec,Fmatrix,OW2,H_EST_LT_LS,H_EST_PS_MMSE); stop = clock();
	// printf("Elapsed time: %f\n",(double) (stop - start));

	// printf("**** Total Elapsed time: %f\n",(double) (stop - start_tot));
	
	// printVect(H_EST_LT_LS,SAMPUTIL,"H_EST_LT_LS");
	// printVect(H_EST_PS_Linear,SAMPUTIL,"H_EST_PS_Linear");
	// printVect(H_EST_PS_Cubic,SAMPUTIL,"H_EST_PS_Cubic");
	// printVect(H_EST_PS_Sinc,SAMPUTIL,"H_EST_PS_Sinc");
	// printVect(H_EST_PS_MMSE,SAMPUTIL,"H_EST_PS_MMSE");
	return 0;
}

void WiFi_channel_estimation_LT_LS(long double complex tx_pre[], long double complex rx_pre[],long double complex H_EST[]){
	long double complex conj1, conj2;
	for(int i=0 ; i<26 ; i++){
		conj1 = creal(tx_pre[i]) - cimag(tx_pre[i]);
		conj2 = creal(tx_pre[i+27]) - cimag(tx_pre[i+27]);
		H_EST[i] = ( conj1*rx_pre[i] ) / ( conj1*tx_pre[i] );
		H_EST[i+27] = ( conj2*rx_pre[i+27] ) / ( conj2*tx_pre[i+27] );
	}
	H_EST[26] = 0.0;
}

void WiFi_channel_estimation_PS_Linear(long double complex tx_symbols[], long double complex rx_symbols[],long double complex H_EST[]){
	long double complex H_PILOTS[4];
	long double complex tx_pilots[4] = {tx_symbols[P0],tx_symbols[P1],tx_symbols[P2],tx_symbols[P3]};
	long double complex rx_pilots[4] = {rx_symbols[P0],rx_symbols[P1],rx_symbols[P2],rx_symbols[P3]};
	long double alpha, delta = P1 - P0;
	for(int i=0 ; i<4 ; i++){
		H_PILOTS[i] = rx_pilots[i] / tx_pilots[i];
	}

	for(int i=0 ; i<SAMPUTIL ; i++){
		if(i<P1){
			alpha = (i-P0)/delta;
			H_EST[i] = H_PILOTS[0]+( (H_PILOTS[1]-H_PILOTS[0] )*alpha );
		} else if ((i>=P1) && (i<P2)){
            alpha = (i-P1)/delta;
            H_EST[i] = H_PILOTS[1]+( (H_PILOTS[2]-H_PILOTS[1] )*alpha );
        } else if ((i>=P2) && (i<P3)){
            alpha = (i-P2)/delta;
            H_EST[i] = H_PILOTS[2]+( (H_PILOTS[3]-H_PILOTS[2] )*alpha );
        } else if (i>=P3) {
            alpha = (i-P2)/delta;
            H_EST[i] = H_PILOTS[2]+( (H_PILOTS[3]-H_PILOTS[2] )*alpha );
        }
	}
}

void WiFi_channel_estimation_PS_Cubic(long double complex tx_symbols[], long double complex rx_symbols[],long double complex H_EST[]){
	long double complex H_PILOTS[4];
	long double complex tx_pilots[4] = {tx_symbols[P0],tx_symbols[P1],tx_symbols[P2],tx_symbols[P3]};
	long double complex rx_pilots[4] = {rx_symbols[P0],rx_symbols[P1],rx_symbols[P2],rx_symbols[P3]};
	long double complex f0, f01, f12, f23, f012, f123, f0123;
	long double delta = P1 - P0;
	for(int i=0 ; i<4 ; i++){
		H_PILOTS[i] = rx_pilots[i] / tx_pilots[i];
	}
	f0 = H_PILOTS[0];
	f01   = (H_PILOTS[1]-H_PILOTS[0]) / delta;
	f12   = (H_PILOTS[2]-H_PILOTS[1]) / delta;
    f23   = (H_PILOTS[3]-H_PILOTS[2]) / delta;
    f012  = (f12-f01) / delta;
    f123  = (f23-f12) / delta;
    f0123 = (f123-f012) / delta;
    for(int k=0 ; k<SAMPUTIL ; k++){
        H_EST[k] = f0 + f01*(k-P0) + f012*(k-P0)*(k-P1) + f0123*(k-P0)*(k-P1)*(k-P2);
    }
}

void WiFi_channel_estimation_PS_Sinc(long double complex tx_symbols[], long double complex rx_symbols[],long double complex H_EST[]){
	double complex H_PILOTS[4], sinc1[SAMPUTIL],sinc2[SAMPUTIL],sinc3[SAMPUTIL],sinc4[SAMPUTIL];
	long double complex tx_pilots[4] = {tx_symbols[P0],tx_symbols[P1],tx_symbols[P2],tx_symbols[P3]};
	long double complex rx_pilots[4] = {rx_symbols[P0],rx_symbols[P1],rx_symbols[P2],rx_symbols[P3]};
	long double delta = P1 - P0;
	double a, b, c, d;

	for(int i=0 ; i<4 ; i++){
		H_PILOTS[i] = rx_pilots[i] / tx_pilots[i];
	}

	for(int k=0 ; k<SAMPUTIL ; k++){
		a = (k-P0) / delta;
		b = (k-P1) / delta;
		c = (k-P2) / delta;
		d = (k-P3) / delta;
		sinc1[k] = H_PILOTS[0]*sinc(a);
		sinc2[k] = H_PILOTS[1]*sinc(b);
		sinc3[k] = H_PILOTS[2]*sinc(c);
		sinc4[k] = H_PILOTS[3]*sinc(d);
		H_EST[k] = sinc1[k] + sinc2[k] + sinc3[k] + sinc4[k];
    }
}

void WiFi_channel_estimation_PS_MMSE(long double complex tx_symbols[], long double complex rx_symbols[], long double complex **F, double ow2, long double complex H_EST_LS[], long double complex H_EST_MMSE[]){

	long double complex **FHermitian; malloc2dLongDoubleComplex(&FHermitian, SAMPUTIL, SAMPUTIL);
	long double complex **X4Hermitian; malloc2dLongDoubleComplex(&X4Hermitian, SAMPUTIL, SAMPUTIL);
	long double complex **X4; malloc2dLongDoubleComplex(&X4, SAMPUTIL, SAMPUTIL);
	long double complex **Rhh; malloc2dLongDoubleComplex(&Rhh, SAMPUTIL, SAMPUTIL);
	long double complex **Rhy; malloc2dLongDoubleComplex(&Rhy, SAMPUTIL, SAMPUTIL);
	long double complex **Ryy; malloc2dLongDoubleComplex(&Ryy, SAMPUTIL, SAMPUTIL);
	long double complex **invRyy; malloc2dLongDoubleComplex(&invRyy, SAMPUTIL, SAMPUTIL);
	long double complex **invF; malloc2dLongDoubleComplex(&invF, SAMPUTIL, SAMPUTIL);
	long double complex **temp1; malloc2dLongDoubleComplex(&temp1, SAMPUTIL, SAMPUTIL);
	long double complex **temp2; malloc2dLongDoubleComplex(&temp2, SAMPUTIL, SAMPUTIL);
	long double complex **temp3; malloc2dLongDoubleComplex(&temp3, SAMPUTIL, SAMPUTIL);
	long double complex **Id; malloc2dLongDoubleComplex(&Id, SAMPUTIL, SAMPUTIL);

	long double complex **rx_symbols1; malloc2dLongDoubleComplex(&rx_symbols1, SAMPUTIL, 1);
	long double complex **H_EST1; malloc2dLongDoubleComplex(&H_EST1, SAMPUTIL, 1);

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
		rx_symbols1[r][0] = rx_symbols[r];
		H_EST1[r][0] = H_EST_LS[r];
	}

	hermitian(F,SAMPUTIL,SAMPUTIL,FHermitian);								// FHermitian = F'
	hermitian(X4,SAMPUTIL,SAMPUTIL,X4Hermitian);							//X4Hermitian = X4'

	inverse(F,SAMPUTIL,invF);												//invF
	multiply(invF,SAMPUTIL,SAMPUTIL,H_EST1,SAMPUTIL,1,temp1);				//temp1 = invF*H_EST
	hermitian(temp1,SAMPUTIL,SAMPUTIL,temp2);								//temp2 = (invF*H_EST)'
	multiplyVxVeqM(temp1,SAMPUTIL,SAMPUTIL,temp2,SAMPUTIL,SAMPUTIL,Rhh);	//Rhh = (invF*H_EST)*(invF*H_EST)'

	multiply(Rhh,SAMPUTIL,SAMPUTIL,FHermitian,SAMPUTIL,SAMPUTIL,temp1);		//temp1 = Rhh*F'
	multiply(temp1,SAMPUTIL,SAMPUTIL,X4,SAMPUTIL,SAMPUTIL,Rhy);				//Rhy

	multiply(temp1,SAMPUTIL,SAMPUTIL,X4Hermitian,SAMPUTIL,SAMPUTIL,temp2);	//temp2 = Rhh*F'*X4'
	multiply(F,SAMPUTIL,SAMPUTIL,temp1,SAMPUTIL,SAMPUTIL,temp1);  			//temp1 = F*Rhh*F'*X4'
	multiply(X4,SAMPUTIL,SAMPUTIL,temp1,SAMPUTIL,SAMPUTIL,temp2);			//temp2 = X4*F*Rhh*F'*X4'

	identity(Id,SAMPUTIL,ow2);
	addition(Id,SAMPUTIL,SAMPUTIL,temp2,SAMPUTIL,SAMPUTIL,Ryy);				//Ryy

	inverse(Ryy,SAMPUTIL,invRyy);											//invRyy

 	multiply(F,SAMPUTIL,SAMPUTIL,Rhy,SAMPUTIL,SAMPUTIL,temp1);				//temp1 = F*Rhy
	multiply(invRyy,SAMPUTIL,SAMPUTIL,rx_symbols1,SAMPUTIL,1,temp3);		//temp3 = invRyy*rx_symbols
	multiply(temp1,SAMPUTIL,SAMPUTIL,temp3,SAMPUTIL,1,temp2);				//H_EST = F*Rhy*invRyy*rx_symbols

	for(int r=0 ; r<SAMPUTIL ; r++)
		H_EST_MMSE[r] = temp2[r][0];

	free(FHermitian);free(X4Hermitian);free(X4);free(Rhh);free(Rhy);free(Ryy);free(invRyy);
	free(invF);free(temp1);free(temp2);free(Id);free(temp3);free(rx_symbols1);free(H_EST1);
}