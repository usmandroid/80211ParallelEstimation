#include "inputs.h"
#include "utils.h"

void WiFi_channel_estimation_LT_LS(double complex tx_pre[], double complex rx_pre[],double complex H_EST[]);
void WiFi_channel_estimation_PS_Linear(double complex tx_symbols[], double complex rx_symbols[],double complex H_EST[]);
void WiFi_channel_estimation_PS_Cubic(double complex tx_symbols[], double complex rx_symbols[],double complex H_EST[]);
void WiFi_channel_estimation_PS_Sinc(double complex tx_symbols[], double complex rx_symbols[],double complex H_EST[]);
void WiFi_channel_estimation_PS_MMSE(double complex tx_symbols[], double complex rx_symbols[], double complex **F, double ow2, double complex H_EST[]);

int main(void) {

	double complex tx_symb_vec[SAMPUTIL],rx_symb_vec[SAMPUTIL];
	double complex H_EST_LT_LS[SAMPUTIL],H_EST_PS_Linear[SAMPUTIL],H_EST_PS_Cubic[SAMPUTIL];
	double complex H_EST_PS_MMSE[SAMPUTIL];

	double complex **Fmatrix = new double complex*[SAMPUTIL];
	for (int i = 0; i < SAMPUTIL; i++) {
		Fmatrix[i] = new double complex[SAMPUTIL];
	}
	for (int f=0; f<SAMPUTIL; f++){
        for (int t=0; t<SAMPUTIL; t++){
            Fmatrix[t][f] = cexp(-2*I*PI*t*f/SAMPUTIL);
		}
    }

	WiFi_channel_estimation_LT_LS(tx_preamble_fft,rx_preamble_fft,H_EST_LT_LS); 

	/*
	for(int c=0 ; c<OFDMBLK ; c++){
		for(int r=0 ; r<SAMPUTIL ; r++){
			tx_symb_vec[r] = tx_symb[SAMPUTIL*c + r];
			rx_symb_vec[r] = rx_symb[SAMPUTIL*c + r];
		}
		WiFi_channel_estimation_PS_Linear(tx_symb_vec,rx_symb_vec,H_EST_PS_Linear);
		WiFi_channel_estimation_PS_Cubic(tx_symb_vec,rx_symb_vec,H_EST_PS_Cubic);
		WiFi_channel_estimation_PS_Sinc(tx_symb_vec,rx_symb_vec,H_EST_PS_Cubic);
		WiFi_channel_estimation_PS_MMSE(tx_symb_vec,rx_symb_vec,Fmatrix,OW2,H_EST_LT_LS);
	}*/
	for(int r=0 ; r<SAMPUTIL ; r++){
		tx_symb_vec[r] = tx_symb[r];
		rx_symb_vec[r] = rx_symb[r];
	}
	WiFi_channel_estimation_PS_MMSE(tx_symb_vec,rx_symb_vec,Fmatrix,OW2,H_EST_LT_LS);
    
	return 0;
}

void WiFi_channel_estimation_LT_LS(double complex tx_pre[], double complex rx_pre[],double complex H_EST[]){
	double complex conj1, conj2;
	for(int i=0 ; i<26 ; i++){
		conj1 = creal(tx_pre[i]) - cimag(tx_pre[i]);
		conj2 = creal(tx_pre[i+27]) - cimag(tx_pre[i+27]);
		H_EST[i] = ( conj1*rx_pre[i] ) / ( conj1*tx_pre[i] );
		H_EST[i+27] = ( conj2*rx_pre[i+27] ) / ( conj2*tx_pre[i+27] );
	}
	H_EST[26] = 0.0;
}

void WiFi_channel_estimation_PS_Linear(double complex tx_symbols[], double complex rx_symbols[],double complex H_EST[]){
	double complex H_PILOTS[4];
	double complex tx_pilots[4] = {tx_symbols[6],tx_symbols[20],tx_symbols[34],tx_symbols[48]};
	double complex rx_pilots[4] = {rx_symbols[6],rx_symbols[20],rx_symbols[34],rx_symbols[48]};
	int alpha, delta = P1 - P0;
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
        } else if (i>=P2) {
            alpha = (i-P2)/delta;
            H_EST[i] = H_PILOTS[2]+( (H_PILOTS[3]-H_PILOTS[2] )*alpha );
        }
	}
}

void WiFi_channel_estimation_PS_Cubic(double complex tx_symbols[], double complex rx_symbols[],double complex H_EST[]){
	double complex H_PILOTS[4];
	double complex tx_pilots[4] = {tx_symbols[6],tx_symbols[20],tx_symbols[34],tx_symbols[48]};
	double complex rx_pilots[4] = {rx_symbols[6],rx_symbols[20],rx_symbols[34],rx_symbols[48]};	double complex f0, f01, f12, f23, f012, f123, f0123;
	int alpha, delta = P1 - P0;
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

void WiFi_channel_estimation_PS_Sinc(double complex tx_symbols[], double complex rx_symbols[],double complex H_EST[]){
	double complex H_PILOTS[4],sinc1[SAMPUTIL],sinc2[SAMPUTIL],sinc3[SAMPUTIL],sinc4[SAMPUTIL];
	double complex tx_pilots[4] = {tx_symbols[6],tx_symbols[20],tx_symbols[34],tx_symbols[48]};
	double complex rx_pilots[4] = {rx_symbols[6],rx_symbols[20],rx_symbols[34],rx_symbols[48]};	double complex f0, f01, f12, f23, f012, f123, f0123;
	int alpha, delta = P1 - P0;
	for(int k=0 ; k<SAMPUTIL ; k++){
		sinc1[k] = H_PILOTS[0]*(sin(PI*(k-P0) / delta)/(PI*(k-P0) / delta));
		sinc2[k] = H_PILOTS[1]*(sin(PI*(k-P1) / delta)/(PI*(k-P1) / delta));
		sinc3[k] = H_PILOTS[2]*(sin(PI*(k-P2) / delta)/(PI*(k-P2) / delta));
		sinc4[k] = H_PILOTS[3]*(sin(PI*(k-P3) / delta)/(PI*(k-P3) / delta));

		H_PILOTS[k] = sinc1[k] + sinc2[k] + sinc3[k] + sinc4[k];
    }
}

void WiFi_channel_estimation_PS_MMSE(double complex tx_symbols[], double complex rx_symbols[], double complex **F, double ow2, double complex H_EST[]){	
	double complex **FHermitian = new double complex*[SAMPUTIL];
	double complex **X4Hermitian = new double complex*[SAMPUTIL];
	double complex **X4 = new double complex*[SAMPUTIL];
	double complex **Rhh = new double complex*[SAMPUTIL];
	double complex **Rhy = new double complex*[SAMPUTIL];
	double complex **Ryy = new double complex*[SAMPUTIL];
	double complex **invRyy = new double complex*[SAMPUTIL];
	double complex **invF = new double complex*[SAMPUTIL];
	double complex **temp1 = new double complex*[SAMPUTIL];
	double complex **temp2 = new double complex*[SAMPUTIL];
	double complex **Id = new double complex*[SAMPUTIL];
	double complex **temp3 = new double complex*[SAMPUTIL];
	double complex **rx_symbols1 = new double complex*[SAMPUTIL];
	double complex **H_EST1 = new double complex*[SAMPUTIL];
	rx_symbols1[0] = new double complex[1];
	H_EST1[0] = new double complex[1];

	for (int i = 0; i < SAMPUTIL; i++) {
		FHermitian[i] = new double complex[SAMPUTIL];
		X4Hermitian[i] = new double complex[SAMPUTIL];
		X4[i] = new double complex[SAMPUTIL];
		Rhh[i] = new double complex[SAMPUTIL];
		Rhy[i] = new double complex[SAMPUTIL];
		Ryy[i] = new double complex[SAMPUTIL];
		invRyy[i] = new double complex[SAMPUTIL];
		invF[i] = new double complex[SAMPUTIL];
		temp1[i] = new double complex[SAMPUTIL];
		temp2[i] = new double complex[SAMPUTIL];
		temp3[i] = new double complex[SAMPUTIL];
		Id[i] = new double complex[SAMPUTIL];
		rx_symbols1[i] = new double complex[1];
	}

	for(int i=0 ; i<SAMPUTIL; i++){
		rx_symbols1[i][0] = rx_symbols[i];
	}

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

	hermitian(F,FHermitian);
	hermitian(X4,X4Hermitian);

	
	int size = 55;
	double complex **Matrix = new double complex*[size];
	double complex **invMatrix = new double complex*[size];
	for(int i=0 ; i<size ; i++){
		Matrix[i] = new double complex[size];
		invMatrix[i] = new double complex[size];
		for(int j=0 ; j<size ; j++){
			Matrix[i][j] = ((int) rand()%20) + I*((int) rand()%20);
		}
	}
	double complex det = (double complex) determinant_impl(Matrix,size);
	printf("%lf + i%lf\n", creal(det), cimag(det));	

	inverse(Matrix,size,invMatrix);		//invF
//	inverse(F,SAMPUTIL,invF);			//invF
//	multiply(invF,H_EST1,temp1);		//temp1 = invF*H_EST
/*	hermitian(temp1,temp2);				//temp2 = (invF*H_EST)'
	multiply(temp1,temp2,Rhh);			//temp1 = invF*H_EST

	multiply(Rhh,FHermitian,temp1);		//temp1 = Rhh*F'
	multiply(temp1,X4,Rhy);				//Rhy

	multiply(temp1,X4Hermitian,temp2);	//temp2 = Rhh*F'*X4'
	multiply(F,temp1,temp1);  			//temp1 = F*Rhh*F'*X4'
	multiply(X4,temp1,temp2);			//temp2 = X4*F*Rhh*F'*X4'

	identity(Id,SAMPUTIL,ow2);
	addition(Id,temp2,Ryy);				//Ryy
	inverse(Ryy,SAMPUTIL,invRyy);		//invRyy

	multiply(F,Rhy,temp1);				//temp1 = F*Rhy
	multiply(invRyy,rx_symbols1,temp3);	//temp2 = invRyy*rx_symbols
	multiply(temp1,temp3,H_EST1);		//H_EST = F*Rhy*invRyy*rx_symbols*/

	free(FHermitian);free(X4Hermitian);free(X4);free(Rhh);free(Rhy);free(Ryy);free(invRyy);
	free(invF);free(temp1);free(temp2);free(Id);free(temp3);free(rx_symbols1);free(H_EST1);
}