## Introduction

In any wireless communication, the impact of the channel on the received signal needs to always be taken into account prior demodulation. Channels can have different distribution according to the scenario where the communication happens. 802.11 is a wireless standard meant to be deployed in indoors, where the multiple reflections and the short delay between them define a slow flat fading channel.

In order to mitigate the effect of the channel, two extra stages on the receiver chain are added: channel estimator and signal equalization. The first one estimates a channel impulse response based on the received known symbols whereas the second one uses that estimation to increase the accuracy at a symbol level.

This project shows the improvement introduced in some of the channel estimators used in the 802.11 standard by exploiting different forms of parallelism. Two forms of parallelism are explored: OpenMP and MPI. A C and Matlab sequential code are built as a base line, and thy are used to compare and justify the improvement on the performance. Substantial improvement is obtained for the most complex estimator, the MMSE, for all parallel methods implemented.

The parallelization is implemented in C and Matlab, using the respective resources available (OpenMP and MPI in C and Parallel Toolbox in Matlab).

## Code in C

- **compile.c**: Program that allow the user to compile the code (Sequential, OpenMP and MPI) in an easier way, by just executing one command on terminal.
- **main.c**: Contains the main program. 5 different channel estimation techniques are implemented and for each, the output generated reflects the channel frequency response.
- **main_openmp.c**: Contains the main program in OpenMP. 5 different channel estimation techniques are implemented and for each, the output generated reflects the channel frequency response.
- **main_openmp.bash**: Bash script of the main_openmp.c file.
- **main_mpi.c**: Contains the main program in MPI. 5 different channel estimation techniques are implemented and for each, the output generated reflects the channel frequency response.
- **main_mpi.bash**: Bash script of the main_mpi.c file.
- **inputs.h**: Contains the inputs of the main program. tx_preamble_fft and rx_preamble_fft represent the transmitted and received preambles in the frequency domain. tx_symb and rx_symb represent a vector of length 53*15, where 53 is the number of information in each OFDM block and 15 is the number of OFDM blocks within a frame. OW2 represents an estimation of the received noise power.
- **utils.h**: Contains the declaration of the main function used to manipulate matrices and operate with them. 
- **utils.c**: Contains the implementation of the main function used to manipulate matrices and operate with them.

## Code in MATLAB

- **WiFi_RX: Main program**. It's inputs are the transmitted and received preamble and frame in the time domain. It performs 5 different channel estimation techniques and outputs the channel frequency response in each. It also includes an equalization stage, where the channel estimation is applied before the symbols are decoded.
- **WiFi_Inputs**: Contains the inputs of the main program. rx_packet and tx_packet are the time domain samples of the received frames without the preamble. tx_lptot and rx_lptot are the transmitted and received preambles, used for channel estimation porpuses only.
- **WiFi_block_extraction**: It transforms the received samples in the time domain into symbols.
- **WiFi_Equelization**: Equalizes the received symbols dividing them by the estimated Channel Frequency Response. 
- **WiFi_channel_estimation_LT_LS**: Performs a Least squares (LS) Channel Estimation using the Training Symbols (TS) and outputs the Channel Frequency Response.
- **WiFi_channel_estimation_PS_Linear**: Performs a Least squares (LS) Channel Estimation using the Pilot Subcarriers (PS) and interpolates the result using a Linear Interpolation. Outputs the channel Frequency Response.
- **WiFi_channel_estimation_PS_Cubic**: Performs a Least squares (LS) Channel Estimation using the Pilot Subcarriers (PS) and interpolates the result using a Cubic Spline Interpolation. Outputs the channel Frequency Response.
- **WiFi_channel_estimation_PS_Sinc**: Performs a Least squares (LS) Channel Estimation using the Pilot Subcarriers (PS) and interpolates the result using a Sinc Interpolation. Outputs the channel Frequency Response.
- **WiFi_channel_estimation_PS_MMSE**: Performs a Minimum Mean Square Error (MMSE) Channel Estimation using the Pilot Subcarriers (PS) and the LT_LS estimation. Outputs the channel Frequency Response.

## Compilation and execution

Compilation of compile.c file:
> g++ -w -o compile compile.c

Compilation and Execution of the Sequential Implementation
> ./compile sequential
> ./main

Compilation and Execution of the OpenMP Implementation:
> ./compile openmp
> ./main_openmp.bash

Compilation and Execution of the MPI Implementation:
> ./compile mpi
> ./main_mpi.bash

## IMPORTANT NOTES: 

1. MPI functions are included in utils.c, forcing the 3 implementations to be compiled using the mpiCC compiler. The user is endorsed to use the compile function due to potential
oversights.
2. MMSE Implementation using MPI balances the process-load automatically, depending on the number of processes that we run the execution on. The same process still has to be done for 
the rest of the functions, which are configured for 20 processes.
