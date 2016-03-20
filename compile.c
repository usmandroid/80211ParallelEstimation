#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

/* --------------------- IMPORTANT NOTE -----------------------
* This code links libraries, compiles C code and generates the 
* executable file for the Sequential, OpenMP and MPI implementation 
* of channel estimation. 
* The code has to be executed with 1 parameter, defining the code being
* compiled. The three options are "sequential", "openmp" and "mpi"
* Example: 
* 	$ g++ -w -o compile compile.c
* 	$ ./compile Sequential
* 
* For the parallel versions, the code may be run using the bash scripts.
* ------------------------------------------------------------- */

int main(int argc, char *argv[]) {
	char usage[] = "/* --------------------- IMPORTANT NOTE ----------------------- \n \
* This code links libraries, compiles C code and generates the \n \
* executable file for the Sequential, OpenMP and MPI implementation \n \
* of channel estimation. \n \
* The code has to be executed with 1 parameter, defining the code being \n \
* compiled. The three options are \"sequential\", \"openmp\" and \"mpi\" \n \
* Example: \n \
* 	$ g++ -w -o compile compile.c \n \
* 	$ ./compile Sequential \n \
* \n \
* For the parallel versions, the code may be run using the bash scripts. \n \
* ------------------------------------------------------------- */ \n \n";

	printf("%s",usage);
	int status;
	char *const parmList1[] = {"/usr/bin/g++", "-w", "-o", "utils.o", "-c", "utils.c",NULL};
	char *const parmList_Seq[] = {"/usr/bin/g++", "-w", "-o", "main", "main.c", "utils.o",NULL};
	char *const parmList_OpenMP[] = {"/usr/bin/g++", "-fopenmp", "-w", "-o", "main_openmp", "main_openmp.c", "utils.o",NULL};
	char *const parmList_MPI[] = {"/opt/ibm/platform_mpi/bin/mpiCC", "-o", "main_mpi", "main_mpi.c", "utils.o",NULL};
	
	if(argc < 2){
		printf("** error: Not enough input arguments\n");
		printf("** error: Argument must be unique and either \"Sequential\", \"OpenMP\" or \"MPI\" \n");
		return -1;
	} else if(strcmp(argv[1],"sequential") && strcmp(argv[1],"openmp") && strcmp(argv[1],"mpi")) {
		printf("Invalid Argument\n");
		printf("Argument must be unique and either \"sequential\", \"openmp\" or \"mpi\" \n");
		return -1;
	}

	pid_t pid = fork();
	if(pid == 0){
		execv("/usr/bin/g++", parmList1);
	}
	else {
		waitpid(0,&status, 0);
		if(!strcmp(argv[1],"sequential")) {
			printf("Compiling Sequential Code... Executable will be stored in \"main\" \n\n");
			execv("/usr/bin/g++", parmList_Seq);
		}
		if(!strcmp(argv[1],"openmp")){
			printf("Compiling OpenMP Code... Executable will be stored in \"main_openmp\" \n\n");
			execv("/usr/bin/g++", parmList_OpenMP);
		}
		if(!strcmp(argv[1],"mpi")) {
			printf("Compiling MPI Code... Executable will be stored in \"main_mpi\" \n\n");
			execv("/opt/ibm/platform_mpi/bin/mpiCC", parmList_MPI);
		}
	}
	return 0;
}