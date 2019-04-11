#include <stdio.h>
#include <omp.h>

/*
 * This source code can be downloaded from supercomputingblog.com
 * The purpose of this code is to provide a basic understanding of OpenMP.
 * */

int main(int argc, char* argv[])
{
//	 This statement should only print once
	 printf("Starting Program!\n");
	
	 int nThreads, tid;
	
	 #pragma omp parallel private(tid)
	 {
	 // This statement will run on each thread.
	 // If there are 4 threads, this will execute 4 times in total

	 tid = omp_get_thread_num();
	 printf("Running on thread %d\n", tid);
	
		 if (tid == 0)
	 	{
	 	nThreads = omp_get_num_threads();
	 	printf("Total number of threads: %d\n", nThreads);
	 	}
         }
	 																										
	 // We're out of the parallelized secion.
	 // Therefor, this should execute only once
	 printf("Finished!\n");
	
	 return 0;
	 																													
}                                                                                                                                                                                                                                                  
