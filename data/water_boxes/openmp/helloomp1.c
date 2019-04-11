#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
int main () {

	int num_thds, tid;
	#pragma omp parallel 
 	{
		num_thds = omp_get_num_threads();
	 	tid = omp_get_thread_num();
 		printf("Hello World from thread number %d out of %d threads\n ",tid,num_thds);

// 	if (tid == 0) 
// 	{
// 	nthreads = omp_get_num_threads();
// 	printf("Number of threads = %d\n", nthreads);
// 	}

 	} 

}
