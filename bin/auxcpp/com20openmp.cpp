#include <iostream>
#include <omp.h>

int main() {
    // Set the number of threads you want to use
    omp_set_num_threads(4);

    // Parallel region starts here
    #pragma omp parallel
    {
        // Get the ID of the current thread
        int thread_id = omp_get_thread_num();
        
        // Get the total number of threads
        int num_threads = omp_get_num_threads();
        
    }
    
    return 0;
}

