#include <iostream>
#include <span>
#include <vector>
#include <mpi.h>


// Function that takes a std::span to view elements without copying them
void printNumbers(std::span<int> numbers) {
    for (int num : numbers) {
        std::cout << num << " ";
    }
    std::cout << std::endl;
}

int main() {

    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    std::vector<int> numbers = {1, 2, 3, 4, 5};

    // Use std::span to pass the vector to the function
    printNumbers(numbers);

    MPI_Finalize();

    return 0;
}

