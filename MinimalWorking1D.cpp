#include <iostream>
#include <random>
#include <math.h>
#include <fstream>
#include <string>

const int N = 100;                  // NUMBER OF ROWS
const int M = 100;                  // NUMBER OF COLUMNS
const int K = M*N;
const double TEMP = 0.5;           // TEMPERATURE
const double INV_TEMP = 1 / TEMP;  // INVERSE OF THE TEMPERATURE
const int STEPS = 20;             // NUMBER OF SIMULATION STEPS
const int J = 1;                   // INTERACTION CONSTANT
const int H = 0;                   // CONSTANT FOR EXTERNAL MAGNETIC FIELD
const std::string OUTFILE = "OUTPUT.jade";

double random_double() {
    
    static std::random_device rd;                                 // Seed generator
    static std::mt19937 gen(rd());                                // Mersenne Twister PRNG
    static std::uniform_real_distribution<double> dist(0.0, 1.0); // Uniform distribution
    return dist(gen);
    
}


void init_array(int (&arr)[K]) {

    // Initializes the array randomly with (i,j) IN {-1,1}
    
    std::random_device rd;  // Random seed
    std::mt19937 gen(rd()); // Mersenne Twister PRNG
    std::uniform_int_distribution<int> dist(0, 1); // Generates 0 or 1
    
    for (int i = 0; i < K; ++i) {
            arr[i] = dist(gen) * 2 - 1; // Convert {0,1} to {-1,1}
    }
};

// Function to write the 1D array to a file and append each step

void write_array_to_file(int (&arr)[K], int step) {
    std::ofstream file(OUTFILE, std::ios::app);  // Open the file in append mode
    if (!file) {
        std::cerr << "Error opening file for writing!" << std::endl;
        return;
    }
    // Write the array as a flat 1D array to the file
    for (int i = 0; i < K; ++i) {
        file << arr[i] << " ";  // Write each element followed by a space
    }
    file << "\n";  // Newline after each step's array

    file.close();  // Close the file
}


void step(int (&arr)[K]){

    int ipp, jpp, inn, jnn, nb, spin, dE;

    for(int i=0;i<N;i++)
    for(int j=0;j<M;j++){

        ipp = ((i+1) < N)? (i+1) : 0;
        jpp = ((j+1) < M)? (j+1) : 0;
        
        inn = ((i-1) >= 0)? (i-1) : (N - 1);
        jnn = ((j-1) >= 0)? (j-1) : (M - 1);

        // CALCULATE NEIGHBORS
        nb = arr[ipp*M+j] + arr[i*M+jpp] + arr[inn*M+j] + arr[i*M+jnn];

        // COMPUTE ENERGY DIFFERENCE
        spin = arr[i*M+j];
        dE = -2*spin*(J * nb + H);

        if(dE > 0){
            arr[i*M+j] = -spin;
        } else {
            
            if(random_double() < exp(dE*INV_TEMP)){
                arr[i*M+j] = -spin;
            } 
        } 
    }
};

int main() {
    int arr[K];       // Declare the 1D array
    init_array(arr);  // Initialize the array with random values

    for (int s = 0; s < STEPS; s++) {
        step(arr);  // Perform one simulation step
        write_array_to_file(arr, s);  // Append the current state of the array to the file
    };

    return 0;
}



