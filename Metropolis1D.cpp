#include <iostream>
#include <random>
#include <math.h>
#include <fstream>
#include <string>

/*
    METROPOLIS ALGORITHM FOR ISING MODEL ON 2D ARRAY REPRESENTED IN 1D VECTOR.
    
*/

const int N = 100;                  // NUMBER OF ROWS
const int M = 100;                  // NUMBER OF COLUMNS
const int K = M*N/2;                // SIZE OF 1D VECTOR FOR *EACH* OF BLACK AND WHITE 
const double TEMP = 0.5;           // TEMPERATURE
const double INV_TEMP = 1 / TEMP;  // INVERSE OF THE TEMPERATURE
const int STEPS = 5000;                                             // NUMBER OF SIMULATION STEPS
const int J = 1;                                                  // INTERACTION CONSTANT
const int H = 0;                                                  // CONSTANT FOR EXTERNAL MAGNETIC FIELD
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

void write_array_to_file(int (&arr)[M*N]) {
    std::ofstream file(OUTFILE, std::ios::app);  // Open the file in append mode
    if (!file) {
        std::cerr << "Error opening file for writing!" << std::endl;
        return;
    }
    // Write the array as a flat 1D array to the file
    for (int i = 0; i < M*N; ++i) {
        file << arr[i] << " ";  // Write each element followed by a space
    }
    file << "\n";  // Newline after each step's array

    file.close();  // Close the file
}

void compute_spins(int (&arr)[K], int (&other_arr)[K], bool FIRST_IS_BLACK){

    int ipp, jpp, inn, jnn, neighbor_spin_sum, spin, dE, joff;

    for(int i=0;i<N;i++)
    for(int j=0;j<M/2;j++){

        ipp = ((i+1) < N)? (i+1) : 0;
        jpp = ((j+1) < M/2)? (j+1) : 0;
        
        inn = ((i-1) >= 0)? (i-1) : (N - 1);
        jnn = ((j-1) >= 0)? (j-1) : (M/2 - 1);

        if (FIRST_IS_BLACK){
            joff = (i%2)? jpp : jnn;
        } else {
            joff = (i%2)? jnn : jpp;
        };

        // CALCULATE NEIGHBORS
        neighbor_spin_sum = other_arr[inn*M+j] + other_arr[i*M+j] + other_arr[ipp*M+j] + other_arr[i*M+joff];

        // COMPUTE ENERGY DIFFERENCE
        spin = arr[i*M/2+j];
        
        dE = -2*spin*(J * neighbor_spin_sum + H);

        if(dE > 0){
            arr[i*M/2+j] = -spin;
        } else {
            
            if(random_double() < exp(dE*INV_TEMP)){
                arr[i*M/2+j] = -spin;
            } 
        } 
    }
};

void join_arrays(int(&arr_black)[K], int(&arr_white)[K], int(&arr_global)[M*N]){

    for(int i=0; i<N; i++)
    for(int j=0; j<M/2; j++){
        
        if(i % 2){
            arr_global[i*M+(2*j+1)] = arr_black[i*M/2+j];
            arr_global[i*M+(2*j)] = arr_white[i*M/2+j];
        } else {
            arr_global[i*M+(2*j)] = arr_black[i*M/2+j];
            arr_global[i*M+(2*j+1)] = arr_white[i*M/2+j];
        };
    }
}



int main(int argc, char* argv[]) {
    int WHITE[K];       // Declare 1D arrays of size of K:= M*(N/2)
    int BLACK[K];
    int RANDOM_VALUES[K];

    int GLOBAL_ARRAY[M*N]; // Declare one array that keeps track of everything (black+white)
    
    init_array(WHITE);  // Initialize the arrays with random values
    init_array(BLACK); 
    init_array(RANDOM_VALUES);


    for (int s = 0; s < STEPS; s++) {
        init_array(RANDOM_VALUES);    // Generate Random Values For Stochastic Spin Flip
        compute_spins(BLACK, WHITE, true);
        init_array(RANDOM_VALUES);    // Generate Random Values For Stochastic Spin Flip
        compute_spins(WHITE, BLACK, false);
        join_arrays(BLACK, WHITE, GLOBAL_ARRAY);
        
        write_array_to_file(GLOBAL_ARRAY);  // Append the current state of the array to the file
    };

    return 0;
}



