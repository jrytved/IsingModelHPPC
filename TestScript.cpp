// include basic libraries for I/O, math and array manipulation

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <random>
#include <chrono>
#include <bitset>


const int N = 128;                 // NUMBER OF ROWS
const int M = N/2;                 // NUMBER OF COLUMNS OF EITHER COLOR
const int K = N*M;                 // SIZE OF 1D VECTOR FOR *EACH* OF BLACK AND WHITE
const int SEED = 1234;             // RANDOM SEED
const int FWRITE = 2;             // WRITE FREQUENCY
const double TEMP = 0.1;           // TEMPERATURE (t)
const int STEPS = 500;             // NUMBER OF SIMULATION STEPS
const int H = 0;                   // CONSTANT FOR EXTERNAL MAGNETIC FIELD (h/J)
const std::string OUTFILE = "OUTPUT.jade";



std::bitset<K>* create_lattice(std::mt19937& generator){
    std::bitset<K>* lattice = new std::bitset<K>;
    
    // fill in random numbers    
    std::uniform_int_distribution<int> distribution(0,1);
    
    for (size_t i = 0; i < K; ++i){
        
        lattice->set(i, distribution(generator)); // Sets the bit at idx. i to what the generator returns
    }

    return lattice;
};
  

int main(void){

    std::mt19937 generator(SEED);
    // create white lattice
    std::bitset<K>* white = create_lattice(generator);

    for(size_t i=0; i<K; i++){
        std::cout << (*white).test(i);
    }
}