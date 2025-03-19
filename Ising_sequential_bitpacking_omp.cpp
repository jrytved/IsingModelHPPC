// include basic libraries for I/O, math and array manipulation

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <random>
#include <chrono>
#include <cstdint>


const int N = 2048;                 // NUMBER OF ROWS
const int M = N/2;                 // NUMBER OF COLUMNS OF EITHER COLOR
const int K = N*M;                 // SIZE OF 1D VECTOR FOR *EACH* OF BLACK AND WHITE
const int SEED = 1234;             // RANDOM SEED
const int FWRITE = 2;             // WRITE FREQUENCY
const double TEMP = 0.1;           // TEMPERATURE (t)
const int STEPS = 500;             // NUMBER OF SIMULATION STEPS
const int H = 0;                   // CONSTANT FOR EXTERNAL MAGNETIC FIELD (h/J)
const std::string OUTFILE = "";

struct Lattice {

    std::vector<uint64_t> bits;

    size_t size;

     
    Lattice(size_t N) : size(N) {
        bits.resize((N + 63) / 64); // We need at *least* N/64 integers to store N bits. Integer division truncates, so add 63.
    }

    inline bool get(size_t i) const {
        // To get bit at index i, find the relevant integer(i/64). 
        // Shift all the bits to the right by (i%64) to isolate the bit and make a logical AND.
        return (bits[i / 64] >> (i % 64)) & 1; 
    }

    inline void set(size_t i, bool value) {
        //Use 1ULL to create a bit mask, that sets bit at (i%64) to 1.
        // If we want tot set to 0, use OR
        if (value)
            bits[i / 64] |= (1ULL << (i % 64)); // Set bit
        // If we want to set to 1, use XOR
        else
            bits[i / 64] &= ~(1ULL << (i % 64)); // Clear bit
    }

    inline void flip(size_t i) {
        // Use XOR to flip the bit. 0 ^ 1 = 1; 1 ^ 1 = 0;
        bits[i / 64] ^= (1ULL << (i % 64)); // Flip bit
    }
};


Lattice create_lattice(std::mt19937& generator, size_t num_bits) {
    
    Lattice lattice(num_bits);

    std::uniform_int_distribution<int> distribution(0, 1);
    
    for (size_t i = 0; i < num_bits; ++i) {
        lattice.set(i, distribution(generator));
    }
    return lattice;
}

void evolve_lattice(std::mt19937& generator, Lattice &white, Lattice &black){
    // update all white
    
    for (size_t i = 0; i < N; ++i){
        for (size_t j = 0; j < M; ++j){
            // i,j are now in local white space
            // transform to global I,J indices
            size_t I = i;
            size_t J = 2*j + (i%2);
            // get neighbors in a single step without conditionals
            size_t leftJ  = (J-1+N) % N;
            size_t rightJ = (J+1)   % N;
            size_t upI    = (I-1+N) % N;
            size_t downI  = (I+1)   % N;

            int left_black  = 2*black.get(I*M     + (leftJ  - 1 + (I%2))/2)-1;
            int right_black = 2*black.get(I*M     + (rightJ - 1 + (I%2))/2)-1;
            int up_black    = 2*black.get(upI*M   + (J - 1 + (upI  %2))/2)-1;
            int down_black  = 2*black.get(downI*M + (J - 1 + (downI%2))/2)-1;
            int spin        = 2*white.get(i*M     + j)-1;
            // energy difference
            double dE = 2.0 * spin * (left_black + right_black + up_black + down_black + H);
            // decide based on a generator with minimal conditionals
            double boltzmann = exp(-dE/TEMP);
            std::uniform_real_distribution<double> distribution(0.0, 1.0);
            if (distribution(generator) < boltzmann){
                white.flip(i*M + j); // Flip that bit!
            }
        } 
    }
    // update all black
    for (size_t i = 0; i < N; ++i){
        for (size_t j = 0; j < M; ++j){
            // i,j are now in local black space
            // transform to global I,J indices
            size_t I = i;
            size_t J = 2*j + 1 - (i%2);
            // get neighbors in a single step without conditionals
            size_t leftJ  = (J-1+N) % N;
            size_t rightJ = (J+1)   % N;
            size_t upI    = (I-1+N) % N;
            size_t downI  = (I+1)   % N;
            //
            int left_white  = 2*white.get(I*M     + (leftJ  - (I%2))/2)-1;
            int right_white = 2*white.get(I*M     + (rightJ - (I%2))/2)-1;
            int up_white    = 2*white.get(upI*M   + (J - (upI  %2))/2)-1;
            int down_white  = 2*white.get(downI*M + (J - (downI%2))/2)-1;
            int spin        = 2*black.get(i*M     + j)-1;
            
            // energy difference
            double dE = 2.0 * spin * (left_white + right_white + up_white + down_white + H);
            // decide based on a generator with minimal conditionals
            double boltzmann = exp(-dE/TEMP);
            std::uniform_real_distribution<double> distribution(0.0, 1.0);
            if (distribution(generator) < boltzmann){
                black.flip(i*M + j); // Flip that bit!
            }
        } 
    }
}

void join_arrays(const Lattice &white, const Lattice &black, Lattice &global) {
    // Here we assume the global lattice is arranged in N rows, each with 2*M columns.
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < M; j++) {
            // Global index: row i, column (2*j + offset)
            global.set(i * (2 * M) + (2 * j + (i % 2)), white.get(i * M + j));
            global.set(i * (2 * M) + (2 * j + 1 - (i % 2)), black.get(i * M + j));
        }
    }
}

void write_to_file(const Lattice &global) {
    std::ofstream file(OUTFILE, std::ios::app);
    if (!file) {
        std::cerr << "Error opening file for writing!" << std::endl;
        return;
    }
    for (size_t i = 0; i < 2 * K; ++i) {
        file << global.get(i) << " ";
    }
    file << "\n";
    file.close();
}

void calculate_quantities(const Lattice &global, double& magnetization, double& energy){
    
    magnetization = 0.0;
    energy = 0.0;
    
    for (size_t i = 0; i < N; ++i){
        for (size_t j = 0; j < N; ++j){
            magnetization += global.get(i*N + j);
            energy += -global.get(i*N + j) * (global.get(((i+1)%N)*N + j) + global.get(((i-1+N)%N)*N + j)
                                        + global.get(i*N + (j+1)%N) + global.get(i*N + (j-1+N)%N) + H);
        }
    }
    magnetization /= N*N;
    magnetization = fabs(magnetization);
    energy /= N*N;
}

void save_quantities(const std::vector<double>& magnetization, const std::vector<double>& energy){
    // save magnetization and energy
    std::string filename_energy = "energy_" + OUTFILE;
    std::string filename_magnetization = "magnetization_" + OUTFILE;
    // save the magnetization
    std::ofstream file_magnetization(filename_magnetization);
    for (size_t s = 0; s < STEPS; ++s){
        file_magnetization << s << " " << magnetization[s] << "\n";
    }
    file_magnetization.close();
    // save the energy
    std::ofstream file_energy(filename_energy);
    for (size_t s = 0; s < STEPS; ++s){
        file_energy << s << " " << energy[s] << "\n";
    }
    file_energy.close();
}

int main(int argc, char* argv[]){
    // set up random generator
    std::mt19937 generator(SEED);

    // Create white lattice
    Lattice white = create_lattice(generator, K);
    // Create black lattice
    Lattice black = create_lattice(generator, K);
    // Create global lattice
    Lattice global(2 * K);
    
    // create containers for magnetization and energy
    std::vector<double> magnetization(STEPS);
    std::vector<double> energy(STEPS);

    // delete files if they exist already
    std::remove(OUTFILE.c_str());
    std::remove(("energy_" + OUTFILE).c_str());
    std::remove(("magnetization_" + OUTFILE).c_str());

    auto begin = std::chrono::steady_clock::now();
    
    // loop over steps
    for (size_t s = 0; s < STEPS; ++s){
        // evolve lattice
        evolve_lattice(generator, white, black);
        // collect in global
        join_arrays(white, black, global);
        // calculate energy and magnetization
        calculate_quantities(global, magnetization[s], energy[s]);        
        // save to file
        if (s % FWRITE == 0 && !OUTFILE.empty()){
            write_to_file(global);
        }
    }
    auto end = std::chrono::steady_clock::now();
    
    // save magnetization and energy
    save_quantities(magnetization, energy);    

    std::cout << "Elapsed Time\t" << (end - begin).count() / 1000000000.0 << " sec" << std::endl;
    
    return 0;
}