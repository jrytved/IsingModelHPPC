// include basic libraries for I/O, math and array manipulation

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <random>
#include <chrono>
#include <bitset>


const int N = 2048;                 // NUMBER OF ROWS
const int M = N/2;                 // NUMBER OF COLUMNS OF EITHER COLOR
const int K = N*M;                 // SIZE OF 1D VECTOR FOR *EACH* OF BLACK AND WHITE
const int SEED = 1234;             // RANDOM SEED
const int FWRITE = 2;             // WRITE FREQUENCY
const double TEMP = 0.1;           // TEMPERATURE (t)
const int STEPS = 500;             // NUMBER OF SIMULATION STEPS
const int H = 0;                   // CONSTANT FOR EXTERNAL MAGNETIC FIELD (h/J)
const std::string OUTFILE = "";


std::bitset<K> create_lattice(std::mt19937& generator){
    
    std::bitset<K> lattice;
    
    // fill in random numbers    
    std::uniform_int_distribution<int> distribution(0,1);
    
    for (size_t i = 0; i < K; ++i){
        
        lattice.set(i, distribution(generator)); // Sets the bit at idx. i to what the generator returns
    }

    return lattice;
};

void evolve_lattice(std::mt19937& generator, std::bitset<K> &white, std::bitset<K> &black){
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
            // Consequense of using test :: bitset.test(i) is 1 if bit i is 1 and zero otherwise.
            int left_black  = 2*black.test(I*M     + (leftJ  - 1 + (I%2))/2)-1;
            int right_black = 2*black.test(I*M     + (rightJ - 1 + (I%2))/2)-1;
            int up_black    = 2*black.test(upI*M   + (J - 1 + (upI  %2))/2)-1;
            int down_black  = 2*black.test(downI*M + (J - 1 + (downI%2))/2)-1;
            int spin        = 2*white.test(i*M     + j)-1;
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
            int left_white  = 2*white.test(I*M     + (leftJ  - (I%2))/2)-1;
            int right_white = 2*white.test(I*M     + (rightJ - (I%2))/2)-1;
            int up_white    = 2*white.test(upI*M   + (J - (upI  %2))/2)-1;
            int down_white  = 2*white.test(downI*M + (J - (downI%2))/2)-1;
            int spin        = 2*black.test(i*M     + j)-1;
            
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

void join_arrays(const std::bitset<K> &white, const std::bitset<K> &black, std::bitset<2*K> &global){

    for(size_t i=0; i < N; i++){
        for(size_t j=0; j < M; j++){
            global.set(i*N + (2*j + (i%2))), white.test(i*M + j);
            global.set(i*N + (2*j + 1 - (i%2)), black.test(i*M + j));
        }
    }
}

void write_to_file(const std::bitset<2*K> &global){
    std::ofstream file(OUTFILE, std::ios::app);
    if (!file) {
        std::cerr << "Error opening file for writing!" << std::endl;
        return;
    }
    for (size_t k = 0; k < 2*K; ++k){
        file << global.test(k) << " ";
    }
    file << "\n";
    file.close();
}

void calculate_quantities(const std::bitset<2*K> &global, double& magnetization, double& energy){
    
    magnetization = 0.0;
    energy = 0.0;
    
    for (size_t i = 0; i < N; ++i){
        for (size_t j = 0; j < N; ++j){
            magnetization += global.test(i*N + j);
            energy += -global.test(i*N + j) * (global.test(((i+1)%N)*N + j) + global.test(((i-1+N)%N)*N + j)
                                        + global.test(i*N + (j+1)%N) + global.test(i*N + (j-1+N)%N) + H);
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
    // create white lattice
    std::bitset<K> white = create_lattice(generator);; 
    // create black lattice
    std::bitset<K> black = create_lattice(generator);
    // create global lattice
    std::bitset<2*K> global;
    
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