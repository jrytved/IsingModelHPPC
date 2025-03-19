// include basic libraries for I/O, math and array manipulation

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <random>
#include <chrono>

class SystemConfig {
    public:
        size_t steps = 500;     // number of steps
        size_t N = 128;      // dimensions
        size_t M = N/2;
        double T = 0.1;     // temperature
        double h = 0.0;     // external magnetic field
        size_t fwrite = 4; // how often to save data
        std::string filename = "";   // name of the output file with trajectory
        int seed = 1234;    // seed for random number generator
        
        SystemConfig(std::vector <std::string> argument){
            for (size_t i = 1; i<argument.size() ; i += 2){
                std::string arg = argument.at(i);
                if(arg=="-h"){ // Write help
                    std::cout << "Ising model --steps <number of steps> --size <square lattice size>"
                              << " --temp <temperature> --extfield <h> --fwrite <io frequency> --out <filename> \n";
                    exit(0);
                    break;
                } else if(arg=="--steps"){
                    steps = std::stoi(argument[i+1]);
                } else if(arg=="--size"){
                    N = std::stoi(argument[i+1]);
                    M = std::stoi(argument[i+1]);
                } else if(arg=="--temp"){
                    T = std::stod(argument[i+1]);
                } else if(arg=="--extfield"){
                    h = std::stod(argument[i+1]);
                } else if(arg=="--fwrite"){
                    fwrite = std::stoi(argument[i+1]);
                } else if(arg=="--out"){
                    filename = argument[i+1];
                } else{
                    std::cout << "---> error: the argument type is not recognized \n";
                }
            }    
        }
    };


int* create_lattice(SystemConfig& config, std::mt19937& generator){
    size_t size = config.N*config.M;
    int* lattice = new int[size];
    // fill in random numbers    
    std::uniform_int_distribution<int> distribution(0,1);
    for (size_t i = 0; i < size; ++i){
        lattice[i] = distribution(generator) * 2 - 1;
    }
    return lattice;
}

void evolve_lattice(SystemConfig& config, std::mt19937& generator, int* white, int* black){
    // update all white
    #pragma omp parallel for
    for (size_t i = 0; i < config.N; ++i){
        for (size_t j = 0; j < config.M; ++j){
            // i,j are now in local white space
            // transform to global I,J indices
            size_t I = i;
            size_t J = 2*j + (i%2);
            // get neighbors in a single step without conditionals
            size_t leftJ  = (J-1+config.N)%config.N;
            size_t rightJ = (J+1)%config.N;
            size_t upI    = (I-1+config.N)%config.N;
            size_t downI  = (I+1)%config.N;
            // 
            int left_black  = black[I*config.M     + (leftJ  - 1 + (I%2))/2];
            int right_black = black[I*config.M     + (rightJ - 1 + (I%2))/2];
            int up_black    = black[upI*config.M   + (J - 1 + (upI  %2))/2];
            int down_black  = black[downI*config.M + (J - 1 + (downI%2))/2];
            int spin        = white[i*config.M     + j];
            // energy difference
            double dE = 2.0 * spin * (left_black + right_black + up_black + down_black + config.h);
            // decide based on a generator with minimal conditionals
            double boltzmann = exp(-dE/config.T);
            std::uniform_real_distribution<double> distribution(0.0, 1.0);
            if (distribution(generator) < boltzmann){
                white[i*config.M + j] = -spin;
            }
        } 
    }
    // update all black

    #pragma omp parallel for
    for (size_t i = 0; i < config.N; ++i){
        for (size_t j = 0; j < config.M; ++j){
            // i,j are now in local black space
            // transform to global I,J indices
            size_t I = i;
            size_t J = 2*j + 1 - (i%2);
            // get neighbors in a single step without conditionals
            size_t leftJ  = (J-1+config.N)%config.N;
            size_t rightJ = (J+1)%config.N;
            size_t upI    = (I-1+config.N)%config.N;
            size_t downI  = (I+1)%config.N;
            //
            int left_white  = white[I*config.M     + (leftJ  - (I%2))/2];
            int right_white = white[I*config.M     + (rightJ - (I%2))/2];
            int up_white    = white[upI*config.M   + (J - (upI  %2))/2];
            int down_white  = white[downI*config.M + (J - (downI%2))/2];
            int spin        = black[i*config.M     + j];
            // energy difference
            double dE = 2.0 * spin * (left_white + right_white + up_white + down_white + config.h);
            // decide based on a generator with minimal conditionals
            double boltzmann = exp(-dE/config.T);
            std::uniform_real_distribution<double> distribution(0.0, 1.0);
            if (distribution(generator) < boltzmann){
                black[i*config.M + j] = -spin;
            }
        } 
    }
}

void join_arrays(SystemConfig& config, const int* white, const int* black, int* global){
    
    // implemented without conditionals

    #pragma omp parallel for
    for(size_t i=0; i<config.N; i++){
        for(size_t j=0; j<config.M; j++){
            global[i*config.N + (2*j + (i%2))] = white[i*config.M + j];
            global[i*config.N + (2*j + 1 - (i%2))] = black[i*config.M + j];
        }
    }
}

void write_to_file(SystemConfig& config, const int* global){
    std::ofstream file(config.filename, std::ios::app);
    if (!file) {
        std::cerr << "Error opening file for writing!" << std::endl;
        return;
    }
    for (size_t i = 0; i < config.N; ++i){
        for (size_t j = 0; j < config.N; ++j){
            file << global[i*config.N + j] << " ";
        }
    }
    file << "\n";
    file.close();
}

void calculate_quantities(SystemConfig& config, const int* global, double& magnetization, double& energy){
    magnetization = 0.0;
    energy = 0.0;

    #pragma omp parallel for reduction(+:magnetization, energy)
    for (size_t i = 0; i < config.N; ++i){
        for (size_t j = 0; j < config.N; ++j){
            magnetization += global[i*config.N + j];
            energy += -global[i*config.N + j] * (global[((i+1)%config.N)*config.N + j] + global[((i-1+config.N)%config.N)*config.N + j]
                                                + global[i*config.N + (j+1)%config.N] + global[i*config.N + (j-1+config.N)%config.N] + config.h);
        }
    }
    magnetization /= config.N*config.N;
    magnetization = fabs(magnetization);
    energy /= config.N*config.N;
}

void save_quantities(SystemConfig& config, const double* magnetization, const double* energy){
    // save magnetization and energy
    std::string filename_energy = "energy_" + config.filename;
    std::string filename_magnetization = "magnetization_" + config.filename;
    // save the magnetization
    std::ofstream file_magnetization(filename_magnetization);
    for (size_t s = 0; s < config.steps; ++s){
        file_magnetization << s << " " << magnetization[s] << "\n";
    }
    file_magnetization.close();
    // save the energy
    std::ofstream file_energy(filename_energy);
    for (size_t s = 0; s < config.steps; ++s){
        file_energy << s << " " << energy[s] << "\n";
    }
    file_energy.close();
}

int main(int argc, char* argv[]){
    SystemConfig config ({argv, argv+argc});
    
    // set up random generator
    std::mt19937 generator(config.seed);
    // create white lattice
    int* white = create_lattice(config, generator); 
    // create black lattice
    int* black = create_lattice(config, generator);
    // create global lattice
    int* global = new int[config.N*config.M*2]; // M is defined as N/2 in SystemConfig
    
    // create containers for magnetization and energy
    double* magnetization = new double[config.steps];
    double* energy = new double[config.steps];

    // delete files if they exist already
    std::remove(config.filename.c_str());
    std::remove(("energy_" + config.filename).c_str());
    std::remove(("magnetization_" + config.filename).c_str());

    auto begin = std::chrono::steady_clock::now();
    
    // loop over steps
    for (size_t s = 0; s < config.steps; ++s){
        // evolve lattice
        evolve_lattice(config, generator, white, black);
        // collect in global
        join_arrays(config, white, black, global);
        // calculate energy and magnetization
        calculate_quantities(config, global, magnetization[s], energy[s]);        
        // save to file
        if (s % config.fwrite == 0 && !config.filename.empty()){
            write_to_file(config, global);
        }
    }

    auto end = std::chrono::steady_clock::now();
    
    // save magnetization and energy
    save_quantities(config, magnetization, energy);

    std::cout << "Elapsed Time\t" << (end - begin).count() / 1000000000.0 << " sec" << std::endl;
    
    // clean up
    delete[] white;
    delete[] black;
    delete[] global;
    delete[] magnetization;
    delete[] energy;

    return 0;

}