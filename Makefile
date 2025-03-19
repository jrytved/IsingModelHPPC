CXX := g++

# Fastest executable (-ffast-math removes checking for NaNs and other things)
OPT=-O3 -ffast-math

# Add thread checking to code
#OPT=-O1 -fsanitize=thread,undefined,bounds -fPIE -pie

# Add memory / address error checking to code
#OPT=-O1 -fsanitize=leak,address,undefined,bounds -fPIE -pie

# Faster compilation time
#OPT=-O1

CXXFLAGS := $(OPT) -Wall -march=native -g -std=c++17

default: IsingSeq IsingOmp IsingBitSet IsingBitPacking IsingBitPackingOmp

IsingSeq: Ising_sequential_cmd.cpp
	$(CXX) Ising_sequential_cmd.cpp $(CXXFLAGS) -o IsingSeq

IsingOmp: Ising_omp_cmd.cpp
	$(CXX) Ising_omp_cmd.cpp $(CXXFLAGS) -fopenmp -o IsingOmp

IsingBitSet: Ising_sequential_bitset.cpp
	$(CXX) Ising_sequential_bitset.cpp $(CXXFLAGS) -o IsingBitSet

IsingBitPacking: Ising_sequential_bitpacking.cpp
	$(CXX) Ising_sequential_bitpacking.cpp $(CXXFLAGS) -o IsingBitPacking

IsingBitPackingOmp: Ising_sequential_bitpacking_omp.cpp
	$(CXX) Ising_sequential_bitpacking_omp.cpp $(CXXFLAGS) -fopenmp -o IsingBitPackingOmp

clean:
	rm -fr seq omp bitset bitpacking