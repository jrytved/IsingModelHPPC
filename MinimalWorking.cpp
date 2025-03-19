#include <iostream>
#include <random>
#include <math.h>

const int N = 10;                  // NUMBER OF ROWS
const int M = 10;                  // NUMBER OF COLUMNS
const double TEMP = 0.5;           // TEMPERATURE
const double INV_TEMP = 1 / TEMP;  // INVERSE OF THE TEMPERATURE
const int STEPS = 100;             // NUMBER OF SIMULATION STEPS
const int J = 1;                   // INTERACTION CONSTANT
const int H = 0;                   // CONSTANT FOR EXTERNAL MAGNETIC FIELD

double random_double() {
    
    static std::random_device rd;                                 // Seed generator
    static std::mt19937 gen(rd());                                // Mersenne Twister PRNG
    static std::uniform_real_distribution<double> dist(0.0, 1.0); // Uniform distribution
    return dist(gen);
}


void init_array(int (&arr)[N][M]) {

    // Initializes the array randomly with (i,j) IN {-1,1}
    
    std::random_device rd;  // Random seed
    std::mt19937 gen(rd()); // Mersenne Twister PRNG
    std::uniform_int_distribution<int> dist(0, 1); // Generates 0 or 1

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            arr[i][j] = dist(gen) * 2 - 1; // Convert {0,1} to {-1,1}
        }
    }
};

void print_arr(int (&arr)[N][M]){

    // Prints the array in nice readable format.

    for(int i=0; i<N; i++){
        if(i>0){printf("\n");};
        for(int j=0; j<M; j++){

            char symbol = (arr[i][j] > 0)? '#' : ' ';
            printf("|%c|", symbol);
        }
    }
    printf("\n\n");
};

void step(int (&arr)[N][M]){

    int ipp, jpp, inn, jnn, nb, spin, dE;

    for(int i=0;i<N;i++)
    for(int j=0;j<M;j++){

        ipp = ((i+1) < N)? (i+1) : 0;
        jpp = ((j+1) < M)? (j+1) : 0;
        
        inn = ((i-1) >= 0)? (i-1) : (N - 1);
        jnn = ((j-1) >= 0)? (j-1) : (M - 1);

        // CALCULATE NEIGHBORS
        nb = arr[ipp][j] + arr[i][jpp] + arr[inn][j] + arr[i][jnn];

        // COMPUTE ENERGY DIFFERENCE
        spin = arr[i][j];
        dE = -2*spin*(J * nb + H);

        if(dE > 0){
            arr[i][j] = -spin;
        } else {
            
            if(random_double() < exp(dE*INV_TEMP)){
                arr[i][j] = -spin;
            }
            
        }
        
    }
    
};


int main(void){
    
    int arr[M][N];
    init_array(arr);
    printf("Initial Conditions:\n");
    print_arr(arr);
    step(arr);
    print_arr(arr);
    step(arr);
    print_arr(arr);

    
}



