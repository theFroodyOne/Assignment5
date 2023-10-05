#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define randomreal0() (rand()/(RAND_MAX+1.0))

int L = 200;
int LS = 40000;
double deltaHsPlus[9];
double deltaHsMinus[9];

//macros for wrapping
#define left(j) (j == 0) ? (L - 1) : (j - 1)
#define right(j) (j == L) ? (0) : (j + 1)
#define up(j) (j < L) ? (LS-L+j): (j - L)
#define down(j) (j >= LS-L) ? (j + L - LS): (j + L)

int main(int argc, char** argv){
    double J = strtod(argv[1], NULL);
    double h = strtod(argv[2], NULL);
    char* filename = malloc(sizeof(char)*101);
    sprintf(filename,"part1data/J %.3f h %.3f L %d.csv", J, h, L);
    FILE *output = fopen(filename, "w");
    fprintf(output, "Sweep, M/L^2, |M|/L^2\n");
    int latticeSize = L*L;
    signed char* lattice = malloc(sizeof(char)*latticeSize);
    memset(lattice, -1, sizeof(char)*latticeSize);
    for(int i = -4; i <= 4; i ++){
        double beta = 1;
        deltaHsPlus[i+4] = beta*exp(-2*(i*J+2*h));
        deltaHsMinus[i+4] = beta*exp(-2*(i*J-2*h));
    }
    int passes = 1000;
    signed char si;
    int sigmaSj;
    double p;
    int M;
    double sumM = 0;
    double sumAbsM = 0;
    double sumM2 = 0;
    for(int i = 0; i < passes; i ++){
        M = 0;
        //printf("Pass %d\n", i);
        for(int j = 0; j < latticeSize; j ++){
            si = lattice[j];
            sigmaSj = 0;
            //check neighbours
            sigmaSj += lattice[left(j)];
            sigmaSj += lattice[right(j)];
            sigmaSj += lattice[up(j)];
            sigmaSj += lattice[down(j)];
            if(si == 1){
                p = deltaHsPlus[si*sigmaSj+4];
                if(randomreal0() < p || p > 1){
                    lattice[j] = -1;
                }
            }else{
                p = deltaHsMinus[si*sigmaSj+4];
                if(randomreal0() < p*RAND_MAX || p > 1){
                    lattice[j] = 1;
                }
            }
        }
        for(int j = 0; j < latticeSize; j ++){
            M += lattice[j];
        }
        fprintf(output, "%d, %f, %f\n", i, (double)M/latticeSize, (double)abs(M)/latticeSize);

        //todo energy stuff
    }
    sumM /= L*L;
    sumAbsM /= L*L;
    sumM2 /= L*L;
    printf("<m>: %f, <|m|>: %f, <m^2>: %f\n", sumM/passes, sumAbsM/passes, sumM2/passes);
}