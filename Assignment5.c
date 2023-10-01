#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define randomreal0() (rand()/(RAND_MAX+1.0))

int L = 10;
int LS = 100;
double deltaHsPlus[9];
double deltaHsMinus[9];

//macros for wrapping
#define left(j) (j == 0) ? (L - 1) : (j - 1)
#define right(j) (j == L) ? (0) : (j + 1)
#define up(j) (j < L) ? (LS-L+j): (j - L)
#define down(j) (j >= LS-L) ? (j + L - LS): (j + L)

int main(int argc, char** argv){
    FILE *output = fopen("data.csv", "w");
    fprintf(output, "Sweep, M, |M|\n");
    int latticeSize = L*L;
    signed char* lattice = malloc(sizeof(char)*latticeSize);
    memset(lattice, -1, sizeof(char)*latticeSize);
    double J = strtod(argv[1], NULL);
    double h = strtod(argv[2], NULL);
    for(int i = -4; i <= 4; i ++){
        double beta = 1;
        deltaHsPlus[i+4] = beta*exp(-2*(i*J+2*h));
        printf("%d, %f, %f\n", i, -2.0*(i*J+2*h), beta*exp(-2*(i*J+2*h)));
        deltaHsMinus[i+4] = beta*exp(-2*(i*J-2*h));
    }
    for(int i = 0; i < 9; i ++){
        printf("%d, %f\n", i, deltaHsMinus[i]);
    }
    int passes = 1000;
    signed char si;
    int sigmaSj;
    double p;
    int M;
    for(int i = 0; i < passes; i ++){
        M = 0;
        printf("Pass %d\n", i);
        for(int j = 0; j < latticeSize; j ++){
            si = lattice[j];
            sigmaSj = 0;
            //check neighbours
            //todo wrapping
            sigmaSj += lattice[left(j)];
            sigmaSj += lattice[right(j)];
            sigmaSj += lattice[up(j)];
            sigmaSj += lattice[down(j)];
            if(si == 1){
                printf("%d\n", si*sigmaSj+4);
                p = deltaHsPlus[si*sigmaSj+4];
                printf("%f\n", p);
                if(randomreal0() < p || p > 1){
                    lattice[j] = -1;
                }
            }else{
                //read backwards?
                //printf("%d\n", si);
                printf("%d\n", si*sigmaSj+4);
                p = deltaHsMinus[si*sigmaSj+4];
                printf("%f\n", p);
                if(randomreal0() < p*RAND_MAX || p > 1){
                    lattice[j] = 1;
                }
            }
        }
        for(int j = 0; j < latticeSize; j ++){
            M += lattice[j];
        }
        fprintf(output, "%d, %d, %d\n", i, M, abs(M));
    }
    //record values
}