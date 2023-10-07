#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define randomreal0() (rand()/(RAND_MAX+1.0))

int L = 200;
int LS = 40000;
double deltaHs[9];

//macros for wrapping
#define left(j) (j == 0) ? (L - 1) : (j - 1)
#define right(j) (j == L) ? (0) : (j + 1)
#define up(j) (j < L) ? (LS-L+j): (j - L)
#define down(j) (j >= LS-L) ? (j + L - LS): (j + L)

int main(int argc, char** argv){
    double h = 0;
    FILE *output = fopen("part2data/data.csv", "w");
    fprintf(output, "J, M/L^2, |M|/L^2\n");
    for(int k = 0; k < 100; k ++) {
        printf("%d\n", k);
        double J = k/100.0;
        int latticeSize = L * L;
        signed char *lattice = malloc(sizeof(char) * latticeSize);
        memset(lattice, 1, sizeof(char) * latticeSize);
        for (int i = -4; i <= 4; i++) {
            double beta = 1;
            deltaHs[i + 4] = beta * exp(-2 * (i * J + 2 * h));
        }
        int passes = 10000;
        int equilibration = passes/10;
        signed char si;
        int sigmaSj;
        double p;
        int M;
        double sumM = 0;
        double sumAbsM = 0;
        for (int i = 0; i < passes; i++) {
            M = 0;
            //printf("Pass %d\n", i);
            for (int j = 0; j < latticeSize; j++) {
                si = lattice[j];
                sigmaSj = 0;
                //check neighbours
                sigmaSj += lattice[left(j)];
                sigmaSj += lattice[right(j)];
                sigmaSj += lattice[up(j)];
                sigmaSj += lattice[down(j)];
                if (si == 1) {
                    p = deltaHs[si * sigmaSj + 4];
                    if (randomreal0() < p || p > 1) {
                        lattice[j] = -1;
                    }
                } else {
                    p = deltaHs[si * sigmaSj + 4];
                    if (randomreal0() < p || p > 1) {
                        lattice[j] = 1;
                    }
                }
            }
            for (int j = 0; j < latticeSize; j++) {
                M += lattice[j];
            }
            if (i >= equilibration) {
                sumM += M;
                sumAbsM += abs(M);
            }
            //todo energy stuff
        }
        sumM /= L * L;
        sumAbsM /= L * L;
        fprintf(output, "%f, %f, %f\n", J, sumM / (passes-equilibration), sumAbsM / (passes-equilibration));
    }
}