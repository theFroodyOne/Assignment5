#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define randomreal0() (rand()/(RAND_MAX+1.0))

int L = 200;
int LS = 40000;
double ps[9];
double dH[9];

//macros for wrapping
#define left(j) (j == 0) ? (L - 1) : (j - 1)
#define right(j) (j == L) ? (0) : (j + 1)
#define up(j) (j < L) ? (LS-L+j): (j - L)
#define down(j) (j >= LS-L) ? (j + L - LS): (j + L)

int main(int argc, char** argv){
    double h = 0;
    double Jc = 0.4407;
    double Tc = 2.269;
    FILE *output = fopen("part3data/data.csv", "w");
    fprintf(output, "t, c/k, chi\n");
    for(double t = -0.1; t < 0.1;) {
        if(t > -0.01 && t < 0.01){
            t += 0.0001;
        } else {
            if (t > -0.1 && t < 0.1) {
                t += 0.001;
            } else {
                t += 0.01;
            }
        }
        printf("t: %f\n", t);
        double J = Jc/(t + 1);
        int latticeSize = L * L;
        signed char *lattice = malloc(sizeof(char) * latticeSize);
        memset(lattice, 1, sizeof(char) * latticeSize);
        for (int i = -4; i <= 4; i++) {
            ps[i + 4] = exp(-2 * (i * J + 2 * h));
            dH[i + 4] = -2 * (i * J + 2 * h);
        }
        int passes = 100000;
        int equilibration = passes/10;
        signed char si;
        int sigmaSj;
        double p;
        int M;
        double E;
        double sumM = 0;
        double sumM2 = 0;
        double sumE = 0; //set initial energy when all spins are pointing up as 0
        double sumE2 = 0;
        for (int i = 0; i < passes; i++) {
            M = 0;
            E = 0;
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
                    p = ps[si * sigmaSj + 4];
                    if (randomreal0() < p || p > 1) {
                        lattice[j] = -1;
                        E += dH[si * sigmaSj + 4];
                    }
                } else {
                    p = ps[si * sigmaSj + 4];
                    if (randomreal0() < p || p > 1) {
                        lattice[j] = 1;
                        E += dH[si * sigmaSj + 4];
                    }
                }
            }
            for (int j = 0; j < latticeSize; j++) {
                M += lattice[j];
            }
            if (i >= equilibration) {
                sumM += M;
                sumM2 += M*M;
                sumE += E;
                sumE2 += E*E;
            }
        }
        //todo replace with proper variables for part 3
        sumM /= (passes - equilibration);
        sumM2 /= (passes - equilibration);
        double beta = 1/(t+1)*Tc;
        sumE /= (passes - equilibration);
        sumE2 /= (passes - equilibration);
        double c = (sumE2 - sumE*sumE)/latticeSize;
        c *= beta*beta;
        double chi = (sumM2 - sumM*sumM)/latticeSize;
        chi *= beta;

        fprintf(output, "%f, %f, %f\n", t, c, chi);
    }
}