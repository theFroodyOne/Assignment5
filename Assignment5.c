#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define randomreal0() (rand()/(RAND_MAX+1.0))

int L = 200;
double ps[9];

//macros for wrapping
#define LS L*L
#define left(j) (j == 0) ? (L - 1) : (j - 1)
#define right(j) (j == L) ? (0) : (j + 1)
#define up(j) (j < L) ? (LS-L+j): (j - L)
#define down(j) (j >= LS-L) ? (j + L - LS): (j + L)

int main(int argc, char** argv){
    double h = 0;
    double Jc = 0.4407;
    FILE *output = fopen("part4data/data.csv", "w");
    fprintf(output, "t, m, G(n),\n");
    for(double t = -0.4; t < 1.0; t += 0.2) {
        printf("t: %f\n", t);
        double J = Jc/(t + 1);
        int latticeSize = L * L;
        signed char *lattice = malloc(sizeof(char) * latticeSize);
        memset(lattice, 1, sizeof(char) * latticeSize);
        for (int i = -4; i <= 4; i++) {
            ps[i + 4] = exp(-2 * (i * J + 2 * h));
        }
        int passes = 400000;
        int equilibration = passes/10;
        signed char si;
        int sigmaSj;
        double p;
        int M;
        double sumM = 0;
        double *Grs = calloc(101, sizeof(double));
        for (int i = 0; i < passes; i++) {
            M = 0;
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
                    }
                } else {
                    p = ps[si * sigmaSj + 4];
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
                for(int k = 0; k < L - 1; k ++){
                    for(int r = 1; r <= 100; r ++){
                        Grs[r] += lattice[k*L]*lattice[k*L + r];
                    }
                }
            }
        }
        sumM /= (passes - equilibration);
        fprintf(output, "%f, %f", t, sumM/latticeSize);
        for(int i = 1; i <= 100; i ++){
            Grs[i] /= (L - 1)*(passes - equilibration);
            fprintf(output, ", %f", Grs[i]);
        }
        fprintf(output, "\n");
    }
}