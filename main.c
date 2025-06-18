#include "lj_potential.h"
#include "xyz_reader.h"
#include <stdio.h>
#include <stdlib.h>
#include "main.h"


int main(int argc, char *argv[]){
  if (argc == 1){
    printf("Program requires one positional argument:\n");
    printf("%s <file.xyz>\n", argv[0]);
    return EXIT_FAILURE;
  }
  int NMAX = get_max_atoms(argv[1]);
  struct Atoms atoms[NMAX];
  main_xyz_reader(argv[1], atoms, NMAX);
  // LJ Epsilon
  double epsilon = 3.1401000;  // eV
  // LJ sigma
  double sigma =  2.4232400;  // Angstroms
  // Cell number
  const int NCLMAX = 10;
  // Cutoff value
  const float R_CUT = 9;
  double Region[3] = {4., 4, 4};
  main_lj(epsilon,
          sigma,
          NMAX,
          NCLMAX,
          R_CUT,
          atoms,
          Region);


  return EXIT_SUCCESS;
}
