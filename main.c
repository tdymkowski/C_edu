#include "lj_potential.h"
#include "xyz_reader.h"
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "read_parameters.h"
#include "verlet.h"
#include "xyz_writer.h"

int main(int argc, char *argv[]){
  if (argc == 1){
    printf("Program requires one positional argument:\n");
    printf("%s <file.xyz>\n", argv[0]);
    return EXIT_FAILURE;
  }
  float R_CUT = 4.;
  int NMAX = get_max_atoms(argv[1]);
  struct Atoms atoms[NMAX];
  main_xyz_reader(argv[1], atoms, NMAX);
  main_params_reader(atoms, NMAX);
  const int NCLMAX = 10;
  // Lx, Ly, Lz
  double Region[3] = {8., 8., 8.};
  start_simulation(atoms,
    NMAX,
    NCLMAX,
    R_CUT,
    Region,
    1,
    100);

  return EXIT_SUCCESS;
}
