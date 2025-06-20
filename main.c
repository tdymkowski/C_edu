#include <stdio.h>
#include <stdlib.h>
#include "readers/xyz_reader.h"
#include "main.h"
#include "verlet.h"


int main(int argc, char *argv[]){
  if (argc == 1){
    printf("Program requires two positional arguments:\n");
    printf("%s <file.xyz> <model_name> <timesteps> [default: 1]>\n", argv[0]);
    return EXIT_FAILURE;
  }
  const char *model_name = argv[2];
  float R_CUT = 4.;
  // Lx, Ly, Lz
  double Region[3] = {20., 20., 21.};

  int NMAX = get_max_atoms(argv[1]);
  struct Atoms atoms[NMAX];
  main_xyz_reader(argv[1], atoms, NMAX);
  double dt = 1;
  start_simulation(model_name,
                   atoms,
                   NMAX,
                   R_CUT,
                   Region,
                   dt,
                   100);

  return EXIT_SUCCESS;
}
