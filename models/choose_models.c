#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "models.h"
#include "lj_potential.h"
#include "../main.h"
#include "../readers/read_parameters.h"


PairModels choose_model(const char *model_name,
                        struct Atoms atoms[],
                        int N_MAX,
                        float R_MAX){
  printf("DEBUG 1 choose_models.c: choosing model lj\n");
  PairModels M;
  if (strcmp(model_name, "lj") == 0){
    printf("DEBUG 2 choose_models.c: reading params\n");
    main_params_reader(atoms, N_MAX);
    printf("DEBUG 3 choose_models.c: parsing energy model\n");
    make_lj_model(&M.energy, R_MAX);
    printf("DEBUG 4 choose_models.c: parsing froce model\n");
    make_lj_force(&M.force, R_MAX);
  }
  else {
    fprintf(stderr, "Unknown model '%s'\n", model_name);
    exit(EXIT_FAILURE);
  }

  return M;
}
