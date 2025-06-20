// xyz_reader.h
#ifndef VERLET_H
#define VERLET_H

#include <stdio.h>
#include <math.h>
#include "main.h"
#include "models/models.h"
#include "models/choose_models.h"

typedef int (*ModelFn)(int NMAX,
                       double R_CUT,
                       struct Atoms atoms[],
                       double Region[3]);

void single_step(struct Atoms atoms[], double dt, int NMAX);
void start_simulation(const char *model_name,
                      struct Atoms atoms[],
                      int NMAX,
                      double R_CUT,
                      double Region[3],
                      double dt,
                      double MAX_T);
void clear_atoms(struct Atoms atoms[], int NMAX);
void acceleration(ForceModel *model,
                  struct Atoms *atoms_i,
                  struct Atoms *atoms_j);


#endif
