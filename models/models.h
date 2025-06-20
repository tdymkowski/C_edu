// physical models
#ifndef MODELS_H
#define MODELS_H

#include <stdio.h>
#include <math.h>
#include "../main.h"


// structs for creating new models
typedef double (*Energy)(void *ctx,
                         struct Atoms *i,
                         struct Atoms *j);

typedef struct{
  void *ctx;
  double R_MAX;
  Energy raw_energy;
} EnergyModel;


typedef double (*RawForceFn)(void *ctx,
                             struct Atoms *i,
                             struct Atoms *j);

typedef struct{
  void *ctx;
  double R_MAX;
  RawForceFn raw_force;
} ForceModel;


typedef struct{
  EnergyModel energy;
  ForceModel force;
} PairModels;

// MODELS

// LJ model
void lj_force(ForceModel *m,
              float R_MAX,
              struct Atoms i,
              struct Atoms j);

void lj_model(EnergyModel *m,
              float R_MAX,
              struct Atoms i,
              struct Atoms j);

// Coulomb model

#endif

