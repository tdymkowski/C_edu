// lj_potential.h

#ifndef LJ_POTENTIAL_H
#define LJ_POTENTIAL_H

#include <stdio.h>
#include <math.h>
#include "../main.h"
#include "models.h"


typedef struct{
  float R_MAX;

} LJParams;

void make_lj_force(ForceModel *m, float R_MAX);
void make_lj_model(EnergyModel *m, float R_MAX);
#endif
