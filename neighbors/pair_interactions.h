// physical models
#ifndef PAIR_INTERACTIONS_H
#define PAIR_INTERACTIONS_H

#include <stdio.h>
#include <math.h>
#include "../models/models.h"
#include "../main.h"


int compute_pair_interactions(int NMAX,
                              double R_CUT,
                              struct Atoms atoms[],
                              double Region[3],
                              EnergyModel *energy_model,
                              ForceModel *force_model);

#endif

