// lj_potential.h

#ifndef LJ_POTENTIAL_H
#define LJ_POTENTIAL_H

#include <stdio.h>
#include <math.h>
#include "main.h"

double lj_potential(float epsilon, float sigma, float r, float R_MAX);
int main_lj(int NMAX, int NCLMAX, double R_CUT, struct Atoms atoms[], double Region[3]);

#endif
