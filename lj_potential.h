// lj_potential.h

#ifndef LJ_POTENTIAL_H
#define LJ_POTENTIAL_H

#include <stdio.h>
#include <math.h>
#include "main.h"

float lj_potential(float epsilon, float sigma, float r);
int main_lj(double epsilon,
            double sigma,
            int NMAX,
            int NCLMAX,
            float R_CUT,
            struct Atoms atoms[],
            double Region[3]);
#endif
