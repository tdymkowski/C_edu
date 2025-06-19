// xyz_reader.h
#ifndef VERLET_H
#define VERLET_H

#include <stdio.h>
#include <math.h>
#include "main.h"

typedef int (*ModelFn)(int,
    int,
    double,
    struct Atoms[],
    double []) ;

void update_velocities(struct Atoms atom, double dt);
void update_positions(struct Atoms atom, double dt);
void half_kick(struct Atoms atoms[], double dt, int NMAX);
void single_step(struct Atoms atoms[], double dt, int NMAX);
void start_simulation(struct Atoms atoms[],
    int NMAX,
    int NCLMAX,
    double R_CUT,
    double Region[3],
    double dt,
    double MAX_T);

#endif
