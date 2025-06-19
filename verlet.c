#include "lj_potential.h"
#include "main.h"
#include "xyz_writer.h"
#include "verlet.h"


void update_velocities(struct Atoms atom, double dt){
  float v_x, v_y, v_z;
  atom.vel_x += dt * atom.a_x;
  atom.vel_y += dt * atom.a_y;
  atom.vel_z += dt * atom.a_z;
}


void update_positions(struct Atoms atom, double dt){
  atom.x += dt * atom.vel_x;
  atom.x += dt * atom.vel_y;
  atom.z += dt * atom.vel_z;
}


void half_kick(struct Atoms atoms[], double dt, int NMAX){
  for (int i = 0; i < NMAX; i++){
    update_velocities(atoms[i], dt / 2);
  }
}


void single_step(struct Atoms atoms[], double dt, int NMAX){
  int i;
  for (i = 0; i < NMAX; i++) half_kick(atoms, dt, NMAX);
  for (i = 0; i < NMAX; i++){
      atoms[i].x = atoms[i].x + dt * atoms[i].vel_x;
      atoms[i].y = atoms[i].y + dt * atoms[i].vel_y;
      atoms[i].z = atoms[i].z + dt * atoms[i].vel_z;
      printf("%f %f %f\n", atoms[i].x, atoms[i].y, atoms[i].z);
  }
  for (i = 0; i < NMAX; i++) half_kick(atoms, dt, NMAX);
}

void propagate_verlet(ModelFn model,
    struct Atoms atoms[],
    double dt,
    double MAX_T,
    int NMAX,
    int NCLMAX,
    double R_CUT,
    double Region[]){
  for (double t = 0; t < MAX_T; t += dt){
    model(NMAX, NCLMAX, R_CUT, atoms, Region);
    single_step(atoms, dt, NMAX);
    write2xyz(atoms, "out.xyz", NMAX);
  }
}

void start_simulation(struct Atoms atoms[],
    int NMAX,
    int NCLMAX,
    double R_CUT,
    double Region[3],
    double dt,
    double MAX_T){
    write2xyz(atoms, "out.xyz", NMAX);
    propagate_verlet(main_lj,
      atoms,
      dt,
      MAX_T,
      NMAX,
      NCLMAX,
      R_CUT,
      Region);
}


