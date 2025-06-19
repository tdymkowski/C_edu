#include "lj_potential.h"
#include "main.h"
#include "xyz_writer.h"
#include <stdio.h>
#include "verlet.h"


void clear_atom(struct Atoms *atom){
  atom->a_x = atom->a_y = atom->a_z = 0;
  atom->vel_x = atom->vel_y = atom->vel_z = 0;
}

void clear_atoms(struct Atoms atoms[], int NMAX) {
    for (int i = 0; i < NMAX; i++) {
        atoms[i].vel_x = atoms[i].vel_y = atoms[i].vel_z = 0.0;
        atoms[i].a_x   = atoms[i].a_y   = atoms[i].a_z   = 0.0;
    }
}


void update_velocities(struct Atoms *atom, double dt){
  float v_x, v_y, v_z;
  atom->vel_x += dt * atom->a_x;
  atom->vel_y += dt * atom->a_y;
  atom->vel_z += dt * atom->a_z;
}


void update_positions(struct Atoms *atom, double dt){
  atom->x += dt * atom->vel_x;
  atom->y += dt * atom->vel_y;
  atom->z += dt * atom->vel_z;
}


void half_kick(struct Atoms atoms[], double dt, int NMAX){
  for (int i = 0; i < NMAX; i++){
    update_velocities(&atoms[i], dt / 2);
  }
}


void single_step(struct Atoms atoms[], double dt, int NMAX){
  int i;
  half_kick(atoms, dt, NMAX);
  for (i = 0; i < NMAX; i++){
      update_positions(&atoms[i], dt);
      printf("%f %f %f\n", atoms[i].x, atoms[i].y, atoms[i].z);
  }
  half_kick(atoms, dt, NMAX);
}

void propagate_verlet(ModelFn model,
                      struct Atoms atoms[],
                      int NMAX,
                      int NCLMAX,
                      double R_CUT,
                      double Region[3],
                      double dt,
                      double MAX_T){
    for (double t = 0; t < MAX_T; t += dt) {
        // compute & store forces into atoms[].a_x/y/z
        model(NMAX, NCLMAX, R_CUT, atoms, Region);

        // advance via velocity‐Verlet
        single_step(atoms, dt, NMAX);

        write2xyz(atoms, "out.xyz", NMAX);
        // Clear forces
        //for (int z=0; z<NMAX; z++) clear_atom(&atoms[z]);
        clear_atoms(atoms, NMAX);
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
    clear_atoms(atoms, NMAX);
    propagate_verlet(main_lj,
                     atoms,
                     NMAX,
                     NCLMAX,
                     R_CUT,
                     Region,
                     dt,
                     MAX_T);
}


