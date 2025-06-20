#include <math.h>
#include <stdio.h>
#include "verlet.h"
#include "main.h"
#include "writers/xyz_writer.h"
#include "models/models.h"
#include "models/choose_models.h"
#include "neighbors/pair_interactions.h"


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


void acceleration(ForceModel *model,
                  struct Atoms *atoms_i,
                  struct Atoms *atoms_j){
  double fx, fy, fz;
  float r_ijx = atoms_i->x - atoms_j->x;
  float r_ijy = atoms_i->y - atoms_j->y;
  float r_ijz = atoms_i->z - atoms_j->z;
  float r_ij = sqrtf(r_ijx*r_ijx + r_ijy*r_ijy + r_ijz * r_ijz);
  if (r_ij == 0.0) return;

  double Fmag = model->raw_force(model->ctx,
                                 atoms_i,
                                 atoms_j);
  fx = Fmag * r_ijx / r_ij;
  fy = Fmag * r_ijy / r_ij;
  fz = Fmag * r_ijz / r_ij;
  
  
  atoms_i->a_x += fx / atoms_i->mass;
  atoms_i->a_y += fy / atoms_i->mass;
  atoms_i->a_z += fz / atoms_i->mass;
  
  atoms_j->a_x -= fx / atoms_j->mass;
  atoms_j->a_y -= fy / atoms_j->mass;
  atoms_j->a_z -= fz / atoms_j->mass;
  }



void update_velocities(struct Atoms *atom, double dt){
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
//      printf("%f %f %f\n", atoms[i].x, atoms[i].y, atoms[i].z);
  }
  half_kick(atoms, dt, NMAX);
}

void propagate_verlet(PairModels M,
                      struct Atoms atoms[],
                      int NMAX,
                      double R_CUT,
                      double Region[3],
                      double dt,
                      double MAX_T){
    int nsteps = (int)ceil(MAX_T / dt);
    fprintf(stderr,
        "DEBUG propagate: dt = %g, MAX_T = %g, ceil(MAX_T/dt) = %d\n",
        dt, MAX_T, nsteps);
    for (int step = 0; step < nsteps; step++) {
        // Compute & store forces into atoms[].a_x/y/z
        fprintf(stdout, "Step: %i\n", step);

        compute_pair_interactions(NMAX,
                                  R_CUT,
                                  atoms,
                                  Region,
                                  &M.energy,
                                  &M.force);
        // Get single step
        single_step(atoms, dt, NMAX);
        // Write to file
        write2xyz(atoms, "out.xyz", NMAX);
        // Clear forces
        clear_atoms(atoms, NMAX);
    }
}


void start_simulation(const char *model_name,
                      struct Atoms atoms[],
                      int NMAX,
                      double R_CUT,
                      double Region[3],
                      double dt,
                      double MAX_T){
//    printf("DEBUG 1 verlet.c: choosing model\n");
    PairModels M = choose_model(model_name, atoms, NMAX, R_CUT);
//    printf("DEBUG 2 verlet.c: writing 1st frame\n");
    write2xyz(atoms, "out.xyz", NMAX);
//    printf("DEBUG 3 verlet.c: clearing forces\n");
    clear_atoms(atoms, NMAX);
//    printf("DEBUG 4 verlet.c: propagating verlet\n");
    propagate_verlet(M,
                     atoms,
                     NMAX,
                     R_CUT,
                     Region,
                     dt,
                     MAX_T);
    printf("Done running %lf steps\n", MAX_T);
}
