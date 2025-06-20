#include <math.h>
#include "lj_potential.h"
#include "../main.h"
#include "models.h"


static double get_R(struct Atoms *i, struct Atoms *j){
  double r_x = i->x - j->x;
  double r_y = i->y - j->y;
  double r_z = i->z - j->z;
  double r2 = r_x * r_x + r_y * r_y + r_z * r_z;
  double r = sqrt(r2);
  return r;

}

static inline double lj_raw(float epsilon, float sigma, float r){
  float s_r = sigma/r;
  float s6 = powf(s_r, 6.0f);
  float V = 4 * epsilon * (s6*s6 - s6);
  return V;
}


static inline double lj_raw_dr(float epislon,
                               float sigma,
                               float r){
  float s_r = sigma/r;
  float s6 = powf(s_r, 6.0f);
  float dVdr = 4 * epislon * (-12 * s6 *s6/r + 6 * s6/r);
  return dVdr;
}


double lj_potential(float epsilon,
                    float sigma,
                    float r,
                    float R_MAX){
  // Calculate potential and slope at the R = R_cutoff
  float Vc = lj_raw(epsilon, sigma, R_MAX);
  float dVdr = lj_raw_dr(epsilon, sigma, R_MAX);

  // Calculate parts of LJ potential
  float V = lj_raw(epsilon, sigma, r);

  if (V > Vc || r >= R_MAX){
    return Vc;
  }

  // Calculate LJ potential
  return V - Vc - (r - R_MAX)*dVdr;
}

static double lj_energy_wrapper(void *vctx,
                                struct Atoms *atom_i,
                                struct Atoms *atom_j){
  LJParams *p = (LJParams*)vctx;
  float r = get_R(atom_i, atom_j);
  return lj_potential(atom_i->epsilon, atom_i->sigma, r, p->R_MAX);
}


static double lj_force_wrapper(void *vctx,
                               struct Atoms *atom_i,
                               struct Atoms *atom_j){
  LJParams *p = (LJParams*)vctx;
  double r = get_R(atom_i, atom_j);
  double dVdrMAX = lj_raw_dr(atom_i->epsilon, atom_i->sigma, p->R_MAX);
  double dVdr = lj_raw_dr(atom_i->epsilon, atom_i->sigma, r);
  if (fabs(dVdr) > fabs(dVdrMAX)){
    return dVdrMAX;
  }
  else {
    return dVdr;
  }

}

void make_lj_force(ForceModel *m,
                   float R_MAX){
  static LJParams params;
  params.R_MAX = R_MAX;
  m->ctx = &params;
  m->R_MAX = R_MAX;
  m->raw_force = lj_force_wrapper;
}

void make_lj_model(EnergyModel *m,
                   float R_MAX){
  static LJParams params;
  params.R_MAX = R_MAX;
  m->ctx = &params;
  m->R_MAX = R_MAX;
  m->raw_energy = lj_energy_wrapper;
}
