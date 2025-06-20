#include <math.h>
#include <stdio.h>
#include <string.h>
#include "lj_potential.h"
#include "models.h"
#include "my_utils.h"
#include "../main.h"


static double kong_s6(struct Atoms *i, struct Atoms *j){
  if (strcmp(i->symbol, j->symbol) == 0){
    return powf(i->sigma, 6.0) * i->epsilon;
  }
  
  double sigma6 = powf(
                      i->sigma * powf(i->sigma, 6.0f) * j->epsilon * powf(j->sigma, 12.0f), 0.50);


  return sigma6;
}

static double kong_s12(struct Atoms *i, struct Atoms *j){
  if (strcmp(i->symbol, j->symbol) == 0){
    return powf(i->sigma, 12.0) * i->epsilon;
  }

  double sigma12 = (
      powf(
      (
      powf(i->epsilon * powf(i->sigma, 12), 1.f/13.f) + 
      powf(j->epsilon * powf(j->sigma, 12), 1.f/13.f)
      ) * 0.5 ,13)
      );
  return sigma12;
}

static inline double lj_raw_kong(float sigma6,
                                 float sigma12,
                                 float r){
  float r6 = powf(r, 6.);
  float s_r6 = sigma6/r6;
  float s_r12 = sigma12/(r6*r6);

  float V = 4 * (s_r12 - s_r6);
  return V;
}


static inline double lj_raw_dr_kong(float sigma6,
                                    float sigma12,
                                    float r){
  if (r <= 0.){
    return 0.0;
  }
  float r13 = powf(r, 13.);
  float r7 = powf(r, 7.);
  float s_r7 = sigma6/r7;
  float s_r13 = sigma12/(r13);
  float dVdr = 48 * ( s_r13 - 0.5 * s_r7);
  return dVdr;
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


double lj_potential_kong(float sigma6,
                         float sigma12,
                         float r,
                         float R_MAX){
  // Calculate potential and slope at the R = R_cutoff
  if (r <= 0.){
    return 0.0;
  }
  float Vc = lj_raw_kong(sigma6, sigma12, R_MAX);
  float dVdr = lj_raw_dr_kong(sigma6, sigma12, R_MAX);

  // Calculate parts of LJ potential
  float V = lj_raw_kong(sigma6, sigma12, r);
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
  double r = get_R(atom_i, atom_j);
  double s6_kong = kong_s6(atom_i, atom_j);
  double s12_kong = kong_s12(atom_i, atom_j);
  //return lj_potential(atom_i->epsilon, atom_i->sigma, r, p->R_MAX);
  return lj_potential_kong(s6_kong, s12_kong, r, p->R_MAX);
}


static double lj_force_wrapper(void *vctx,
                               struct Atoms *atom_i,
                               struct Atoms *atom_j){
  LJParams *p = (LJParams*)vctx;
  double r = get_R(atom_i, atom_j);
  double s6_kong = kong_s6(atom_i, atom_j);
  double s12_kong = kong_s12(atom_i, atom_j);
  double dVdrMAX = lj_raw_dr_kong(s6_kong, s12_kong, p->R_MAX);
  double dVdr = lj_raw_dr_kong(s6_kong, s12_kong, r);
//  double dVdrMAX = lj_raw_dr(atom_i->epsilon, atom_i->sigma, p->R_MAX);
//  double dVdr = lj_raw_dr(atom_i->epsilon, atom_i->sigma, r);
  if (fabs(dVdr) > fabs(dVdrMAX)){
    return -dVdrMAX;
  }
  else {
    return -dVdr;
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
