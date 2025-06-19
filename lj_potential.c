#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lj_potential.h"
#include "main.h"


static inline double lj_raw(float epsilon, float sigma, float r){
  float s_r = sigma/r;
  float s6 = powf(s_r, 6.0f);
  float V = 4 * epsilon * (s6*s6 - s6);
  return V;
}


static inline double lj_raw_dr(float epislon, float sigma, float r){
  float s_r = sigma/r;
  float s6 = powf(s_r, 6.0f);
  float dVdr = 4 * epislon * (-12 * s6 *s6/r + 6 * s6/r);
  return dVdr;
}


void lj_i_acc(float epsilon,
    float sigma,
    struct Atoms atoms_i,
    struct Atoms atoms_j){

  float r_ij_x = atoms_i.x - atoms_j.x;
  float r_ij_y = atoms_i.y - atoms_j.y;
  float r_ij_z = atoms_i.z - atoms_j.z;
  float r_ij = sqrtf(r_ij_x*r_ij_x + r_ij_y*r_ij_y + r_ij_z * r_ij_z);

  double r1 = sqrt(r_ij/sigma);  // double chek it
  double r2 = 1.0/r_ij;
  double r6 = r2 * r2 * r2;
  double fcVal = epsilon * 48.0 * r2 * r6 * (r6 - 0.5);

  double f_x = fcVal * r_ij_x;
  double f_y = fcVal * r_ij_y;
  double f_z = fcVal * r_ij_z;
  
  atoms_i.a_x += f_x / atoms_i.mass;
  atoms_i.a_y += f_y / atoms_i.mass;
  atoms_i.a_z += f_z / atoms_i.mass;
  
  atoms_j.a_x -= f_x / atoms_j.mass;
  atoms_j.a_y -= f_y / atoms_i.mass;
  atoms_j.a_z -= f_z / atoms_i.mass;
  }


void clear_acc(struct Atoms atom){
   atom.a_x = 0;
   atom.a_y = 0;
   atom.a_z = 0;
}

double lj_potential(float epsilon,
    float sigma,
    float r,
    float R_MAX){
  if (r >= R_MAX){
    return 0.0f;
  }
  // Calculate potential and slope at the R = R_cutoff
  float Vc = lj_raw(epsilon, sigma, R_MAX);
  float dVdr = lj_raw_dr(epsilon, sigma, R_MAX);

  // Calculate parts of LJ potential
  float V = lj_raw(epsilon, sigma, r);

  // Calculate LJ potential
  return V - Vc - (r - R_MAX)*dVdr;
}


int main_lj(int NMAX,
            int NCLMAX,
            double R_CUT,
            struct Atoms atoms[],
            double Region[3]){
  // Create NN list with linked-lists
  const int EMPTY = -1;
  int lc[3] = {
    floorl(Region[0]/R_CUT),
    floorl(Region[1]/R_CUT),
    floorl(Region[2]/R_CUT),
  };
  double rc[3] = {
    Region[0]/lc[0],
    Region[1]/lc[1],
    Region[2]/lc[2],
  };
  int lscl[NMAX];
  int head[NCLMAX];
  int lcyz = lc[1] * lc[2];
  int lcxyz = lcyz * lc[0];
  
  int mc[3];
  for (int c=0; c<lcxyz; c++) head[c] = EMPTY;

  for (int i=0; i<NMAX; i++){
    float r[3] = {atoms[i].x, atoms[i].y, atoms[i].z};
    for (int a=0; a<3; a++) {
      mc[a] = r[a]/rc[a];
    }
    int c = mc[0] * lcyz + mc[1] * lc[2] + mc[2];
    lscl[i] = head[c];
    head[c] = i;
    }

  // Clear forces
  for (int z=0; z<NMAX; z++) clear_acc(atoms[z]);

  // Calculate LJ
  double E_tot = 0;
  double F_tot = 0;
  int mcl[3]; 
  double rshift[3];

  for (mc[0]=0; mc[0]<lc[0]; (mc[0])++)
  for (mc[1]=0; mc[1]<lc[1]; (mc[1])++)
  for (mc[2]=0; mc[2]<lc[2]; (mc[2])++){  // Scan inner cells
    int c = mc[0] * lcyz + mc[1] * lc[2] + mc[2];  // Calculate a scalar cell index
    for (mcl[0]=mc[0]-1; mcl[0]<=mc[0]+1; (mcl[0])++)
    for (mcl[1]=mc[1]-1; mcl[1]<=mc[1]+1; (mcl[1])++)
    for (mcl[2]=mc[2]-1; mcl[2]<=mc[2]+1; (mcl[2])++){ // Scan neighbor cells
      for (int a=0; a<3; a++){ // Unwrapping the pbc
          if (mcl[a] < 0) {
            rshift[a] = -Region[a];
          }
          else if (mcl[a] >= lc[a]) {
            rshift[a] = Region[a];
          }
          else {
          rshift[a] = 0.0;
          }
      }

      int cl = ((mcl[0] + lc[0]) % lc[0]) * lcyz
        + ((mcl[1] + lc[1]) % lc[1]) * lc[2]
        + ((mcl[2] + lc[2]) % lc[2]);  // Scalar cell index of the neighor cell
      int i = head[c];  // Scan atom i in cell c
      while (i != EMPTY){
        int j = head[cl]; // Scan atom j in cell cl
        while (j != EMPTY){
          if (i < j){
            float r_ijx = atoms[i].x - atoms[j].x + rshift[0];
            float r_ijy = atoms[i].y - atoms[j].y + rshift[1];
            float r_ijz = atoms[i].z - atoms[j].z + rshift[2];
            float r_ij = r_ijx*r_ijx + r_ijy*r_ijy + r_ijz*r_ijz;
            if (r_ij  < R_CUT * R_CUT){  // Calculate LJ potential
              float r = sqrtf(r_ij);
              double E_i = lj_potential(atoms[i].epsilon, atoms[i].sigma, r, R_CUT);
              E_tot += E_i;
              double F_i = lj_raw_dr(atoms[i].epsilon, atoms[i].sigma, r);
              F_tot += F_i;
              lj_i_acc(atoms[i].epsilon,
                  atoms[i].sigma,
                  atoms[i],
                  atoms[j]);
            }
          }
          j = lscl[j];
        }
        i = lscl[i];
      }
    }
  }
  printf("LJ Potential: %3lf eV\n", E_tot);
  printf("LJ Forces:    %3lf eV/Ang\n", F_tot);
  return EXIT_SUCCESS;
}
