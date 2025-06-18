#include <stdio.h>
#include <math.h>
#include "lj_potential.h"
#include "xyz_reader.h"

float lj_potential(float epsilon,
    float sigma, float r){
  float V;
  float A_1;
  float A_2;
  // Calculate parts of LJ potential
  A_1 = powf(sigma/r, 12.0f);
  A_2 = powf(sigma/r, 6.0f);
  // Calculate LJ potential
  V = 4 * epsilon * (A_1 - A_2);
  return V;
}


int main_lj(double epsilon,
            double sigma,
            int NMAX,
            int NCLMAX,
            float R_CUT,
            //float r[][3,
            struct Atoms atoms[],
            double Region[3]){
  /*
  // LJ Epsilon
  double epsilon = 3.1401000;  // eV
  // LJ sigma
  double sigma =  2.4232400;  // Angstroms
  // Atom number
  const int NMAX = 5;
  // Cell number
  const int NCLMAX = 4;
  // Cutoff value
  const float R_CUT = 9.6929800;
  // Positions
  float r[][3] = {
    {0., 0., 1.},
    {1., 1., 0.},
    {2., 0., 4.},
    {3., 0., 1.},
    {4., 0., 0.}
  };
  // Cell size
  double Region[3] = {10., 10., 10.};
  */

  // Create NN list with linked-lists
  const int EMPTY = -1;
  int lc[3] = {
    floorl(Region[0]/R_CUT),
    floorl(Region[1]/R_CUT),
    floorl(Region[2]/R_CUT),
  };
  double rc[3] = {
    floorl(Region[0]/lc[0]),
    floorl(Region[1]/lc[1]),
    floorl(Region[2]/lc[2]),
  };
  double lscl[NMAX];
  double head[NCLMAX];
  double lcyz = lc[1] * lc[2];
  double lcxyz = lcyz * lc[0];
  
  double mc[3];
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
  // Calculate LJ
  double E_tot = 0;
  int mcl[3]; 
  double rshift[3];
  for (mc[0]=0; mc[0]<lc[0]; (mc[0])++)
  for (mc[1]=0; mc[1]<lc[0]; (mc[1])++)
  for (mc[2]=0; mc[2]<lc[0]; (mc[2])++){  // Scan inner cells
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
              E_tot += lj_potential(epsilon, sigma, r_ij);
            }
          }
          j = lscl[j];
        }
        i = lscl[i];
      }
    }
  }
  printf("LJ Potential: %lf eV\n", E_tot);
  return 0;
}
