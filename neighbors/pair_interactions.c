#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "pair_interactions.h"
#include "../main.h"
#include "../models/models.h"
#include "../verlet.h"


int compute_pair_interactions(int NMAX,
                              double R_CUT,
                              struct Atoms atoms[],
                              double Region[3],
                              EnergyModel *energy_model,
                              ForceModel *force_model){
  // Create NN list with linked-lists
//  printf("DEBUG pair_interactions.c: creating lc\n");
  const int EMPTY = -1;
  int lc[3] = {
    floorl(Region[0]/R_CUT),
    floorl(Region[1]/R_CUT),
    floorl(Region[2]/R_CUT),
  };
//  printf("DEBUG pair_interactions.c: creating rc\n");
  double rc[3] = {
    Region[0]/lc[0],
    Region[1]/lc[1],
    Region[2]/lc[2],
  };
  int NCLMAX = lc[0] * lc[1] * lc[2];
  int lscl[NMAX];
  int head[NCLMAX];
  int lcyz = lc[1] * lc[2];
  int lcxyz = lcyz * lc[0];
  
  int mc[3];

//  printf("DEBUG pair_interactions.c: clearing head array\n");
  for (int c=0; c<lcxyz; c++) head[c] = EMPTY;

//  printf("DEBUG pair_interactions.c: creating cutoff radius\n");
  for (int i=0; i<NMAX; i++){
    float r[3] = {atoms[i].x, atoms[i].y, atoms[i].z};
    for (int a=0; a<3; a++) {
      mc[a] = r[a]/rc[a];
    }
    int c = mc[0] * lcyz + mc[1] * lc[2] + mc[2];
    lscl[i] = head[c];
    head[c] = i;
    }

  // Calculate
  double E_tot = 0;
  double F_tot = 0;
  int mcl[3]; 
  double rshift[3];

//  printf("DEBUG pair_interactions.c: scanning inner cells\n");
  for (mc[0]=0; mc[0]<lc[0]; (mc[0])++)
  for (mc[1]=0; mc[1]<lc[1]; (mc[1])++)
  for (mc[2]=0; mc[2]<lc[2]; (mc[2])++){  // Scan inner cells
    int c = mc[0] * lcyz + mc[1] * lc[2] + mc[2];  // Calculate a scalar cell index
//    printf("DEBUG pair_interactions.c: scanning neighbor cells\n");
    for (mcl[0]=mc[0]-1; mcl[0]<=mc[0]+1; (mcl[0])++)
    for (mcl[1]=mc[1]-1; mcl[1]<=mc[1]+1; (mcl[1])++)
    for (mcl[2]=mc[2]-1; mcl[2]<=mc[2]+1; (mcl[2])++){ // Scan neighbor cells
//      printf("DEBUG pair_interactions.c: unwrapping pbc\n");
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

//      printf("DEBUG pair_interactions.c: scanning i atoms\n");
      while (i != EMPTY){
        int j = head[cl]; // Scan atom j in cell cl
//        printf("DEBUG pair_interactions.c: scanning j atoms\n");
        while (j != EMPTY){
          if (i < j){
            float r_ijx = atoms[i].x - atoms[j].x + rshift[0];
            float r_ijy = atoms[i].y - atoms[j].y + rshift[1];
            float r_ijz = atoms[i].z - atoms[j].z + rshift[2];
            float r_ij = r_ijx*r_ijx + r_ijy*r_ijy + r_ijz*r_ijz;
//            printf("DEBUG pair_interactions.c: calculating neighbours\n");
            if (r_ij  < R_CUT * R_CUT){  // Calculate LJ potential
                E_tot += energy_model->raw_energy(
                energy_model->ctx,
                &atoms[i],
                &atoms[j]);
                F_tot += force_model->raw_force(
                force_model->ctx,
                &atoms[i],
                &atoms[j]);
                acceleration(force_model, &atoms[i], &atoms[j]);
            }
          }
          j = lscl[j];
        }
        i = lscl[i];
      }
    }
  }
  printf("Total E:         %3g eV    \n", E_tot);
  printf("Total Forces:    %3g eV/Ang\n", F_tot);
  return EXIT_SUCCESS;
}
