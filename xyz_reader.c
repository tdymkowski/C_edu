#include <stdio.h>
#include <stdlib.h>
#include "xyz_reader.h"
#include "main.h"


int goto_line(FILE *fptr, int targer_line){
  if (targer_line <= 1){
    rewind(fptr);
    return EXIT_SUCCESS;
  }
  rewind(fptr);
  char buf[256];
  for (int i = 1; i < targer_line; i++){
    if (!fgets(buf, sizeof buf, fptr)){
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}


int get_max_atoms(char filename[]){
  int NMAX;
  FILE *fptr = fopen(filename, "r");
  if (fscanf(fptr, "%d", &NMAX) != 1){
    fprintf(stderr, "Error: could not read atom count from first line\n");
    return EXIT_FAILURE;
  }
  // Reset the pointer
  goto_line(fptr, 1);
  return NMAX;
}


int read_coords(FILE *fptr, struct Atoms coords[], int max_num){
  int num = 0;
  if (goto_line(fptr, 3) != 0){
    return EXIT_FAILURE;
  }
   
  printf("%d\n", max_num);
  while (num < max_num && fscanf(fptr, "%3s %f %f %f",
                                 coords[num].symbol,
                                 &coords[num].x,
                                 &coords[num].y,
                                 &coords[num].z) == 4){
    ++num;
  }
  return num;
}


void print_coords(FILE *fptr, struct Atoms coords[], int num){
  fprintf(fptr, "x y z\n");
  for (int i = 0; i < num; ++i){
    fprintf(fptr, "%3s %f %f %f\n", coords[i].symbol, coords[i].x, coords[i].y, coords[i].z);
  }
}

int main_xyz_reader(char filename[], struct Atoms coords[], int NMAX){
  FILE *fp = fopen(filename, "r");
  if (!fp){
    fprintf(stderr, "Error: could not open '%s'\n", filename);
    return EXIT_FAILURE;
  }
  // Get number of atoms
//  int NMAX = get_max_atoms(fp);
  if (NMAX < 0){
    fclose(fp);
    return EXIT_FAILURE;
  }

//  struct Atoms coords[NMAX];

  int num = read_coords(fp, coords, NMAX);
  if (num < 0){
    fprintf(stderr, "Error: '%s' is missing the header or has an invalid header\n", filename);
    fclose(fp);
    return EXIT_FAILURE;
  }
  print_coords(stdout, coords, num);
  fclose(fp);
  return EXIT_SUCCESS;
}
