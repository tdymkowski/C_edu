// xyz_reader.h
#ifndef XYZ_READER_H
#define XYZ_READER_H

#include <stdio.h>
#include <math.h>
#include "main.h"

int goto_line(FILE *fptr, int target_line);
int get_max_atoms(char filename[]);
int read_coords(FILE *fptr, struct Atoms coords[], int max_num);
void print_coords(FILE *out, struct Atoms coords[], int num);
int main_xyz_reader(char filename[], struct Atoms coords[], int NMAX);

#endif
