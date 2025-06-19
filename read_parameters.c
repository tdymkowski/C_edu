#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include "xyz_reader.h"
#include "main.h"

#define MAX_LINE 128
#define MAX_TYPE 16


typedef struct {
    char species_i[MAX_TYPE];
    double cutoff;
    double epsilon;
    double sigma;
} LJEntry;


static void trim(char *s){
  char *end;
  while (isspace((unsigned char)*s)) s++;
  if (*s == 0) return;
  end = s + strlen(s) - 1;
  while (end > s && isspace((unsigned char)*end)) end--;
  *(end+1) = '\0';
}

int main_params_reader(struct Atoms atoms[], int N_MAX){
  char filename[] = "LJ.params";
  FILE *fptr = fopen(filename, "r");
  char line[MAX_LINE];
  LJEntry entry;
  goto_line(fptr, 7);
  int line_no = 0;
  while (fgets(line, sizeof(line), fptr)){
    line_no++;
    char *p = line;
    
    while (isspace((unsigned char)*p)) p++;
    if (*p == '#' || *p == '\0') continue;

    // Check species i 
    char *tok = strtok(line, ",");
    if (!tok){
      fprintf(stderr, "Malformed line %d\n", line_no);
      continue;
    }
    strncpy(entry.species_i, tok ? tok : "", MAX_TYPE - 1);
    trim(entry.species_i);
   
    // Skip species j
    strtok(NULL, ",");

    // Get cutoff
    tok = strtok(NULL,",");
    entry.cutoff = tok ? strtod(tok, NULL) : 0.0;

    // get epsilon
    tok = strtok(NULL,",");
    entry.epsilon = tok ? strtod(tok, NULL) : 0.0;
    
    // get sigma
    tok = strtok(NULL,",");
    entry.sigma = tok ? strtod(tok, NULL) : 0.0;

    // Compare species_i with atom[i].symbol
    bool printed = false;
    for (int i = 0; i < N_MAX; i++){
      if (strcmp(entry.species_i, atoms[i].symbol) == 0){
        if (!printed && strcmp(entry.species_i, atoms[i].symbol) == 0){
        fprintf(stdout, "Match for atom %s:\ncutoff=%.3f Ang\nε=%.6f eV\nσ=%.6f Ang\n",
            atoms[i].symbol, entry.cutoff, entry.epsilon, entry.sigma);
        printed = true;
        }
        atoms[i].r_cut = entry.cutoff;
        atoms[i].epsilon = entry.epsilon;
        atoms[i].sigma = entry.sigma;
      }
    }

  }

  return EXIT_SUCCESS;
}
