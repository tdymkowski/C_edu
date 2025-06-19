// main.h
#ifndef MAIN_H
#define MAIN_H

#include <stdio.h>
#include <math.h>

struct Atoms {
    char symbol[3];
    float x;
    float y;
    float z;
    float epsilon;
    float sigma;
    float r_cut;
};

struct Params {
    char symbol[8];
    float epsilon;
    float sigma;
    float r_cut;

};

#endif
