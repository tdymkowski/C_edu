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
    float vel_x;
    float vel_y;
    float vel_z;
    float a_x;
    float a_y;
    float a_z;
    float mass;
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
