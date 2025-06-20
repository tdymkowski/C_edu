// physical models
#ifndef MODEL_PARSER_H
#define MODEL_PARSER_H

#include <stdio.h>
#include <math.h>
#include "models.h"
#include "../main.h"

// parser
PairModels choose_model(const char *model_name,
                        struct Atoms atoms[],
                        int N_MAX,
                        float R_MAX);


#endif

