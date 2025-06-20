#include <math.h>
#include "../main.h"


double get_R(struct Atoms *i, struct Atoms *j){
  double r_x = i->x - j->x;
  double r_y = i->y - j->y;
  double r_z = i->z - j->z;
  double r2 = r_x * r_x + r_y * r_y + r_z * r_z;
  double r = sqrt(r2);
  return r;

}
