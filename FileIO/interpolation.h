#ifndef __INTERPOLATION_H
#define __INTERPOLATION_H

#include <math.h>
#include "ifm_common.h"
#include "grid.h"

class Interpolation {
  public:
    Interpolation() { };
    ~Interpolation() { };

    bool   getXY(Grid *latGrid, Grid *lonGrid, real_t lat, real_t lon, uint64_t &sx, uint64_t &ex, uint64_t &sy, uint64_t &ey);
    real_t bilinear(Grid *data, BoundingBox_t *dataBBox, Grid *latGrid, Grid *lonGrid, BoundingBox_t *latlonBBox, real_t lat, real_t lon);
};

#endif
