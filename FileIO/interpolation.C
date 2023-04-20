#include <limits.h>
#include "interpolation.h"
#include "engine.h"

#ifndef ULLONG_MAX
#define ULLONG_MAX 18446744073709551615ULL
#endif

bool Interpolation::getXY(Grid *latGrid, Grid *lonGrid, real_t lat, real_t lon, 
  uint64_t &sx, uint64_t &ex, uint64_t &sy, uint64_t &ey)
{
  if (lat > GRID(latGrid, 0, 0) || lat < GRID(latGrid, 0, latGrid->rows-1) ||
      lon < GRID(lonGrid, 0, 0) || lon > GRID(lonGrid, lonGrid->cols-1, 0))
    return false;

  sy = ULLONG_MAX;
  ey = ULLONG_MAX; 
  for (uint64_t i=0; i<latGrid->rows; ++i) {
    if (GRID(latGrid, 0, i) <= lat && sy == ULLONG_MAX)
      sy = i == 0 ? 0 : i-1;
    if (GRID(latGrid, 0, i) < lat && ey == ULLONG_MAX && sy != ULLONG_MAX) {
      ey = i;
      break;
    }
  }
  if (ey == ULLONG_MAX)
    ey = sy;

  sx = ULLONG_MAX;
  ex = ULLONG_MAX;
  for (int64_t i=lonGrid->cols-1; i>=0; --i) {
    if (GRID(lonGrid, i, 0) <= lon && ex == ULLONG_MAX)
      ex = (uint64_t) i == lonGrid->cols-1 ? lonGrid->cols-1 : i+1;
    if (GRID(lonGrid, i, 0) < lon && sx == ULLONG_MAX && ex != ULLONG_MAX) {
      sx = i;
      break;
    }
  }
  if (sx == ULLONG_MAX)
    sx = ex;

  if (ex == ULLONG_MAX || sy == ULLONG_MAX)
    return false;

  return true;
}
    
real_t Interpolation::bilinear(Grid *data, BoundingBox_t *dataBBox, 
  Grid *latGrid, Grid *lonGrid, BoundingBox_t *latlonBBox, real_t lat, real_t lon)
{
  assert(data);
  assert(dataBBox);
  assert(latGrid);
  assert(lonGrid);
  assert(latlonBBox);

  uint64_t x1 = latlonBBox->startx;
  uint64_t x2 = latlonBBox->endx;
  uint64_t y1 = latlonBBox->starty;
  uint64_t y2 = latlonBBox->endy;
  
  real_t x1val = GRID(lonGrid, x1, 0);
  real_t x2val = GRID(lonGrid, x2, 0);
  real_t y1val = GRID(latGrid, 0, y1);
  real_t y2val = GRID(latGrid, 0, y2);
  real_t offsetx = fabs(lon - x1val) / fabs(x2val - x1val);
  real_t offsety = fabs(lat - y1val) / fabs(y2val - y1val);

  // Linear interpolation between the values at columns x1 and x2, at rows y1 and y2
  x1 = dataBBox->startx;
  x2 = dataBBox->endx;
  y1 = dataBBox->starty;
  y2 = dataBBox->endy;
  real_t x1y1 = GRID(data, x1, y1) != data->nodata ? GRID(data, x1, y1) : 0;
  real_t x1y2 = GRID(data, x1, y2) != data->nodata ? GRID(data, x1, y2) : 0;
  real_t x2y1 = GRID(data, x2, y1) != data->nodata ? GRID(data, x2, y1) : 0;
  real_t x2y2 = GRID(data, x2, y2) != data->nodata ? GRID(data, x2, y2) : 0;
  real_t interpx1 = x1y1 * (1-offsetx) + x2y1 * offsetx;
  real_t interpx2 = x1y2 * (1-offsetx) + x2y2 * offsetx;

  // Bilinear interpolation
  //printf("%f %f, %f %f\n", GRID(data, x1, y1), GRID(data, x2, y1), GRID(data, x1, y2), GRID(data, x2, y2));
  return (1 - offsety) * interpx1 + offsety * interpx2;
}
