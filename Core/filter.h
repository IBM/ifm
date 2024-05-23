#ifndef __FILTER_H
#define __FILTER_H

#include <stdio.h>
#include "ifm_common.h"
#include "netcdf.h"

class Grid;

int gaussian_filter(Grid *grid, int64_t, real_t);
int fix_boundaries(real_t*, short*, uint64_t, uint64_t, real_t, real_t);

#endif
