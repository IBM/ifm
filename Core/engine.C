/*
 * @brief   Defines top level simulation driver.
 * @note    StarPU powers shared-memory parallelism.
 *
 * @author <main author>
 * @email  <main author's email>
 * @author  Maksims Abalenkovs
 * @email   maksims.abalenkovs@stfc.ac.uk
 * @date    Oct 8, 2024
 * @version 1.4
 */

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <starpu.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include <zlib.h>

#include "engine.h"
#include "filter.h"
#include "ifm_common.h"
#include "netcdf.h"
#include "rain_table.h"

#define MANNINGS \
  0.246,  0.41,  0.235, 0.184, 0.02, 0.15, 0.090, 0.24,  0.400, 0.450, \
  0.192,  0.012, 0.035, 0.17,  0.15, 0.05, 0.025, 0.13,  0.029, 0.320, \
  0.0137, 0.035, 0.24,  0.192, 0.17, 0.15, 0.184, 0.192, 0.40,  0.05, \
  0.235

#define RET_COEF \
  0.0,    0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0, \
  0.0,    0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0, \
  0.0,    0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0, \
  0.0

const real_t Engine::_mmhr = 3600.0*1000.0; // mm/hr -> m/s
const real_t Engine::_mm2m = 1000.0;        // mm->m
const real_t Engine::_cmhr = 3600.0*100.0;  // cm/hr -> m/s
const real_t Engine::_cm2m = 100.0;     

void Engine::setup()
{
  assert(_dem);
  assert(_mask);
  assert(_landuse);
  assert(_out_x);
  assert(_out_y);
  assert(_out_s);
  assert(_pt);

  /* Manning's N, no conversion needed, index start @ 0 !! */
  real_t default_manning[] = { MANNINGS }; // 20 values, start at 0
  if (! _manning) {
    _manning = default_manning;
    _manning_size = sizeof(default_manning)/sizeof(default_manning[0]);
  }
 
  /* Interception, depends on land-use, in mm, should be the same size as manning's N */
  real_t default_retention_coef[] = { RET_COEF };
  if (! _retention_coef) {
    _retention_coef = default_retention_coef;
    _retention_coef_size = sizeof(default_retention_coef)/sizeof(default_retention_coef[0]);
  }

  /*
     We now read the values from file, but keep them here as reference

     Hydraulic conductivity, cm/hr 
     real_t infilt_coef1[] = {0.0, 0.336, 0.3072, 0.3552, 0.3648, 0.3456, 0.432, 0.384}; // start@1,

     Pressure head @ wetting front, in cm 
     real_t infilt_coef2[] = {0.0, 22, 14, 17, 22, 18, 22, 15 };  // start @ 1, pressure head

     Solid moisture deficit, no conversion needed. Consider splitting into effective
     porosity and initial moisture content 
     real_t infilt_coef3[] = {0.0, 0.29,0.29,0.29,0.29,0.29,0.29,0.29};  // start @ 1
   */

  // convert the input coefficients into MKS
  for (uint64_t ii = 0; ii < _retention_coef_size; ii++) _retention_coef[ii] /= _mm2m;

  /* this is no longer needed, but we keep them here
     for (uint64_t ii=0; ii<(uint64_t)(sizeof(infilt_coef1)/sizeof(real_t)); ii++) infilt_coef1[ii] /= _cmhr;
     for (uint64_t ii=0; ii<(uint64_t)(sizeof(infilt_coef2)/sizeof(real_t)); ii++) infilt_coef2[ii] /= _cm2m;
   */
 
  fix_boundaries(_dem->data, _mask, _dem->rows, _dem->cols, _dem->nodata, -8.0);
#ifdef PRE_FILTER
  gaussian_filter(_dem, 2, 1.4);
#endif

  // Setup
  _ws.Setup(_dem->rows, _dem->cols, (real_t)_dem->cellsize, 0.0, 0.0);  // not providing llx and lly values
  _ws.Init(_dem->data, _mask, NULL, NULL);  // no lakes, initial value depth set to zero
  _ws.SetN(_landuse, _manning, _retention_coef, _manning_size);
  _ws.SetOutlets(_num_outlets, _out_x, _out_y, _out_s);
  if (_soil_hc && _soil_ph && _soil_ep && _soilMoisture)
    _ws.SetSoil(_soil_hc->data, 1.0/_cmhr, _soil_ph->data, 1.0/_cm2m, _soil_ep->data, 1.0, _soilMoisture->data);
}

void Engine::run(uint64_t kk, real_t dt, bool infiltrate_flag)
{
  // uint64_t pre_cntr = (uint64_t)(pre_window/dt);

  // process the precipitation (specified in mm/hr)
  if (_pt->isUniform()) {
    real_t pre_rate = _pt->Get( dt*kk );
    // printf("T=%.1f Use precipitation rate %.1f mm/hr\n", dt*kk,pre_rate);
    _ws.SetPrecipitation(pre_rate/_mmhr);
  } else {
    Grid *pp = _pt->Get_Grid( dt*kk );
    _ws.SetPrecipitation(pp->data, _mmhr);
  }
#ifndef ROUTING_ONLY
  // _ws.CompIntercept(dt);
  _ws.comp_intercept_starpu(dt);
#endif
  // _ws.CompOverlandDepth(dt);
  _ws.comp_overland_depth_starpu(dt);
#ifndef ROUTING_ONLY
  if (infiltrate_flag && _soil_hc && _soil_ph && _soil_ep) {
    // _ws.CompInfiltration(dt);
    _ws.comp_infiltration_starpu(dt);
  }
#endif
  _ws.CompDiffusiveRouting(dt);
  // _ws.comp_diffusive_routing_starpu(dt);
  // _ws.CompOutlet(dt);
  _ws.comp_outlet_starpu(dt);
  _ws.comp_storm_starpu(dt);
}

void Engine::profile_run(uint64_t kk, real_t dt, bool infiltrate_flag,
    real_t *t_intercept, real_t *t_overland, real_t *t_infiltration,
    real_t *t_diffusive, real_t *t_outlet)
{
  // uint64_t pre_cntr = (uint64_t)(pre_window/dt);

  // process the precipitation (specified in mm/hr)
  if (_pt->isUniform()) {
    real_t pre_rate = _pt->Get( dt*kk );
    // printf("T=%.1f Use precipitation rate %.1f mm/hr\n", dt*kk,pre_rate);
    _ws.SetPrecipitation(pre_rate/_mmhr);
  } else {
    Grid *pp = _pt->Get_Grid( dt*kk );
    _ws.SetPrecipitation(pp->data, _mmhr);
  }
#ifndef ROUTING_ONLY
  double t_intercept_0 = starpu_timing_now();
  // _ws.CompIntercept(dt);
  _ws.comp_intercept_starpu(dt);
  *t_intercept += (starpu_timing_now() - t_intercept_0);
#endif
  double t_overland_0  = starpu_timing_now();
  // _ws.CompOverlandDepth(dt);
  _ws.comp_overland_depth_starpu(dt);
  *t_overland += (starpu_timing_now() - t_overland_0);
#ifndef ROUTING_ONLY
  if (infiltrate_flag && _soil_hc && _soil_ph && _soil_ep) {
    double t_infiltration_0 = starpu_timing_now();
    // _ws.CompInfiltration(dt);
    _ws.comp_infiltration_starpu(dt);
    *t_infiltration += (starpu_timing_now() - t_infiltration_0);
  }
#endif
  double t_diffusive_0 = starpu_timing_now();
  _ws.CompDiffusiveRouting(dt);
  // _ws.comp_diffusive_routing_starpu(dt);
  *t_diffusive += (starpu_timing_now() - t_diffusive_0);

  double t_outlet_0 = starpu_timing_now();
  // _ws.CompOutlet(dt);
  _ws.comp_outlet_starpu(dt);
  _ws.comp_storm_starpu(dt);
  *t_outlet += (starpu_timing_now() - t_outlet_0);
}

// @eof engine.C
