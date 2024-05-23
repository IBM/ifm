#ifndef __ENGINE_H
#define __ENGINE_H

#include <stdio.h>
#include <string>
#include <zlib.h>
#include "ifm_common.h"
#include "watershed.h"
#include "grid.h"

using namespace std;
class precip_table;

#define GZDUMP(fname,func) do { \
  gzFile outfp = gzopen((fname), "w9"); \
  (func)(outfp, -9); \
  gzclose(outfp); \
} while(0)

// Copies private data from real_t* members of the Watershed class
#define FRIENDLY_GET_METHOD(funcname,member) \
  Grid *funcname() { \
    Grid *g = new Grid; \
    g->cols = _dem->cols; \
    g->rows = _dem->rows; \
    g->cellsize = _dem->cellsize; \
    g->nodata = _dem->nodata; \
    g->data  = (real_t *) malloc(sizeof(real_t) * (g->cols * g->rows)); \
    memcpy(g->data, _ws.member, sizeof(real_t) * (g->cols * g->rows)); \
    return g; \
  }

#define FRIENDLY_SET_METHOD(funcname,member) \
  void funcname(Grid *in) { \
    memcpy(_ws.member, in->data, sizeof(real_t) * (in->cols * in->rows)); \
  }

class Engine {
  private:
    Grid       *_soil_hc;
    Grid       *_soil_ph;
    Grid       *_soil_ep;
    Grid       *_soilMoisture;
    Grid       *_dem;
    short      *_mask;
    short      *_landuse;
    uint64_t   *_out_x;
    uint64_t   *_out_y;
    real_t     *_out_s;
    real_t     *_manning;        // should be const.. we do not manage this pointer's allocation
    real_t     *_retention_coef; // should be const.. we do not manage this pointer's allocation
    uint64_t    _manning_size;
    uint64_t    _retention_coef_size;
    uint64_t    _num_outlets;
    WaterShed   _ws;
    precip_table *_pt;

    // conversion factors
    static const real_t _mmhr;
    static const real_t _mm2m;
    static const real_t _cmhr;
    static const real_t _cm2m;

  public:
    Engine():_soil_hc(NULL),_soil_ph(NULL),_soil_ep(NULL),_soilMoisture(NULL),
             _dem(NULL),_mask(NULL),_landuse(NULL),
             _out_x(NULL),_out_y(NULL),_out_s(NULL),
             _manning(NULL), _retention_coef(NULL), _manning_size(0), _retention_coef_size(0),
             _num_outlets(0),_pt(NULL) { }

    ~Engine() { }

    void setSoil(Grid *soil_hc, Grid *soil_ph, Grid *soil_ep, Grid *soilMoisture) {
      _soil_hc = soil_hc;
      _soil_ph = soil_ph;
      _soil_ep = soil_ep;
      _soilMoisture = soilMoisture;
    }

    void setLand(short *landuse, real_t *manning, real_t *retention_coef, uint64_t size) {
      assert ((manning && retention_coef && size>0) || (manning == NULL && retention_coef == NULL));
      _landuse = landuse;
      _manning = manning;
      _manning_size = size;
      _retention_coef = retention_coef;
      _retention_coef_size = size;
    }

    void setOutlets(uint64_t *out_x, uint64_t *out_y, real_t *out_s, uint64_t num) {
      _out_x = out_x;
      _out_y = out_y;
      _out_s = out_s;
      _num_outlets = num;
    }

    void setDomain(Grid *dem, short *mask) {
      _dem = dem;
      _mask = mask;
    }

    void setPrecipitation(precip_table *pt) {
      _pt = pt;
    }

    real_t getHeight(uint64_t x, uint64_t y) {
      return _ws.GetHeight(x, y);
    }

    void dumpHeight(char *fname) {
      GZDUMP(fname, _ws.DumpH);
    }

    void dumpVSAT(char *fname) {
      GZDUMP(fname, _ws.DumpVSAT);
    }

    void dumpMaxH(char *fname) {
      GZDUMP(fname, _ws.DumpMaxH);
    }

    void setup();
    void run(uint64_t kk, real_t dt);
    Grid* getDEM() {
        return _dem;
    }
    // The methods below have the following signature:
    // Grid *funcname();
    // The returned pointer needs to be deleted by the caller after used.
    FRIENDLY_GET_METHOD(getHeight,      _H);
    FRIENDLY_GET_METHOD(getOLR,         _OLR);
    FRIENDLY_GET_METHOD(getOLRDIM0_OLD, _OLRDIM0_OLD);
    FRIENDLY_GET_METHOD(getOLRDIM1_OLD, _OLRDIM1_OLD);
    FRIENDLY_GET_METHOD(getPRE,         _PRE);
    FRIENDLY_GET_METHOD(getRET,         _RET);
    FRIENDLY_GET_METHOD(getMaxH,        _MAXH);
    FRIENDLY_GET_METHOD(getIntH,        _INTH);
    FRIENDLY_GET_METHOD(getVSAT,        _VSAT);
    FRIENDLY_GET_METHOD(getVolume,      _VOL);
    
    // The methods below have the following signature:
    // void funcname(Grid *);
    FRIENDLY_SET_METHOD(setHeight,      _H);
    FRIENDLY_SET_METHOD(setOLR,         _OLR);
    FRIENDLY_SET_METHOD(setOLRDIM0_OLD, _OLRDIM0_OLD);
    FRIENDLY_SET_METHOD(setOLRDIM1_OLD, _OLRDIM1_OLD);
    FRIENDLY_SET_METHOD(setPRE,         _PRE);
    FRIENDLY_SET_METHOD(setRET,         _RET);
    FRIENDLY_SET_METHOD(setMaxH,        _MAXH);
    FRIENDLY_SET_METHOD(setIntH,        _INTH);
    FRIENDLY_SET_METHOD(setVSAT,        _VSAT);
};

#endif
