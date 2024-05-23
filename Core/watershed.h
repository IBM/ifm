/* definition in one water shed */

#ifndef _WATERSHED_H
#define _WATERSHED_H

#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "ifm_common.h"

// Macros for branch prediction hints
#ifdef __GNUC__
#define likely(x)    __builtin_expect((x),1)
#define unlikely(x)  __builtin_expect((x),0)
#else
#define likely(x)    (x)
#define unlikely(x)  (x)
#endif

#define DUMP_TO_FILE(F,nodata_tag,grid,print_func) do { \
  uint64_t jj, kk; \
  assert((F)); \
  (print_func)(F, "NCOLS %ld\n", _ncol); \
  (print_func)(F, "NROWS %ld\n", _nrow); \
  (print_func)(F, "XLLCORNER 0\nYLLCORNER 0\nCELLSIZE 90\nNODATA_VALUE %d\n", nodata_tag); \
  for (jj=0; jj<_nrow; jj++) { \
    for (kk=0; kk<_ncol; kk++) { \
      if ( _MASK[jj*_ncol+kk] )  (print_func)(F,"%9.3e ", (grid)[jj*_ncol+kk]); \
      else                       (print_func)(F,"%d ", nodata_tag); \
    } \
    print_func(F,"\n"); \
  } \
} while(0)

// Forward declaration
class Engine;

class WaterShed {
 private:
  real_t    *_ELE;        // elevation
  real_t    *_H;          // depth
  real_t    *_VOL;        // volume
  real_t    *_OLR;        // overland routing
  real_t    *_PRE;        // precipitation 
  real_t    *_N;          // manning's N
  real_t    *_STORE;      // storage ponds, or lakes
  real_t    *_RET;        // retention capability, for intercept

  real_t    *_MAXH;       // max depth
  real_t    *_INTH;       // integral of depth

                          // these three are needed by the infiltration model
  real_t    *_VSAT;       // saturation volume, 
  real_t    *_HCON;       // conductivity
  real_t    *_P2;         // second term in GA model
  //real_t    *_tmpinf1;  // temporary workspace 1 for infiltration
  //real_t    *_tmpinf2;  // temporary workspace 2 for infiltration

  short     *_MASK;       // mask to limit where the boundary of watershed
  short     *_IDX_SOIL;   // index to soil type, only support a few discrete values
  short     *_IDX_N;      // index to manning's N, only support a few discrete values

  /* for outlets */
  uint64_t  *_OUTLETS;
  uint64_t   _N_OUT;
  real_t    *_OUT_SLOPES;

  /* for storm drains */
  uint64_t  *_SD;         // locations of storm drains
  uint64_t   _N_STORM;    // number of storm drains
  real_t    *_threshold;  // threshold depth
  real_t    *_rate;       // coefficinet of rate, model as Q = rate * y^(8/3);
  real_t    *_saturation; // saturation point
  
  /* temp workspace */
  //real_t    *_tmpsf;
  //real_t    *_tmpn;
  //real_t    *_tmph;
  //real_t    *_tmpp;
  //real_t    *_tmpt;

  //real_t    *_OLRDIM0;
  //real_t    *_OLRDIM1;

  real_t    *_OLRDIM0_OLD;
  real_t    *_OLRDIM1_OLD;

  /* for time step control */
  real_t    *_minTx;
  real_t    *_minTy;

  /* it might be worthwhile to allocate the values into arrays, wiht duplications. The
     setup to the index array sucks when used with openMP */

  real_t    _gsz;     // grid size, assume equi-distance for X and Y
  real_t    _llx;     // llx coordinates
  real_t    _lly;     // lly coordinates;
  uint64_t  _nrow;    // dimensions, maximal size would be 32,766
  uint64_t  _ncol;
  uint64_t  _store_size;  // size of storage

  /* debugging control */
  bool      _printed_outlet_underflow;
  bool      _printed_depth_underflow;

 public:
  WaterShed() {
    _ELE = 0;
    _H = 0;
    _VOL = 0;
    _OLR = 0;
    _PRE = 0;
    _N = 0;
    _STORE = 0;
    _RET = 0;
    _VSAT = 0;
    _HCON = 0;
    _P2 = 0;
    _MAXH = 0;
    _INTH = 0;

    //_tmpinf1 = 0;
    //_tmpinf2 = 0;

    //_tmpsf = 0;
    //_tmpn = 0;
    //_tmph = 0;
    //_tmpp = 0;
    //_tmpt = 0;

    //_OLRDIM0 = 0;
    //_OLRDIM1 = 0;
    _OLRDIM0_OLD = 0;
    _OLRDIM1_OLD = 0;

    _OUTLETS = 0;
    _OUT_SLOPES = 0;
    _N_OUT = 0;

    _SD = 0;
    _N_STORM = 0;
    _threshold = 0;
    _rate = 0;
    _saturation = 0;

    _MASK = 0;
    _IDX_SOIL = 0;
    _IDX_N = 0;
 
    _gsz = _llx = _lly = 0.0;
    _nrow = _ncol = 0;

    _minTx = 0;
    _minTy = 0;
    // need to add more

    _printed_outlet_underflow = false;
    _printed_depth_underflow = false;
  }

  ~WaterShed() {
    if ( _ELE ) { free(_ELE); _ELE=0; }
    if ( _H ) { free(_H); _H=0; }
    if ( _MAXH) {free(_MAXH); _MAXH=0; }
    if ( _INTH) {free(_INTH); _INTH=0; }

    if ( _VOL ) { free(_VOL); _VOL=0; }
    if ( _OLR ) { free(_OLR); _OLR=0; }
    if ( _PRE ) { free(_PRE); _PRE = 0; }
    if ( _N ) { free(_N); _N = 0; }
    if ( _STORE ) { free(_STORE); _STORE = 0; }
    if ( _RET ) { free(_RET); _RET = 0; }
    if ( _VSAT ) { free(_VSAT); _VSAT = 0; }
    if ( _HCON ) { free(_HCON); _HCON = 0; }
    if ( _P2 ) { free(_P2); _P2 = 0; }
    //if ( _tmpinf1 ) { free(_tmpinf1); _tmpinf1 = 0; }
    //if ( _tmpinf2 ) { free(_tmpinf2); _tmpinf2 = 0; }

    if ( _OUTLETS ) { free(_OUTLETS); _OUTLETS = 0; }
    if ( _OUT_SLOPES) { free(_OUT_SLOPES); _OUT_SLOPES=0;}

    if ( _SD ) { free(_SD); _SD = 0; }
    if ( _threshold ) { free(_threshold); _threshold = 0; }
    if ( _rate ) { free(_rate); _rate = 0; }
    if ( _saturation) { free(_saturation); _saturation=0;}

    //if ( _OLRDIM0 ) { free(_OLRDIM0); _OLRDIM0=0; }
    //if ( _OLRDIM1 ) { free(_OLRDIM1); _OLRDIM1=0; }

    if ( _OLRDIM0_OLD ) { free(_OLRDIM0_OLD); _OLRDIM0_OLD=0; }
    if ( _OLRDIM1_OLD ) { free(_OLRDIM1_OLD); _OLRDIM1_OLD=0; }
    
    //if ( _tmpsf ) { free(_tmpsf); _tmpsf=0; }
    //if ( _tmpn ) { free(_tmpn); _tmpn=0; }
    //if ( _tmph ) { free(_tmph); _tmph=0; }
    //if ( _tmpp ) { free(_tmpp); _tmpp=0; }
    //if ( _tmpt ) { free(_tmpt); _tmpt=0; }

    if ( _IDX_SOIL ) { free(_IDX_SOIL); _IDX_SOIL=0; }
    if ( _MASK ) { free(_MASK); _MASK=0; }
    if ( _IDX_N ) { free(_IDX_N); _IDX_N=0; }

    if ( _minTx ) { free(_minTx); _minTx = 0; }
    if ( _minTy ) { free(_minTy); _minTy = 0; }
  }

  real_t GetHeight(uint64_t x, uint64_t y) {
    uint64_t idx = y*_ncol+x;
    assert(_H);
    assert(idx < _store_size);
    return _H[idx];
  }

  void Setup(uint64_t nrow, uint64_t ncol, real_t gsz, real_t llx, real_t lly);
  void Init(real_t *ele, short *mask, real_t *lakes, real_t *init_h=NULL);
  void SetN(short *idx_n, real_t *nvals, real_t *rvals, uint64_t nvals_size);
  void SetSoil(short *idx_s, real_t *p1, real_t *p2, real_t *p3, real_t *p4);
  void SetSoil(real_t *soil_hc, real_t f1, real_t *soil_ph, real_t f2, real_t *soil_ep, real_t f3, real_t *soilMoisture);
  void SetOutlets(uint64_t nout, uint64_t *xx, uint64_t *yy, real_t *slopes);
  void SetStormDrain(uint64_t nout, uint64_t *xx, uint64_t *yy, real_t *offset, real_t *rates, real_t *saturation);

  void SetPrecipitation(real_t *pa, real_t div) { 
    uint64_t jj;
    assert(pa); 
#ifdef PARA			//yh@June 27th
#pragma omp parallel for default(shared) private(jj)
#endif
    for(jj=0;jj<_store_size; jj++) 
      _PRE[jj]=pa[jj]/div; 
  }

  void SetPrecipitation(real_t prep) {
    uint64_t jj;
#ifdef PARA			//yh@June 27th
#pragma omp parallel for default(shared) private(jj)
#endif
    for(jj=0;jj<_store_size; jj++) 
      _PRE[jj]=prep;
  }

  void CompIntercept(real_t dt);
  int  CompOverlandDepth(real_t dt);
  int  CompInfiltration(real_t dt);
  int  CompDiffusiveRouting(real_t dt);
  int  CompOutlet(real_t dt);

  void DumpH(FILE *F, int nodata) {
    DUMP_TO_FILE(F, nodata, _H, fprintf);
  }
  void DumpMaxH(FILE *F, int nodata) {
    DUMP_TO_FILE(F, nodata, _MAXH, fprintf);
  }
  void DumpIntH(FILE *F, int nodata) {
    DUMP_TO_FILE(F, nodata, _INTH, fprintf);
  }
  void DumpOLR(FILE *F, int nodata) {
    DUMP_TO_FILE(F, nodata, _OLR, fprintf);
  }
  void DumpVSAT(FILE *F, int nodata) {
    DUMP_TO_FILE(F, nodata, _VSAT, fprintf);
  }
  void DumpPRE(FILE *F, int nodata) {
    DUMP_TO_FILE(F, nodata, _PRE, fprintf);
  }
  void DumpRET(FILE *F, int nodata) {
    DUMP_TO_FILE(F, nodata, _RET, fprintf);
  }
  // overloaded versions to deal with gzipped files
  void DumpH(gzFile F, int nodata) {
    DUMP_TO_FILE(F, nodata, _H, gzprintf);
  }
  void DumpMaxH(gzFile F, int nodata) {
    DUMP_TO_FILE(F, nodata, _MAXH, gzprintf);
  }
  void DumpIntH(gzFile F, int nodata) {
    DUMP_TO_FILE(F, nodata, _INTH, gzprintf);
  }
  void DumpOLR(gzFile F, int nodata) {
    DUMP_TO_FILE(F, nodata, _OLR, gzprintf);
  }
  void DumpVSAT(gzFile F, int nodata) {
    DUMP_TO_FILE(F, nodata, _VSAT, gzprintf);
  }
  void DumpPRE(gzFile F, int nodata) {
    DUMP_TO_FILE(F, nodata, _PRE, gzprintf);
  }
  void DumpRET(gzFile F, int nodata) {
    DUMP_TO_FILE(F, nodata, _RET, gzprintf);
  }
  // need to add channel routing

  // Let the Engine class access our private members
  friend class Engine;
};

#endif

/* end */
