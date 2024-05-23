/* methods for watershed */
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <omp.h>

#include "watershed.h"

#define MY_SIGN(A) ( (A) >=0 ? 1 : -1)
#define MY_ISNAN(A) ( (A) != (A) )

#define EIGHT3RD 2.66666667

void WaterShed::Setup(uint64_t nrow, uint64_t ncol, real_t gsize, real_t llx, real_t lly) {
  //uint64_t mysize;

  assert( nrow>0 && ncol>0);
  assert( gsize > 0.0 );

  _nrow = nrow;
  _ncol = ncol;

  _store_size = nrow * ncol;
  _gsz = gsize;
  _llx = llx;
  _lly = lly;

  _ELE = (real_t*) malloc(sizeof(real_t) * _store_size);
  _H = (real_t*) malloc(sizeof(real_t) * _store_size);
  _MAXH = (real_t*) malloc(sizeof(real_t) * _store_size);
  _INTH = (real_t*) malloc(sizeof(real_t) * _store_size);

  _VOL = (real_t*) malloc(sizeof(real_t) * _store_size);
  _OLR = (real_t*) malloc(sizeof(real_t) * _store_size);
  _PRE = (real_t*) malloc(sizeof(real_t) * _store_size);
  _N = (real_t*) malloc(sizeof(real_t) * _store_size);
  _STORE = (real_t*) malloc(sizeof(real_t) * _store_size);
  _RET = (real_t*) malloc(sizeof(real_t) * _store_size);
  _VSAT = (real_t*) malloc(sizeof(real_t) * _store_size);
  _HCON = (real_t*) malloc(sizeof(real_t) * _store_size);
  _P2 = (real_t*) malloc(sizeof(real_t) * _store_size);

  //_tmpinf1 = (real_t*) malloc(sizeof(real_t) * _store_size);
  //_tmpinf2 = (real_t*) malloc(sizeof(real_t) * _store_size);

  //_OLRDIM0 = (real_t*) malloc(sizeof(real_t) * _store_size);
  //_OLRDIM1 = (real_t*) malloc(sizeof(real_t) * _store_size);

  _OLRDIM0_OLD = (real_t*) malloc(sizeof(real_t) * _store_size);
  _OLRDIM1_OLD = (real_t*) malloc(sizeof(real_t) * _store_size);

  _minTx = (real_t*)malloc(sizeof(real_t)* _store_size);
  _minTy = (real_t*)malloc(sizeof(real_t)* _store_size);

  _MASK = (short*) malloc(sizeof(short) * _store_size);
  _IDX_SOIL = (short*) malloc(sizeof(short) * _store_size);
  _IDX_N = (short*) malloc(sizeof(short) * _store_size);


  //mysize = nrow >= ncol ? nrow : ncol; // max
  //_tmpsf = (real_t*)malloc(sizeof(real_t) * mysize);
  //_tmpn = (real_t*)malloc(sizeof(real_t) * mysize);
  //_tmph = (real_t*)malloc(sizeof(real_t) * mysize);
  //_tmpp = (real_t*)malloc(sizeof(real_t) * mysize);
  //_tmpt = (real_t*)malloc(sizeof(real_t) * mysize);
}

void WaterShed::Init(real_t *elevation, short *mask, real_t *lakes, real_t *init_h) {
  // limited check performed, use at your risk
  assert(elevation);
  assert(mask);

  memcpy(_ELE, elevation, _store_size*sizeof(real_t));
  memcpy(_MASK, mask,      _store_size*sizeof(short));

  if ( init_h )  memcpy(_H, init_h,    _store_size*sizeof(real_t));  
  else           memset(_H, 0,         _store_size*sizeof(real_t));

  memset(_MAXH, 0, _store_size*sizeof(real_t));
  memset(_INTH, 0, _store_size*sizeof(real_t));
  memset(_VOL,  0, _store_size*sizeof(real_t));

  if ( lakes )   memcpy(_STORE, lakes,   _store_size*sizeof(real_t));  
  else           memset(_STORE, 0,      _store_size*sizeof(real_t));

  memset(_OLR, 0, _store_size*sizeof(real_t));
  //memset(_OLRDIM0, 0, _store_size*sizeof(real_t));
  //memset(_OLRDIM1, 0, _store_size*sizeof(real_t));

  memset(_minTx, 0, _store_size*sizeof(real_t));
  memset(_minTy, 0, _store_size*sizeof(real_t));

  uint64_t jj;
#ifdef PARA //yh@June 27th
  #pragma omp parallel for default(shared) private(jj)
#endif
  for (jj=0; jj<_store_size; jj++) {
    // quickly go over the mask, set values to NAN if not within mask
    if ( !_MASK[jj] ) _ELE[jj] = NAN;
    _OLRDIM0_OLD[jj]=1e3;   // hope they are large enough
    _OLRDIM1_OLD[jj]=1e3;
  }
}

// copy and set up manning's n
void WaterShed::SetN(short *idx_n, real_t *nvals, real_t *retcoef, uint64_t nvals_size) {
  uint64_t jj;

  assert( idx_n );
  assert( nvals );
  assert( retcoef);

  memcpy(_IDX_N, idx_n, _store_size*sizeof(short));

  // this could be slow, and uses memory
  // but we want to trade in for the runtime efficiency
#ifdef PARA //yh@June 27th
  #pragma omp parallel for default(shared) private(jj)
#endif
  for (jj=0; jj<_store_size; jj++) {
    assert( _IDX_N[jj] >= 0 );
    if ((uint64_t) _IDX_N[jj] >= nvals_size) {
      fprintf(stderr, "Error: LandUse data contains values greater than %ld (ManningsN and RetCoef array sizes)\n", nvals_size);
      abort();
    }
    _N[jj] = nvals[ _IDX_N[jj] ];
    _RET[jj] = retcoef[ _IDX_N[jj] ];
  }
}

// copy and set up G-A soil model
void WaterShed::SetSoil(short *idx_s, real_t *p1, real_t *p2, real_t *p3, real_t *p4) {
  uint64_t jj;

  assert( idx_s );
  assert( p1 );
  assert( p2 );
  assert( p3 );
  assert( p4 );

  memcpy(_IDX_SOIL, idx_s, _store_size*sizeof(short));

  for (jj=0; jj<_store_size; jj++) {
    assert( _IDX_SOIL[jj] >= 0 );
    _HCON[jj] = p1[ _IDX_SOIL[jj] ];
    _P2[jj] = p2[ _IDX_SOIL[jj] ] * p3[ _IDX_SOIL[jj] ];
  }

  // initializ VSAT
  memcpy(_VSAT, p4, _store_size*sizeof(real_t));
}


// overloaded function for soil setup
void WaterShed::SetSoil(real_t *soil_hc, real_t f1, real_t *soil_ph, real_t f2, real_t *soil_ep, real_t f3, real_t *soilMoisture) {
  uint64_t jj;

  assert(soil_hc);
  assert(soil_ph);
  assert(soil_ep);
  assert(soilMoisture);

  // we are not going to use this array
  memset(_IDX_SOIL, 0, _store_size*sizeof(short));

#ifdef PARA //yh@June 27th
  #pragma omp parallel for default(shared) private(jj)
#endif
  for (jj=0; jj<_store_size; jj++) {
    _HCON[jj] = soil_hc[jj]*f1;  // scaling and copying only 
    _P2[jj] = soil_ph[jj]*soil_ep[jj]*f2*f3; // product of the two, plus the scaling factor 
  }

  // initialize VSAT
  memcpy(_VSAT, soilMoisture, _store_size*sizeof(real_t));
}
 
// set up the outlets
void WaterShed::SetOutlets(uint64_t nout, uint64_t *xx, uint64_t *yy, real_t *slopes) {
  uint64_t jj;

  assert(xx);
  assert(yy);
  assert(slopes);
  assert(nout>0);

  _OUTLETS = (uint64_t*)malloc(sizeof(uint64_t)*nout);
  _OUT_SLOPES = (real_t*)malloc(sizeof(real_t)*nout);
  _N_OUT = nout;

  for (jj=0; jj<nout; jj++) {
    _OUTLETS[jj] = xx[jj]*_ncol+yy[jj];
    _OUT_SLOPES[jj] = slopes[jj];
  }
}

// storm drains
void WaterShed::SetStormDrain(uint64_t nout, uint64_t *xx, uint64_t *yy, real_t *offset, real_t *rates, real_t *saturation) {
  uint64_t jj;

  assert( xx );
  assert( yy );
  assert( offset );
  assert( rates );
  assert( saturation );
  
  _SD = (uint64_t*)malloc(sizeof(uint64_t)*nout);
  _threshold = (real_t*)malloc(sizeof(real_t)*nout);
  _rate = (real_t*)malloc(sizeof(real_t)*nout);
  _saturation = (real_t*)malloc(sizeof(real_t)*nout);
  _N_STORM = nout;

  for (jj=0; jj<nout; jj++) {
    _SD[jj] = xx[jj]*_ncol + yy[jj];
    _threshold[jj] = offset[jj];
    _rate[jj] = rates[jj];
    _saturation[jj] = saturation[jj];
  }
}

// intercept
void WaterShed::CompIntercept(real_t dt) {
  uint64_t jj;

#ifdef PARA //yh@June 27th
  #pragma omp parallel for default(shared) private(jj) 
#endif
  for (jj=0; jj<_store_size; jj++) {
    if ( _PRE[jj]*dt >= _RET[jj]) { 
      _PRE[jj] -= _RET[jj]/dt;
      _RET[jj] = 0.0;
    } else {
      _RET[jj] -= _PRE[jj]*dt;
      _PRE[jj] = 0.0;
    }
   }
}

// calculate overland depth
int WaterShed::CompOverlandDepth(real_t dt) {
  uint64_t jj;
  real_t dtdx2 = dt/(_gsz*_gsz);

#ifdef PARA //yh@June 27th
  #pragma omp parallel for default(shared) private(jj) 
#endif
  for (jj=0; jj<_store_size; jj++) {
    _H[jj] += _OLR[jj]*dtdx2 + _PRE[jj]*dt; // should we worry about stability?

    if ( _MASK[jj]==-1) _H[jj]=0.001;       // we take the water away at those boundary pixels

    if ( unlikely(_H[jj] < 0 )) {  // chicken
      if (! _printed_depth_underflow) {
        printf("Possible numerical instability: %10ld: %.5e out of %5e\n", jj, _H[jj], _OLR[jj]);
        _printed_depth_underflow = true;
      }
      _H[jj] = real_sqrt(REAL_EPSILON);
    }

    // find the maximal depth and store it
    _MAXH[jj]  = _MAXH[jj] > _H[jj] ? _MAXH[jj] : _H[jj];
    _VOL[jj]   = _MAXH[jj] * _gsz*_gsz;
    _INTH[jj] += _H[jj]*dt;
  }

  // Should call the infiltration routine next
  return 0;
}

int WaterShed::CompInfiltration(real_t dt) {
  uint64_t jj;
  real_t eight_dt = 8.0*dt;
  real_t two_dt = 2.0*dt;
  real_t tmpinf;

#ifdef PARA //yh@June 27th
  #pragma omp parallel for default(shared) private(jj,tmpinf) 
#endif
  for (jj=0; jj<_store_size; jj++) {
    tmpinf = _HCON[jj]*dt-2*_VSAT[jj];  // tmpinf = hcon * dt - 2*vsat
    // tmpinf = ( real_sqrt(8dt*hcon*(vsat+p2) + tmpinf*tmpinf) +tmpinf)/(2*dt);
    tmpinf = ( real_sqrt((_VSAT[jj]+_P2[jj])*_HCON[jj]*eight_dt + tmpinf*tmpinf) +tmpinf)/two_dt;

    if ( _H[jj]/dt <= tmpinf ) {
      tmpinf = _H[jj]/dt;
      _H[jj] = 0.0;
    } else {
      _H[jj] -= tmpinf*dt;
    }
    _VSAT[jj] += tmpinf*dt;

  }

  return 0;
}


// 2D diffusive routing
int WaterShed::CompDiffusiveRouting(real_t dt) {
  uint64_t cur,top,rgt;
  real_t cellsize = _gsz;   // we might need more for openMP
  const real_t chicken = -real_sqrt(REAL_EPSILON);
  const real_t cfl = 0.7*cellsize*cellsize/dt;

  memset(_OLR, 0, _store_size*sizeof(real_t));  // reset overland

  // Note: need to take care of local storage and channel routing !!!!
  real_t tmpsf,tmpn,tmph,tmpp;
  real_t OLRDIM0,OLRDIM1;
  real_t curfabs,oldfabs;
  int64_t mysign;

#ifdef PARA //yh@June 27th
  #pragma omp parallel for default(shared) private(cur,top,rgt,tmpsf,tmpn,tmph,tmpp,OLRDIM0,OLRDIM1,mysign,curfabs,oldfabs) 
#endif
  for (cur=0; cur<(_nrow)*(_ncol); cur++) {
    if (cur/_ncol < _nrow-1) {
      top = cur+_ncol;
      // initialize 
      tmpsf=0.0;
      tmpn=1.0;  // tmpn first use on denominator, avoid NaN
      tmph=0.0;
      tmpp=0.0;
      if ( _MASK[cur] && _MASK[top] ) { // stupid checks, should
        // do something smarter
        tmpsf = ( _ELE[cur]-_ELE[top]+_H[cur]-_H[top] )/cellsize + REAL_EPSILON;

        if (tmpsf>=0.0) {
          tmph = _H[cur];
          tmpn = _N[cur];
          tmpp = _STORE[cur];
        } else {
          tmph = _H[top];
          tmpn = _N[top];
          tmpp = _STORE[top];
        }

        //tmpt[cur] = _H[cur] < _H[top] ? _H[cur] : _H[top];

        if (unlikely( tmph < chicken )) {
          assert( tmph > chicken );  // this is a chicken switch, consider removing it?
        }
      }

      tmpn = real_sqrt(real_fabs(tmpsf))/tmpn;
      tmph -= tmpp;    // substract pond/lake storage

      OLRDIM0 = MY_SIGN(tmpsf) * cellsize * tmpn * tmph * cbrt( tmph * tmph );

      mysign = MY_SIGN(OLRDIM0);
      curfabs = real_fabs(OLRDIM0);
      oldfabs = real_fabs(_OLRDIM0_OLD[cur]);
      // bounding
      if ( oldfabs < 1e-6 ) {
        // no-op
      } else if ( mysign != MY_SIGN(_OLRDIM0_OLD[cur]) ) {
        if ( curfabs > oldfabs ) {
          OLRDIM0 = mysign*oldfabs;
          curfabs = real_fabs(OLRDIM0);
        }
      } else {
        if ( curfabs > 10*oldfabs ) {
          OLRDIM0 = mysign*10*oldfabs;
          curfabs = real_fabs(OLRDIM0);
        }
      }

      if ( curfabs > cfl ) {
        OLRDIM0 = mysign*cfl;
      }

      _OLRDIM0_OLD[cur] = OLRDIM0;
      _OLR[cur] -= OLRDIM0; //combine the routing in x
      _OLR[top] += OLRDIM0;
    }

    if (cur%_ncol < _ncol-1) {
      rgt = cur+1;
      // initialize 
      tmpsf=0.0;
      tmpn=1.0;
      tmph=0.0;
      tmpp=0.0;
      if ( _MASK[cur] && _MASK[rgt] ) { // stupid checks, should do more efficient

        tmpsf = ( _ELE[cur]-_ELE[rgt]+_H[cur]-_H[rgt])/cellsize + REAL_EPSILON;
        if (tmpsf>=0.0) {
          tmph = _H[cur];
          tmpn = _N[cur];
          tmpp = _STORE[cur];
        } else {
          tmph = _H[rgt];
          tmpn = _N[rgt];
          tmpp = _STORE[rgt];
        }

        //tmpt[cur] = _H[cur] < _H[rgt] ? _H[cur] : _H[rgt];
        if (unlikely( tmph < chicken )) {
          assert(tmph > chicken); // chicken switch
        }
      }

      tmpn = real_sqrt(real_fabs(tmpsf))/tmpn;
      tmph -= tmpp;

      OLRDIM1 = MY_SIGN(tmpsf) * cellsize * tmpn * tmph * cbrt( tmph * tmph );

      mysign = MY_SIGN(OLRDIM1);
      curfabs = real_fabs(OLRDIM1);
      oldfabs = real_fabs(_OLRDIM1_OLD[cur]);
      // bounding
      if ( real_fabs( _OLRDIM1_OLD[cur]) < 1e-6 ) {
        // no-op
      } else if ( mysign != MY_SIGN(_OLRDIM1_OLD[cur]) ) {
        if ( curfabs > oldfabs ) {
          OLRDIM1 = mysign*oldfabs;
          curfabs = real_fabs(OLRDIM1);
        }
      } else {
        if ( curfabs > 10*oldfabs ) {
          OLRDIM1 = mysign*10*oldfabs;
          curfabs = real_fabs(OLRDIM1);
        }
      }

      if ( curfabs > cfl ) {
        OLRDIM1 = mysign*cfl;
      }

      _OLRDIM1_OLD[cur] = OLRDIM1;
      _OLR[cur] -= OLRDIM1; //combine the routing in y
      _OLR[rgt] += OLRDIM1; 
    }
  }

  return 0;
}

// outlet flow, this is will not work well in openMP since we are not expecting many outlets
int WaterShed::CompOutlet(real_t dt) {
  uint64_t jj;
  uint64_t idx;
  real_t qout;
  real_t dtdx2 = dt/(_gsz*_gsz);
  real_t tt;

  if ( _N_OUT <= 0 && _N_STORM <= 0) return 0;

  // we only have a few of them
  for (jj=0; jj<_N_OUT; jj++) {
    idx = _OUTLETS[jj];
    if ( _H[ idx ] > 0.0 ) {
      tt = _H[idx]-_STORE[idx];
      qout = _gsz * real_sqrt( _OUT_SLOPES[jj])/_N[idx] * tt * real_cbrt(tt * tt);

      _H[ idx ] -= qout * dtdx2;

      if (unlikely( _H[idx] < 0 )) {  // check
        if (! _printed_outlet_underflow) {
          printf("too much water draw at outlet %ld\n", idx);
          _printed_outlet_underflow = true;
        }
        _H[idx]=0.0;
      }
    }
  }
  
  // go through the process for the storm drains
  // assuming there is no overlap between outlets and storm drains
  for (jj=0; jj<_N_STORM; jj++) {
    idx = _SD[jj];
    if ( _H[ idx ] > _saturation[jj] ) {
      qout = _rate[jj] * real_pow(_saturation[jj]-_threshold[jj]-_STORE[idx], EIGHT3RD);
    } else if ( _H[ idx] > _threshold[jj]) {
      qout = _rate[jj] * real_pow(_H[idx] - _threshold[jj]-_STORE[idx], EIGHT3RD);
    } else {
      qout = 0;
    }
    _H[idx ] -= qout *dtdx2;
  }

  // more bookkeeping might be needed here

  return 0;
}

/* end */
