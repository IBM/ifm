/*
 * @brief   Defines methods for watershed.
 * @note    StarPU powers shared-memory parallelism.
 *
 * @author <main author>
 * @email  <main author's email>
 * @author  Maksims Abalenkovs
 * @email   maksims.abalenkovs@stfc.ac.uk
 * @date    Oct 8, 2024
 * @version 1.4
 */

#include <assert.h>
#include <float.h>
#include <math.h>
#include <starpu.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

    // calculate no. of blocks
    nt = _store_size/NB;

    // register data arrays with StarPU (intercept)
    starpu_vector_data_register(&pre_h, 0, (uintptr_t)_PRE, _store_size, sizeof(_PRE[0]));
    starpu_vector_data_register(&ret_h, 0, (uintptr_t)_RET, _store_size, sizeof(_RET[0]));

    // register data arrays with StarPU (overland depth)
    starpu_vector_data_register(&h_h,    0, (uintptr_t)_H,    _store_size, sizeof(_H[0]));
    starpu_vector_data_register(&olr_h,  0, (uintptr_t)_OLR,  _store_size, sizeof(_OLR[0]));
    starpu_vector_data_register(&mask_h, 0, (uintptr_t)_MASK, _store_size, sizeof(_MASK[0]));
    starpu_vector_data_register(&maxh_h, 0, (uintptr_t)_MAXH, _store_size, sizeof(_MAXH[0]));
    starpu_vector_data_register(&vol_h,  0, (uintptr_t)_VOL,  _store_size, sizeof(_VOL[0]));
    starpu_vector_data_register(&inth_h, 0, (uintptr_t)_INTH, _store_size, sizeof(_INTH[0]));

    // create block filter
    block_filter = {
        .filter_func = starpu_vector_filter_block,
        .nchildren   = nt,
    };

    // partition data arrays (intercept)
    starpu_data_partition(pre_h, &block_filter);
    starpu_data_partition(ret_h, &block_filter);

    // partition data arrays (overland routing)
    starpu_data_partition(h_h,    &block_filter);
    starpu_data_partition(olr_h,  &block_filter);
    starpu_data_partition(mask_h, &block_filter);
    starpu_data_partition(maxh_h, &block_filter);
    starpu_data_partition(vol_h,  &block_filter);
    starpu_data_partition(inth_h, &block_filter);
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
  // @todo (omp->xpu) parallelise with StarPU
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
  // @todo (omp->xpu) parallelise with StarPU
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

  // @todo (omp->xpu) parallelise with StarPU
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

// @brief StarPU intercept kernel function for CPU execution
// @param[in] buffers pointer to array of StarPU vector interfaces
// @param[in] cl_args array of StarPU inline arguments
void intercept_cpu_func(void *buffers[], void *cl_args) {

    printf("Entering 'intercept_cpu_func'...\n");

    // retrive precipitation and retention vector handles
    struct starpu_vector_interface *pre_h = (starpu_vector_interface*) buffers[0];
    struct starpu_vector_interface *ret_h = (starpu_vector_interface*) buffers[1];

    // obtain no. of elements and base pointers
    int _store_size = STARPU_VECTOR_GET_NX(pre_h);
    real_t *_PRE = (real_t*) STARPU_VECTOR_GET_PTR(pre_h);
    real_t *_RET = (real_t*) STARPU_VECTOR_GET_PTR(ret_h);

    // obtain time increment as inline argument
    real_t dt;
    starpu_codelet_unpack_args(cl_args, &dt);
    
    // kernel body
    for (int jj = 0; jj < _store_size; jj++) {

        if (_PRE[jj]*dt >= _RET[jj]) {
            _PRE[jj] -= _RET[jj]/dt;
            _RET[jj]  = 0.0;
        }
        else {
            _PRE[jj]  = 0.0;
            _RET[jj] -= _PRE[jj]*dt;
        }
    }

    printf("Exiting 'intercept_cpu_func'...\n");
}

// StarPU codelet for computing intercept
struct starpu_codelet intercept_cl {
    .cpu_func = {intercept_cpu_func},
    .nbuffers = 2,
    .modes    = {STARPU_RW, STARPU_RW},
};

// @brief Computes intercept
// @param[in] nb StarPU block size
// @param[in] dt time increment
void WaterShed::comp_intercept_starpu(real_t dt) {

    printf("Entering 'WaterShed::comp_intercept_starpu'...\n");

    /*
    // no. of StarPU blocks
    uint32_t nt = _store_size/nb;

    // define StarPU handles for precipitation and retention data arrays
    // starpu_data_handle_t pre_h, ret_h;

    // register data arrays with StarPU
    // starpu_vector_data_register(&pre_h, 0, (uintptr_t)_PRE, _store_size, sizeof(_PRE[0]));
    // starpu_vector_data_register(&ret_h, 0, (uintptr_t)_RET, _store_size, sizeof(_RET[0]));

    // divide data arrays into blocks
    struct starpu_data_filter block_filter = {
        .filter_func = starpu_vector_filter_block,
        .nchildren   = nt,
    };

    starpu_data_partition(pre_h, &block_filter);
    starpu_data_partition(ret_h, &block_filter);
    */

    // for each block do: submit StarPU tasks non-blockingly
    int status = 0;

    for (int b = 0; b < nt; b++) {

        // obtain handles for blocks
        starpu_data_handle_t pre_nb_h = starpu_data_get_sub_data(pre_h, 1, b);
        starpu_data_handle_t ret_nb_h = starpu_data_get_sub_data(ret_h, 1, b);

        printf("Submitting StarPU task for block %d...\n", b);

        // submit StarPU task
        status = starpu_task_insert(
            &intercept_cl,
            STARPU_RW,     pre_nb_h,
            STARPU_RW,     ret_nb_h,
            STARPU_VALUE, &dt, sizeof(dt),
            0);

        STARPU_CHECK_RETURN_VALUE(status, "starpu_task_insert");

        printf("StarPU task for block %d was submitted successfully.\n", b);
    }

    // wait for all tasks submitted so far
    starpu_task_wait_for_all();

    /*
    // unpartition data
    starpu_data_unpartition(pre_h, 0);
    starpu_data_unpartition(ret_h, 0);

    // unregister data arrays
    starpu_data_unregister(pre_h);
    starpu_data_unregister(ret_h);
    */

    printf("Exiting 'WaterShed::comp_intercept_starpu'...\n");
}

// Computes intercept
void WaterShed::CompIntercept(real_t dt) {
  uint64_t jj;

  for (jj=0; jj<_store_size; jj++) {
    if (_PRE[jj]*dt >= _RET[jj]) {
      _PRE[jj] -= _RET[jj]/dt;
      _RET[jj] = 0.0;
    } else {
      _PRE[jj] = 0.0;
      _RET[jj] -= _PRE[jj]*dt;
    }
   }
}

// @brief StarPU overland depth kernel function for CPU execution
// @param[in] buffers pointer to array of StarPU vector interfaces
// @param[in] cl_args array of StarPU inline arguments
void overland_depth_cpu_func(void *buffers[], void *cl_args) {

    printf("Entering 'overland_depth_cpu_func'...\n");

    // retrieve vector handles
    struct starpu_vector_interface *h_h    = (starpu_vector_interface*) buffers[0];
    struct starpu_vector_interface *olr_h  = (starpu_vector_interface*) buffers[1];
    struct starpu_vector_interface *pre_h  = (starpu_vector_interface*) buffers[2];
    struct starpu_vector_interface *mask_h = (starpu_vector_interface*) buffers[3];
    struct starpu_vector_interface *maxh_h = (starpu_vector_interface*) buffers[4];
    struct starpu_vector_interface *vol_h  = (starpu_vector_interface*) buffers[5];
    struct starpu_vector_interface *inth_h = (starpu_vector_interface*) buffers[6];

    // obtain no. of elements and base pointers
    int _store_size = STARPU_VECTOR_GET_NX(h_h);
    real_t *_H    = (real_t*) STARPU_VECTOR_GET_PTR(h_h);
    real_t *_OLR  = (real_t*) STARPU_VECTOR_GET_PTR(olr_h);
    real_t *_PRE  = (real_t*) STARPU_VECTOR_GET_PTR(pre_h);
    real_t *_MASK = (real_t*) STARPU_VECTOR_GET_PTR(mask_h);
    real_t *_MAXH = (real_t*) STARPU_VECTOR_GET_PTR(maxh_h);
    real_t *_VOL  = (real_t*) STARPU_VECTOR_GET_PTR(vol_h);
    real_t *_INTH = (real_t*) STARPU_VECTOR_GET_PTR(inth_h);

    // obtain inline arguments
    bool   printed_depth_underflow;
    real_t dt, dtdx2, gsz2;
    starpu_codelet_unpack_args(cl_args, &printed_depth_underflow, &dt, &dtdx2, &gsz2);

    // kernel body
    // @todo call original kernel function 'WaterShed::CompOverlandDepth' to
    //       avoid extra mistakes
    for (int jj = 0; jj < _store_size; jj++) {

        _H[jj] += _OLR[jj]*dtdx2 + _PRE[jj]*dt;

        // remove water @ boundary cells
        if (_MASK[jj] == -1) _H[jj] = 0.001;

        // check for underflow
        if (unlikely(_H[jj] < 0)) {
            if (!printed_depth_underflow) {
                printf("Possible numerical instability: %10ld: %.5e out of %5e\n", jj, _H[jj], _OLR[jj]);
                printed_depth_underflow = true;
            }
            _H[jj] = real_sqrt(REAL_EPSILON);
        }

        // find maximum depth, store it
        _MAXH[jj]  = _MAXH[jj] > _H[jj] ? _MAXH[jj] : _H[jj];
        _VOL[jj]   = _MAXH[jj] * gsz2;
        _INTH[jj] += _H[jj] * dt;
    }

    printf("Exiting 'overland_depth_cpu_func'...\n");
}

// StarPU codelet for computing overland depth
struct starpu_codelet overland_depth_cl {
    .cpu_func = {overland_depth_cpu_func},
    .nbuffers = 7,
    .modes    = {STARPU_RW, STARPU_R, STARPU_R, STARPU_R, STARPU_RW, STARPU_W, STARPU_RW},
};

// Computes overland depth
int WaterShed::comp_overland_depth_starpu(real_t dt) {

    printf("Entering 'WaterShed::comp_overland_depth_starpu'...\n");

    /*
    // define StarPU data handles for depth, overland routing, precipitation, mask, volume and integral of depth data arrays
    starpu_data_handle_t h_h, olr_h, pre_h, mask_h, maxh_h, vol_h, inth_h;

    // register arrays with StarPU
    starpu_vector_data_register(&h_h,    0, (uintptr_t)_H,    _store_size, sizeof(_H[0]));
    starpu_vector_data_register(&olr_h,  0, (uintptr_t)_OLR,  _store_size, sizeof(_OLR[0]));
    starpu_vector_data_register(&pre_h,  0, (uintptr_t)_PRE,  _store_size, sizeof(_PRE[0]));
    starpu_vector_data_register(&mask_h, 0, (uintptr_t)_MASK, _store_size, sizeof(_MASK[0]));
    starpu_vector_data_register(&maxh_h, 0, (uintptr_t)_MAXH, _store_size, sizeof(_MAXH[0]));
    starpu_vector_data_register(&vol_h,  0, (uintptr_t)_VOL,  _store_size, sizeof(_VOL[0]));
    starpu_vector_data_register(&inth_h, 0, (uintptr_t)_INTH,  _store_size, sizeof(_INTH[0]));

    // divide arrays into blocks
    struct starpu_data_filter block_filter = {
        .filter_func = starpu_vector_filter_block,
        .nchildren   = nt,
    };

    starpu_data_partition(h_h,    &block_filter);
    starpu_data_partition(olr_h,  &block_filter);
    starpu_data_partition(pre_h,  &block_filter);
    starpu_data_partition(mask_h, &block_filter);
    starpu_data_partition(maxh_h, &block_filter);
    starpu_data_partition(vol_h,  &block_filter);
    starpu_data_partition(inth_h, &block_filter);
    */

    real_t gsz2  = _gsz * _gsz;
    real_t dtdx2 = dt / gsz2;

    // for each block do: submit StarPU tasks non-blockingly
    int status = 0;

    for (int b = 0; b < nt; b++) {

        // obtain handles for blocks
        starpu_data_handle_t h_nb_h    = starpu_data_get_sub_data(h_h,    1, b);
        starpu_data_handle_t olr_nb_h  = starpu_data_get_sub_data(olr_h,  1, b);
        starpu_data_handle_t pre_nb_h  = starpu_data_get_sub_data(pre_h,  1, b);
        starpu_data_handle_t mask_nb_h = starpu_data_get_sub_data(mask_h, 1, b);
        starpu_data_handle_t maxh_nb_h = starpu_data_get_sub_data(maxh_h, 1, b);
        starpu_data_handle_t vol_nb_h  = starpu_data_get_sub_data(vol_h,  1, b);
        starpu_data_handle_t inth_nb_h = starpu_data_get_sub_data(inth_h, 1, b);

        printf("Submitting StarPU task for block %d...\n", b);

        // submit StarPU task
        status = starpu_task_insert(
            &overland_depth_cl,
            STARPU_RW, h_nb_h,
            STARPU_R,  olr_nb_h,
            STARPU_R,  pre_nb_h,
            STARPU_R,  mask_nb_h,
            STARPU_RW, maxh_nb_h,
            STARPU_W,  vol_nb_h,
            STARPU_RW, inth_nb_h,
            STARPU_VALUE, &_printed_depth_underflow, sizeof(_printed_depth_underflow),
            STARPU_VALUE, &dt,    sizeof(dt),
            STARPU_VALUE, &dtdx2, sizeof(dtdx2),
            STARPU_VALUE, &gsz2,  sizeof(gsz2),
            0);

        STARPU_CHECK_RETURN_VALUE(status, "starpu_task_insert");

        printf("StarPU task for block %d was submitted successfully.\n", b);
    }

    // wait for all tasks submitted so far
    starpu_task_wait_for_all();

    /*
    // unpartition data
    starpu_data_unpartition(h_h,    0);
    starpu_data_unpartition(olr_h,  0);
    starpu_data_unpartition(pre_h,  0);
    starpu_data_unpartition(mask_h, 0);
    starpu_data_unpartition(maxh_h, 0);
    starpu_data_unpartition(vol_h,  0);
    starpu_data_unpartition(inth_h, 0);

    // unregister data arrays
    starpu_data_unregister(h_h);
    starpu_data_unregister(olr_h);
    starpu_data_unregister(pre_h);
    starpu_data_unregister(mask_h);
    starpu_data_unregister(maxh_h);
    starpu_data_unregister(vol_h);
    starpu_data_unregister(inth_h);
    */

    printf("Exiting 'WaterShed::comp_overland_depth_starpu'...\n");

    return 0;
}

// Computes overland depth
int WaterShed::CompOverlandDepth(real_t dt) {
  uint64_t jj;
  real_t dtdx2 = dt/(_gsz*_gsz);

  // @todo (omp->xpu) parallelise with StarPU
  for (jj=0; jj<_store_size; jj++) {
    _H[jj] += _OLR[jj]*dtdx2 + _PRE[jj]*dt; // should we worry about stability?

    if (_MASK[jj]==-1) _H[jj]=0.001;       // we take the water away at those boundary pixels

    if (unlikely(_H[jj] < 0 )) {  // chicken
      if (!_printed_depth_underflow) {
        printf("Possible numerical instability: %10ld: %.5e out of %5e\n", jj, _H[jj], _OLR[jj]);
        _printed_depth_underflow = true;
      }
      _H[jj] = real_sqrt(REAL_EPSILON);
    }

    // find the maximal depth and store it
    _MAXH[jj]  = _MAXH[jj] > _H[jj] ? _MAXH[jj] : _H[jj];
    _VOL[jj]   = _MAXH[jj] * _gsz * _gsz;
    _INTH[jj] += _H[jj] * dt;
  }

  // Should call the infiltration routine next
  return 0;
}

// StarPU infiltration kernel function for CPU execution
void infiltrate_cpu_func(void *buffers[], void *cl_args) {

    // retrive vector handles
    struct starpu_vector_interface *hcon_h = (starpu_vector_interface*) buffers[0];
    struct starpu_vector_interface *vsat_h = (starpu_vector_interface*) buffers[1];
    struct starpu_vector_interface *p2_h   = (starpu_vector_interface*) buffers[2];
    struct starpu_vector_interface *h_h    = (starpu_vector_interface*) buffers[3];

    // obtain no. of elements and base pointers
    int n = STARPU_VECTOR_GET_NX(hcon_h);
    real_t *hcon = (real_t*) STARPU_VECTOR_GET_PTR(hcon_h);
    real_t *vsat = (real_t*) STARPU_VECTOR_GET_PTR(vsat_h);
    real_t *p2   = (real_t*) STARPU_VECTOR_GET_PTR(p2_h);
    real_t *h    = (real_t*) STARPU_VECTOR_GET_PTR(h_h);

    // obtain inline arguments
    real_t dt, two_dt, eight_dt;
    starpu_codelet_unpack_args(cl_args, &dt, &two_dt, &eight_dt);

    real_t inf;

    // kernel body
    for (int i = 0; i < n; i++) {

        // inf = hcon * dt - 2*vsat
        inf = hcon[i]*dt - 2.0*vsat[i];
    
        // inf = (real_sqrt(8dt*hcon*(vsat+p2) + inf*inf) + inf)/(2*dt);
        inf = (real_sqrt((vsat[i]+p2[i])*hcon[i]*eight_dt + inf*inf) + inf)/two_dt;
    
        if (h[i]/dt <= inf) {
          inf = h[i]/dt;
          h[i]   = 0.0;
        } else {
          h[i]  -= inf*dt;
        }
        vsat[i] += inf*dt;
    }
}

// StarPU codelet for computing infiltration
struct starpu_codelet infiltrate_cl {
    .cpu_func = {infiltrate_cpu_func},
    .nbuffers = 4,
    .modes    = {STARPU_R, STARPU_RW, STARPU_R, STARPU_RW},
};

// Computes infiltration
// @note Uses StarPU for shared-memory parallelism
int WaterShed::comp_infiltration_starpu(real_t dt) {

    // Number of StarPU blocks
    // @todo Change value in future
    int const NBLOCKS = 8;

    // define StarPU data handles for conductivity, saturation volume, second term in GA model and depth data arrays
    starpu_data_handle_t hcon_h, vsat_h, p2_h, h_h;

    // register data arrays with StarPU
    // @todo confirm size of data arrays
    starpu_vector_data_register(&hcon_h, 0, (uintptr_t)_HCON, _store_size, sizeof(_HCON[0]));
    starpu_vector_data_register(&vsat_h, 0, (uintptr_t)_VSAT, _store_size, sizeof(_VSAT[0]));
    starpu_vector_data_register(&p2_h,   0, (uintptr_t)_P2,   _store_size, sizeof(_P2[0]));
    starpu_vector_data_register(&h_h,    0, (uintptr_t)_H,    _store_size, sizeof(_H[0]));

    // divide data arrays into blocks
    struct starpu_data_filter block_filter = {
        .filter_func = starpu_vector_filter_block,
        .nchildren   = NBLOCKS,
    };

    starpu_data_partition(hcon_h, &block_filter);
    starpu_data_partition(vsat_h, &block_filter);
    starpu_data_partition(p2_h,   &block_filter);
    starpu_data_partition(h_h,    &block_filter);

    real_t two_dt   = 2.0*dt;
    real_t eight_dt = 8.0*dt;

    // for each block do: submit StarPU tasks non-blockingly
    for (int b = 0; b < NBLOCKS; b++) {

        // obtain handles for blocks
        starpu_data_handle_t hcon_nb_h = starpu_data_get_sub_data(hcon_h, 1, b);
        starpu_data_handle_t vsat_nb_h = starpu_data_get_sub_data(vsat_h, 1, b);
        starpu_data_handle_t p2_nb_h   = starpu_data_get_sub_data(p2_h,   1, b);
        starpu_data_handle_t h_nb_h    = starpu_data_get_sub_data(h_h,    1, b);

        // submit StarPU task
        starpu_task_insert(
            &infiltrate_cl,
            STARPU_R,      hcon_nb_h,
            STARPU_RW,     vsat_nb_h,
            STARPU_R,      p2_nb_h,
            STARPU_RW,     h_nb_h,
            STARPU_VALUE, &dt,       sizeof(dt),
            STARPU_VALUE, &two_dt,   sizeof(two_dt),
            STARPU_VALUE, &eight_dt, sizeof(eight_dt),
            0);
    }

    // wait for all tasks submitted so far
    starpu_task_wait_for_all();

    // unpartition data
    starpu_data_unpartition(hcon_h, 0);
    starpu_data_unpartition(vsat_h, 0);
    starpu_data_unpartition(p2_h,   0);
    starpu_data_unpartition(h_h,    0);

    // unregister data arrays
    starpu_data_unregister(hcon_h);
    starpu_data_unregister(vsat_h);
    starpu_data_unregister(p2_h);
    starpu_data_unregister(h_h);

    return 0;
}

// Computes infiltration
int WaterShed::CompInfiltration(real_t dt) {

    printf("Entering 'WaterShed::CompInfiltration'...\n");

  uint64_t jj;
  real_t eight_dt = 8.0*dt;
  real_t two_dt = 2.0*dt;
  real_t tmpinf;

  // @todo (omp->xpu) parallelise with StarPU
  for (jj=0; jj<_store_size; jj++) {

    // tmpinf = hcon * dt - 2*vsat
    tmpinf = _HCON[jj]*dt-2*_VSAT[jj];

    // tmpinf = (real_sqrt(8dt*hcon*(vsat+p2) + tmpinf*tmpinf) + tmpinf)/(2*dt);
    tmpinf = ( real_sqrt((_VSAT[jj]+_P2[jj])*_HCON[jj]*eight_dt + tmpinf*tmpinf) +tmpinf)/two_dt;

    if ( _H[jj]/dt <= tmpinf ) {
      tmpinf = _H[jj]/dt;
      _H[jj] = 0.0;
    } else {
      _H[jj] -= tmpinf*dt;
    }
    _VSAT[jj] += tmpinf*dt;

  }

    printf("Exiting 'WaterShed::CompInfiltration'...\n");

  return 0;
}

// StarPU diffusive routing kernel function for CPU execution
// @todo Apply two dimensional filter to divide data arrays
// @todo Think how to omit calculation in boundary regions (last row, last column)
void diffusive_routing_cpu_func(void *buffers[], void *cl_args) {

    // retrieve data array vector handles
    struct starpu_vector_interface *mask_h        = (starpu_vector_interface*) buffers[0];
    struct starpu_vector_interface *ele_h         = (starpu_vector_interface*) buffers[1];
    struct starpu_vector_interface *h_h           = (starpu_vector_interface*) buffers[2];
    struct starpu_vector_interface *n_h           = (starpu_vector_interface*) buffers[3];
    struct starpu_vector_interface *store_h       = (starpu_vector_interface*) buffers[4];
    struct starpu_vector_interface *olrdim0_old_h = (starpu_vector_interface*) buffers[5];
    struct starpu_vector_interface *olrdim1_old_h = (starpu_vector_interface*) buffers[6];
    struct starpu_vector_interface *olr_h         = (starpu_vector_interface*) buffers[7];

    // obtain no. of elements and base pointers
    int n = STARPU_VECTOR_GET_NX(mask_h);
    real_t *_MASK        = (real_t*) STARPU_VECTOR_GET_PTR(mask_h);
    real_t *_ELE         = (real_t*) STARPU_VECTOR_GET_PTR(ele_h);
    real_t *_H           = (real_t*) STARPU_VECTOR_GET_PTR(h_h);
    real_t *_N           = (real_t*) STARPU_VECTOR_GET_PTR(n_h);
    real_t *_STORE       = (real_t*) STARPU_VECTOR_GET_PTR(store_h);
    real_t *_OLRDIM0_OLD = (real_t*) STARPU_VECTOR_GET_PTR(olrdim0_old_h);
    real_t *_OLRDIM1_OLD = (real_t*) STARPU_VECTOR_GET_PTR(olrdim1_old_h);
    real_t *_OLR         = (real_t*) STARPU_VECTOR_GET_PTR(olr_h);

    real_t   _gsz;
    real_t    dt;
    uint64_t _nrow;
    uint64_t _ncol;
    // cellsize
    // REAL_EPSILON

    // obtain no. of rows and columns as inline arguments
    starpu_codelet_unpack_args(cl_args, &_gsz, &dt, &_nrow, &_ncol);

    uint64_t cur, top, rgt;
    real_t cellsize      = _gsz;
    const real_t chicken = -real_sqrt(REAL_EPSILON);
    const real_t cfl     = 0.7*cellsize*cellsize/dt;

    // reset overland
    memset(_OLR, 0, n*sizeof(real_t));

    real_t  tmpsf, tmpn, tmph, tmpp;
    real_t  OLRDIM0, OLRDIM1;
    real_t  curfabs, oldfabs;
    int64_t mysign;

    // kernel body
    // @todo Rename 'top' variable into 'bottom', since it denotes _bottom_ element
    //       relative to _current_
    // @note Original loop range changed to omit processing last row of elements
    //       This eliminated need for if statement
    // for (cur = 0; cur < (_nrow)*(_ncol); cur++) {
    //     if (cur/_ncol < _nrow-1) {

    for (uint64_t i = 0; i < _nrow-1; i++) {
        for (uint64_t j = 0; j < _ncol; j++) {

            cur = i*(_ncol)+j;
            top = cur+_ncol;
    
            // initialize 
            tmpsf = 0.0;
            tmpn  = 1.0;  // tmpn first use on denominator, avoid NaN
            tmph  = 0.0;
            tmpp  = 0.0;
    
            if (_MASK[cur] && _MASK[top]) {  // stupid checks, should do something smarter
    
                tmpsf = (_ELE[cur]-_ELE[top]+_H[cur]-_H[top])/cellsize + REAL_EPSILON;
        
                if (tmpsf >= 0.0) {
                    tmph = _H[cur];
                    tmpn = _N[cur];
                    tmpp = _STORE[cur];
                }
                else {
                    tmph = _H[top];
                    tmpn = _N[top];
                    tmpp = _STORE[top];
                }
        
                // tmpt[cur] = _H[cur] < _H[top] ? _H[cur] : _H[top];
    
                // this is a chicken switch, consider removing it?
                if (unlikely(tmph < chicken)) {
                    assert(tmph > chicken);
                }
            }
      
            tmpn  = real_sqrt(real_fabs(tmpsf))/tmpn;
            tmph -= tmpp;  // substract pond/lake storage
      
            OLRDIM0 = MY_SIGN(tmpsf) * cellsize * tmpn * tmph * cbrt(tmph * tmph);
      
            mysign  = MY_SIGN(OLRDIM0);
            curfabs = real_fabs(OLRDIM0);
            oldfabs = real_fabs(_OLRDIM0_OLD[cur]);
    
            // bounding
            if (oldfabs < 1e-6) {
                // no-op
            }
            else if (mysign != MY_SIGN(_OLRDIM0_OLD[cur])) {
    
                if (curfabs > oldfabs) {
                    OLRDIM0 = mysign*oldfabs;
                }
    
            }
            else {
                if (curfabs > 10*oldfabs) {
                    OLRDIM0 = mysign*10*oldfabs;
                }
            }
      
            if (curfabs > cfl) {
                OLRDIM0 = mysign*cfl;
            }
      
            _OLRDIM0_OLD[cur] = OLRDIM0;
            _OLR[cur] -= OLRDIM0;  // combine the routing in x
            _OLR[top] += OLRDIM0;
        }
    }

    // @note Original loop range changed to omit processing last column of elements
    //       This eliminated need for if statement
    // for (cur = 0; cur < (_nrow)*(_ncol); cur++) {
    //     if (cur%_ncol < _ncol-1) {

    for (uint64_t i = 0; i < _nrow; i++) {
        for (uint64_t j = 0; j < _ncol-1; j++) {

            cur = i*(_ncol)+j;
            rgt = cur+1;
    
            // initialize 
            tmpsf = 0.0;
            tmpn  = 1.0;
            tmph  = 0.0;
            tmpp  = 0.0;
    
            // stupid checks, should do more efficient
            if (_MASK[cur] && _MASK[rgt]) {
    
                tmpsf = (_ELE[cur]-_ELE[rgt]+_H[cur]-_H[rgt])/cellsize + REAL_EPSILON;
                if (tmpsf >= 0.0) {
                    tmph = _H[cur];
                    tmpn = _N[cur];
                    tmpp = _STORE[cur];
                }
                else {
                    tmph = _H[rgt];
                    tmpn = _N[rgt];
                    tmpp = _STORE[rgt];
                }
        
                // tmpt[cur] = _H[cur] < _H[rgt] ? _H[cur] : _H[rgt];
    
                // chicken switch
                if (unlikely(tmph < chicken)) {
                    assert(tmph > chicken);
                }
            }
    
            tmpn  = real_sqrt(real_fabs(tmpsf))/tmpn;
            tmph -= tmpp;
    
            OLRDIM1 = MY_SIGN(tmpsf) * cellsize * tmpn * tmph * cbrt(tmph * tmph);
    
            mysign  = MY_SIGN(OLRDIM1);
            curfabs = real_fabs(OLRDIM1);
            oldfabs = real_fabs(_OLRDIM1_OLD[cur]);
    
            // bounding
            if (real_fabs(_OLRDIM1_OLD[cur]) < 1e-6) {
                // no-op
            }
            else if (mysign != MY_SIGN(_OLRDIM1_OLD[cur])) {
    
                if (curfabs > oldfabs) {
                    OLRDIM1 = mysign*oldfabs;
                }
            }
            else {
                if (curfabs > 10*oldfabs) {
                    OLRDIM1 = mysign*10*oldfabs;
                }
            }
    
            if (curfabs > cfl) {
                OLRDIM1 = mysign*cfl;
            }
    
            _OLRDIM1_OLD[cur] = OLRDIM1;
            _OLR[cur] -= OLRDIM1;  // combine the routing in y
            _OLR[rgt] += OLRDIM1; 
        }
    }
}

// StarPU codelet for computing diffusive routing
struct starpu_codelet diffusive_routing_cl {
    .cpu_func = {diffusive_routing_cpu_func},
    .nbuffers = 8,
    .modes    = {STARPU_R, STARPU_R, STARPU_R, STARPU_R, STARPU_R, STARPU_RW, STARPU_RW, STARPU_RW}
};

// Computes diffusive routing
// @note Uses StarPU for shared-memory parallelism
int WaterShed::comp_diffusive_routing_starpu(real_t dt) {

    // number of StarPU blocks
    // @todo change value in future
    int const NBLOCKS = 8;

    // define StarPU handles for data arrays
    starpu_data_handle_t mask_h, ele_h, h_h, n_h, store_h, olrdim0_old_h, olrdim1_old_h, olr_h;

    // register data arrays with StarPU
    // @todo check array sizes
    starpu_vector_data_register(&mask_h,        0, (uintptr_t)_MASK,        _store_size, sizeof(_MASK[0]));
    starpu_vector_data_register(&ele_h,         0, (uintptr_t)_ELE,         _store_size, sizeof(_ELE[0]));
    starpu_vector_data_register(&h_h,           0, (uintptr_t)_H,           _store_size, sizeof(_H[0]));
    starpu_vector_data_register(&n_h,           0, (uintptr_t)_N,           _store_size, sizeof(_N[0]));
    starpu_vector_data_register(&store_h,       0, (uintptr_t)_STORE,       _store_size, sizeof(_STORE[0]));
    starpu_vector_data_register(&olrdim0_old_h, 0, (uintptr_t)_OLRDIM0_OLD, _store_size, sizeof(_OLRDIM0_OLD[0]));
    starpu_vector_data_register(&olrdim1_old_h, 0, (uintptr_t)_OLRDIM1_OLD, _store_size, sizeof(_OLRDIM1_OLD[0]));
    starpu_vector_data_register(&olr_h,         0, (uintptr_t)_OLR,         _store_size, sizeof(_OLR[0]));

    // divide data arrays into blocks
    struct starpu_data_filter block_filter = {
        .filter_func = starpu_vector_filter_block,
        .nchildren   = NBLOCKS,
    };

    starpu_data_partition(mask_h,        &block_filter);
    starpu_data_partition(ele_h,         &block_filter);
    starpu_data_partition(h_h,           &block_filter);
    starpu_data_partition(n_h,           &block_filter);
    starpu_data_partition(store_h,       &block_filter);
    starpu_data_partition(olrdim0_old_h, &block_filter);
    starpu_data_partition(olrdim1_old_h, &block_filter);
    starpu_data_partition(olr_h,         &block_filter);

    // for each block do: submit StarPU tasks non-blockingly
    for (int b = 0; b < NBLOCKS; b++) {

        // obtain handles for blocks
        starpu_data_handle_t mask_nb_h, ele_nb_h, h_nb_h, n_nb_h, store_nb_h, olrdim0_old_nb_h, olrdim1_old_nb_h, olr_nb_h;

        mask_nb_h        = starpu_data_get_sub_data(mask_h,        1, b);
        ele_nb_h         = starpu_data_get_sub_data(ele_h,         1, b);
        h_nb_h           = starpu_data_get_sub_data(h_h,           1, b);
        n_nb_h           = starpu_data_get_sub_data(n_h,           1, b);
        store_nb_h       = starpu_data_get_sub_data(store_h,       1, b);
        olrdim0_old_nb_h = starpu_data_get_sub_data(olrdim0_old_h, 1, b);
        olrdim1_old_nb_h = starpu_data_get_sub_data(olrdim1_old_h, 1, b);
        olr_nb_h         = starpu_data_get_sub_data(olr_h,         1, b);

        // submit StarPU task
        starpu_task_insert(
            &diffusive_routing_cl,
            STARPU_R,  mask_nb_h,
            STARPU_R,  ele_nb_h,
            STARPU_R,  h_nb_h,
            STARPU_R,  n_nb_h,
            STARPU_R,  store_nb_h,
            STARPU_RW, olrdim0_old_nb_h,
            STARPU_RW, olrdim1_old_nb_h,
            STARPU_RW, olr_nb_h,
            STARPU_VALUE, &_gsz,  sizeof(_gsz),
            STARPU_VALUE, &dt,    sizeof(dt),
            STARPU_VALUE, &_nrow, sizeof(_nrow),
            STARPU_VALUE, &_ncol, sizeof(_ncol),
            0);
    }

    // wait for all tasks submitted so far
    starpu_task_wait_for_all();

    // unpartition data
    starpu_data_unpartition(mask_h,        0);
    starpu_data_unpartition(ele_h,         0);
    starpu_data_unpartition(h_h,           0);
    starpu_data_unpartition(n_h,           0);
    starpu_data_unpartition(store_h,       0);
    starpu_data_unpartition(olrdim0_old_h, 0);
    starpu_data_unpartition(olrdim1_old_h, 0);
    starpu_data_unpartition(olr_h,         0);

    // unregister data arrays
    starpu_data_unregister(mask_h);
    starpu_data_unregister(ele_h);
    starpu_data_unregister(h_h);
    starpu_data_unregister(n_h);
    starpu_data_unregister(store_h);
    starpu_data_unregister(olrdim0_old_h);
    starpu_data_unregister(olrdim1_old_h);
    starpu_data_unregister(olr_h);

    return 0;
}

// Computes diffusive routing
int WaterShed::CompDiffusiveRouting(real_t dt) {

    printf("Entering 'WaterShed::CompDiffusiveRouting'...\n");

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

  // @todo (omp->xpu) parallelise with StarPU
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
        }
      } else {
        if ( curfabs > 10*oldfabs ) {
          OLRDIM0 = mysign*10*oldfabs;
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
        }
      } else {
        if ( curfabs > 10*oldfabs ) {
          OLRDIM1 = mysign*10*oldfabs;
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

    printf("Exiting 'WaterShed::CompDiffusiveRouting'...\n");

  return 0;
}

// StarPU outlet kernel function for CPU execution
void outlet_cpu_func(void *buffers[], void *cl_args) {

    uint64_t jj;
    uint64_t idx;

    real_t qout;
    real_t tt;

    // retrieve vector handles
    struct starpu_vector_interface *outlets_h    = (starpu_vector_interface*) buffers[0];
    struct starpu_vector_interface *h_h          = (starpu_vector_interface*) buffers[1];
    struct starpu_vector_interface *store_h      = (starpu_vector_interface*) buffers[2];
    struct starpu_vector_interface *out_slopes_h = (starpu_vector_interface*) buffers[3];
    struct starpu_vector_interface *n_h          = (starpu_vector_interface*) buffers[4];

    // obtain no. of elements and base pointers
    int _N_OUT = STARPU_VECTOR_GET_NX(outlets_h);

    uint64_t *_OUTLETS  = (uint64_t*) STARPU_VECTOR_GET_PTR(outlets_h);     // _n_out
    real_t *_H          = (real_t*)   STARPU_VECTOR_GET_PTR(h_h);           // _store_size
    real_t *_STORE      = (real_t*)   STARPU_VECTOR_GET_PTR(store_h);       // _store_size
    real_t *_OUT_SLOPES = (real_t*)   STARPU_VECTOR_GET_PTR(out_slopes_h);  // _n_out
    real_t *_N          = (real_t*)   STARPU_VECTOR_GET_PTR(n_h);           // _store_size

    // obtain _gsz, dtdx2, _printed_outlet_underflow as inline arguments
    real_t _gsz, dtdx2;
    bool _printed_outlet_underflow;
    starpu_codelet_unpack_args(cl_args, &_gsz, &dtdx2, &_printed_outlet_underflow);

    // kernel body
    // for each outlet do
    for (jj = 0; jj < _N_OUT; jj++) {

        idx = _OUTLETS[jj];

        if (_H[idx] > 0.0) {

            tt   = _H[idx] - _STORE[idx];
            qout = _gsz * real_sqrt(_OUT_SLOPES[jj])/_N[idx] * tt * real_cbrt(tt * tt);
      
            _H[idx] -= qout * dtdx2;
      
            if (unlikely(_H[idx] < 0)) {  // check
                if (! _printed_outlet_underflow) {
                    printf("too much water draw at outlet %ld\n", idx);
                    _printed_outlet_underflow = true;
                }

                _H[idx] = 0.0;
            }
        }
    }
}

// StarPU codelet for computing outlet
struct starpu_codelet outlet_cl {
    .cpu_func = {outlet_cpu_func},
    .nbuffers = 5,
    .modes    = {STARPU_R, STARPU_RW, STARPU_R, STARPU_R, STARPU_R},
};

// Computes outlet flow
// @note Uses StarPU for shared-memory parallelism
// @note h, store, n data arrays are _not_ divided into nblocks since their block size is different from outlets and out_slopes block size
int WaterShed::comp_outlet_starpu(real_t dt) {

    if (_N_OUT > 0) {

        // number of StarPU blocks
        // @todo change value in future
        int const NBLOCKS = 1;
    
        // define StarPU handles for data arrays
        starpu_data_handle_t outlets_h, h_h, store_h, out_slopes_h, n_h;
    
        // register data arrays with StarPU
        starpu_vector_data_register(&outlets_h,    0, (uintptr_t)_OUTLETS,    _N_OUT,        sizeof(_OUTLETS[0]));
        starpu_vector_data_register(&h_h,          0, (uintptr_t)_H,          _store_size, sizeof(_H[0]));
        starpu_vector_data_register(&store_h,      0, (uintptr_t)_STORE,      _store_size, sizeof(_STORE[0]));
        starpu_vector_data_register(&out_slopes_h, 0, (uintptr_t)_OUT_SLOPES, _N_OUT,        sizeof(_OUT_SLOPES[0]));
        starpu_vector_data_register(&n_h,          0, (uintptr_t)_N,          _store_size, sizeof(_N[0]));
    
        // divide arrays into blocks
        struct starpu_data_filter block_filter = {
            .filter_func = starpu_vector_filter_block,
            .nchildren   = NBLOCKS,
        };
    
        starpu_data_partition(outlets_h,    &block_filter);
        // starpu_data_partition(h_h,          &block_filter);
        // starpu_data_partition(store_h,      &block_filter);
        starpu_data_partition(out_slopes_h, &block_filter);
        // starpu_data_partition(n_h,          &block_filter);
    
        real_t dtdx2 = dt/(_gsz*_gsz);
    
        // for each _outlet_ block do: submit StarPU tasks non-blockingly
        for (int b = 0; b < NBLOCKS; b++) {
    
            // obtain handles for _outlet_ blocks
            starpu_data_handle_t outlets_nb_h    = starpu_data_get_sub_data(outlets_h,    1, b);
            // starpu_data_handle_t h_nb_h          = starpu_data_get_sub_data(h_h,          1, b);
            // starpu_data_handle_t store_nb_h      = starpu_data_get_sub_data(store_h,      1, b);
            starpu_data_handle_t out_slopes_nb_h = starpu_data_get_sub_data(out_slopes_h, 1, b);
            // starpu_data_handle_t n_nb_h          = starpu_data_get_sub_data(n_h,          1, b);
    
            starpu_task_insert(
                &outlet_cl,
                STARPU_R,  outlets_nb_h,
                STARPU_RW, h_h,
                STARPU_R,  store_h,
                STARPU_R,  out_slopes_nb_h,
                STARPU_R,  n_h,
                STARPU_VALUE, &_gsz,  sizeof(_gsz),
                STARPU_VALUE, &dtdx2, sizeof(dtdx2),
                STARPU_VALUE, &_printed_outlet_underflow, sizeof(_printed_outlet_underflow),
                0);
        }
    
        // wait for all tasks submitted so far
        starpu_task_wait_for_all();
    
        // unpartition data
        starpu_data_unpartition(outlets_h,    0);
        // starpu_data_unpartition(h_h,          0);
        // starpu_data_unpartition(store_h,      0);
        starpu_data_unpartition(out_slopes_h, 0);
        // starpu_data_unpartition(n_h,          0);
    
        // unregister data arrays
        starpu_data_unregister(outlets_h);
        starpu_data_unregister(h_h);
        starpu_data_unregister(store_h);
        starpu_data_unregister(out_slopes_h);
        starpu_data_unregister(n_h);
    }

    return 0;
}

// StarPU storm kernel function for CPU execution
void storm_cpu_func(void *buffers[], void *cl_args) {

    uint64_t jj;
    uint64_t idx;

    real_t qout;

    // retrieve vector handles
    struct starpu_vector_interface *sd_h         = (starpu_vector_interface*) buffers[0];
    struct starpu_vector_interface *h_h          = (starpu_vector_interface*) buffers[1];
    struct starpu_vector_interface *saturation_h = (starpu_vector_interface*) buffers[2];
    struct starpu_vector_interface *rate_h       = (starpu_vector_interface*) buffers[3];
    struct starpu_vector_interface *threshold_h  = (starpu_vector_interface*) buffers[4];
    struct starpu_vector_interface *store_h      = (starpu_vector_interface*) buffers[5];

    // obtain no. of elements and base pointers
    int _N_STORM = STARPU_VECTOR_GET_NX(sd_h);

    uint64_t *_SD       = (uint64_t*) STARPU_VECTOR_GET_PTR(sd_h);          // _n_out
    real_t *_H          = (real_t*)   STARPU_VECTOR_GET_PTR(h_h);           // _store_size
    real_t *_saturation = (real_t*)   STARPU_VECTOR_GET_PTR(saturation_h);  // _store_size
    real_t *_rate       = (real_t*)   STARPU_VECTOR_GET_PTR(rate_h);        // _n_out
    real_t *_threshold  = (real_t*)   STARPU_VECTOR_GET_PTR(threshold_h);   // _n_out
    real_t *_STORE      = (real_t*)   STARPU_VECTOR_GET_PTR(store_h);       // _store_size

    // obtain dtdx2 as inline argument
    real_t dtdx2;
    starpu_codelet_unpack_args(cl_args, &dtdx2);

    // kernel body
    // for each storm do
    for (jj = 0; jj < _N_STORM; jj++) {

        idx = _SD[jj];

        if (_H[idx] > _saturation[jj]) {
            qout = _rate[jj] * real_pow(_saturation[jj] - _threshold[jj] - _STORE[idx], EIGHT3RD);
        }
        else if (_H[ idx] > _threshold[jj]) {
            qout = _rate[jj] * real_pow(_H[idx] - _threshold[jj] - _STORE[idx], EIGHT3RD);
        }
        else {
          qout = 0;
        }

        _H[idx] -= qout * dtdx2;
    }
}

// StarPU codelet for computing storm
struct starpu_codelet storm_cl {
    .cpu_func = {storm_cpu_func},
    .nbuffers = 6,
    .modes    = {STARPU_R, STARPU_RW, STARPU_R, STARPU_R, STARPU_R, STARPU_R},
};

// Computes storm flow
// @note Uses StarPU for shared-memory parallelism
// @note h, threshold, store data arrays are _not_ divided into nblocks since their block size is different from outlets and out_slopes block size
int WaterShed::comp_storm_starpu(real_t dt) {

    if (_N_STORM > 0) {

        // number of StarPU blocks
        // @todo change value in future
        int const NBLOCKS = 1;
    
        // define StarPU handles for data arrays
        starpu_data_handle_t sd_h, h_h, saturation_h, threshold_h, store_h, rate_h;
    
        // register data arrays with StarPU
        starpu_vector_data_register(&sd_h,         0, (uintptr_t)_SD,         _N_OUT,      sizeof(_SD[0]));
        starpu_vector_data_register(&h_h,          0, (uintptr_t)_H,          _store_size, sizeof(_H[0]));
        starpu_vector_data_register(&saturation_h, 0, (uintptr_t)_saturation, _N_OUT,      sizeof(_saturation[0]));
        starpu_vector_data_register(&threshold_h,  0, (uintptr_t)_threshold,  _store_size, sizeof(_threshold[0]));
        starpu_vector_data_register(&store_h,      0, (uintptr_t)_STORE,      _store_size, sizeof(_STORE[0]));
        starpu_vector_data_register(&rate_h,       0, (uintptr_t)_rate,       _N_OUT,      sizeof(_rate[0]));
    
        // divide arrays into blocks
        struct starpu_data_filter block_filter = {
            .filter_func = starpu_vector_filter_block,
            .nchildren   = NBLOCKS,
        };
    
        starpu_data_partition(sd_h,         &block_filter);
        // starpu_data_partition(h_h,          &block_filter);
        starpu_data_partition(saturation_h, &block_filter);
        // starpu_data_partition(threshold_h,  &block_filter);
        // starpu_data_partition(store_h,      &block_filter);
        starpu_data_partition(rate_h,       &block_filter);
    
        real_t dtdx2 = dt/(_gsz*_gsz);
    
        // for each _storm_ block do: submit StarPU tasks non-blockingly
        for (int b = 0; b < NBLOCKS; b++) {
    
            // obtain handles for _storm_ blocks
            starpu_data_handle_t sd_nb_h         = starpu_data_get_sub_data(sd_h,         1, b);
            starpu_data_handle_t saturation_nb_h = starpu_data_get_sub_data(saturation_h, 1, b);
            starpu_data_handle_t rate_nb_h       = starpu_data_get_sub_data(rate_h,       1, b);
    
            starpu_task_insert(
                &storm_cl,
                STARPU_R,  sd_nb_h,
                STARPU_RW, h_h,
                STARPU_R,  saturation_nb_h,
                STARPU_R,  rate_nb_h,
                STARPU_R,  threshold_h,
                STARPU_R,  store_h,
                STARPU_VALUE, &dtdx2, sizeof(dtdx2),
                0);
        }
    
        // wait for all tasks submitted so far
        starpu_task_wait_for_all();
    
        // unpartition data
        starpu_data_unpartition(sd_h,         0);
        // starpu_data_unpartition(h_h,          0);
        starpu_data_unpartition(saturation_h, 0);
        starpu_data_unpartition(rate_h,       0);
        // starpu_data_unpartition(threshold_h,  0);
        // starpu_data_unpartition(store_h,      0);
    
        // unregister data arrays
        starpu_data_unregister(sd_h);
        starpu_data_unregister(h_h);
        starpu_data_unregister(saturation_h);
        starpu_data_unregister(rate_h);
        starpu_data_unregister(threshold_h);
        starpu_data_unregister(store_h);
    }

    return 0;
}

// outlet flow, this is will not work well in openMP since we are not expecting many outlets
int WaterShed::CompOutlet(real_t dt) {

    printf("Entering 'WaterShed::CompOutlet'...\n");

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
      qout = _gsz * real_sqrt(_OUT_SLOPES[jj])/_N[idx] * tt * real_cbrt(tt * tt);

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
    _H[idx] -= qout * dtdx2;
  }

  // more bookkeeping might be needed here

    printf("Exiting 'WaterShed::CompOutlet'...\n");

  return 0;
}

// @eof watershed.C
