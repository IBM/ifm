/* define the precipitation table */

#ifndef _PRECIP_TABLH_H
#define _PRECIP_TABLH_H

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <vector>
#include "ifm_common.h"

// simple design, not intended to hold a lot of data
class precip_table {
  private:
    uint64_t _size;        // storage size
    uint64_t _num;         // number of elements
    uint64_t _current;
    real_t  *_time;
    real_t  *_rate;
    Grid   **_netcdfGrids;

    void _grow(uint64_t nsz) {
      if ( nsz <= _size) return;
      if (_time) free (_time);
      if (_rate) free (_rate);
      if (_netcdfGrids) {
        for (uint64_t i=0; i<_size; ++i) if (_netcdfGrids[i]) delete _netcdfGrids[i];
        free (_netcdfGrids);
      }
      _size = nsz;
      _time = (real_t*)malloc(_size*sizeof(real_t));
      _rate = (real_t*)malloc(_size*sizeof(real_t));
      _netcdfGrids = (Grid**)calloc(_size,sizeof(Grid*));
    }

    // assuming values in t are in ascending order
    void _assign(real_t *t, real_t *r, Grid **netcdfGrids, uint64_t num) {
      _grow(num);
      _num = num;
      memcpy(_time, t, num*sizeof(real_t));
      memcpy(_rate, r, num*sizeof(real_t));
      memcpy(_netcdfGrids, netcdfGrids, num*sizeof(Grid*));
    }

    uint64_t _getIndex(real_t t) {
      uint64_t jj;
      if ( t <= _time[0] ) return 0;
      else if ( t>=_time[_num-1] ) return _num-1;
      else if ( t >= _time[_current] && t < _time[_current+1] ) return _current;
      else {
        for (jj=1; jj<_num; jj++) if ( t >= _time[jj-1] && t < _time[jj] ) break;
        _current = jj-1;
        return _current;
      }
    }

    float getMultiplicationFactor(string wrfname1, string wrfname2);

  public:
    precip_table(uint64_t sz=127):_current(0),_time(NULL),_rate(NULL),_netcdfGrids(NULL) {
      assert(sz>0);
      _size = sz;
      _num = 0;
      _time = (real_t*)malloc(sizeof(real_t)*_size);
      _rate = (real_t*)malloc(sizeof(real_t)*_size);
      _netcdfGrids = (Grid**)calloc(_size,sizeof(Grid*));
    }

    ~precip_table() {
      if (_time) free (_time);
      if (_rate) free (_rate);
      if (_netcdfGrids) {
        for (uint64_t i=0; i<_size; ++i) if (_netcdfGrids[i]) delete _netcdfGrids[i];
        free(_netcdfGrids);
      }
    }

    bool isUniform() {
      return (_netcdfGrids && _netcdfGrids[0] == NULL);
    }

    // returns the current reading, using piece-wise-linear assumption
    // we assume the access pattern is sequential, from _time[0] to _time[end]
    real_t Get(real_t t) {
      return _rate[_getIndex(t)];
    }

    Grid *Get_Grid(real_t t) {
      uint64_t index = _getIndex(t);
      assert(_netcdfGrids[index]);
      return _netcdfGrids[index];
    }

    int parse(string filename, Grid *dem);
};

#endif

/* end */
