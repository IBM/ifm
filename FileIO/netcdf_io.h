#ifndef _IFM_NETCDF_H
#define _IFM_NETCDF_H

#include <list>
#include <string>
#include <stdio.h>
#include <netcdfcpp.h>
#include <assert.h>
#include "ifm_common.h"
#include "grid.h"

using namespace std;

class NetCDF {
  public:
    NetCDF(int _cellsize) : cellsize(_cellsize), m_dumpFileInfo(true) { }
    ~NetCDF() { }

    string getAttribute(string filename, string attrname);
    Grid *parse(string filename, string varname="", BoundingBox_t *cropArea=NULL);
    Grid *parseWrfout(string filename, Grid *cropArea=NULL);
    bool  writeHistogram(string filename, Grid *grid, string varname);
    bool  write(string filename, Grid *grid, string varname, short *mask=NULL);
    bool  write(string filename, list<Grid*> gridList, list<string> varnames);
    bool  copyVar(string fromFilename, string varname, string toFilename);

    int   cellsize;

  private:
    bool  m_dumpFileInfo;

    Grid *parseVar(NcFile *nc, string filename, string varname, BoundingBox_t *cropArea=NULL, bool autoFlip=true);
    bool  guessVarnames(NcFile *nc, string *precipitation, string *lat, string *lon);
    bool  hasVar(NcFile *nc, string varname);
    bool  writeGrid(NcFile *nc, Grid *grid, NcVar *data, int timeIndex, bool verticalFlip=false);
    Grid *resize(Grid *precipitation, Grid *lat, Grid *lon, Grid *cropArea);
    NcFile *getHandle(string filename, bool *is_new_file);

    inline bool needsFlip(NcFile *nc, NcVar *var) {
      string varname = var->name();
      if (this->hasVar(nc, "lat")) {
        if (varname.compare("lon") == 0) {
          // Never flip the uint64_titude
          return false;

        } else if (varname.compare("lat") == 0) {
          // Cannot read-ahead..
          return false;

        } else {
          // Sniff the latitude grid to figure if it's upside down or not
          Grid *lat = this->parseVar(nc, "dummy", "lat", NULL, false);
          if (lat) {
            bool upsideDown = (GRID(lat, 0, 0) < GRID(lat, 0, lat->rows-1));
            delete lat;
            return upsideDown;
          } else {
            // Flip this variable
            fprintf(stderr, "Failed to sniff the 'lat' variable, assuming %s is upside down\n", varname.c_str());
            return true;
          }
        }
      // } else if (this->hasVar(nc, "XLAT")) {
      //   if (varname.compare("XLONG") == 0) {
      //     // Never flip the uint64_titude
      //     return false;
      //
      //   } else if (varname.compare("XLAT") == 0) {
      //     // Cannot read-ahead..
      //     return false;
      //
      //   } else {
      //     // Sniff the latitude grid to figure if it's upside down or not
      //     Grid *lat = this->parseVar(nc, "dummy", "XLAT", NULL, false);
      //     if (lat) {
      //       bool upsideDown = (GRID(lat, 0, 0) < GRID(lat, 0, lat->rows-1));
      //       delete lat;
      //       return upsideDown;
      //     } else {
      //       // Flip this variable
      //       fprintf(stderr, "Failed to sniff the 'XLAT' variable, assuming %s is upside down\n", varname.c_str());
      //       return true;
      //     }
      //   }
      } else {
        // When in doubt, flip the variable
        return true;
      }
    }

    inline real_t *alloc(uint64_t &cols, uint64_t &rows, BoundingBox_t *cropArea=NULL) {
      if (cropArea) {
        uint64_t startx = cropArea ? cropArea->startx : 0;
        uint64_t starty = cropArea ? cropArea->starty : 0;
        uint64_t endx   = cropArea ? cropArea->endx+1 : cols;
        uint64_t endy   = cropArea ? cropArea->endy+1 : rows;
        cols = endx - startx;
        rows = endy - starty;
      }
      return (real_t *) malloc(sizeof(real_t) * cols * rows);
    }

    inline bool read(NcFile *nc, NcVar *var, real_t *data, uint64_t a, BoundingBox_t *cropArea=NULL, bool autoFlip=true) {
      if (cropArea) {
        uint64_t startx = cropArea ? cropArea->startx : 0;
        uint64_t starty = cropArea ? cropArea->starty : 0;
        uint64_t offset = startx > starty ? startx : starty;
        uint64_t swappedOffset = cropArea->endx > cropArea->endy ? cropArea->endx-offset : cropArea->endy-offset;
        var->set_cur(swappedOffset);
      }
      var->get(data, a);
      if (autoFlip && string(var->name()).compare("lat") == 0) {
        bool upsideDown = data[0] < data[a-1];
        if (upsideDown) {
          // Flip
          for (uint64_t i=0; i<a/2; ++i) {
            real_t old = data[a-i-1];
            data[a-i-1] = data[i];
            data[i] = old;
          }
        }
      }
      return true;
    }

    inline bool read(NcFile *nc, NcVar *var, real_t *data, uint64_t a, uint64_t b, BoundingBox_t *cropArea=NULL, bool autoFlip=true) {
      uint64_t endCol   = cropArea ? cropArea->endx+1 : b;
      uint64_t endRow   = cropArea ? cropArea->endy+1 : a;
      uint64_t startCol = cropArea ? cropArea->startx : 0;
      uint64_t startRow = cropArea ? cropArea->starty : 0;
      bool flip = autoFlip ? this->needsFlip(nc, var) : false;
      if (flip) {
        // Swap lines thanks to NetCDF's damn BottomUp feature.
        endRow = cropArea ? a - cropArea->starty : a;
        startRow = cropArea ? endRow - (cropArea->endy+1 - cropArea->starty) : 0;
      }
      for (uint64_t i=startRow; i<endRow; i++) {
        NcBool ok = flip ? var->set_cur((endRow-1-i)+startRow, startCol) : var->set_cur(i, startCol);
        if (! ok) {
          fprintf(stderr, "[a,b] NetCDF::read: failed to seek\n");
          fprintf(stderr, "Original rows:    %ld .. %ld\n", cropArea ? cropArea->starty : 0, cropArea ? cropArea->endy+1 : (int) a);
          fprintf(stderr, "Transformed rows: %ld .. %ld\n", startRow, endRow);
          fprintf(stderr, "Seek offset: (row=%ld, col=%ld)\n", flip ? (endRow-1-i)+startRow : i, startCol);
          return false;
        }
        real_t *ptr = &data[(i-startRow) * (endCol-startCol)];
        if (! var->get(ptr, 1, endCol-startCol)) {
          fprintf(stderr, "NetCDF::read: failed to read from offset %ld after seek to (%ld, %ld)\n",
              (i-startRow) * (endCol-startCol), flip ? (endRow-1-i)+startRow : i, startCol);
          if (cropArea)
            fprintf(stderr, "NetCDF::read: cropArea: x=%ld..%ld y=%ld..%ld\n",
                cropArea->startx, cropArea->endx, cropArea->starty, cropArea->endy);
          return false;
        }
      }
      return true;
    }

    inline bool read(NcFile *nc, NcVar *var, real_t *data, uint64_t a, uint64_t b, uint64_t c, BoundingBox_t *cropArea=NULL, bool autoFlip=true) {
      uint64_t endCol   = cropArea ? cropArea->endx+1 : c;
      uint64_t endRow   = cropArea ? cropArea->endy+1 : b;
      uint64_t startCol = cropArea ? cropArea->startx : 0;
      uint64_t startRow = cropArea ? cropArea->starty : 0;
      bool flip = autoFlip ? this->needsFlip(nc, var) : false;
      if (flip) {
        // Swap lines thanks to NetCDF's damn BottomUp feature.
        endRow = cropArea ? b - cropArea->starty : b;
        startRow = cropArea ? endRow - (cropArea->endy+1 - cropArea->starty) : 0;
      }
      for (uint64_t i=startRow; i<endRow; i++) {
        NcBool ok = flip ? var->set_cur(a-1, (endRow-1-i)+startRow, startCol) : var->set_cur(a-1, i, startCol);
        if (! ok) {
          fprintf(stderr, "[a,b] NetCDF::read: failed to seek\n");
          fprintf(stderr, "Original rows:    %ld .. %ld\n", cropArea ? cropArea->starty : 0, cropArea ? cropArea->endy+1 : (int) b);
          fprintf(stderr, "Transformed rows: %ld .. %ld\n", startRow, endRow);
          fprintf(stderr, "Seek offset: (row=%ld, col=%ld)\n", flip ? (endRow-1-i)+startRow : i, startCol);
          return false;
        }
        real_t *ptr = &data[(i-startRow) * (endCol-startCol)];
        if (! var->get(ptr, 1, 1, endCol-startCol)) {
          fprintf(stderr, "NetCDF::read: failed to read from offset %ld after seek to (%ld, %ld, %ld)\n",
              (i-startRow) * (endCol-startCol), a-1, flip ? (endRow-1-i)+startRow : i, startCol);
          if (cropArea)
            fprintf(stderr, "NetCDF::read: cropArea: x=%ld..%ld y=%ld..%ld\n",
                cropArea->startx, cropArea->endx, cropArea->starty, cropArea->endy);
          return false;
        }
      }
      return true;
    }
};

#endif /* _IFM_NETCDF_H */
