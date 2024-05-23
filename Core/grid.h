#ifndef __GRID_H
#define __GRID_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string>
#include "ifm_common.h"

#define GRID(raster,col,row) (raster)->data[(row) * (raster)->cols + (col)]

using namespace std;

typedef struct BoundingBox {
  int64_t startx;
  int64_t starty;
  int64_t endx;
  int64_t endy;
} BoundingBox_t;

class Grid {
  public:
    Grid() : data(NULL), cols(0), rows(0), cellsize(90), nodata(-9999), filename("") { }
    
    Grid(uint64_t cols, uint64_t rows) : cols(cols), rows(rows), cellsize(90), nodata(-9999), filename("") {
      data = (real_t *) malloc(sizeof(real_t) * cols * rows);
    }

    ~Grid()
    {
      if (data)
        free(data);
    }

    void reset(real_t value=0.0) {
      uint64_t size = cols * rows;
      for (uint64_t i=0; i<size; ++i)
        data[i] = value;
    }

    void dump(string name)
    {
      printf("%s:\n", name.c_str());
      for (uint64_t jj=0; jj<rows; jj++) {
        for (uint64_t kk=0; kk<cols; kk++)
          printf("%.9f ", data[jj*cols+kk]);
        printf("\n");
      }
    }

    void saveAsInt(string filename)
    {
      FILE *F = fopen(filename.c_str(), "w");
      for (uint64_t jj=0; jj<rows; jj++) {
        for (uint64_t kk=0; kk<cols; kk++) {
          fprintf(F, " %1d", (int) data[jj*cols+kk]);
        }
        fprintf(F,"\n");
      }
      fclose(F);
    }

    void saveAsDouble(string filename)
    {
      FILE *F = fopen(filename.c_str(), "w");
      for (uint64_t jj=0; jj<rows; jj++) {
        for (uint64_t kk=0; kk<cols; kk++) {
          fprintf(F, " %.3f", data[jj*cols+kk]);
        }
        fprintf(F,"\n");
      }
      fclose(F);
    }

    void flip(void) {
      uint64_t lastIndex = rows-1;
      for (uint64_t i=0; i<cols; ++i)
        for (uint64_t j=0; j<rows/2; ++j) {
          real_t tmp = GRID(this, i, j);
          GRID(this, i, j) = GRID(this, i, lastIndex-j);
          GRID(this, i, lastIndex-j) = tmp;
        }
    }

    real_t  *data;
    uint64_t cols;
    uint64_t rows;
    int cellsize;
    int  nodata;
    string   filename;
};

#endif
