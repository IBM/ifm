#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <libgen.h>
#include "grid.h"
#include "engine.h"
#include "rain_table.h"
#include "FileIO/netcdf_io.h"
    
bool isNumber(char *buf) {
  char *ptr = buf;
  while (*ptr) {
    if (! isdigit(*ptr) && (*ptr != '.' && *ptr != '-' && *ptr != 'e'))
      return false;
    ptr++;
  }
  return true;
}

// returns the size of the array
int precip_table::parse(string filename, Grid *dem) {
  const char lim[]=" \n\t,";  // delimitor allowed
  char buffer[1024];
  char *ptr;
  real_t *time=0, *rates=0;
  Grid **netcdfGrids=NULL;
  int sz=0, jj, field;
  FILE *F;
  NetCDF netcdf(dem->cellsize);

  F = fopen(filename.c_str(),"r");
  if (!F) {
    printf("Bummer: unable to open file \"%s\"\n", filename.c_str());
    return -1;
  }

  while ( (fgets(buffer,1023,F)) != NULL ) {
    ptr = strtok(buffer,lim);
    if ( ptr != NULL && ptr[0] != '#') sz++;
  }

  if (sz <= 0 ) {
    printf("Bummer: buggy precipitation file: no valid values found.\n");
    fclose(F);
    return -1;
  }

  rewind(F);
  time = (real_t*)malloc(sz*sizeof(real_t));
  rates = (real_t*)calloc(sz,sizeof(real_t));
  netcdfGrids = (Grid**)calloc(sz,sizeof(Grid*));

  int firstFileTime = -1;
  int multiplicationFactor = -1;

  jj=0;
  while ( (fgets(buffer,1023,F)) != NULL ) {
    ptr = strtok(buffer,lim);
    if ( !ptr || ptr[0] == '#') continue;
    field = 0;
    while ( ptr != NULL ){
      field++;
      if ( field == 1 ) time[jj]=atof(ptr);
      else if ( field == 2 && isNumber(ptr)) rates[jj]=atof(ptr);
      else if ( field == 2) {
        // Guess multiplicationFactor based on file names
        // Name example: wrfout_d03_2013-09-10_17:30:00
        string wrfname = basename(ptr);
        int namelen = wrfname.size();
        if (namelen == 30 && wrfname.find("wrfout") == 0) {
          int min = ((wrfname[namelen-5] - '0') * 10) + wrfname[namelen-4] - '0';
          int hour = ((wrfname[namelen-8] - '0') * 10) + wrfname[namelen-7] - '0';
       //   printf("hour: %02d, min=%02d\n", hour, min);
          if (firstFileTime == -1) {
            firstFileTime = hour * 60 + min;
          } else if (multiplicationFactor == -1) {
            int diff = hour * 60 + min - firstFileTime;
            switch (diff) {
              case 10: // one wrf file every 10min
                multiplicationFactor = 6; // 10*6 = 1h
                break;
              case 60: // one wrf file every hour
                multiplicationFactor = 1; // we already have mm/h
                break;
              default:
                fprintf(stderr, "Warning: unexpected file name sequence given\n");
                multiplicationFactor = 1; // we already have mm/h
            }
            printf("Multiplication factor: %d\n", multiplicationFactor);
          }
        }
        // Parse the file
        netcdfGrids[jj] = netcdf.parseWrfout(ptr, dem);
        if (! netcdfGrids[jj]) {
          fclose(F);
          return -1;
        }
      } else printf("More than 2 columns in some rows, excessive columns ignored.\n");
      ptr = strtok(NULL,lim);
    }
    if ( jj>=1 ) {
      if ( time[jj] <= time[jj-1] ) {
        printf("Bummer: time sequence not in ascending order.\n");
        fclose(F);
        return (-1);
      }
    }
    jj++;
  }
  // The precipitation values in the NetCDF files represent the accumulation since
  // the beginning of the forecast. Therefore, we must convert them into precipitation
  // rate by subtracting netcdfGrids[i] by netcdfGrids[i-1] for i=1..jj.
  while (--jj>=0) {
    uint64_t gridSize = netcdfGrids[jj] ? netcdfGrids[jj]->cols * netcdfGrids[jj]->rows : 0;
    if (jj>0 && netcdfGrids[jj] && netcdfGrids[jj-1]) {
      string this_fname = netcdfGrids[jj]->filename;
      string this_forecast = netcdf.getAttribute(this_fname, "SIMULATION_START_DATE");
      string last_fname = netcdfGrids[jj-1]->filename;
      string last_forecast = netcdf.getAttribute(last_fname, "SIMULATION_START_DATE");
      if (this_forecast.compare(last_forecast) == 0) {
        for (uint64_t idx=0; idx<gridSize; ++idx) {
          netcdfGrids[jj]->data[idx] -= netcdfGrids[jj-1]->data[idx];
        }
      }
    }
    if (netcdfGrids[jj] && multiplicationFactor != 1) {
      for (uint64_t idx=0; idx<gridSize; ++idx)
        netcdfGrids[jj]->data[idx] *= multiplicationFactor;
    }
  }

  // The members of the netcdfGrids array are passed over to a private
  // structure in _assign(). The members allocated with netcdf.parse()
  // are deleted when that private structure is freed.
  _assign(time, rates, netcdfGrids, (uint64_t) sz);

  if (time) free(time);
  if (rates) free(rates);
  if (netcdfGrids) free(netcdfGrids);

  fclose(F);
  return sz;
}
