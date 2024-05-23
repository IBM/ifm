/* top level driver to make it work */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <list>
#include <libgen.h>

#include "watch.h"
#include "engine.h"
#include "ifm_common.h"
#include "rain_table.h"
#include "interpolation.h"
#include "option_parser.h"
#include "FileIO/mapfile.h"
#include "FileIO/netcdf_io.h"

// generate the rotating wheel
void rotating_wheel() {
  static int pos=0;
  char cursors[4]={'-','\\','|','/'};
  printf("%c\b", cursors[pos]);
  fflush(stdout);
  pos = (pos+1) % 4;
}

short *float2short(Grid *grid) {
  uint64_t len = grid->cols * grid->rows;
  short *out = (short*)malloc(sizeof(short) * len);
  for (uint64_t jj=0; jj<len; jj++)
    out[jj] = (short) grid->data[jj];
  free(grid->data);
  grid->data = NULL;
  return out;
}

int read_watch_file(Grid *latGrid, Grid *lonGrid, WatchPoint **watchList, string fname) {
  char buffer[1024];
  int w = 0, numlines = 0;
  WatchPoint *_watchList = NULL;

  FILE *F = fopen(fname.c_str(), "r");
  if (! F) {
    printf("Bummer: unable to open file \"%s\"\n", fname.c_str());
    return -1;
  }
  while ( (fgets(buffer, 1023, F)) != NULL ) {
    if (buffer[0] != '#')
      numlines++;
  }
  rewind(F);

  if (numlines) {
    _watchList = new WatchPoint[numlines];
    while ( (fgets(buffer, 1023, F)) != NULL ) {
      if (buffer[0] != '#') {
        char *copy = strdup(buffer);
        char *lat = strtok(copy, " ");
        char *lon = strtok(NULL, " ");
        char *name = strstr(buffer, ":");
        if (! lat || ! lon || ! name) {
          printf("Bummer: invalid syntax found while parsing file \"%s\"\n", fname.c_str());
          delete [] _watchList;
          free(copy);
          fclose(F);
          return -1;
        }
        while (name[0] == ' ' || name[0] == ':')
          name++;
        if (strlen(name) > 2 && name[strlen(name)-1] == '\n')
          name[strlen(name)-1] = '\0';

        Interpolation i;
        uint64_t sx=0, ex=0, sy=0, ey=0;
        if (i.getXY(latGrid, lonGrid, atof(lat), atof(lon), sx, ex, sy, ey) == true)
          _watchList[w++].init((sx+ex)/2, (sy+ey)/2, name);
        else
          printf("Warning: '%s' is outside the DEM boundaries\n", name);
        free(copy);
      }
    }
    *watchList = _watchList;
  }
  fclose(F);
  return w;
}

Grid* parse_as_double(string fname, int cellsize, bool print_error=true) {
  NetCDF netcdf(cellsize);
  Grid *grid = netcdf.parse(fname);
  if (! grid) {
    if (print_error)
      printf("Bummer: unable to parse file \"%s\"\n", fname.c_str());
    return NULL;
  }
  return grid;
}

Grid* parse_as_short(string fname, short **out, int cellsize, bool print_error=true) {
  Grid *grid = parse_as_double(fname, cellsize, print_error);
  *out = NULL;
  if (grid) {
    // TODO: right now the output will be put on an array of doubles. Here we basically
    // copy that into an array of shorts. We should update the internals of the Grid
    // class so that it stores the data in the appropriate array elegantly.
    *out = float2short(grid);
  }
  return grid;
}

void set_nodata_cells(Grid *dem, vector<real_t> nodata) {
  uint64_t jj, len = dem->rows * dem->cols;
  size_t n;
  for (jj=0; jj < len; jj++) {
    for (n=0; n<nodata.size(); n++) {
      if (FLOAT_EQ(dem->data[jj], nodata[n])) {
        dem->data[jj] = dem->nodata;
        break;
      }
    }
  }
}

/* generate mask from dem data, anywhere is not NODATA, we generate the mask */
Grid* gen_mask_from_dem(Grid *dem, short **out) {
  uint64_t jj, len;
  Grid *grid;

  grid = new Grid;
  grid->cols = dem->cols;
  grid->rows = dem->rows;
  grid->nodata = dem->nodata;
  grid->cellsize = dem->cellsize;

  len = dem->rows * dem->cols;
  *out = (short*)malloc(sizeof(short) * len);

  for (jj=0; jj < len; jj++) {
    /* valid pixel=1, NODATA=0 */
    if ( !(FLOAT_EQ(dem->data[jj], dem->nodata)) ) (*out)[jj]=1;
    else (*out)[jj]=0;
  }
  return grid;
}

/* read the outlet file, each line represents one outlet, assuming 128 max!!
   returns the number of outlets found */
int read_outlet(FILE *F, uint64_t **xcoord, uint64_t **ycoord, real_t **slopes) {
  int idx=0, numlines=0;
  char buffer[1024];
  uint64_t xx, yy;
  double s;

  while ( (fgets(buffer, 1023, F)) != NULL )
    numlines++;
  rewind(F);

  *xcoord = (uint64_t*) malloc(sizeof(uint64_t)*numlines);
  *ycoord = (uint64_t*) malloc(sizeof(uint64_t)*numlines);
  *slopes = (real_t*) malloc(sizeof(real_t)*numlines);

  while ( (fgets(buffer, 1023, F)) != NULL ) {
    sscanf(buffer, "%ld %ld %lf", &xx, &yy, &s);
    (*xcoord)[idx]= xx;
    (*ycoord)[idx]= yy;
    if ( s < 1e-5) {
      fprintf(stdout,"  Dangerously small slope found at location (%ld,%ld) of value %.3e\n",
          xx,yy,s);
    }
    (*slopes)[idx++] = (real_t) s;
  }

  return idx;
}

/* find the first instance of ".", and truncate everything afterwards, remove all  */
char* find_basename(char *input) {
  char *ptr = basename(input);
  char *rc = ptr;
  while ( *ptr != '.' ) ptr++;
  *ptr = '\0';
  return rc;
}


void georeference(NetCDF *netcdf, string targetFile, string referenceFile)
{
  netcdf->copyVar(referenceFile, "", targetFile);
  netcdf->copyVar(referenceFile, "lat", targetFile);
  netcdf->copyVar(referenceFile, "lon", targetFile);
  netcdf->copyVar(referenceFile, "crs", targetFile);
}

bool saveState(Engine *engine, const char *stateFile) {
  struct stat statbuf;
  if (stat(stateFile, &statbuf) == 0)
    unlink(stateFile);

  Grid *h           = engine->getHeight();
  Grid *olr         = engine->getOLR();
  Grid *olrdim0_old = engine->getOLRDIM0_OLD();
  Grid *olrdim1_old = engine->getOLRDIM1_OLD();
  Grid *pre         = engine->getPRE();
  Grid *ret         = engine->getRET();
  Grid *inth        = engine->getIntH();
  Grid *vsat        = engine->getVSAT();
  Grid *grids[] = { h, olr, olrdim0_old, olrdim1_old, pre, ret, inth, vsat, NULL };
  string varnames[] = { "Band1", "Band2", "Band3", "Band4", "Band5", "Band6", "Band7", "Band8" };

  list<Grid*> gridList;
  list<string> varnameList;
  for (int i=0; grids[i]; ++i) {
    gridList.push_back(grids[i]);
    varnameList.push_back(varnames[i]);
  }

  NetCDF netcdf(h->cellsize);
  netcdf.write(stateFile, gridList, varnameList);

  for (int i=0; grids[i]; ++i)
    delete grids[i];
  return true;
}

bool restoreState(Engine *engine, const char *stateFile) {
  struct stat statbuf;
  if (stat(stateFile, &statbuf) != 0) {
    printf("Warning: simulation state file '%s' does not exist, doing a cold start\n", stateFile);
    return false;
  }

  NetCDF netcdf(90);
  Grid *h           = netcdf.parse(stateFile, "Band1");
  Grid *olr         = netcdf.parse(stateFile, "Band2");
  Grid *olrdim0_old = netcdf.parse(stateFile, "Band3");
  Grid *olrdim1_old = netcdf.parse(stateFile, "Band4");
  Grid *pre         = netcdf.parse(stateFile, "Band5");
  Grid *ret         = netcdf.parse(stateFile, "Band6");
  Grid *inth        = netcdf.parse(stateFile, "Band7");
  Grid *vsat        = netcdf.parse(stateFile, "Band8");

  engine->setHeight(h);
  engine->setOLR(olr);
  engine->setOLRDIM0_OLD(olrdim0_old);
  engine->setOLRDIM1_OLD(olrdim1_old);
  engine->setPRE(pre);
  engine->setRET(ret);
  engine->setIntH(inth);
  engine->setVSAT(vsat);

  delete h;
  delete olr;
  delete olrdim0_old;
  delete olrdim1_old;
  delete pre;
  delete ret;
  delete inth;
  delete vsat;
  return true;
}

void dumpInformation(const char *fname, const char *base, int ismaskgen, OptionParser *options)
{
  FILE *F = fopen(fname,"w");
  if (!F) {
    printf("Warning: unable to write summary file \"%s\"\n", fname);
  } else {
    time_t now;
    struct tm *timeinfo;
    now = time(0);
    timeinfo = localtime( &now);
    fprintf(F, "===== IFM run of %s =====\n", base);
    fprintf(F, "Simulation performed on          \t\t%s\n", asctime( timeinfo));
    fprintf(F, "Floating point type:             \t\t%s\n", sizeof(real_t) == sizeof(float) ? "FLOAT" : "DOUBLE");
    fprintf(F, "\n");
    fprintf(F, "DEM input file:                  \t\t%s\n", options->dem.c_str());
    fprintf(F, "Mask file:                       \t\t%s\n", ismaskgen==0?options->mask.c_str():"GENERATED FROM DEM");
    fprintf(F, "LandUse file:                    \t\t%s\n", options->landuse.c_str());
#ifndef ROUTING_ONLY
    fprintf(F, "Soil hydraulic conductivity:     \t\t%s\n", options->soil_hc.c_str());
    fprintf(F, "Soil pressure head:              \t\t%s\n", options->soil_ph.c_str());
    fprintf(F, "Soil effective porosity:         \t\t%s\n", options->soil_ep.c_str());
    fprintf(F, "Soil moisture:                   \t\t%s\n", options->soilMoisture.c_str());
#endif
    fprintf(F, "Outlet file:                     \t\t%s\n", options->outlet.c_str());
    fprintf(F, "Precipitation file:              \t\t%s\n", options->precipitation.c_str());
    fprintf(F, "Output flow rate:                \t\t%s\n", options->saveFlowRate ? "YES" : "NO");
    fprintf(F, "Output depth:                    \t\t%s\n", options->saveDepth ? "YES" : "NO");
    fprintf(F, "Output max volume:               \t\t%s\n", options->saveVolume ? "YES" : "NO");
    fprintf(F, "Output OLR:                      \t\t%s\n", options->saveOLR ? "YES" : "NO");
    fprintf(F, "Output VSAT:                     \t\t%s\n", options->saveVSAT ? "YES" : "NO");
    fprintf(F, "\n");
    fprintf(F, "Simulation time step:            \t\t%.1f sec\n", options->tstep);
    fprintf(F, "Simulation Tstop:                \t\t%d sec\n", options->tfinal);
    fprintf(F, "Printing interval:               \t\t%d sec\n", options->outputRate);
    fclose(F);
  }
}

#define DELETE_OBJS() do { \
  if (dem) delete dem; \
  if (mask) delete mask; \
  if (landuse) delete landuse; \
  if (landuseManningsN) delete landuseManningsN; \
  if (landuseRetention) delete landuseRetention; \
  if (soil_hc) delete soil_hc; \
  if (soil_ph) delete soil_ph; \
  if (soil_ep) delete soil_ep; \
  if (soilMoisture) delete soilMoisture; \
  if (maskData) free(maskData); \
  if (landuseData) free(landuseData); \
  if (watchList) delete [] watchList; \
  if (out_x) free(out_x); \
  if (out_y) free(out_y); \
  if (out_s) free(out_s); \
  if (pt) delete pt; \
} while(0)

#define SANITY_CHECK(obj,ref) do { \
  if (! (obj)) { \
    DELETE_OBJS(); \
    return 1; \
  } \
  if ( (obj)->cellsize != (ref)->cellsize ) { \
    printf("Warning: inconsistent grid size found in file \"%s\". Ignored\n", (obj)->filename.c_str()); \
  } \
  if ( (obj)->rows != (ref)->rows || (obj)->cols != (ref)->cols ) { \
    printf("Bummer: inconsistent row/col sizes found in file \"%s\". Found (%ld,%ld), expecting (%ld,%ld)\n", \
        (obj)->filename.c_str(), (obj)->rows, (obj)->cols, (ref)->rows, (ref)->cols); \
    DELETE_OBJS(); \
    return 1; \
  } \
} while(0)

#define UPDATE_NODATA(obj,objdata) do { \
  uint64_t len = (uint64_t) (obj)->rows * (uint64_t) (obj)->cols; \
  for (uint64_t jj=0; jj<len; jj++) \
    (objdata)[jj] = ((objdata)[jj] == (obj)->nodata ? 0 : (objdata)[jj]); \
} while(0)

#define HARDCODE_VALUE(obj,objdata,val) do { \
  uint64_t len = (uint64_t) (obj)->rows * (uint64_t) (obj)->cols; \
  for (uint64_t jj=0; jj<len; jj++) \
    (objdata)[jj] = val; \
} while(0)

#define MULTIPLY_BY(obj,objdata,val) do { \
  uint64_t len = (uint64_t) (obj)->rows * (uint64_t) (obj)->cols; \
  for (uint64_t jj=0; jj<len; jj++) \
    (objdata)[jj] *= val; \
} while(0)

int main(int argc, char* argv[]) {
  // File I/O
  precip_table *pt=NULL;
  FILE         *F=NULL;
  char          fname[512];
  WatchPoint   *watchList = NULL;
  int           num_watches=0;
  struct stat   statbuf;

  // DEM, landuse, soil
  Grid         *soil_hc=NULL, *soil_ph=NULL, *soil_ep=NULL, *soilMoisture=NULL;
  Grid         *dem=NULL, *mask=NULL, *landuse=NULL;
  short        *maskData=NULL, *landuseData=NULL;
  int           ismaskgen;
  real_t       *landuseManningsN=NULL;
  real_t       *landuseRetention=NULL;
  uint64_t      landuseMapSize=0;

  // Outlets
  uint64_t     *out_x=NULL, *out_y=NULL;
  real_t       *out_s=NULL;
  uint64_t      num_outlets=0;

  // Command line argument parsing
  OptionParser options(argc, argv);
  options.parse();

  dem = parse_as_double(options.dem, options.cellsize);
  if (! dem)
    return 1;

  assert( dem->rows > 0 && dem->cols > 0);

  if (options.nodata.size())
    set_nodata_cells(dem, options.nodata);

  printf("Read DEM data from file \"%s\"\n", options.dem.c_str());

  if (options.mask.size()) {
    mask = parse_as_short(options.mask, &maskData, options.cellsize, false);
    ismaskgen = 0;
    SANITY_CHECK(mask, dem);
    // Make sure we're using a uniform mask
    for (uint64_t jj=0; jj<mask->rows * mask->cols; jj++)
      maskData[jj] = maskData[jj] == mask->nodata ? 0 : 1;
    printf("Read mask data from file \"%s\"\n", options.mask.c_str());
  } else {
    mask = gen_mask_from_dem(dem, &maskData);
    ismaskgen = 1;
    printf("Generated mask from dem data.\n");
  }

  NetCDF netcdf(options.cellsize);

  if (options.watch.size()) {
    Grid *lat = netcdf.parse(options.dem, "lat");
    Grid *lon = netcdf.parse(options.dem, "lon");
    if (! lat || ! lon) {
      printf("Bummer: cannot get the lat/lon from the DEM file \"%s\"\n", options.dem.c_str());
      return 1;
    }
    num_watches = read_watch_file(lat, lon, &watchList, options.watch);
    if (num_watches < 0) {
      return 1;
    }
    delete lat;
    delete lon;
  }

  if (stat(options.landuse.c_str(), &statbuf) == 0) {
    landuse = parse_as_short(options.landuse, &landuseData, options.cellsize);
    SANITY_CHECK(landuse, dem);
    UPDATE_NODATA(landuse, landuseData);
    if (options.landuseMap.size()) {
      MapFile mapfile(options.landuseMap);
      landuseManningsN = mapfile.getList(1, &landuseMapSize);
      assert(landuseManningsN);
      landuseRetention = mapfile.getList(2);
      assert(landuseRetention);
    }
    printf("Read land use data from file \"%s\"\n", options.landuse.c_str());
  } else {
    landuse = new Grid(dem->cols, dem->rows);
    landuseData = float2short(landuse);
    HARDCODE_VALUE(landuse, landuseData, 0);
    landuseMapSize = 1;
    landuseManningsN = new real_t[landuseMapSize];
    landuseRetention = new real_t[landuseMapSize];
    landuseManningsN[0] = atof(options.landuse.c_str());
    landuseRetention[0] = 0.0;
    printf("Assigned uniform value %s to land use\n", options.landuse.c_str());
  }

  if (options.soil_hc.size() && options.soil_ph.size() && options.soil_ep.size()) {
    if (stat(options.soil_hc.c_str(), &statbuf) == 0) {
      soil_hc = parse_as_double(options.soil_hc, options.cellsize);
      SANITY_CHECK(soil_hc, dem);
      UPDATE_NODATA(soil_hc, soil_hc->data);
      MULTIPLY_BY(soil_hc, soil_hc->data, options.soil_hc_multiplier);
      printf("Read soil hydraulic conductivity data from file \"%s\" and multiplied by %f\n",
        options.soil_hc.c_str(), options.soil_hc_multiplier);
    } else {
      soil_hc = new Grid(dem->cols, dem->rows);
      HARDCODE_VALUE(soil_hc, soil_hc->data, atof(options.soil_hc.c_str()));
      printf("Assigned uniform value %s to soil hydraulic conductivity\n", options.soil_hc.c_str());
    }

    if (stat(options.soil_ph.c_str(), &statbuf) == 0) {
      soil_ph = parse_as_double(options.soil_ph, options.cellsize);
      SANITY_CHECK(soil_ph, dem);
      UPDATE_NODATA(soil_ph, soil_ph->data);
      MULTIPLY_BY(soil_ph, soil_ph->data, options.soil_ph_multiplier);
      printf("Read soil pressure head data from file \"%s\"\n", options.soil_ph.c_str());
    } else {
      soil_ph = new Grid(dem->cols, dem->rows);
      HARDCODE_VALUE(soil_ph, soil_ph->data, atof(options.soil_ph.c_str()));
      printf("Assigned uniform value %s to soil pressure head\n", options.soil_ph.c_str());
    }

    if (stat(options.soil_ep.c_str(), &statbuf) == 0) {
      soil_ep = parse_as_double(options.soil_ep, options.cellsize);
      SANITY_CHECK(soil_ep, dem);
      UPDATE_NODATA(soil_ep, soil_ep->data);
      MULTIPLY_BY(soil_ep, soil_ep->data, options.soil_ep_multiplier);
      printf("Read soil effective porosity data from file \"%s\"\n", options.soil_ep.c_str());
    } else {
      soil_ep = new Grid(dem->cols, dem->rows);
      HARDCODE_VALUE(soil_ep, soil_ep->data, atof(options.soil_ep.c_str()));
      printf("Assigned uniform value %s to soil effective porosity\n", options.soil_ep.c_str());
    }
  }

  if (stat(options.soilMoisture.c_str(), &statbuf) == 0) {
    soilMoisture = parse_as_double(options.soilMoisture, options.cellsize);
    SANITY_CHECK(soilMoisture, dem);
    UPDATE_NODATA(soilMoisture, soilMoisture->data);
    printf("Read soil moisture data from file \"%s\"\n", options.soilMoisture.c_str());
  } else {
    soilMoisture = new Grid(dem->cols, dem->rows);
    HARDCODE_VALUE(soilMoisture, soilMoisture->data, atof(options.soilMoisture.c_str()));
    printf("Assigned uniform value %s to soil moisture\n", options.soilMoisture.c_str());
  }

  F = fopen(options.outlet.c_str(),"r");  // outlet file
  if (F) {
    num_outlets = read_outlet(F, &out_x, &out_y, &out_s);
    fclose(F);
    printf("Read outlet from file \"%s\", found %ld outlets to the watershed\n", options.outlet.c_str(), num_outlets);
  } else {
    out_x = (uint64_t*) malloc(sizeof(uint64_t));
    out_y = (uint64_t*) malloc(sizeof(uint64_t));
    out_s = (real_t *) malloc(sizeof(real_t));
    *out_x = *out_y = 1;
    *out_s = 0.003539;
    num_outlets = 1;
    printf("Assigned a dummy outlet at 1,1\n");
  }
  pt = new precip_table();
  if (pt->parse(options.precipitation, dem) < 0) {
    return 1;
  }
  printf("Read precipitation from %s.\n", options.precipitation.c_str());

  /* print statistics to the summary fine */
  char *base, basebuffer[512];
  strncpy(basebuffer, options.dem.c_str(), 511);
  basebuffer[511]='\0';
  base = find_basename(basebuffer);
  sprintf(fname,"%s/%s.summary", options.outdir.c_str(), base);
  dumpInformation(fname, base, ismaskgen, &options);
  dumpInformation("/dev/stdout", base, ismaskgen, &options);

  Engine engine;
  engine.setSoil(soil_hc, soil_ph, soil_ep, soilMoisture);
  engine.setLand(landuseData, landuseManningsN, landuseRetention, landuseMapSize);
  engine.setOutlets(out_x, out_y, out_s, num_outlets);
  engine.setDomain(dem, maskData);
  engine.setPrecipitation(pt);
  engine.setup();
  if (options.simulationState.size())
    restoreState(&engine, options.simulationState.c_str());

  bool draining = true;
  bool reported = false;
  real_t dt = options.tstep;
  int magic = (int)((real_t)options.outputRate/dt);
  int start_pt = (int)((real_t)options.tbegin/dt);
  int num_pt =(int)((real_t)options.tfinal/dt);

  for (int kk=start_pt; kk<num_pt || draining == false; kk++) {
    engine.run(kk, dt);
    rotating_wheel();

    /* selective printing */
    if ( kk%magic==(magic-1)) {
      int num_seconds = (int)((kk+1)*dt);
      char fname[512];

      // Check if all watch points are draining
      draining = true;
      for (int w=0; w<num_watches; ++w) {
        uint64_t x = watchList[w].x;
        uint64_t y = watchList[w].y;
        watchList[w].push(engine.getHeight(x, y));
        if (watchList[w].isDraining())
          watchList[w].report(num_seconds);
        else
          draining = false;
      }
      if (kk >= num_pt && ! reported) {
        printf("Extending simulation time until all watch points start to drain.\n");
        reported = true;
      }

      if (options.saveDepth) {
        // Dump the depth
        sprintf(fname, "%s/depth_%07d.nc", options.outdir.c_str(), num_seconds);
        Grid *h = engine.getHeight();
        netcdf.write(fname, h, "H", maskData);
        georeference(&netcdf, fname, options.dem);
        delete h;
      }
      if (options.saveOLR) {
        // Dump the OLR
        sprintf(fname, "%s/olr_%07d.nc", options.outdir.c_str(), num_seconds);
        Grid *g = engine.getOLR();
        netcdf.write(fname, g, "OLR", maskData);
        georeference(&netcdf, fname, options.dem);
        delete g;
      }
      if (options.saveVSAT) {
        // Dump the VSAT
        sprintf(fname, "%s/vsat_%07d.nc", options.outdir.c_str(), num_seconds);
        Grid *g = engine.getVSAT();
        netcdf.write(fname, g, "VSAT", maskData);
        georeference(&netcdf, fname, options.dem);
        delete g;
      }
      if (options.saveVolume) {
        // Dump the volume
        sprintf(fname, "%s/maxvolume_%07d.nc", options.outdir.c_str(), num_seconds);
        Grid *v = engine.getVolume();
        netcdf.write(fname, v, "VOL", maskData);
        georeference(&netcdf, fname, options.dem);
        delete v;
      }
      if (options.saveFlowRate) {
        // Dump the flow rate (Q)
        // This quantity says how much water flows at different locations.
        // Q is computed in x-direction, then y-direction. The Pythagorean theorem
        // is used to compute the compound Q: Q = sqrt( q_x*q_x + q_y*q_y).
        sprintf(fname, "%s/flowrate_%07d.nc", options.outdir.c_str(), num_seconds);
        Grid *olr_x = engine.getOLRDIM0_OLD();
        Grid *olr_y = engine.getOLRDIM1_OLD();
        Grid q(olr_x->cols, olr_x->rows);
        uint64_t len = olr_x->rows * olr_x->cols;
        for (uint64_t qi=0; qi<len; ++qi)
          q.data[qi] = sqrt(olr_x->data[qi]*olr_x->data[qi] + olr_y->data[qi]*olr_y->data[qi]);
        netcdf.write(fname, &q, "Q", maskData);
        georeference(&netcdf, fname, options.dem);
        delete olr_y;
        delete olr_x;
      }
    }
  }
  printf("Done\n");

  if (options.simulationState.size()) {
    printf("Saving simulation state to %s\n", options.simulationState.c_str());
    saveState(&engine, options.simulationState.c_str());
  }

  DELETE_OBJS();
  return 0;
}

/* end */
