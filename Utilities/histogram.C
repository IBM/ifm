#include "FileIO/netcdf_io.h"
#include "Core/engine.h"
#include "FileIO/interpolation.h"
#include "histogram.h"

void Histogram::usage(const char *appname, int errorCode)
{
  cout << "Usage: " << basename((char *) appname) << " <options>" << endl << endl;
  cout << "The valid arguments are:" << endl;
  cout << "   -b, --bin-size=NUM      Bin size (default: " << m_binSize << ")" << endl;
  cout << "   -d, --dem=FILE          Georeferenced DEM" << endl;
  cout << "   -h, --help              This help" << endl;
  cout << "   -D, --output-dem=FILE   Output file in DEM format" << endl;
  cout << "   -P, --output-png=FILE   Output file in PNG format" << endl;
  cout << "   -r, --radius=NUM        Distance in pixels between the input coordinate and the bbox corner";
  cout <<                             " (default: " << m_radius << ")" << endl;
  cout << "   -v, --verbose           Verbose mode" << endl;
  cout << "   -L, --latlon=<lat,lon>  Input coords, given as latitude,longitude" << endl;
  cout << endl;
  exit(errorCode);
}

int Histogram::parseOptions(int argc, char **argv)
{
  int c, option_index = 0;
  const char * const short_options = "b:d:hD:P:r:vL:";
  const struct option long_options[] = {
    { "bin-size", 1, 0, 'b' },
    { "dem", 1, 0, 'd' },
    { "help", 0, 0, 'h' },
    { "output-dem", 1, 0, 'D' },
    { "output-png", 1, 0, 'P' },
    { "radius", 1, 0, 'r' },
    { "verbose", 0, 0, 'v' },
    { "latlon", 1, 0, 'L' },
    { NULL, 0, 0, 0 },
  };

  while (true) {
    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1)
      break;
    switch (c) {
      case 'b':
        m_binSize = atoi(optarg);
        break;
      case 'd':
        m_dem = optarg;
        break;
      case 'h':
        usage(argv[0], 0);
        break;
      case 'D':
        m_outputDEM.assign(optarg);
        break;
      case 'P':
        m_outputPNG.assign(optarg);
        break;
      case 'r':
        m_radius = atoi(optarg);
        break;
      case 'v':
        m_verbose = true;
        break;
      case 'L':
        this->parseCoords(optarg);
        break;
      case '?':
      default:
        break;
    }
  }

  if (! m_dem.size()) {
    fprintf(stderr, "Error: --dem argument must be given\n");
    this->usage(argv[0], 1);
  } else if (! m_outputPNG.size()) {
    fprintf(stderr, "Error: --output-png argument must be given\n");
    this->usage(argv[0], 1);
  } else if (! m_outputDEM.size()) {
    fprintf(stderr, "Error: --output-dem argument must be given\n");
    this->usage(argv[0], 1);
  }

  return 0;
}

void Histogram::parseCoords(const char *coords)
{
  char *dup = strdup(coords);
  if (! dup) {
    fprintf(stderr, "Out of memory\n");
    exit(1);
  }
  char *sep = strstr(dup, ",");
  if (! sep) {
    fprintf(stderr, "Failed to find separator (',') in string '%s'\n", dup);
    free(dup);
    exit(1);
  }
  *sep = '\0';
  sep++;
  m_lat = atof(dup);
  m_lon = atof(sep);
  free(dup);
}

BoundingBox_t *Histogram::calculatePOI(Grid *lat, Grid *lon, Grid *dem, double *poi)
{
  assert(poi);

  Interpolation i;
  BoundingBox_t *bbox = new BoundingBox_t;
  if (i.getXY(lat, lon, m_lat, m_lon, bbox->startx, bbox->endx, bbox->starty, bbox->endy) == false) {
    fprintf(stderr, "Error: the provided lat/lon fall outside the area covered by the input file(s)\n");
    fprintf(stderr, "Lat=%f, Lon=%f\n", m_lat, m_lon);
    fprintf(stderr, "DEM Lat=%f..%f Lon=%f..%f\n", GRID(lat,0,0), GRID(lat,0,lat->rows-1), GRID(lon,0,0), GRID(lon,lon->cols-1,0));
    delete bbox;
    return NULL;
  }
  if (bbox->startx == bbox->endx) { if (bbox->startx>0) bbox->startx--; else bbox->endx++; }
  if (bbox->starty == bbox->endy) { if (bbox->starty>0) bbox->starty--; else bbox->endy++; }
  *poi = i.bilinear(dem, bbox, lat, lon, bbox, m_lat, m_lon);
  return bbox;
}

int *Histogram::calculateIndex(Grid *dem, BoundingBox_t *bbox, double poi, double min, double max, int *poi_index)
{
  assert(poi_index);

  // Calculate the index of the point of interest and create the histogram  
  int *bin = new int[m_binSize];
  *poi_index = m_binSize-1;

  double step = (max - min) / m_binSize;
  memset(bin, 0, sizeof(int) * m_binSize);
  for (int bi=0; bi<m_binSize; ++bi) {
    float binMin = min + (step * bi);
    float binMax = min + (step * (bi+1));
    if (poi >= binMin && poi < binMax) {
      bin[bi] += 1;
      *poi_index = bi;
      break;
    }
  }
  for (int i=0; i<m_radius*2; ++i) {
    int x = (bbox->startx-m_radius) + i;
    if (x < 0 || x >= dem->cols)
      continue;

    for (int j=0; j<m_radius*2; ++j) {
      int y = (bbox->starty-m_radius) + j;
      if (y < 0 || y >= dem->rows)
        continue;

      double ele = GRID(dem, x, y);
      if ((int) ele == dem->nodata)
        continue;

      for (int bi=0; bi<m_binSize; ++bi) {
        double binMin = min + (step * bi);
        double binMax = min + (step * (bi+1));
        if (ele >= binMin && ele < binMax) {
          bin[bi] += 1;
          break;
        }
      }
    }
  }

  if (m_verbose)
    printf("POI: %f at index %d and offset %d,%d (min=%f, max=%f, step=%f)\n", 
      poi, *poi_index, bbox->startx, bbox->starty, min, max, step);

  return bin;
}

bool Histogram::plot(int *bin, int poi_index, double min, double max)
{
  FILE *fp = fopen("gnuplot.data", "w");
  if (! fp) {
    perror("gnuplot.data");
    return false;
  }
  double step = (max - min) / m_binSize;
  for (int i=0; i<m_binSize; ++i) {
    if (m_verbose)
      printf("%.2f %d %d\n", i*step, bin[i], i == poi_index ? 1 : 0);
    fprintf(fp, "%.2f %d %d\n", i*step, bin[i], i == poi_index ? 1 : 0);
  }
  fclose(fp);

  fp = fopen("gnuplot.cmd", "w");
  if (! fp) {
    perror("gnuplot.cmd");
    unlink("gnuplot.data");
    return false;
  }
  fprintf(fp, "reset\n");
  fprintf(fp, "set terminal png font '/usr/share/fonts/dejavu/DejaVuSans.ttf' 10\n");
  fprintf(fp, "set output \"%s\"\n", m_outputPNG.c_str());
  fprintf(fp, "set style line 1 lt 1 lc rgb \"#ff0000\"\n");
  fprintf(fp, "set style line 2 lt 1 lc rgb \"#008b8b\"\n");
  fprintf(fp, "set style fill solid 1.00 border 0\n");
  fprintf(fp, "set label \"*\" at %.2f,%d\n", poi_index-0.05, bin[poi_index]+2);
  fprintf(fp, "set style histogram\n");
  fprintf(fp, "set style data histogram\n");
  fprintf(fp, "set xlabel \"Elevation difference to lowest point\" \n");
  fprintf(fp, "set ylabel \"Number of points\"\n");
  fprintf(fp, "plot 'gnuplot.data' u (column(0)):2:(0.5):($3>0?1:2):xtic(1) ti 'Histogram' with boxes lc variable\n");
  fclose(fp);

  system("gnuplot gnuplot.cmd");
  unlink("gnuplot.data");
  unlink("gnuplot.cmd");
  return true;
}

void Histogram::calculateMinMax(const Grid *dem, const BoundingBox_t *bbox, Grid *newDem, double *min, double *max)
{
  int newCol = 0, newRow = 0;

  assert(min);
  assert(max);

  *max = -99999999.0;
  *min =  99999999.0;
  if (m_verbose) printf("Grid around the point of interest:\n");
  for (int x=(bbox->startx-m_radius); x<bbox->endx+m_radius; ++x) {
    if (x >= 0 && x < dem->cols) {
      for (int y=(bbox->starty-m_radius); y<bbox->endy+m_radius; ++y) {
        if (y >= 0 && y < dem->rows) {
          double ele = GRID(dem, x, y);
          if (newDem) GRID(newDem, newCol++, newRow) = ele;
          if (m_verbose) printf("%.3f ", ele);
          if ((int) ele != dem->nodata) {
            if (ele < *min)
              *min = ele;
            if (ele > *max)
              *max = ele;
          }
        }
      }
      if (m_verbose) printf("\n");
      if (newDem) { newCol=0; newRow++; }
    }
  }
  if (m_verbose) printf("\n");
}

bool Histogram::run()
{
  NetCDF netcdf;

  Grid *lat = netcdf.parse(m_dem, "lat");
  if (! lat) {
    return false;
  }
  Grid *lon = netcdf.parse(m_dem, "lon");
  if (! lon) {
    delete lat;
    return false;
  }
  Grid *dem = netcdf.parse(m_dem);
  if (! dem) {
    delete lat;
    delete lon;
    return false;
  }

  double poi = 0.0;
  BoundingBox_t *bbox = this->calculatePOI(lat, lon, dem, &poi);
  if (! bbox) {
    delete lat;
    delete lon;
    delete dem;
    return false;
  }
  
  // Create the new DEM as we iterate the cells
  int cols = (bbox->endx+m_radius) - (bbox->startx-m_radius);
  int rows = (bbox->endy+m_radius) - (bbox->starty-m_radius);
  Grid *newDem = new Grid(cols, rows);
  newDem->gridsize = dem->gridsize;
  newDem->nodata = dem->nodata;
  newDem->filename = m_outputDEM;

  // Get the minimum, maximum elevation values of the neighborhood and the
  // bin step size.
  double max, min;
  this->calculateMinMax(dem, bbox, newDem, &min, &max);
  
  // Calculate the index of the point of interest and create the histogram
  int poi_index = 0;
  int *bin = this->calculateIndex(dem, bbox, poi, min, max, &poi_index);
  bool retval = this->plot(bin, poi_index, min, max);

  // Save the new DEM to disk
  netcdf.write(newDem->filename, newDem, "Band1");

  delete [] bin;
  delete newDem;
  delete bbox;
  delete dem;
  delete lon;
  delete lat;
  return retval;
}
