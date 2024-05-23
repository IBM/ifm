#include <string>
#include <cstdlib>

#include <string.h>
#include <getopt.h>
#include <libgen.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "FileIO/netcdf_io.h"
#include "FileIO/interpolation.h"
#include "Core/ifm_common.h"

class GetValue {
  public:
    GetValue() : m_dem(""), m_netcdfLayer(""), m_verbose(false) {
      m_filenames.clear();
      m_coords.clear();
      m_pixels.clear();
    };

    void parsePixels(const char *coords);
    void parseCoords(const char *coords);
    int  parseOptions(int argc, char **argv);
    bool parseInputFile(const char *filename);
    bool interpolate(Grid *lat, Grid *lon, real_t latitude, real_t longitude);
    bool getpixel(uint64_t x, uint64_t y);
    void usage(const char *appname, int errorCode);
    bool run();

  private:
    string m_dem;
    string m_netcdfLayer;
    bool m_verbose;
    list<string> m_filenames;
    list<pair<real_t,real_t> > m_coords;
    list<pair<uint64_t,uint64_t> > m_pixels;
};

void GetValue::usage(const char *appname, int errorCode)
{
  cout << "Usage: " << basename((char *) appname) << " [options] <file(s)>" << endl << endl;
  cout << "Optional arguments are:" << endl;
  cout << "   -h, --help              This help" << endl;
  cout << "   -d, --dem=FILE          Georeferenced DEM" << endl;
  cout << "   -v, --verbose           Verbose mode" << endl;
  cout << "   -f, --input=file=FILE   Take input coords from FILE (one lat,lon or x,y pair per line)" << endl;
  cout << "   -l, --latlon=<lat,lon>  Geographic coordinates. May be entered multiple times." << endl;
  cout << "   -n, --netcdf-layer=NAME NetCDF layer to process, if given files are in NetCDF format." << endl;
  cout << "   -p, --pixel=<x,y>       Pixel coordinates. May be entered multiple times." << endl;
  cout << endl;
  cout << "Options --latlon and --pixel are mutually exclusive." << endl;
  cout << endl;
  exit(errorCode);
}

int GetValue::parseOptions(int argc, char **argv)
{
  int c, option_index = 0;
  const char * const short_options = "d:f:hl:n:p:v";
  const struct option long_options[] = {
    { "dem",          required_argument, 0, 'd' },
    { "file",         required_argument, 0, 'f' },
    { "help",         no_argument,       0, 'h' },
    { "latlon",       required_argument, 0, 'l' },
    { "netcdf-layer", required_argument, 0, 'n' },
    { "pixel",        required_argument, 0, 'p' },
    { "verbose",      no_argument,       0, 'v' },
    { NULL, 0, 0, 0 },
  };

  while (true) {
    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1)
      break;
    switch (c) {
      case 'd':
        m_dem = optarg;
        break;
      case 'f':
        this->parseInputFile(optarg);
        break;
      case 'h':
        usage(argv[0], 0);
        break;
      case 'l':
        this->parseCoords(optarg);
        break;
      case 'n':
        m_netcdfLayer = optarg;
        break;
      case 'p':
        this->parsePixels(optarg);
        break;
      case 'v':
        m_verbose = true;
        break;
      case '?':
      default:
        break;
    }
  }

  while (optind < argc)
    m_filenames.push_back(argv[optind++]);

  if (! m_filenames.size()) {
    fprintf(stderr, "Error: missing file list\n");
    this->usage(argv[0], 1);
  } else if (! m_dem.size()) {
    NetCDF nc(90);
    Grid *g = nc.parse(argv[optind-1], "XLAT");
    if (! g) {
      fprintf(stderr, "Error: --dem argument must be given\n");
	  this->usage(argv[0], 1);
    }
    delete g;
  } else if (m_coords.size() && m_pixels.size()) {
    fprintf(stderr, "Error: --latlon and --pixel are mutually exclusive\n");
    this->usage(argv[0], 1);
  }

  return 0;
}

#define GET_PAIR(coords) \
  char *coord1 = strdup(coords); assert(coord1); \
  char *coord2 = strstr(coord1, ","); \
  if (! coord2) { \
    fprintf(stderr, "Failed to find separator (',') in string '%s'\n", coord1); \
    exit(1); \
  } \
  *coord2++ = '\0'

void GetValue::parseCoords(const char *coords)
{
  GET_PAIR(coords);
  real_t lat = atof(coord1);
  real_t lon = atof(coord2);
  free(coord1);
  m_coords.push_back(make_pair(lat,lon));
}

void GetValue::parsePixels(const char *coords)
{
  GET_PAIR(coords);
  uint64_t x = atol(coord1);
  uint64_t y = atol(coord2);
  free(coord1);
  m_pixels.push_back(make_pair(x,y));
}

bool GetValue::parseInputFile(const char *filename)
{
  char *nl, buf[256];
  FILE *fp = fopen(filename, "r");
  if (! fp) {
    fprintf(stderr, "Failed to open %s: %s\n", filename, strerror(errno));
    return false;
  }
  while (! feof(fp)) {
    memset(buf, 0, sizeof(buf));
    fgets(buf, sizeof(buf), fp);
    nl = strchr(buf, '\n');
    if (nl)
      *nl = '\0';
    if (strlen(buf) > 0 && strstr(buf, "."))
      // Assume buf contains a pair of floating point numbers
      this->parseCoords(buf);
    else if (strlen(buf) > 0)
      // Assume buf contains a pair of integers
      this->parsePixels(buf);
  }
  fclose(fp);
  return true;
}

bool GetValue::run()
{
  bool ret = true;
  Grid *lat, *lon;

  if (m_coords.size()) {
    NetCDF netcdf(90);
    if (m_dem.size()) {
      lat = netcdf.parse(m_dem, "lat");
      if (! lat) return false;
      lon = netcdf.parse(m_dem, "lon");
      if (! lon) { delete lat; return false; }
    } else {
      list<string>::iterator it = m_filenames.begin();
      lat = netcdf.parse(*it, "XLAT");
      if (! lat) return false;
      lon = netcdf.parse(*it, "XLONG");
      if (! lon) { delete lat; return false; }
    }

    list<pair<real_t,real_t> >::iterator coord_it;
    for (coord_it=m_coords.begin(); coord_it != m_coords.end(); coord_it++) {
      if (this->interpolate(lat, lon, (*coord_it).first, (*coord_it).second) == false) {
        ret = false;
        break;
      }
    }
    delete lon;
    delete lat;
  }

  if (m_pixels.size() && ret == true) {
    list<pair<uint64_t,uint64_t> >::iterator pixel_it;
    for (pixel_it=m_pixels.begin(); pixel_it != m_pixels.end(); pixel_it++) {
      if (this->getpixel((*pixel_it).first, (*pixel_it).second) == false) {
        ret = false;
        break;
      }
    }
  }

  return ret;
}

bool GetValue::interpolate(Grid *lat, Grid *lon, real_t latitude, real_t longitude)
{
  Interpolation i;
  BoundingBox_t bbox, dataBBox;
  uint64_t sx=0, ex=0, sy=0, ey=0;
  if (i.getXY(lat, lon, latitude, longitude, sx, ex, sy, ey) == false) {
    fprintf(stderr, "Error: the provided lat/lon fall outside the area covered by the input file(s)\n");
    return false;
  }
  bbox.startx = sx;
  bbox.endx   = ex;
  bbox.starty = sy;
  bbox.endy   = ey;
  dataBBox.startx = dataBBox.starty = 0;
  dataBBox.endx = bbox.endx - bbox.startx;
  dataBBox.endy = bbox.endy - bbox.starty;

  if (m_verbose)
    printf("Lat,Lon=%f,%f -> X=(%ld,%ld), Y=(%ld,%ld)\n", 
        latitude, longitude, bbox.startx, bbox.endx, bbox.starty, bbox.endy);

  int seq = 0;
  NetCDF netcdf(90);
  list<string>::iterator it;
  for (it=m_filenames.begin(); it != m_filenames.end(); it++) {
    Grid *pp = netcdf.parse(*it, m_netcdfLayer, &bbox);
    if (! pp)
      return false;
    printf("%d %f\n", seq++, i.bilinear(pp, &dataBBox, lat, lon, &bbox, latitude, longitude));
    delete pp;
  }
  return true;
}

bool GetValue::getpixel(uint64_t x, uint64_t y)
{
  BoundingBox_t bbox;
  bbox.startx = x; bbox.endx = x;
  bbox.starty = y; bbox.endy = y;

  if (m_verbose)
    printf("Pixel value at X=%ld, Y=%ld\n", x, y);

  int seq = 0;
  NetCDF netcdf(90);
  list<string>::iterator it;
  for (it=m_filenames.begin(); it != m_filenames.end(); it++) {
    Grid *pp = netcdf.parse(*it, m_netcdfLayer, &bbox);
    if (! pp)
      return false;
    printf("%d %f\n", seq++, GRID(pp,0,0));
    delete pp;
  }
  return true;
}

int main(int argc, char **argv)
{
  GetValue gv;
  int ret;

  ret = gv.parseOptions(argc, argv);
  if (ret < 0)
    return ret;

  return gv.run() == true ? 0 : 1;
}
