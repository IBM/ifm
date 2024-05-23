#include "FileIO/netcdf_io.h"
#include "Core/engine.h"
#include <getopt.h>
#include <vector>
#include <libgen.h>

#define MIN(x,y) ((x)<(y)?(x):(y))

class RiskEstimator {
  public:
    RiskEstimator() : m_dem(""), m_layer("Band1"), m_output("risk.nc"), m_verbose(false) {
      m_filenames.clear();
    };

    int  parseOptions(int argc, char **argv);
    void usage(const char *appname, int errorCode);
    bool run();

  private:
    string m_dem;
    string m_layer;
    string m_output;
    bool m_verbose;
    vector<string> m_filenames;
};

void RiskEstimator::usage(const char *appname, int errorCode)
{
  cout << "Usage: " << basename((char *) appname) << " [options] <file(s)>" << endl << endl;
  cout << "Optional arguments are:" << endl;
  cout << "   -h, --help                This help" << endl;
  cout << "   -d, --dem=FILE            DEM file" << endl;
  cout << "   -o, --output=FILE         Output file (default: risk.nc)" << endl;
  cout << "   -l, --layer=NAME          Layer name to process (default: Band1)" << endl;
  cout << "   -v, --verbose             Verbose mode" << endl;
  cout << endl;
  exit(errorCode);
}

int RiskEstimator::parseOptions(int argc, char **argv)
{
  int c, option_index = 0;
  const char * const short_options = "d:hl:o:v";
  const struct option long_options[] = {
    { "help",    no_argument,       0, 'h' },
    { "dem",     required_argument, 0, 'd' },
    { "layer",   required_argument, 0, 'l' },
    { "output",  required_argument, 0, 'o' },
    { "verbose", no_argument,       0, 'v' },
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
      case 'h':
        usage(argv[0], 0);
        break;
      case 'l':
        m_layer = optarg;
        break;
      case 'o':
        m_output = optarg;
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
    fprintf(stderr, "Error: --dem argument must be given\n");
    this->usage(argv[0], 1);
  }

  return 0;
}

bool RiskEstimator::run()
{
  NetCDF netcdf(30);
  Grid *dem = netcdf.parse(m_dem);
  if (! dem)
    return false;

  // Compute the maximum depth observed on each grid cell
  Grid maxDepth(dem->cols, dem->rows);
  maxDepth.reset(maxDepth.nodata);

  for (size_t n=0; n<m_filenames.size(); ++n) {
    if (m_verbose) printf("Processing file %s\n", m_filenames[n].c_str());
    Grid *depth = netcdf.parse(m_filenames[n], m_layer);
    if (! depth)
      return false;
    for (uint64_t i=0; i<(dem->cols * dem->rows); ++i) {
      if (dem->data[i] == dem->nodata || depth->data[i] != depth->data[i])
        maxDepth.data[i] = 0;
	  else if (n == 0 || depth->data[i] > maxDepth.data[i])
        maxDepth.data[i] = depth->data[i];
    }
    delete depth;
  }

  // Compute a flood risk histogram using 10 buckets.
  // 1 meter is considered a drastic enough event.
  const int numBuckets = 10;
  const float highWatermark = 1.0 / numBuckets;

  Grid riskIndex(numBuckets, 1);
  riskIndex.reset(0);
  for (uint64_t i=0; i<(dem->cols * dem->rows); ++i) {
    int index = MIN((int) (maxDepth.data[i] / highWatermark), numBuckets-1);
    riskIndex.data[index] += 1;
  }

  // Write grids to disk
  unlink(m_output.c_str());
  netcdf.write(m_output.c_str(), &maxDepth, "MaxDepth", NULL);
  netcdf.writeHistogram(m_output.c_str(), &riskIndex, "Histogram");

  delete dem;
  return true;
}

int main(int argc, char **argv)
{
  RiskEstimator risk;
  int ret = risk.parseOptions(argc, argv);
  if (ret < 0)
    return ret;
  return risk.run() == true ? 0 : 1;
}
