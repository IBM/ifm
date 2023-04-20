#include <stdio.h>
#include <getopt.h>
#include <errno.h>
#include <sys/stat.h>
#include <gd.h>
#include <list>
#include <vector>
#include <string>
#include <fstream>
#include "FileIO/netcdf_io.h"
#include "Core/engine.h"
#include "colormaps.h"

//#define NODATA_COLOR     0x0080ff
#define NODATA_COLOR     0xffffff
#define MAX_COLORS       sizeof(colormap_mild)/sizeof(colormap_mild[0])
#define MAX_RIVER_COLORS sizeof(colormap_river)/sizeof(colormap_river[0])

class OptionParser {
  public:
    OptionParser() : colors(NULL), colorMap(""), scale(1.0), min(0.0), river(false)
    {
      inputFiles.clear();
      riverFiles.clear();
    }

    void usage(char *appName, int error)
    {
      printf("Usage: %s <options> <input file(s)>\n"
             "Valid options are:\n"
             "    -c, --colormap=NAME      One of 'drastic', 'mild', or 'grayscale'.\n"
             "    -m, --minimum=NUM        Minimum value to render (default: 0.0).\n"
             "    -s, --scale=NUM          Multiply height values by NUM (default: 1.0).\n"
             "    -h, --help               This help.\n"
             "    -r, --river              Include river data.\n\n",
             appName);
      printf("River data is automatically detected. The river data files must be in\n"
             "the same directory as the input files and named after them, with the\n"
             "extension .spt. (e.g.: depth_00300.nc and depth_00300.spt)\n\n");
      printf("The default downscale factor (1) will cause water heights higher than 500mm\n"
             "to be painted in extreme colors. A factor or 2 will shift the extreme color\n"
             "boundary to 1000mm (1 meter), and so on.\n\n");
      exit(error);
    }

    void parseArgs(int argc, char **argv)
    {
      const char *short_options = "c:hm:s:r";
      const struct option long_options[] = {
        { "colormap",  required_argument, 0, 'c' },
        { "scale",     required_argument, 0, 'd' },
        { "help",      no_argument,       0, 'h' },
        { "minimum",   required_argument, 0, 'm' },
        { "river",     no_argument,       0, 'r' },
        { NULL,        0,                 0,  0  }
      };

      while (2) {
        int option_index = 0;
        int c = getopt_long(argc, argv, short_options, long_options, &option_index);
        if (c == -1)
          break;
        switch (c) {
          case 'c':
            colorMap.assign(optarg);
            break;
          case 'h':
            usage(argv[0], 0);
            break;
          case 'm':
            min = atof(optarg);
            break;
          case 'r':
            river = true;
            break;
          case 's':
            scale = atof(optarg);
            break;
          case '?':
          default:
            break;
        }
      }
      while (optind < argc) {
        const char *sep = strrchr(argv[optind], '.');
        struct stat statbuf;
        if (sep && river == true) {
          string riverFile((const char *) argv[optind], sep-argv[optind]);
          riverFile.append(".spt");
          if (stat(riverFile.c_str(), &statbuf) == 0)
            riverFiles.push_back(riverFile);
        }
        inputFiles.push_back(argv[optind++]);
      }

      // Colormap check
      if (colorMap.compare("mild") == 0)
        colors = colormap_mild;
      else if (colorMap.compare("drastic") == 0)
        colors = colormap_drastic;
      else if (colorMap.compare("grayscale") == 0 || colorMap.compare("greyscale") == 0)
        colors = colormap_grayscale;
      else {
        fprintf(stderr, "Error: please specify a valid colormap.\n\n");
        usage(argv[0], 1);
      }

      // Input file(s) check
      if (inputFiles.size() == 0) {
        fprintf(stderr, "Error: please provide at least one input file to process.\n\n");
        usage(argv[0], 1);
      }

      // Scale factor check
      if (scale == 0) {
        fprintf(stderr, "Error: the scale factor cannot be zero.\n\n");
        usage(argv[0], 1);
      }

      if (riverFiles.size() > 0 && riverFiles.size() != inputFiles.size()) {
        fprintf(stderr, "Warning: there are %ld input files versus %ld river files\n", inputFiles.size(), riverFiles.size());
        fprintf(stderr, "Ignoring river files.\n");
        riverFiles.clear();
        river = false;
      }
    }

    // public variables
    const int *colors;
    string appName;
    string colorMap;
    float scale;
    float min;
    bool river;
    vector<string> inputFiles;
    vector<string> riverFiles;
};

void initColors(gdImagePtr img, int *color, const int len, const int *referenceColors, bool setBgColor)
{
  int idx = 0;
  if (setBgColor) {
    // Background color (first allocated)
    color[0] = gdImageColorAllocateAlpha(img, 0, 0, 0, 0);
    gdImageColorTransparent(img, 0);
    idx++;
  }
  for (int i=idx; i<len; ++i) {
    int r = (referenceColors[i] >> 16) & 0xff;
    int g = (referenceColors[i] >> 8) & 0xff;
    int b = (referenceColors[i] & 0xff);
    color[i] = gdImageColorAllocateAlpha(img, r, g, b, 0);
  }
}

struct painterData {
  const OptionParser *opt;
  pthread_mutex_t *netcdf_lock;
  float scale;
  string riverFile;
  string inputFile;
  string outputFile;
};

void *painter(void *inputData)
{
  struct painterData *data = (struct painterData *) inputData;
  NetCDF netcdf(90);

  // Accumulated precipitation file
  pthread_mutex_lock(data->netcdf_lock);
  Grid *input = netcdf.parse(data->inputFile);
  pthread_mutex_unlock(data->netcdf_lock);
  if (! input) {
    fprintf(stderr, "Failed to parse input file %s\n", data->inputFile.c_str());
    delete data;
    pthread_exit(NULL);
  }

  // Initialize the colors array
  gdImagePtr img = gdImageCreateTrueColor(input->cols, input->rows);
  int color, colors[MAX_COLORS], riverColors[MAX_RIVER_COLORS], nodataColor[1];
  int nodataRefColor[1] = { NODATA_COLOR };
  initColors(img, colors, MAX_COLORS, data->opt->colors, true);
  initColors(img, riverColors, MAX_RIVER_COLORS, colormap_river, false);
  initColors(img, nodataColor, 1, nodataRefColor, false);

  // Paint the image
  for (uint64_t i=0; i<input->cols; ++i) {
    for (uint64_t j=0; j<input->rows; ++j) {
      float depth = GRID(input, i, j);
      if ((int64_t) depth == input->nodata || depth <= data->opt->min)
        gdImageSetPixel(img, i, j, nodataColor[0]);
      else if (depth > 0) {
        unsigned int depth_in_mm = (unsigned int) ((depth * 1000.0) * data->scale);
        color = depth_in_mm >= MAX_COLORS ? colors[MAX_COLORS-1] : colors[depth_in_mm];
        gdImageSetPixel(img, i, j, color);
      }
    }
  }

  // Include river points, if set
  if (data->riverFile.size() > 0) {
    ifstream infile(data->riverFile.c_str());
    uint64_t x=0, y=0;
    double baseSpeed=0.0, maxSpeed=0.0, curSpeed=0.0, baseDepth=0.0, maxDepth=0.0, curDepth=0.0;
    string header;
    std::getline(infile, header);
    // 'speed' is river flow rate
    while (infile >> x >> y >> baseSpeed >> maxSpeed >> curSpeed >> baseDepth >> maxDepth >> curDepth) {
      if (curSpeed > maxSpeed) curSpeed = maxSpeed;
      if (curSpeed < baseSpeed) curSpeed = baseSpeed;
      if (curDepth > maxDepth) curDepth = maxDepth;
      if (curDepth < baseDepth) curDepth = baseDepth;
      double speedStep = (maxSpeed - baseSpeed) / MAX_RIVER_COLORS;
      double depthStep = (maxDepth - baseDepth) / MAX_RIVER_COLORS;
      double radius = 2;
      for (unsigned int i=0; i<MAX_RIVER_COLORS; ++i)
        if (curDepth < (baseDepth+((i+1)*depthStep))) {
          radius = i/2.0;
          break;
        }

      color = riverColors[0];
      for (unsigned int i=0; i<MAX_RIVER_COLORS; ++i)
        if (curSpeed < (baseSpeed+((i+1)*speedStep))) {
          color = riverColors[i];
          break;
        }

      gdImageFilledEllipse(img, x, y, (int) radius*2, (int) radius*2, color);
    }
  }

  // Save as PNG
  FILE *output = fopen(data->outputFile.c_str(), "wb+");
  gdImagePng(img, output);
  fclose(output);
  gdImageDestroy(img);

  delete input;
  delete data;
  pthread_exit(NULL);
}

int main(int argc, char **argv)
{
  OptionParser opt;
  opt.parseArgs(argc, argv);

  pthread_mutex_t netcdf_lock = PTHREAD_MUTEX_INITIALIZER;
  pthread_t painters[opt.inputFiles.size()];
  list<string>::iterator it;
  unsigned int i;

  for (i=0; i<opt.inputFiles.size(); ++i) {
    string riverFile = opt.riverFiles.size() > 0UL ? opt.riverFiles[i] : "";
    string inputFile = opt.inputFiles[i];
    string outputFile = inputFile;
    size_t sep = outputFile.find_last_of('.');
    if (sep != string::npos && outputFile.substr(sep).compare(".gz") == 0) {
      size_t sep2 = outputFile.substr(0, sep-1).find_last_of('.');
      if (sep2 != string::npos)
        sep = sep2;
    }
    if (sep != string::npos) {
      outputFile.replace(sep, outputFile.size()-sep, ".png");
    } else {
      outputFile.append(".png");
    }

    struct painterData *data = new struct painterData;
    data->opt = (const OptionParser *) &opt;
    data->netcdf_lock = &netcdf_lock;
    data->scale = opt.scale;
    data->riverFile = riverFile;
    data->inputFile = inputFile;
    data->outputFile = outputFile;
    pthread_create(&painters[i], NULL, painter, data);
  }
  while (i > 0)
    pthread_join(painters[--i], NULL);
  return 0;
}
