#ifndef __HISTOGRAM_H
#define __HISTOGRAM_H

#include <string>
#include <cstdlib>

#include <string.h>
#include <getopt.h>
#include <libgen.h>
#include <errno.h>
#include <float.h>
#include <sys/stat.h>
#include <sys/types.h>

class Histogram {
  public:
    Histogram() : m_dem(""), m_verbose(false), m_radius(10), m_binSize(10), 
                  m_lat(0.0), m_lon(0.0), m_outputDEM(""), m_outputPNG("") { };

    BoundingBox_t *calculatePOI(Grid *lat, Grid *lon, Grid *dem, double *poi);
    void           calculateMinMax(const Grid *dem, const BoundingBox_t *bbox, Grid *newDem, double *min, double *max);
    int           *calculateIndex(Grid *dem, BoundingBox_t *bbox, double poi, double min, double max, int *poi_index);
    bool           plot(int *bin, int poi_index, double min, double max);
    void           parseCoords(const char *coords);
    int            parseOptions(int argc, char **argv);
    void           usage(const char *appname, int errorCode);
    bool           run();

    // all members are public
    string m_dem;
    bool   m_verbose;
    int    m_radius;
    int    m_binSize;
    double m_lat;
    double m_lon;
    string m_outputDEM;
    string m_outputPNG;
};

#endif /* __HISTOGRAM_H */
