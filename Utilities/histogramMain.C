#include "FileIO/netcdf_io.h"
#include "Core/engine.h"
#include "histogram.h"

int main(int argc, char **argv)
{
  Histogram hist;
  int ret = hist.parseOptions(argc, argv);
  if (ret < 0)
    return ret;
  return hist.run() == true ? 0 : 1;
}
