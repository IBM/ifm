#include <getopt.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

class OptionParser {
  private:
    int _argc;
    char **_argv;
  public:
    string dem;
    string mask;
    string landuse;
    string landuseMap;
    string soil_hc; // Hydraulic Conductivity
    string soil_ph; // Pressure Head
    string soil_ep; // Effective Porosity
    string soilMoisture;
    string outlet;
    string outdir;
    string watch;   // Points of interest to watch
    string precipitation;
    string simulationState; // File holding the simulation state
    string outputVars;
    vector<real_t> nodata;
    int    cellsize;        // Grid cell size, in meters
    int    tbegin;
    int    tfinal;
    int    outputRate;      // How often to output files (in seconds)
    double soil_hc_multiplier;
    double soil_ph_multiplier;
    double soil_ep_multiplier;
    double tstep;
    bool   saveFlowRate;    // Output flow rate?
    bool   saveDepth;       // Output water height?
    bool   saveVolume;      // Output accumulated water volume?
    bool   saveOLR;         // Output OLR?
    bool   saveVSAT;        // Output VSAT?

    OptionParser(int argc, char **argv)
    {
      _argc = argc;
      _argv = argv;
      dem = "";
      mask = "";
      landuse = "";
      landuseMap = "";
      nodata.clear();
      soil_hc = "";
      soil_ph = "";
      soil_ep = "";
      soil_ep_multiplier = 1.0;
      soil_hc_multiplier = 1.0;
      soil_ph_multiplier = 1.0;
      soilMoisture = "0";
      outlet = "";
      outdir = ".";
      watch = "";
      precipitation = "";
      simulationState = "";
      cellsize = 90;
      tbegin = 0;
      tfinal = 1800;
      tstep = 1.0;
      outputRate = 300;
      saveFlowRate = false;
      saveVolume = false;
      saveOLR = false;
      saveVSAT = false;
      saveDepth  = true;
      outputVars = "depth";
    }

    ~OptionParser()
    {
    }

    void parse()
    {
      const char *short_options = "b:c:d:f:hl:L:m:M:n:o:O:p:r:s:w:1:2:3:H:E:P:S:V:";
      const struct option long_options[] = {
        { "cellsize",      required_argument, 0, 'c' },
        { "dem",           required_argument, 0, 'd' },
        { "help",          no_argument,       0, 'h' },
        { "landuse",       required_argument, 0, 'l' },
        { "landuse-map",   required_argument, 0, 'L' },
        { "mask",          required_argument, 0, 'm' },
        { "nodata",        required_argument, 0, 'n' },
        { "outdir",        required_argument, 0, 'O' },
        { "outlet",        required_argument, 0, 'o' },
        { "output-rate",   required_argument, 0, 'r' },
        { "precipitation", required_argument, 0, 'p' },
        { "state",         required_argument, 0, 'S' },
        { "soil-ep",       required_argument, 0, 'E' },
        { "soil-hc",       required_argument, 0, 'H' },
        { "soil-ph",       required_argument, 0, 'P' },
        { "soil-ep-multiplier", required_argument, 0, '1' },
        { "soil-hc-multiplier", required_argument, 0, '2' },
        { "soil-ph-multiplier", required_argument, 0, '3' },
        { "soil-moisture", required_argument, 0, 'M' },
        { "tfinal",        required_argument, 0, 'f' },
        { "tbegin",        required_argument, 0, 'b' },
        { "tstep",         required_argument, 0, 's' },
        { "watch",         required_argument, 0, 'w' },
        { "output-vars",   required_argument, 0, 'V' },
        { NULL,            0,                 0,  0  }
      };

      while (2) {
        int option_index = 0;
        int c = getopt_long(_argc, _argv, short_options, long_options, &option_index);
        if (c == -1)
          break;
        switch (c) {
          case 'b':
            tbegin = atoi(optarg);
            break;
          case 'c':
            cellsize = atoi(optarg);
            break;
          case 'd':
            dem.assign(optarg);
            break;
          case 'h':
            usage(0);
            break;
          case 'l':
            landuse.assign(optarg);
            break;
          case 'L':
            landuseMap.assign(optarg);
            break;
          case 'm':
            mask.assign(optarg);
            break;
          case 'n':
            nodata.push_back((real_t) atof(optarg));
            break;
          case 'o':
            outlet.assign(optarg);
            break;
          case 'O':
            outdir.assign(optarg);
            break;
          case 'p':
            precipitation.assign(optarg);
            break;
          case 'f':
            tfinal = atoi(optarg);
            break;
          case 'r':
            outputRate = atoi(optarg);
            break;
          case 's':
            tstep = atof(optarg);
            break;
          case 'w':
            watch.assign(optarg);
            break;
          case '1':
            soil_ep_multiplier = atof(optarg);
            break;
          case '2':
            soil_hc_multiplier = atof(optarg);
            break;
          case '3':
            soil_ph_multiplier = atof(optarg);
            break;
          case 'H':
            soil_hc.assign(optarg);
            break;
          case 'E':
            soil_ep.assign(optarg);
            break;
          case 'P':
            soil_ph.assign(optarg);
            break;
          case 'M':
            soilMoisture.assign(optarg);
            break;
          case 'S':
            simulationState.assign(optarg);
            break;
          case 'V':
            outputVars.assign(optarg);
            break;
          case '?':
          default:
            break;
        }
      }

      if (!dem.size() || !landuse.size() || !precipitation.size()) {
        cerr << "Error: missing arguments" << endl << endl;
        usage(1);
      }

      if (! ((!soil_hc.size() && !soil_ep.size() && !soil_ph.size()) ||
             (soil_hc.size() && soil_ep.size() && soil_ph.size()))) {
        cerr << "Error: missing arguments" << endl << endl;
        cerr << "You must either provide options for all soil properties or no options at all" << endl << endl;
        usage(1);
      }

      saveOLR      = outputVars.find("olr") != string::npos;
      saveVSAT     = outputVars.find("vsat") != string::npos;
      saveDepth    = outputVars.find("depth") != string::npos;
      saveFlowRate = outputVars.find("flowrate") != string::npos;
      saveVolume   = outputVars.find("maxvolume") != string::npos;
      if (! saveDepth && ! saveFlowRate && ! saveVolume) {
        cerr << "Error: neither 'depth', 'flowrate', 'olr' nor 'maxvolume' outputs have been selected" << endl;
        usage(1);
      }
    }

    void usage(int retval)
    {
      cout << "Usage: " << _argv[0] << " <options>" << endl;
      cout << "GENERAL options:" << endl;
      cout << "    -h, --help                    This help." << endl;
      cout << endl;

      cout << "INPUT options:" << endl;
      cout << "    -c, --cellsize=SIZE           DEM cell size, in meters (default: " << cellsize << ")" << endl;
      cout << "    -d, --dem=FILE                DEM file" << endl;
      cout << "    -l, --landuse=<FILE|VALUE>    Landuse file or uniform value" << endl;
      cout << "    -L, --landuse-map=FILE        Landuse file (optional)" << endl;
      cout << "    -m, --mask=FILE               Mask file" << endl;
      cout << "    -n, --nodata=VALUE            Extra NoData value. Can be entered multiple times (default: take from DEM)" << endl;
      cout << "    -o, --outlet=FILE             Outlet file (optional)" << endl;
      cout << "    -O, --outdir=DIR              Output directory (default: " << outdir << "/)" << endl;
      cout << "    -p, --precipitation=FILE      Precipitation file" << endl;
      cout << "    -w, --watch=FILE              Optional list of coordinates to watch" << endl;
#ifndef ROUTING_ONLY
      cout << "    -H, --soil-hc=<FILE|VALUE>    Soil hydraulic conductivity file or uniform value" << endl;
      cout << "    -P, --soil-ph=<FILE|VALUE>    Soil pressure head file or uniform value" << endl;
      cout << "    -E, --soil-ep=<FILE|VALUE>    Soil effective porosity file or uniform value" << endl;
      cout << "    -M, --soil-moisture=<FILE|VALUE> Soil moisture file or uniform value (default: 0.0)" << endl;
#endif
      cout << endl;
#ifndef ROUTING_ONLY
      cout << "MODIFIERS:" << endl;
      cout << "    -1, --soil-ep-multiplier=VALUE   Multiply soil-ep grid points by VALUE (default: 1.0)" << endl;
      cout << "    -2, --soil-hc-multiplier=VALUE   Multiply soil-hc grid points by VALUE (default: 1.0)" << endl;
      cout << "    -3, --soil-ph-multiplier=VALUE   Multiply soil-ph grid points by VALUE (default: 1.0)" << endl;
      cout << endl;
#endif

      cout << "OUTPUT options:" << endl;
      cout << "    -r, --output-rate=NUM_SECONDS How often to output depth files (default: " << outputRate << ")" << endl;
      cout << "    -V, --output-vars=VARIABLES   Comma-separated list of variables to save. Accepted arguments are:" << endl;
      cout << "                                  depth, flowrate, olr, vsat, maxvolume. (default: " << outputVars << ")" << endl;
      cout << endl;

      cout << "SIMULATION options:" << endl;
      cout << "    -b, --tbegin=NUM_SECONDS      Simulation start time (default: " << tbegin << ")" << endl;
      cout << "    -f, --tfinal=NUM_SECONDS      Simulation end time (default: " << tfinal << ")" << endl;
      cout << "    -s, --tstep=NUM_SECONDS       Time step in seconds (default: " << tstep << ")" << endl;
      cout << "    -S, --state=FILE              Save/restore simulation state from the given file" << endl;
      cout << endl;

      cout << "UNIFORM precipitation file example:" << endl;
      cout << "0    100" << endl;
      cout << "3599 100" << endl;
      cout << "3600 0" << endl;
      cout << endl;

      cout << "NETCDF precipitation file example:" << endl;
      cout << "0       wrfout_d04_2010-04-07_00:00:00" << endl;
      cout << "600     wrfout_d04_2010-04-07_00:10:00" << endl;
      cout << "1200    wrfout_d04_2010-04-07_00:20:00" << endl;
      cout << endl;

      cout << "OUTLET file example:" << endl;
      cout << "1 1 0.003539" << endl;
      cout << endl;

      cout << "WATCH file example:" << endl;
      cout << "# Lat Lon : Location name" << endl;
      cout << "-81.584575 28.413437 : Walt Disney World" << endl;
      cout << endl;
      exit(retval);
    }
};
