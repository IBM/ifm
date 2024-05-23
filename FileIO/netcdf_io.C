#include <math.h>
#include "Core/engine.h"
#include "FileIO/netcdf_io.h"

Grid *NetCDF::parse(string filename, string varname, BoundingBox_t *cropArea)
{
  assert(filename.size() > 0);

  Grid *grid = NULL;
  NcError err(NcError::silent_nonfatal);
  NcFile *nc = new NcFile(filename.c_str(), NcFile::ReadOnly);
  if (nc->is_valid()) {
    if (varname.size() == 0) {
      // Guess variable name
      NcVar *var = nc->get_var("Band1");
      if (var == NULL)
        var = nc->get_var(0);
      if (var)
        varname.assign((char *) var->name());
    }
    grid = this->parseVar(nc, filename, varname, cropArea);
  }
  if (! grid)
    fprintf(stderr, "Failed to parse file '%s'\n", filename.c_str());
  delete nc;
  return grid;
}

string NetCDF::getAttribute(string filename, string attrname)
{
  string value = "";
  NcFile *nc = new NcFile(filename.c_str(), NcFile::ReadOnly);
  if (nc->is_valid()) {
    NcAtt *att = nc->get_att(attrname.c_str());
    if (att)
      value = string(att->as_string(0));
  }
  delete nc;
  return value;
}

Grid *NetCDF::parseWrfout(string filename, Grid *cropArea)
{
  assert(filename.size() > 0);

  // Parse the main precipitation variable
  Grid *grid = NULL;
  string ppVar(""), latVar(""), lonVar("");
  NcFile *nc = new NcFile(filename.c_str(), NcFile::ReadOnly);
  if (nc->is_valid()) {
    if (this->guessVarnames(nc, &ppVar, &latVar, &lonVar) == true) {
      // Check if we have a single variable (e.g.: APCP) or a combination (e.g.: RAINC+RAINNC)
      if (ppVar.find("+") >= 0) {
        size_t plus = ppVar.find("+");
        string var1 = ppVar.substr(0, plus);
        string var2 = ppVar.substr(plus+1);
        Grid *p1 = this->parseVar(nc, filename, var1, NULL);
        Grid *p2 = this->parseVar(nc, filename, var2, NULL);
        if (p1 && p2) {
          uint64_t size = p1->cols * p1->rows;
          for (uint64_t i=0; i<size; ++i)
            p1->data[i] += p2->data[i];
          grid = p1;
          delete p2;
        } else {
          if (p1) delete p1;
          if (p2) delete p2;
          // fall through
        }
      } else
        // TODO: bounding box
        grid = this->parseVar(nc, filename, ppVar, NULL);
    }
  }
  // Parse the latitude and longitude variables
  Grid *lat = this->parse(filename, latVar);
  Grid *lon = this->parse(filename, lonVar);

  if (!grid || !lat || !lon) {
    fprintf(stderr, "Failed to parse file '%s'\n", filename.c_str());
    if (grid) delete grid;
    if (lat) delete lat;
    if (lon) delete lon;
    delete nc;
    return NULL;
  } else {
    uint64_t size = grid->cols * grid->rows;
    for (unsigned long i=0; i<size; ++i)
      if (grid->data[i] < 0)
        grid->data[i] = 0.0;
  }
#if 0
  if (lonVar.compare("ELON") == 0) {
    uint64_t size = lon->cols * lon->rows;
    for (unsigned long i=0; i<size; ++i)
      lon->data[i] = lon->data[i] > 180 ? lon->data[i]-360 : 360-lon->data[i];
  }
#endif
  if (cropArea) {
    // Crop the precipitation grid to the extents of the bounding box
    Grid *pp = this->resize(grid, lat, lon, cropArea);
    delete grid;
    grid = pp;
  }

  delete lat;
  delete lon;
  delete nc;
  return grid;
}

Grid *NetCDF::parseVar(NcFile *nc, string filename, string varname, BoundingBox_t *cropArea, bool autoFlip)
{
  NcVar *var = nc->get_var(varname.c_str());
  if (! var) {
    fprintf(stderr, "Failed to parse variable '%s' from file '%s'\n", varname.c_str(), filename.c_str());
    return NULL;
  }

  int nodata = -9999;
  NcError err(NcError::silent_nonfatal);
  NcAtt *att = var->get_att("_FillValue");
  if (att) {
    nodata = att->as_int(0);
    delete att;
  }
  err.set_err(NcError::verbose_nonfatal);

  bool success = true;
  bool fixOffset = false;
  real_t *data = NULL;
  uint64_t cols=1, rows=1;
  long *edges = var->edges();
  assert(edges);

  switch (var->num_dims()) {
    case 1:
      if (varname.compare("lat") == 0) {
        fixOffset = true;
        rows = edges[0];
      } else if (varname.compare("lon") == 0) {
        fixOffset = true;
        cols = edges[0];
      } else {
        cols = edges[0];
      }
      data = this->alloc(cols, rows, cropArea);
      success = this->read(nc, var, data, edges[0], cropArea, autoFlip);
      break;
    case 2:
      cols = edges[1];
      rows = edges[0];
      data = this->alloc(cols, rows, cropArea);
      success = this->read(nc, var, data, edges[0], edges[1], cropArea, autoFlip);
      break;
    case 3:
      cols = edges[2];
      rows = edges[1];
      data = this->alloc(cols, rows, cropArea);
      success = this->read(nc, var, data, 1, edges[1], edges[2], cropArea, autoFlip);
      break;
    default:
      fprintf(stderr, "Cannot parse NetCDF variables with dimensions other than 2 and 3.\n");
      fprintf(stderr, "Variable %s from '%s' has %d dimensions.\n", 
        varname.c_str(), filename.c_str(), var->num_dims());
      success = false;
  }
  if (! success) {
    fprintf(stderr, "Failed to parse variable %s from file %s\n", varname.c_str(), filename.c_str());
    delete [] edges;
    return NULL;
  }

  Grid *grid = new Grid;
  grid->data = data;
  grid->cols = cols;
  grid->rows = rows;
  grid->cellsize = cellsize;
  grid->nodata = nodata;
  grid->filename = filename;

  if (fixOffset) {
    // Update the offset. See the operations with adfTempGeoTransform[] on
    // frmts/netcdf/netcdfdataset.cpp (GDAL) for more information on how
    // and why this is performed.
    uint64_t n = edges[0];
    real_t offset = (data[n-1] - data[0]) / (n-1);
    for (uint64_t i=0; i<n; ++i)
      data[i] -= offset / 2.0;
  }
  delete [] edges;
  return grid;
}

bool NetCDF::write(string filename, Grid *grid, string varname, short *mask)
{
  assert(grid != NULL);

  if (mask == NULL) {
    list<Grid*> gridList(1, grid);
    list<string> varnames(1, varname);
    return this->write(filename, gridList, varnames);
  }

  Grid tmpGrid(grid->cols, grid->rows);
  tmpGrid.filename = grid->filename;
  tmpGrid.cellsize = grid->cellsize;
  tmpGrid.nodata = grid->nodata;
  uint64_t size = grid->cols * grid->rows;
  for (uint64_t i=0; i<size; ++i)
    tmpGrid.data[i] = mask[i] ? grid->data[i] : tmpGrid.nodata;

  list<Grid*> gridList(1, &tmpGrid);
  list<string> varnames(1, varname);
  return this->write(filename, gridList, varnames);
}

bool NetCDF::write(string filename, list<Grid*> gridList, list<string> varnames)
{
  NcDim *rowsDim = NULL, *colsDim = NULL, *timeDim = NULL;
  NcVar *data = NULL;
  NcBool retval = true;
  bool isNewFile = false;
  int timeIndex = 0;

  assert(filename.size() > 0);
  assert(gridList.size() > 0);
  assert(gridList.size() == varnames.size());

  NcFile *nc = getHandle(filename, &isNewFile);
  if (! nc)
    return false;

  list<Grid*>::iterator it;
  for (it=gridList.begin(); it != gridList.end(); it++) {
    Grid *grid = *it;
    string varname = varnames.front();
    varnames.pop_front();

    if (isNewFile) {
      timeDim = timeDim ? timeDim : nc->add_dim("Time");
      rowsDim = rowsDim ? rowsDim : nc->add_dim("lat", grid->rows);
      colsDim = colsDim ? colsDim : nc->add_dim("lon", grid->cols);
      data = nc->add_var(varname.c_str(), ncReal, timeDim, rowsDim, colsDim);
      data->add_att("_FillValue", (real_t) grid->nodata);
    } else {
      rowsDim = nc->get_dim("lat");
      colsDim = nc->get_dim("lon");
      timeDim = nc->get_dim("Time");
      data = nc->get_var(varname.c_str());
      if (data) {
        // Replace
        timeIndex = timeDim->size();
      } else {
        // Append
        data = nc->add_var(varname.c_str(), ncReal, timeDim, rowsDim, colsDim);
        data->add_att("_FillValue", (real_t) grid->nodata);
      }
    }

    retval = (NcBool) this->writeGrid(nc, grid, data, timeIndex, false);
    if (! retval) {
      fprintf(stderr, "Failed to write data to the NcFile %s\n", filename.c_str());
      break;
    }
  }
  nc->close();
  delete nc;
  return retval;
}

bool NetCDF::writeGrid(NcFile *nc, Grid *grid, NcVar *data, int timeIndex, bool verticalFlip)
{
  NcError err(NcError::verbose_nonfatal);

  if (! verticalFlip)
    return data->put_rec(grid->data);

  for (int64_t row=(int64_t)grid->rows-1; row>=0; --row) {
    if (! data->set_cur(timeIndex, grid->rows-1-row, 0))
      return false;
    if (! data->put(&grid->data[row*grid->cols], 1, 1, grid->cols))
      return false;
  }
  return true;
}

bool NetCDF::writeHistogram(string filename, Grid *grid, string varname)
{
  assert(filename.size() > 0);
  assert(grid != NULL);

  bool isNewFile = false;
  NcFile *nc = getHandle(filename, &isNewFile);
  if (! nc)
    return false;

  NcDim *buckets = nc->add_dim("buckets", grid->cols);
  NcVar *data = nc->add_var(varname.c_str(), ncReal, buckets);
  data->put(grid->data, grid->cols);
  data->add_att("BinSize", "100");
  data->add_att("BinUnits", "mm");
  nc->close();
  delete nc;
  return true;
}

NcFile *NetCDF::getHandle(string filename, bool *is_new_file)
{
  NcFile *nc = new NcFile(filename.c_str(), NcFile::Write);
  if (! nc->is_valid()) {
    delete nc;
    *is_new_file = true;
    return new NcFile(filename.c_str(), NcFile::Replace);
  }
  if (! nc->is_valid()) {
    fprintf(stderr, "Failed to create a NcFile handle for %s in write mode.\n", filename.c_str());
    delete nc;
    return NULL;
  }
  return nc;
}

bool NetCDF::copyVar(string fromFile, string varname, string toFile)
{
  NcError err(NcError::verbose_nonfatal);
  int in_nc, out_nc;
  int in_err = nc_open(fromFile.c_str(), NC_NOWRITE, &in_nc);
  int out_err = nc_open(toFile.c_str(), NC_WRITE, &out_nc);

  if (in_err != NC_NOERR) {
    fprintf(stderr, "Failed to open file '%s'\n", fromFile.c_str());
    if (out_nc != NC_NOERR)
      nc_close(out_nc);
    return false;
  }
  if (out_err != NC_NOERR) {
    fprintf(stderr, "Failed to open file '%s'\n", toFile.c_str());
    nc_close(in_nc);
    return false;
  }

  bool retval = false;
  int in_var = NC_GLOBAL;
  int out_var = NC_GLOBAL;
  int num_atts = 0;

  // Copy variable
  string ignore_att = "";
  if (varname.size()) {
    in_err = nc_inq_varid(in_nc, varname.c_str(), &in_var);
    if (in_err != NC_NOERR) {
      fprintf(stderr, "Variable '%s' does not exist in '%s'\n", varname.c_str(), fromFile.c_str());
      goto out;
    }

    // Don't copy the history of changes from the source file to the destination
    int num_dims = 0;
    in_err = nc_inq_varndims(in_nc, in_var, &num_dims);
    if (in_err != NC_NOERR) {
      fprintf(stderr, "Failed to identify number of dimensions of '%s' in '%s': %d\n",
        varname.c_str(), fromFile.c_str(), in_err);
      goto out;
    }
    if (num_dims == 0)
      ignore_att = "history";

    int err = nc_copy_var(in_nc, in_var, out_nc);
    if (err != NC_NOERR) {
      fprintf(stderr, "Failed to copy variable '%s' from '%s' to '%s': %d\n", varname.c_str(), fromFile.c_str(), toFile.c_str(), err);
      goto out;
    }
    out_err = nc_inq_varid(out_nc, varname.c_str(), &out_var);
    if (out_err != NC_NOERR) {
      fprintf(stderr, "Variable '%s' does not exist in '%s'\n", varname.c_str(), fromFile.c_str());
      goto out;
    }
  }

  // Copy attributes of the selected variable
  if (nc_inq_natts(in_var, &num_atts) != NC_NOERR && in_var != NC_GLOBAL) {
    if (nc_inq_natts(in_nc, &num_atts) != NC_NOERR) {
      fprintf(stderr, "Failed to retrieve number of attributes from '%s'\n", fromFile.c_str());
      goto out;
    }
  }

  nc_redef(out_nc);
  for (int attnum=0; attnum<num_atts; attnum++) {
    char aname[NC_MAX_NAME];
    if (nc_inq_attname(in_nc, in_var, attnum, aname) != NC_NOERR && in_var != NC_GLOBAL) {
      if (nc_inq_attname(in_nc, NC_GLOBAL, attnum, aname) != NC_NOERR) {
        fprintf(stderr, "Failed to retrieve attribute %d from '%s'\n", attnum, fromFile.c_str());
        goto out_enddef;
      }
      // Update attribute source
      in_var = NC_GLOBAL;
    }

    if (ignore_att.size() == 0 || ignore_att.compare(aname) != 0) {
      int err = nc_copy_att(in_nc, in_var, aname, out_nc, out_var);
      if (err != NC_NOERR) {
        fprintf(stderr, "Failed to copy attribute '%s:%s' from '%s' to '%s': %d\n",
          varname.c_str(), aname, fromFile.c_str(), toFile.c_str(), err);
        goto out_enddef;
      }
    }
  }

  retval = true;

out_enddef:
  nc_enddef(out_nc);
out:
  nc_close(out_nc);
  nc_close(in_nc);
  return retval;
}

bool NetCDF::hasVar(NcFile *nc, string varname)
{
  for (int i=0; i<nc->num_vars(); ++i) {
    NcVar *var = nc->get_var(i);
    if (! strcmp(var->name(), varname.c_str()))
      return true;
  }
  return false;
}

bool NetCDF::guessVarnames(NcFile *nc, string *precipitation, string *lat, string *lon)
{
  // There are two flavors of NetCDF files that we are interested on. The first of them will have the
  // precipitation data coming in two variables called RAINC and RAINNC (and we will want the sum of
  // these two, as they account for the precipitation present in two different kind of clouds).
  // The second will have the precipitation data coming in one variable called APCP.
  // Despite those two, Deep Thunder also used some other different variable names in the past, but
  // including support for those is not necessary nowadays.
  //
  // The latitude and longitude matrices may also come encoded either in variables called XLAT/XLONG or
  // NLAT/ELON. The way that the former is parsed is not the same as the latter.

  list<string> *varnames = new list<string>(nc->num_vars());
  for (int i=0; i<nc->num_vars(); ++i) {
    NcVar *var = nc->get_var(i);
    varnames->push_back(var->name());
  }

  bool has_APCP = false, has_RAINC = false, has_RAINNC = false;
  bool has_SFROFF = false, has_UDROFF = false, has_NLAT = false, has_ELON = false;
  for (list<string>::iterator it = varnames->begin(); it != varnames->end(); ++it) {
    if (! (*it).compare("APCP"))
      has_APCP = true;
    else if (! (*it).compare("RAINC"))
      has_RAINC = true;
    else if (! (*it).compare("RAINNC"))
      has_RAINNC = true;
    else if (! (*it).compare("SFROFF"))
      has_SFROFF = true;
    else if (! (*it).compare("UDROFF"))
      has_UDROFF = true;
    else if (! (*it).compare("NLAT"))
      has_NLAT = true;
    else if (! (*it).compare("ELON"))
      has_ELON = true;
  }
  delete varnames;

  if (has_NLAT && has_ELON) {
    *precipitation = "APCP";
    *lat = "NLAT";
    *lon = "ELON";
  } else if (has_RAINC && has_RAINNC) {
    *precipitation = "RAINC+RAINNC";
    *lat = "XLAT";
    *lon = "XLONG";
  } else if (has_APCP) {
    *precipitation = "APCP";
    *lat = "XLAT";
    *lon = "XLONG";
  } else {
    fprintf(stderr, "Error: could not find any known variables in the provided WRFOUT NetCDF file\n");
    return false;
  }
#ifdef ROUTING_ONLY
  if (has_UDROFF && has_SFROFF) {
    // Take sub-surface runoff as calculated from WRF
    *precipitation = "UDROFF+SFROFF";
  }
#endif
  if (m_dumpFileInfo) {
    fprintf(stderr, "Precipitation=%s, lat=%s, lon=%s\n", precipitation->c_str(), lat->c_str(), lon->c_str());
    m_dumpFileInfo = false;
  }
  return true;
}

Grid *NetCDF::resize(Grid *pre, Grid *preLatitude, Grid *preLongitude, Grid *cropArea)
{
  assert(pre);
  assert(preLatitude);
  assert(preLongitude);
  assert(cropArea);

  // We assume that the DEM NetCDF grid has been produced with gdal_translate,
  // for it creates two auxiliary variables named 'lat' and 'lon' to store that
  // grid's geo-reference.
  Grid *demLatitude  = this->parse(cropArea->filename, "lat");
  Grid *demLongitude = this->parse(cropArea->filename, "lon");
  uint64_t demEndx = demLongitude->cols-1;
  uint64_t demEndy = demLatitude->rows-1;
  assert(demEndx > 0);
  assert(demEndy > 0);

  // Interpolate points between first and last indexes, overwriting original values.
  // This is to enforce the same output between the Serial and MPI implementations.
  real_t latDelta = (GRID(demLatitude, 0, demEndy) - GRID(demLatitude, 0, 0)) / demEndy;
  real_t lonDelta = (GRID(demLongitude, demEndx, 0) - GRID(demLongitude, 0, 0)) / demEndx;
  for (uint64_t i=1; i<cropArea->rows; ++i)
    GRID(demLatitude, 0, i) = GRID(demLatitude, 0, 0) + (latDelta * i);
  for (uint64_t i=1; i<cropArea->cols; ++i)
    GRID(demLongitude, i, 0) = GRID(demLongitude, 0, 0) + (lonDelta * i);

  int latDirection = (GRID(preLatitude, 0, 0) < GRID(demLatitude, 0, 0)) ? -1 : 1;
  int lonDirection = (GRID(preLongitude, 0, 0) < GRID(demLongitude, 0, 0)) ? 1 : -1;
  if (latDirection != 1 || lonDirection != 1) {
    fprintf(stderr, "Error: please verify that the precipitation file covers the extent of the DEM\n");
    delete demLongitude;
    delete demLatitude;
    return NULL;
  }

  int64_t cropStartx = -1;
  int64_t cropStarty = -1;
  int64_t cropEndx   = -1;
  int64_t cropEndy   = -1;

  // Define the range to which the precipitation latitude array (Y axis) will be cropped
  if (latDirection == 1) {
    real_t demStart = GRID(demLatitude, 0, 0);
    real_t demEnd   = GRID(demLatitude, 0, demEndy);
    for (uint64_t i=0; i<preLatitude->rows; ++i) {
      if (GRID(preLatitude, 0, i) <= demStart && cropStarty == -1)
        cropStarty = i == 0 ? 0 : i-1;
      if (GRID(preLatitude, 0, i) <= demEnd && cropEndy == -1) {
        cropEndy = i;
        break;
      }
    }
  }

  // Define the range to which the precipitation longitude array (X axis) will be cropped
  if (lonDirection == 1) {
    real_t demStart = GRID(demLongitude, 0, 0);
    real_t demEnd   = GRID(demLongitude, demEndx, 0);
    for (uint64_t i=preLongitude->cols-1; i>=0; --i) {
      if (GRID(preLongitude, i, 0) <= demEnd && cropEndx == -1)
        cropEndx = i == preLongitude->cols-1 ? preLongitude->cols-1 : i+1;
      if (GRID(preLongitude, i, 0) <= demStart && cropStartx == -1) {
        cropStartx = i;
        break;
      }
    }
  }

  if (cropStartx == -1 || cropEndx == -1 || cropStarty == -1 || cropEndy == -1) {
    fprintf(stderr, "Error: the DEM must be completely within the area covered by the precipitation file(s)\n");
    delete demLongitude;
    delete demLatitude;
    return NULL;
  }

//  printf("Search bounding box: x=%d..%d, y=%d..%d\n", cropStartx, cropEndx, cropStarty, cropEndy);
//  printf("Search bounding box: latitude=%.20f..%.20f\n", GRID(preLatitude, 0, cropStarty), GRID(preLatitude, 0, cropEndy));
//  printf("Search bounding box: longitude=%.20f..%.20f\n", GRID(preLongitude, cropStartx, 0), GRID(preLongitude, cropEndx, 0));

  uint64_t x1 = cropStartx, x2 = cropStartx + 1, y1 = 0, y2 = 0;
  real_t x1val, x2val, y1val, y2val;
  real_t latitude, longitude, offsetx, offsety;

  Grid *out     = new Grid;
  out->cols     = cropArea->cols;
  out->rows     = cropArea->rows;
  out->cellsize = cropArea->cellsize;
  out->nodata   = cropArea->nodata;
  out->filename = pre->filename;
  out->data     = (real_t *) malloc(sizeof(real_t) * cropArea->cols * cropArea->rows);

  for (uint64_t i=0; i<out->cols; ++i) {
    x1val = GRID(preLongitude, x1, 0);
    x2val = GRID(preLongitude, x2, 0);
    longitude = GRID(demLongitude, i, 0);
    if (longitude > x2val) {
      assert(x2+1 < pre->cols);
      x1val = GRID(preLongitude, ++x1, 0);
      x2val = GRID(preLongitude, ++x2, 0);
    }
    if (longitude < x1val || longitude > x2val) {
      fprintf(stderr, "Bad longitude interpolation!\n");
      fprintf(stderr, "DEM: %f (%ld), PRE1: %f (%ld), PRE2: %f (%ld)\n", longitude, i, x1val, x1, x2val, x2);
      delete demLongitude;
      delete demLatitude;
      delete out;
      return NULL;
    }
    offsetx = fabs(longitude - x1val) / fabs(x2val - x1val);

    y1 = cropStarty;
    y2 = cropStarty + 1;
    for (uint64_t j=0; j<out->rows; ++j) {
      y1val = GRID(preLatitude, 0, y1);
      y2val = GRID(preLatitude, 0, y2);
      latitude = GRID(demLatitude, 0, j);
      if (latitude < y2val) {
        assert(y2+1 < pre->rows);
        y1val = GRID(preLatitude, 0, ++y1);
        y2val = GRID(preLatitude, 0, ++y2);
      }
      if (latitude < y2val || latitude > y1val) {
        fprintf(stderr, "Bad latitude interpolation!\n");
        fprintf(stderr, "DEM: %f (%ld), PRE1: %f (%ld), PRE2: %f (%ld)\n", latitude, j, y1val, y1, y2val, y2);
        delete demLongitude;
        delete demLatitude;
        delete out;
        return NULL;
      }
      offsety = fabs(latitude - y1val) / fabs(y2val - y1val);

      // Linear interpolation between the values at columns x1 and x2, at rows y1 and y2
      real_t interpx1 = GRID(pre, x1, y1) * (1-offsetx) + GRID(pre, x2, y1) * offsetx;
      real_t interpx2 = GRID(pre, x1, y2) * (1-offsetx) + GRID(pre, x2, y2) * offsetx;

      // Bilinear interpolation
      real_t interp = (1 - offsety) * interpx1 + offsety * interpx2;
      GRID(out, i, j) = interp;
    }
  }

  delete demLongitude;
  delete demLatitude;
  return out;
}
