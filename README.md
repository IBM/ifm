# Integrated Flood Model (IFM)

# Introduction

The Integrated Flood Model (IFM) is a 2D diffusiion wave-based flood model that calculates water depth across a 2D structured grid as a result of rainfall events. More details are given in [Singhal et al., 2014](https://doi.org/10.1007/978-3-319-09873-9_58). The IFM consists of two components: a soil and an overland routing model. The soil model estimates the surface-runoff based on the incoming precipitation, soil, and land use properties. These runoff estimates are then input to the overland flood routing engine, which calculates the water in-flows and out-flows on a two-dimensional grid based on topological characteristics. The remnant water-flow from a simulation step is then fed back to the soil model to more accurately determine the water height in the next simulation step. This release contains an example set of test input files which are described below.

This repository branch contains shared-memory version of the Integrated Flood Model (IFM) code. IFM was parallelised for shared-memory machines with OpenMP instructions.

# Dependencies

Shared-memory version of IFM depends on

* [OpenMP](https://www.openmp.org)
* [NetCDF C library](https://downloads.unidata.ucar.edu/netcdf)
* [Legacy NetCDF C++ library](https://downloads.unidata.ucar.edu/netcdf) ([netcdf-cxx 4.2.17](https://downloads.unidata.ucar.edu/netcdf-cxx/4.2/netcdf-cxx-4.2.tar.gz))


 Compiler Toolchain

For better performance always use hardware vendor's C/C++ compiler and OpenMP,
i.e. on a machine with Intel CPUs use Intel compiler toolchain and Intel
implementation of OpenMP.

# Custom Flags (Build Rules)

To specify custom compiler and linker flags (`USER_CXXFLAGS`, `USER_LDFLAGS`)
create a local text file with your user and hostname under the directory
`BuildRules`.

For example, for a user `user` on the machine `machine` create a file
called `user@machine.mk`.

Define `USER` and `HOSTNAME` variables in the `Makefile` to pick up your custom
flags. For example set

```
USER     = user
HOSTNAME = machine
```

Add custom definitions to `USER_CXXFLAGS`, `USER_LDFLAGS`. For example, set

```
USER_CXXFLAGS = -I/user/local/apps/intel/netcdf/4.9.2/include \
                -I${HCBASE}/software/netcdf-cxx-4.2/include

USER_LDFLAGS  = -L/user/local/apps/intel/netcdf/4.9.2/lib \
                -L${HCBASE}/software/netcdf-cxx-4.2/lib
```

in `BuildRules/user@machine.mk`.


# Compilation

Compile the IFM source code with:

```
make clean && make
```

# Testing

```
./ifm -Oexp/output/ -dexp/input/DEM.nc -pexp/input/Precipitation.csv -f10800 -s60 -r120 -c50 -lexp/input/LandCover.nc -Lexp/input/LandCover.map -Hexp/input/SoilHydraulicConductivity.nc -Pexp/input/SoilCapillaryHead.nc  -Eexp/input/SoilEffectivePorosity.nc -M0.3 -n0 -Vdepth
```

Input Data:
- DEM.nc: Digital elevation model grid, m
- Landcover.nc: A landcover grid with each value corresponding to the Darcy flow coeffiecint and water retention described in landcover.map
- Landcover.map: A text file with the colums: LandCover value, Darcy flow coeffiecint, water retention
- Precipitation.csv: csv file with the colums: time (in seconds), precipitation grid for that time, mm
- SoilCapillaryHead.nc: wetting front suction head (capillary head), cm
- SoilHydraulicConductivity.nc: saturated hydraulic conductivity, cm/h
- SoilEffectivePorosity.nc: effective porosity/saturation, cm3/cm3

All grid are to supplied in the netcdf format with the fields
- lat: (of size image height)
- lon: (of size image width)
- Band1: (input variable of size lat x lon)
- crs: with variables: longitude_of_prime_meridian, semi_major_axis, inverse_flattening

## Inputs

- c: DEM cell size, in meters
- d: DEM file
- h: help
- l: Landuse file or uniform value
- L: Landuse file (optional)
- m: Mask file
- n: Extra NoData value. Can be entered multiple times (default: take from DEM)
- O: Output directory
- o: Outlet file (optional)"
- r: How often to output depth files
- p: Precipitation file
- S: Save/restore simulation state from the given file
- E: Soil effective porosity file or uniform value
- H: Soil hydraulic conductivity file or uniform value
- P: Soil pressure head file or uniform value
- 1: Multiply soil-ep grid points by VALUE (default: 1.0)
- 2: Multiply soil-hc grid points by VALUE (default: 1.0)
- 3: Multiply soil-ph grid points by VALUE (default: 1.0)
- M: Soil moisture file or uniform value (default: 0.0)
- f: Simulation end time
- b: Simulation start time
- s: Time step in seconds
- w: Optional list of coordinates to watch
- V: Comma-separated list of variables to save. Accepted arguments are: depth, flowrate, olr, vsat, maxvolume
