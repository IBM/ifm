```
    ____      __                        __           __
   /  _/___  / /____  ____ __________ _/ /____  ____/ /
   / // __ \/ __/ _ \/ __ `/ ___/ __ `/ __/ _ \/ __  /
 _/ // / / / /_/  __/ /_/ / /  / /_/ / /_/  __/ /_/ /
/___/_/ /_/\__/\___/\__, /_/   \__,_/\__/\___/\__,_/
                   /____/
    ________                __   __  ___          __     __
   / ____/ /___  ____  ____/ /  /  |/  /___  ____/ /__  / /
  / /_  / / __ \/ __ \/ __  /  / /|_/ / __ \/ __  / _ \/ /
 / __/ / / /_/ / /_/ / /_/ /  / /  / / /_/ / /_/ /  __/ /
/_/   /_/\____/\____/\__,_/  /_/  /_/\____/\__,_/\___/_/
```
# Introduction

This repository branch contains shared-memory version of the Integrated
Flood Model (IFM) code. IFM was parallelised for shared-memory machines
with OpenMP instructions.

# Dependencies

Shared-memory version of IFM depends on

* OpenMP
* NetCDF C library
* Modern NetCDF C++ library (libnetcdf_c++) v4.9.2 or later

# Compiler Toolchain

For better performance always use hardware vendor's C/C++ compiler and OpenMP,
i.e. on a machine with Intel CPUs use Intel compiler toolchain and Intel
implementation of OpenMP.

# Custom Flags (Build Rules)

To specify custom compiler and linker flags (`USER_CXXFLAGS`, `USER_LDFLAGS`)
create a local text file with your user and hostname under the directory
`BuildRules`.

For example, for a user `mabalenk` on the machine `scafell` create a file
called `mabalenk@scafell.mk`.

Define `USER` and `HOSTNAME` variables in the `Makefile` to pick up your custom
flags. For example set

```
USER     = mabalenk
HOSTNAME = scafell
```

Add custom definitions to `USER_CXXFLAGS`, `USER_LDFLAGS`. For example, set

```
USER_CXXFLAGS = -I/lustre/scafellpike/local/apps/intel/netcdf/4.9.2/include \
                -I${HCBASE}/software/netcdf-cxx-4.2/include

USER_LDFLAGS  = -L/lustre/scafellpike/local/apps/intel/netcdf/4.9.2/lib \
                -L${HCBASE}/software/netcdf-cxx-4.2/lib
```

in `BuildRules/mabalenk@scafell.mk`.

# Modules

Load the necessary modules, if you intend to run the code on a supercomputer.
For example, on Scafell Pike

```
load use.dev gcc12 netcdf

```

# Modern NetCDF C++ library

Currently Scafell Pike has no modern NetCDF C++ library module. This software
needs to be downloaded from an official [NetCDF
website](https://downloads.unidata.ucar.edu/netcdf) and installed from source
into your home directory.

Add the location of modern NetCDF C++ header and library files to custom
`CXXFLAGS` and `LDFLAGS` variable in the custom make rules file, e.g.

```
USER_CXXFLAGS = -I/lustre/scafellpike/local/apps/intel/netcdf/4.9.2/include \
                -I${HCBASE}/software/netcdf-cxx-4.2/include

USER_LDFLAGS  = -L/lustre/scafellpike/local/apps/intel/netcdf/4.9.2/lib \
                -L${HCBASE}/software/netcdf-cxx-4.2/lib
```

# Compilation

Compile the IFM source code with:

```
make clean && make
```
