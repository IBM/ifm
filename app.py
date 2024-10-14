import os
import json
import glob
import sys
# import time
# import stat
# import logging
# import zipfile
import datetime
# import requests
# import rasterio
# import subprocess
import numpy as np
# import concurrent.futures
# from celery import Celery
# from shapely.geometry import shape
# from concurrent.futures import Executor, ThreadPoolExecutor
# from dispatcher.rabbitmq import RabbitMQ, RabbitMessage
# from dispatcher.configparser import Config
from netCDF4 import Dataset
import utm

########################################################
# Flood model orchestration
# For illustration purposes - how to set up the model
########################################################

def launch_ifm(inputdir, outputdir):

    # SRTM resolution, in meters
    # spatial_resolution = 30

    # Extract the SRTM resolution from the input DEM file, in meters
    f = Dataset(inputdir + '/DEM.nc','r')

    lat0     = f.variables['lat'][:].data[0]; lat1 = f.variables['lat'][:].data[1]
    lon0     = f.variables['lon'][:].data[0]; lon1 = f.variables['lon'][:].data[1]
    llc00    = utm.from_latlon(lat0, lon0); llc01 = utm.from_latlon(lat0, lon1); llc10 = utm.from_latlon(lat1, lon0)
    ll_dists = [abs(llc10[1]-llc00[1]), abs(llc01[0]-llc00[0])]

    spatial_resolution = sum(ll_dists)/len(ll_dists)

    output_rate           = int(os.environ.get('OUTRATE', 3600))  # Output rate, in seconds. We don't want too many output files.
    simulation_timestep   = float(os.environ.get('SIMRATE', 1.0)) # Simulation timestep, default 1.0 second
    initial_soil_moisture = os.environ.get('SMINPUT', '0.0')      # Initial soil moisture value or input filename
    outputLayers          = os.environ.get('OUTLAYERS','depth')   # Comma separated list of output layers - depth,flowrate,olr,vsat,maxvolume

    outputLayers = outputLayers.replace('-',',')
    print(outputLayers, file=sys.stderr)

    # if os.environ.get('OUTRATE') is not None:
    #     output_rate = int(os.environ['OUTRATE'])
    # else:
    #     output_rate = 3600

    # Simulation timestep, in seconds
    # simulation_timestep = 1
    # if os.environ.get('SIMRATE') is not None:
    #     simulation_timestep = int(os.environ['SIMRATE'])
    # else:
    #     simulation_timestep = 1

    # Define uniform values for properties that we don't retrieve
    # dynamically. Please see PARAMETERS.md in IFM for detais.

    # initial_soil_moisture = str(0.0)
    # if os.environ.get('SMINPUT') is not None:
    #     if os.environ.get('SMINPUT') != '':
    #         initial_soil_moisture = os.environ['SMINPUT']

    # Each instance of IFM will attempt to make the best use of CPU
    # resources, so we don't want to launch parallel instances of the
    # model. Here we serialize the calls to IFM by creating a wrapper
    # script that's ultimately executed through subprocess.Popen.
    precipitation_file = inputdir + "/Precipitation.csv"
    # for precipitation_file in pp_files:
    with open(precipitation_file) as f:
        # Check how many simulation seconds will be required. The
        # format of the file is quite simple, with each line being
        # comprised of two elements: '<timestamp> <precipitation_rate>'.
        pp_data = [x for x in f.readlines() if len(x.strip("\n"))]
        t1 = int(pp_data[-1].split(" ")[0])
        if t1 == 0:
            # A single precipitation rate was given. Let IFM input
            # that amount of rainfall for 1 simulation hour.
            simulation_time_in_sec = output_rate
        else:
            # We have an actual time series of precipitation events
            t2 = int(pp_data[-2].split(" ")[0])
            simulation_time_in_sec = t1 + (t1-t2)

    print("Simulation time in seconds     = " + str(simulation_time_in_sec), file=sys.stderr)
    print("Simulation timestep in seconds = " + str(simulation_timestep), file=sys.stderr)
    print("Simulation output rate         = " + str(output_rate), file=sys.stderr)
    print("Initial soil moisture          = " + str(initial_soil_moisture), file=sys.stderr)
    print("Simulation output layers       = " + str(outputLayers), file=sys.stderr)

    wrapper = ""
    ifm_command = [
        "./ifm",
        "-O", f"{outputdir}",
        "-d", f"{inputdir}/DEM.nc",
        "-p", f"{inputdir}/Precipitation.csv",
        # "-f", str(10),  # @note set for debug purposes only
        "-f", str(simulation_time_in_sec),
        "-s", str(simulation_timestep),
        "-r", str(output_rate),
        "-c", str(spatial_resolution),
        "-i", str(1),   # infiltrate flag
        "-t", str(1),   # profile    flag (measure elapsed time)
        "-l", f"{inputdir}/LandCover.nc",
        "-L", f"{inputdir}/LandCover.map",
        "-H", f"{inputdir}/SoilHydraulicConductivity.nc",
        "-P", f"{inputdir}/SoilCapillaryHead.nc",
        "-E", f"{inputdir}/SoilEffectivePorosity.nc",
        "-M", str(initial_soil_moisture),
        "-n", str(0),
        "-V", outputLayers
    ]
    wrapper += " ".join(ifm_command) #+ "\n"

    print(wrapper, file=sys.stderr)
    os.system(wrapper)
    print(wrapper, file=sys.stderr)

def combineCompress(outputPath):
    outputLayers = os.environ.get('OUTLAYERS','depth') # Comma separated list of output layers - depth,flowrate,olr,vsat,maxvolume
    outputLayers = outputLayers.split('-')
    for p in outputLayers:
        print('Combining outputs from ' + p)
        os.system('/combine_outputs.sh ' + p)

if __name__ == '__main__':
    launch_ifm('./input','./output')

    # if os.environ.get('COMBCOMP','Y') == 'Y':
        # combineCompress('/')
