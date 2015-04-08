#!/usr/bin/env python
# Copyright (C) 2015 Andy Aschwanden
#



"""This converts a matlab file into a gtiff.
"""

from argparse import ArgumentParser
import numpy as np
from scipy.io import loadmat
from pyproj import Proj


def write_file(fname, array, rasterOrigin, pixelWidth, pixelHeight, _FillValue, proj4str):

    import gdal
    import osr

    cols = array.shape[1]
    rows = array.shape[0]

    originX, originY = rasterOrigin

    driver = gdal.GetDriverByName('netCDF')
    outRaster = driver.Create(fname, cols, rows, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, -pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.SetNoDataValue(_FillValue)
    outband.WriteArray(np.flipud(array))
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromProj4(proj4str)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()
    
if __name__ == "__main__":
    # Set up the option parser
    description = '''A script to extract data along (possibly multiple) profile using
    piece-wise constant or bilinear interpolation.
    The profile must be given as a ESRI shape file.'''
    parser = ArgumentParser()
    parser.description = description
    parser.add_argument("--log",dest="log", action="store_true",
                      help="Write log10 of data", default=False)

    parser.add_argument("INPUTFILE", nargs=1, help="input NetCDF file name")
    parser.add_argument("OUTPUTFILE", nargs=1, help="output NetCDF file name", default="out.nc")

    options = parser.parse_args()
    log = options.log
    indata = loadmat(options.INPUTFILE[0])

    LAT = indata['Apdem100_lat']
    LON = indata['Apdem100_lon']
    DATA = indata['Flux']
    mask = np.zeros_like(DATA)
    mask[DATA==1] = 1
    fill_value = -2e9
    data = np.ma.array(data=DATA, mask=mask, fill_value=fill_value) 
    if log:
        data = np.log(data)
    
    p = Proj(init='EPSG:3976')

    origin = p(LAT[-1,-1], LON[-1,-1])
    origin = p(LON[0,0], LAT[0,0])
    # guessing
    origin = [-2559236, 1695825]
    proj4str = '+init=epsg:3031'
    write_file(options.OUTPUTFILE[0], np.flipud(data), origin, 100, 100, 1, proj4str)
