#!/usr/bin/env python
import os
from pathlib import Path

import pytest
import numpy as np

import subprocess

from ..process_meteo_drivers import retrieve_pixel_value
from ..process_meteo_drivers import humidity_from_dewpoint
from ..process_meteo_drivers import grab_meteo_data


DATA_PATH = os.path.dirname(__file__)

def test_retrieve_dem():
    lon = -0.1340
    lat = 51.5246
    url = "/vsicurl/http://www2.geog.ucl.ac.uk/~ucfafyi/eles/global_dem.vrt"
    reference = subprocess.check_output(["gdallocationinfo",  "-valonly", 
                "-wgs84", url, f"{lon:g}", f"{lat:g}"]).decode("utf-8")
    retval = retrieve_pixel_value(lon, lat, url)
    assert abs(float(reference)-float(retval)< 0.1)

def test_humidity_conversion():
    tdew = 298.0
    reference = 3.14 # say
    retval = humidity_from_dewpoint(tdew)
    assert (abs(reference-retval) < 0.1)


def test_meteo_creation():
    ref_cabo = Path(DATA_PATH)/"data/test_cabo_file"
    test_cabo = np.loadtxt(str(ref_cabo), skiprows=20)
    lon = -0.1340
    lat = 51.5246
    year = 2017

    data_dir="tests/data/"
    retval = grab_meteo_data(lat, lon, year, "050_10d_2017.nc", 
                             data_dir=data_dir)

    this_file = np.loadtxt(retval.as_posix(), skiprows=20)
    assert np.allclose(test_cabo[:, 4], this_file[:, 4])
    
