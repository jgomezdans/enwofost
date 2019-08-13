#!/usr/bin/env python


import pytest
import numpy as np

import subprocess

from ..process_meteo_drivers import retrieve_pixel_value
from ..process_meteo_drivers import humidity_from_dewpoint

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