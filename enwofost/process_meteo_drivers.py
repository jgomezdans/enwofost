#!/usr/bin/env python
"""Some convenience files to grab meteorological data and convert
to CABO format to use within WOFOST. So far, using ERA5
"""
import struct
import logging 

import datetime as dt

from textwrap import dedent

from collections import namedtuple

from pathlib import Path

import numpy as np

from osgeo import gdal

from netCDF4 import Dataset, date2index

import cdsapi


ERAPARAMS = namedtuple("ERAPARAMS", ["ssrd", "mx2t", "mn2t", "tp", 
                                     "u10", "v10", "d2m"])
LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)
if not LOG.handlers:
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - ' +
                                  '%(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    LOG.addHandler(ch)
LOG.propagate = False


def humidity_from_dewpoint(tdew):
    """Calculates humidity from dewpoint temperature
    
    Parameters
    ----------
    tdew : float
        Dewpoint temperature in degrees Kelvin
    
    Returns
    -------
    Relative humidity
    """
    tdew = tdew - 273.15
    tmp = (17.27 * tdew) / (tdew + 237.3)
    ea = 0.6108 * np.exp(tmp)
    return ea


def retrieve_pixel_value(lon, lat, data_source):
    """Retrieve pixel value from a GDAL-friendly dataset.

    We assume the data type of the raster here!!!!
    
    Parameters
    ----------
    lon : float
        Longitude in decimal degrees
    lat : float
        Latitude in decimal degrees
    data_source : str
        An existing GDAL-readable dataset. Can be remote.
    
    Returns
    -------
    int
       The value of the pixel.
    """
    dataset = gdal.Open(data_source)

    gt = dataset.GetGeoTransform()
    the_band = dataset.GetRasterBand(1)
    px = int((lon - gt[0]) / gt[1]) #x pixel
    py = int((lat - gt[3]) / gt[5]) #y pixel

    buf = the_band.ReadRaster(px, py, 1, 1, 
                              buf_type=gdal.GDT_Int16)
    elev = struct.unpack('h' , buf)

    return elev[0]

def grab_era5(output_fname, year, mylat, mylon):
    """A function to downlaod the ERA5 data for one year for a given location.
    Note that this takes a while! Also, you need to have the Copernicus Data
    Service API configured for this to work.
    
    Parameters
    ----------
    output_fname : str
        The output file name. This will be a netCDF file
    year : int
        The year of interest. Must lie within the ERA5 years!
    mylat : float
        Latitude in decimal degrees.
    mylon : float
        Longitude in decimal degrees.
    """
    #[60, -10, 50, 2], # North, West, South, East
    area = f"{int(mylat)+1:d}/{int(mylon)-1:d}/" + \
           f"{int(mylat)-1:d}/{int(mylon)+1:d}"
    c = cdsapi.Client()                
    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'variable':[
                'surface_solar_radiation_downwards',
                'maximum_2m_temperature_since_previous_post_processing','minimum_2m_temperature_since_previous_post_processing',
                'total_precipitation',
                '10m_u_component_of_wind','10m_v_component_of_wind',
                '2m_dewpoint_temperature'                            
            ],
            'product_type':'reanalysis',
            'year':year,
            'month':[
                '01','02','03',
                '04','05','06',
                '07','08','09',
                '10','11','12'
            ],
            'day':[
                '01','02','03',
                '04','05','06',
                '07','08','09',
                '10','11','12',
                '13','14','15',
                '16','17','18',
                '19','20','21',
                '22','23','24',
                '25','26','27',
                '28','29','30',
                '31'
            ],
            'time':[
                '00:00','01:00','02:00',
                '03:00','04:00','05:00',
                '06:00','07:00','08:00',
                '09:00','10:00','11:00',
                '12:00','13:00','14:00',
                '15:00','16:00','17:00',
                '18:00','19:00','20:00',
                '21:00','22:00','23:00'
            ],
            #'area' : "%d/%d/%d/%d"%
            # (int(mylat/10.+1.)*10,
            # int(mylon/10.)*10,
            # int(mylat/10.)*10,
            # int(mylon/10.+1.)*10,),
            'area' : area,
            'format':'netcdf'
        },
        output_fname)


def grab_meteo_data(lat, lon, start_year, end_year, data_dir="./",
                    size=0.25, c1=-0.18, c2=-0.55, station_number=1,
                    dem_file="/vsicurl/http://www2.geog.ucl.ac.uk/" + 
                             "~ucfafyi/eles/global_dem.vrt"):
    """Grab meteorological data and set it up to use with WOFOST.
    At present, we download the data from ERA5, but other sources may
    be considered.

    Parameters
    ----------
    lat : float
        Latitude in decimal degrees.
    lon : float
        Longitude in decimal degrees.
    start_year : int
        The starting year
    end_year : int
        Tne end year (inclusive)
    data_dir : str
        The location where the ERA5 file will be downloaded to, and 
        also where the CABO file will be written to
    size : float
        The grid size in degrees
    c1 : float
        I think this is a parameter related to the Armstrong exponent
    c2 : float
        I think this is a parameter related to the Armstrong exponent
    station_number : int
        A random number needed for CABO file
    dem_file : str
        A GDAL-readable file with a DEM to fish out the elevation of the site.

    Returns
    -------
        dict
        A dictionary with the relevant CABO files indexed by year.
    """    
    # This is the site name. Use the longitude/latitude to make it unique
    site_name = "%05.2f_%05.2f" % (int((lon+size/2.)/size)*size,
                                 int((lat+size/2.)/size)*size)
    # Grab the elevation
    elevation = retrieve_pixel_value(lon, lat, dem_file)
    # These are the parameters
    parnames = ["ssrd", "mx2t", "mn2t", "tp", "u10", "v10", "d2m"]
    return_files = {}
    
    for year in range(start_year, end_year+1):
        cabo_file = Path(data_dir)/f"{site_name:s}.{year:d}"
        if not cabo_file.exists():
            LOG.info(f"No CABO file for {year:d}...")
            nc_file = Path(data_dir)/(f"{int(lon/10.)*10:d}"+
                        f"{int(lat/10.)*10:d}_10d_{year:d}.nc")
            if not nc_file.exists():
                LOG.info(f"No netCDF file for {year:d}...")
                LOG.info(f"Starting downloading...")
                grab_era5(str(nc_file), year, lat, lon)
                LOG.info(f"Done downloading...")
            LOG.info(f"Converting {str(nc_file):s} to CABO")
            LOG.info("Converting units to daily etc.")
            # Open netCDF file, and stuff parameters into useful 
            # data structure
            ds = Dataset(str(nc_file))
            variables = (ds.variables[var][:] for var in parnames)
            pars = ERAPARAMS(*variables)
            # Check corners
            uplat = ds.variables["latitude"][:].max()
            dnlat = ds.variables["latitude"][:].min()
            uplon = ds.variables["longitude"][:].max()
            dnlon = ds.variables["longitude"][:].min()
            x = int((lon-dnlon+size/2)/size)
            y = int((lat-uplat-size/2)/-size)
            times = ds.variables["time"]
            
            # Preprocess data: calculate daily means/aggregates
            # Get the right units.
            rad = np.sum(pars.ssrd.reshape(-1, 24,
                         pars.ssrd.shape[1], 
                         pars.ssrd.shape[2]), axis=1)/1000.
            tmax = np.max(pars.mx2t.reshape(-1, 24,
                         pars.mx2t.shape[1],
                         pars.mx2t.shape[2]), axis=1)-273.15
            tmin = np.min(pars.mn2t.reshape(-1, 24,
                         pars.mn2t.shape[1],
                         pars.mn2t.shape[2]), axis=1)-273.15
            prec = np.sum(pars.tp.reshape(-1, 24,
                          pars.tp.shape[1],
                          pars.tp.shape[2]), axis=1)*1000.
            prec[prec < 0.01] = 0
            wind_u = np.mean(pars.u10.reshape(-1,24,
                             pars.u10.shape[1],
                             pars.u10.shape[2]),axis=1)
            wind_v = np.mean(pars.v10.reshape(-1,24,
                             pars.v10.shape[1],
                             pars.v10.shape[2]),axis=1)
            wind = np.sqrt(np.square(wind_u)+np.square(wind_v))
            hum = humidity_from_dewpoint(np.mean(pars.d2m.reshape(-1, 24,
                                        pars.d2m.shape[1],
                                        pars.d2m.shape[2]), axis=1))
            hdr_chunk = f"""\
                *---------------------------------------------------
                * Station: {site_name:s}
                * Year: {year:d}
                * Origin: ERA5-Reanalysis
                * Columns & units
                * ===================
                * 1. station number
                * 2. year
                * 3. Day of Year
                * 4. Irradiance   (kJ·m-2·d-1)
                * 5. Daily minimum temperature (degC)
                * 6. Daily maximum temperature (degC)
                * 7. Vapour pressure (kPa)
                * 8. Mean wind speed (m·s-1)
                * 9. Precipitation (mm·d-1)
                ** WCCDESCRIPTION={site_name:s}
                ** WCCFORMAT=2
                ** WCCYEARNR={year:d}
                """
            hdr_chunk = dedent(hdr_chunk)
            year_s = date2index(dt.datetime.strptime
                                (f'{year:d}-01-01 00:00:00',
                                '%Y-%m-%d %H:%M:%S'), times)
            year_e = date2index(dt.datetime.strptime
                                (f'{year:d}-12-31 00:00:00',
                                '%Y-%m-%d %H:%M:%S'), times)
            # Dump data file...
            with cabo_file.open("w") as fp:
                fp.write(hdr_chunk)
                for d in range(rad.shape[0]):
                    fp.write(f"{station_number:d}\t{year:d}\t{d+1:d}\t" + 
                            f"{round(rad[d,y,x]):5.1f}\t" + 
                            f"{round(tmin[d,y,x]*10/10):5.1f}\t" +
                            f"{round(tmax[d,y,x]*10/10):5.1f}\t" +
                            f"{round(hum[d,y,x]*1000/1000):5.3f}\t" +
                            f"{round(wind[d,y,x]*10/10):4.1f}\t" +
                            f"{round(prec[d,y,x]*10/10):4.1f}\n"
                            )
            LOG.info(f"Saved CABO file {str(cabo_file):s}.")
        return_files[year] = cabo_file
    return return_files

if __name__ == "__main__":
    lon = -0.1340
    lat = 51.5246
    start_year = 2018
    end_year = 2018
    data_dir="/home/ucfajlg/temp/"
    retval = grab_meteo_data(lat, lon, start_year, end_year, 
                data_dir)
    print(retval)