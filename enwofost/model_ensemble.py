#!/usr/bin/env python
"""Creating ensembles of a crop growth model.
"""
import logging
import datetime as dt
import copy


from pathlib import Path
from .process_meteo_drivers import grab_meteo_data

import numpy as np

import scipy.stats as ss

from bisect import bisect


from pcse.fileinput import CABOFileReader
from pcse.base_classes import ParameterProvider
from pcse.util import WOFOST71SiteDataProvider
from pcse.fileinput import YAMLAgroManagementReader
from pcse.fileinput import CABOWeatherDataProvider
from pcse.models import Wofost71_PP

from pcse.db import NASAPowerWeatherDataProvider

LOG = logging.getLogger(__name__)



def define_prior_distributions(chunk,tsum1=None,tsum2=None):
    """A function to interpret the prior distributions from an Excel c&p job.
    Returns a dictionary indexed by parameter name and a pointer to a
    scipy.stats function that provides a logpdf method ;-). It also returns a
    list that can be used to map from vector to dictionary."""
    prior_dist = {}
    param_list = []
    param_xvalue = {}
    param_type = {}
    if Path(chunk).exists():
        # we have a file
        with Path(chunk).open('r') as fp:
            chunk_data = fp.readlines()
    else:
        chunk_data = chunk.split("\n")
        
    for line in chunk_data:
        if not line.strip().startswith("#"):
            if line.find("Uniform") >= 0:
                param,ty, xp, yp, xmin, xmax, sigma, _ = line.split(",")
                xmin = float(xmin)
                xmax = float(xmax)               
                dist = ss.uniform(xmin, (xmax - xmin))
            elif line.find("Gaussian") >= 0:
                param,ty, xp, yp, xmin, xmax, sigma, _ = line.split(",")
                if param =="TSUM1" and tsum1 is not None:  yp = tsum1
                if param =="TSUM2" and tsum2 is not None:  yp = tsum2
                if ty[0] == "X":
                    xp,yp = yp,xp
                lower = float(xmin)
                upper = float(xmax)
                mu = float(yp)
                sigma = float(sigma)
                dist = ss.truncnorm((lower - mu) / sigma,
                                    (upper - mu) / sigma, loc=mu, scale=sigma)            
            prior_dist[param] = dist
            param_list.append(param)
            param_xvalue[param] = float(xp)
            param_type[param] = ty
    return prior_dist, param_list, param_xvalue,param_type


def set_ensemble_parameters(en_size, prior_list, 
                            param_type, param_xvalue,z_start, crop):
    for i in range(en_size):
        theta_dict = dict(zip(prior_list, z_start[:,i]))
        cropdata = copy.deepcopy(crop)
        tb_x = {}
        tb_y = {}
        tb_t = {}
        tmp_dict={}
        for par in theta_dict.keys():
            try:
                if param_type[par] != 'S':
                    tb_index = par.find("TB")
                    if tb_index < 0:
                        #print(param_xvparam_typealue[par])
                        raise Exception("Are you sure %s is a table value?"%par)
                    tb_name = par[:tb_index+2]
                    tmp_list = [param_xvalue[par],theta_dict[par]]
                    if not tb_name in tb_x:
                        tb_x[tb_name] = np.array([param_xvalue[par]])
                        tb_y[tb_name] = np.array([theta_dict[par]])
                        tb_t[tb_name] = param_type[par]
                    else:
                        tb_x[tb_name] = np.append(tb_x[tb_name],param_xvalue[par])
                        tb_y[tb_name] = np.append(tb_y[tb_name],theta_dict[par])
            except KeyError:
                raise Exception("There's something wrong with %s, please check it."%par)
        tmp_dict={}
        for par in tb_x.keys():  # Table parameters
            s_i = np.argsort(tb_x[par])
            s_x = tb_x[par][s_i]
            s_v = tb_y[par][s_i]
            par_tb = []
    #         print(par,tb_t[par],cropdata[par],s_x,s_v)
            if tb_t[par][1] == 'P':   
                for i in range(len(tb_x[par])):
                    if tb_t[par][0] == 'Y':               # Partly change table Y values
                        if s_x[i] in cropdata[par][::2]:  # change old value
                            c_i = cropdata[par][::2].index(s_x[i])
                            cropdata[par][c_i*2] = s_v[i]
                        else:                             # insert new value
                            array_X = cropdata[par][::2]
                            array_Y = cropdata[par][1:][::2]
                            ins_i = bisect(array_X, s_x[i])
                            cropdata[par].insert( ins_i*2, s_x[i])
                            cropdata[par].insert( ins_i*2+1, s_v[i])
                        #print(cropdata[par])
                    else:                                 # Partly change table X values
                        if s_x[i] in cropdata[par][1:][::2]:  # change old value
                            c_i = cropdata[par][1:][::2].index(s_x[i])
                            cropdata[par][c_i*2] = s_v[i]
                        else:                             # insert new value
                            array_X = cropdata[par][::2]
                            array_Y = cropdata[par][1:][::2]
                            ins_i = bisect(array_X, s_x[i])
                            cropdata[par].insert( ins_i*2, s_x[i])
                            cropdata[par].insert( ins_i*2+1, s_v[i])
                        #print(cropdata[par])
            elif tb_t[par][1] == 'A':                     
                if tb_t[par][0] == 'Y':                  # Totally change table Y values
                    for i in range(len(tb_x[par])):
                        par_tb.append(s_x[i])
                        par_tb.append( s_v[i])
                else:                                    # Totally change table X values
                    for i in range(len(tb_x[par])):
                        par_tb.append(s_v[i])
                        par_tb.append(s_x[i])
                tmp_dict[par] = par_tb
                #print(tmp_dict[par])        
                theta_dict.update(tmp_dict)  
            else:
                raise Exception("There's something wrong with %s, please check it."%par)
        ##########################################################################
    return theta_dict

def create_weather(lat, lon, start_time, end_time, meteo_src):
    """Returns a PCSE weather object based on location, times and data origin.
    
    Parameters
    ----------
    lat : float
        Latitude in decimal degrees
    lon : float
        Longitude in decimal degrees
    start_time : datetime
        Start time
    end_time : datetime
        End time
    meteo_src : Meteo source
        Source of meteo data. Can be either "ERA5", "NASA" or a path to
        a CABO file.
    
    Returns
    -------
    pcse.weather
        Weather object
    
    Raises
    ------
    IOError
        If CABO file doesn't exist.
    """
    if meteo_src == "NASA":
        LOG.info("Using NASA Power meteo data")
        weather = NASAPowerWeatherDataProvider(latitude=lat, longitude=lon)
    elif meteo_src == "ERA5":
        LOG.info("Using ERA5 data!")
        meteo_file = grab_meteo_data(lat, lon,
                        start_time.year, end_time.year,
                        data_dir=data_dir)
        weather = CABOWeatherDataProvider(meteo_file.name,
                        fpath=str(meteo_file))
    else:
        meteo_file = Path(meteo_src)
        if not meteo_file.exists():
            raise IOError(f"You want to use your own CABO file {meteo_file:s}"
                        + f", but it doesn't exist!")
        weather = CABOWeatherDataProvider(meteo_file.name,
                 fpath=str(meteo_file))
    return weather

def ensemble_wofost(lon, lat,  start_time, parameter_priors,
                    crop_param_file,
                    end_time = None, ensemble_size = 3,
                    meteo_src = "NASA", weather_path = None, 
                    out_en_file = None, data_dir=None):
    # Sanity check tests
    assert -180 <= lon <= 180, f"Longitude {lon:g} has to be sensible!"
    assert -90 <= lat <= 90, f"Latitude {lat:g} has to be sensible!"
    if end_time is not None:
        assert end_time > start_time, (f"End time {end_time:s} " + 
                    f"has to be later than start {start_time:s}")
    weather = create_weather(lat, lon, start_time, end_time, meteo_src)    
    crop_file = CABOFileReader(crop_param_file)
