#!/usr/bin/env python
import os
from pathlib import Path

import pytest
import numpy as np

import subprocess

from ..model_ensemble import define_prior_distributions
from ..model_ensemble import set_ensemble_parameters


DATA_PATH = os.path.dirname(__file__)


def test_define_prior_gaussian():
    chunk = "\n".join([
        "#PARAM_CODE,Variation,PARAM_XVALUE,PARAM_YVALUE,Min,Max,StdDev,Distribution",
        "AMAXTB00,YA,0,35.83,10,70,12.5,Gaussian"])
    prior_dist, param_list, param_xvalue, param_type = \
            define_prior_distributions(chunk)
    assert param_list[0] == 'AMAXTB00'
    assert param_type['AMAXTB00'] == "YA"
    assert np.allclose(param_xvalue['AMAXTB00'], 0)

# def test_define_prior_uniform():
#         chunk = "\n".join([
#             "#PARAM_CODE,Variation,PARAM_XVALUE,PARAM_YVALUE,Min,Max,StdDev,Distribution",
#             ])
#         prior_dist, param_list, param_xvalue, param_type = \
#                 define_prior_distributions(chunk)
#         assert param_list[0] == 'AMAXTB00'
#         assert param_type['AMAXTB00'] == "YA"
#         assert np.allclose(param_xvalue['AMAXTB00'], 0)