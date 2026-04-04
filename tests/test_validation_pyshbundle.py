#!/usr/bin/env python
# coding: utf-8

# Unit test script file for pyshbundle package
# This scripts follows the notebook and `04_TWS_time_series.ipynb` `validation_pyshbundle.ipynb` and is used to validate the
# pyshbundle package. Please read them for more details.
# Author: Vivek Kumar Yadav, IISc Bengaluru, 2024-09-02
import numpy as np
import os
import scipy.io
from datetime import datetime
from collections import OrderedDict
from pyshbundle.hydro import TWSCalc
from pyshbundle.io import (
    extract_SH_data,
    extract_deg1_coeff_tn13,
    extract_deg2_3_coeff_tn14,
)

from importlib.resources import files, as_file


# Rest of the code...
def validation_pyshbundle(
    path_sh: str, path_tn13: str, path_tn14: str, source: str = "jpl"
):
    """
    Load spherical harmonic data, apply corrections, and compute TWS fields.

    Parameters
    ----------
    path_sh   : path to directory containing .gz GRACE files
    path_tn13 : path to TN-13 degree-1 coefficient file
    path_tn14 : path to TN-14 C20/C30 coefficient file
    source    : data centre identifier, e.g. 'jpl', 'itsg'

    Returns
    -------
    tws_fields : np.ndarray, float32, shape (n_months, nlat, nlon)
    """
    long_mean_file_path = files("pyshbundle").joinpath(
        "data/long_mean/SH_long_mean_jpl.npy"
    )

    sh_files = os.listdir(path_sh)
    file_paths = [
        os.path.join(path_sh, file)
        for file in sh_files
        if os.path.splitext(file)[1] == ".gz"
    ]
    extracted_data = {}

    for file_path in file_paths:
        file_data = extract_SH_data(file_path, source=source)
        if file_data["time_coverage_start"]:
            # Convert time_coverage_start to a datetime object and then format it as yyyy-mm
            if source == "itsg":
                start_date = datetime.strptime(
                    file_data["time_coverage_start"][-7:], "%Y-%m"
                ).strftime("%Y-%m")
            else:
                start_date = datetime.strptime(
                    file_data["time_coverage_start"], "%Y-%m-%dT%H:%M:%S.%f"
                ).strftime("%Y-%m")
            # Use the formatted date as the key
            extracted_data[start_date] = file_data["coefficients"]

    # Time Sort the dictionary by keys (dates)
    sorted_data = OrderedDict(sorted(extracted_data.items()))
    temp_tn14 = extract_deg2_3_coeff_tn14(path_tn14)
    for date_key in temp_tn14.keys():
        if date_key in sorted_data.keys():
            sorted_data[date_key][(2, 0)]["Clm"] = temp_tn14[date_key]["c20"]
            if temp_tn14[date_key]["c30"] is not None:
                sorted_data[date_key][(3, 0)]["Clm"] = temp_tn14[date_key]["c30"]
    temp_tn13 = extract_deg1_coeff_tn13(path_tn13)
    if str.upper(source) == "JPL":
        for key in sorted_data:
            sorted_data[key][(0, 0)] = {
                "Clm": 0.0,
                "Slm": 0.0,
                "Clm_sdev": 0.0,
                "Slm_sdev": 0.0,
            }
            sorted_data[key][(1, 0)] = {
                "Clm": 0.0,
                "Slm": 0.0,
                "Clm_sdev": 0.0,
                "Slm_sdev": 0.0,
            }
            sorted_data[key][(1, 1)] = {
                "Clm": 0.0,
                "Slm": 0.0,
                "Clm_sdev": 0.0,
                "Slm_sdev": 0.0,
            }
    else:
        pass

    for date_key in temp_tn13.keys():
        if date_key[0] in sorted_data.keys():
            sorted_data[date_key[0]][(date_key[1], date_key[2])]["Clm"] = temp_tn13[
                date_key
            ]["Clm"]
            sorted_data[date_key[0]][(date_key[1], date_key[2])]["Slm"] = temp_tn13[
                date_key
            ]["Slm"]
            sorted_data[date_key[0]][(date_key[1], date_key[2])]["Clm_sdev"] = (
                temp_tn13[date_key]["Clm_sdev"]
            )
            sorted_data[date_key[0]][(date_key[1], date_key[2])]["Slm_sdev"] = (
                temp_tn13[date_key]["Slm_sdev"]
            )
    max_degree = np.max(
        [
            degree
            for date in sorted_data.keys()
            for degree, order in sorted_data[date].keys()
        ]
    )
    number_of_months = len(sorted_data.keys())
    sc_mat = np.zeros(
        [number_of_months, max_degree + 1, 2 * (max_degree + 1)], dtype=np.double
    )

    for index, key in enumerate(sorted_data.keys()):
        temp = sorted_data[key]
        for l in range(0, max_degree + 1):
            for m in range(0, l + 1):
                sc_mat[index, l, max_degree + m + 1] = temp[(l, m)]["Clm"]
                sc_mat[index, l, max_degree - m] = temp[(l, m)]["Slm"]
        del temp
    sc_mat = np.delete(sc_mat, max_degree, axis=2)
    # Load the long-term mean SH coefficients for the specified data source
    with as_file(long_mean_file_path) as long_mean_file:
        SH_long_mean_jpl = np.load(long_mean_file)

    delta_sc = np.ones_like(sc_mat) * np.nan
    delta_sc = sc_mat - SH_long_mean_jpl
    lmax, gs, half_rad_gf = 96, 1, 500
    tws_fields = TWSCalc(delta_sc, lmax, gs, half_rad_gf, number_of_months)
    tws_fields = np.float32(tws_fields)

    return tws_fields


def load_matlab_reference():
    """Load the packaged MATLAB reference TWS field."""
    matlab_file_path = files("pyshbundle").joinpath("data/validation_data/tws_sh.mat")
    with as_file(matlab_file_path) as f:
        data = scipy.io.loadmat(f)
    return data["tws_m"]
