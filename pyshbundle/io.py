# new implementationhek for reading files - combining some asspects from existing reader_replacer codes
# Author: Abhishek Mhamane, MS-Research Geoinformatics, IIT Kanpur
# 2024-06-10, cleaned, updated: Vivek Kumar Yadav, IISc Bengaluru

# - - - - - - - - - - - - - -
# License:
#    This file is part of PySHbundle.
#    PySHbundle is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Acknowledgement Statement:
#    Please note that PySHbundle has adapted the following code packages,
#    both licensed under GNU General Public License
#       1. SHbundle: https://www.gis.uni-stuttgart.de/en/research/downloads/shbundle/

#       2. Downscaling GRACE Total Water Storage Change using Partial Least Squares Regression
#          https://springernature.figshare.com/collections/Downscaling_GRACE_Total_Water_Storage_Change_using_Partial_Least_Squares_Regression/5054564

# Key Papers Referred:
#    1. Vishwakarma, B. D., Horwath, M., Devaraju, B., Groh, A., & Sneeuw, N. (2017).
#       A data‐driven approach for repairing the hydrological catchment signal damage
#       due to filtering of GRACE products. Water Resources Research,
#       53(11), 9824-9844. https://doi.org/10.1002/2017WR021150

#    2. Vishwakarma, B. D., Zhang, J., & Sneeuw, N. (2021).
#       Downscaling GRACE total water storage change using
#       partial least squares regression. Scientific data, 8(1), 95.
#       https://doi.org/10.1038/s41597-021-00862-6
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

from __future__ import annotations

from datetime import date, datetime, timedelta
from typing import Any

# from pyshbundle.reshape_SH_coefficients import sc2cs
import julian
import gzip
import numpy as np
import re


def extract_SH_data(file_path: str, source: str) -> dict:
    """Extracts the spherical harmonic coefficients from all the given files.

    Currently supports JPL, CSR, and ITSG data sources ONLY. Extracts the spherical harmonic
    coefficients from the given file and returns them in a dictionary. Uses the degree and
    order of a coefficient as the key and the coefficient values as the value.

    Args:
        file_path (str): Absolute path to the file.
        source (str): Source of the data (JPL, CSR, or ITSG).

    Returns:
        (dict): Dictionary containing the coefficients and time coverage start and end dates.
    """
    # Initialize an empty dictionary to store the coefficients and dates
    data: dict[str, Any] = {"coefficients": {}, "time_coverage_start": None, "time_coverage_end": None}

    # Regular expression pattern to match the lines with coefficients
    coeff_pattern_csr = re.compile(
        r"^GRCOF2\s+(\d+)\s+(\d+)\s+([-+]?\d*\.\d+E[-+]?\d+)\s+([-+]?\d*\.\d+E[-+]?\d+)\s+([-+]?\d*\.\d+E[-+]?\d+)\s+([-+]?\d*\.\d+E[-+]?\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+"
    )
    coeff_pattern_jpl = re.compile(
        r"^GRCOF2\s+(\d+)\s+(\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)"
    )
    coeff_pattern_itsg = re.compile(
        r"^gfc\s+(\d+)\s+(\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)$"
    )

    if source == "jpl":
        coeff_pattern = coeff_pattern_jpl
    elif source == "csr":
        coeff_pattern = coeff_pattern_csr
    elif source == "itsg":
        coeff_pattern = coeff_pattern_itsg
    else:
        raise ValueError("Invalid source, pyshbundle only supports JPL, CSR, and ITSG")

    # Regular expression patterns to match the time coverage start and end lines
    start_pattern = re.compile(r"time_coverage_start\s*:\s*([\d\-T:.]+)")
    end_pattern = re.compile(r"time_coverage_end\s*:\s*([\d\-T:.]+)")
    timeindex_itsg = re.compile(r"^modelname\s+(.+)$")

    # Open and read the gzipped file to extract the time coverage start and end dates
    if source == "itsg":
        with open(file_path, "rt") as file:
            for line in file:
                # Strip any leading/trailing whitespace characters
                line = line.strip()

                # Search for time coverage start
                start_match = timeindex_itsg.search(line)
                if start_match:
                    data["time_coverage_start"] = start_match.group(1)

                # Break the loop if both dates are found
                if data["time_coverage_start"]:
                    break
            # File is automatically closed here due to the 'with' statement
        with open(file_path, "rt") as file:
            for line in file:
                # Strip any leading/trailing whitespace characters
                line = line.strip()
                # print(line)

                # Search for the coefficient pattern in the line
                coeff_match = coeff_pattern.search(line)
                if coeff_match:
                    # Extract degree, order, Clm, and Slm
                    degree = int(coeff_match.group(1))
                    order = int(coeff_match.group(2))
                    clm = np.double(coeff_match.group(3))
                    slm = np.double(coeff_match.group(4))
                    clm_sdev = np.double(coeff_match.group(5))
                    slm_sdev = np.double(coeff_match.group(6))

                    # Store the coefficients in the dictionary
                    data["coefficients"][(degree, order)] = {
                        "Clm": clm,
                        "Slm": slm,
                        "Clm_sdev": clm_sdev,
                        "Slm_sdev": slm_sdev,
                    }

    elif source == "csr" or source == "jpl":
        with gzip.open(file_path, "rt") as file:  # gzip.open
            for line in file:
                # Strip any leading/trailing whitespace characters
                line = line.strip()

                # Search for time coverage start
                start_match = start_pattern.search(line)
                if start_match:
                    data["time_coverage_start"] = start_match.group(1)

                # Search for time coverage end
                end_match = end_pattern.search(line)
                if end_match:
                    data["time_coverage_end"] = end_match.group(1)

                # Break the loop if both dates are found
                if data["time_coverage_start"] and data["time_coverage_end"]:
                    break
            # File is automatically closed here due to the 'with' statement

        # Open and read the gzipped file again to extract the coefficients
        with gzip.open(file_path, "rt") as file:
            for line in file:
                # Strip any leading/trailing whitespace characters
                line = line.strip()
                # print(line)

                # Search for the coefficient pattern in the line
                coeff_match = coeff_pattern.search(line)
                if coeff_match:
                    # Extract degree, order, Clm, and Slm
                    degree = int(coeff_match.group(1))
                    order = int(coeff_match.group(2))
                    clm = np.double(coeff_match.group(3))
                    slm = np.double(coeff_match.group(4))
                    clm_sdev = np.double(coeff_match.group(5))
                    slm_sdev = np.double(coeff_match.group(6))

                    # Store the coefficients in the dictionary
                    data["coefficients"][(degree, order)] = {
                        "Clm": clm,
                        "Slm": slm,
                        "Clm_sdev": clm_sdev,
                        "Slm_sdev": slm_sdev,
                    }
    return data


def extract_deg1_coeff_tn13(file_path: str) -> dict:
    """Extracts the degree 1 coefficients from the given TN-13 file.

    Ensure the TN-13 file used is the one recommended by respective data centres (JPL, CSR, or ITSG).
    Similar to extract_SH_data, but specifically for TN-13 files.
    Returns degree 1 replacement coefficients as a dictionary.
    Uses the degree and order of a coefficient as the key and the coefficient values as the value.

    Args:
        file_path (str): Absolute path to the file.

    Returns:
        (dict): Dictionary containing the degree 1 (order 1) coefficients and time coverage start and end dates.
    """
    data_dict = {}

    with open(file_path, "rt") as file:
        lines = file.readlines()
        for line in lines:
            # Extract data using regex
            pattern = re.compile(
                r"^GRCOF2\s+(\d+)\s+(\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)"
            )
            match = pattern.match(line)

            if match:
                degree = int(match.group(1))
                order = int(match.group(2))
                Clm = float(match.group(3))
                Slm = float(match.group(4))
                Clm_sdev = np.double(match.group(5))
                Slm_sdev = np.double(match.group(6))
                epoch_begin = match.group(7)
                epoch_end = match.group(8)

                # Use epoch start as key but in yyyy-mm-dd format
                epoch_key = datetime.strptime(epoch_begin, "%Y%m%d.%H%M%S").strftime(
                    "%Y-%m"
                )
                data_dict[epoch_key, degree, order] = {
                    "degree": degree,
                    "order": order,
                    "Clm": Clm,
                    "Slm": Slm,
                    "Clm_sdev": Clm_sdev,
                    "Slm_sdev": Slm_sdev,
                    "epoch_begin": epoch_begin,
                    "epoch_end": epoch_end,
                }
    # Print a sample of the data to check if it's parsed correctly
    # for key in sorted(data_dict.keys())[:5]:  # print first 5 entries
    #     print(f"{key}: {data_dict[key]}")
    return data_dict


def extract_deg2_3_coeff_tn14(file_path: str) -> dict:
    """Extracts the degree 2 and 3 coefficients from the given file.

    Ensure the TN-14 file used is the one recommended by respective data centres (JPL, CSR, or ITSG).
    Similar to extract_SH_data, but specifically for TN-14 files.
    Returns degree 2, 3 replacement coefficients as a dictionary.
    Uses the degree and order of a coefficient as the key and the coefficient values as the value.

    Args:
        file_path (str): Absolute path to the file.

    Returns:
        (dict): Dictionary containing the degree 2, 3 (order 0) coefficients and time coverage start and end dates.
    """
    data_dict = {}

    with open(file_path, "rt") as file:
        lines = file.readlines()
        for line in lines:
            # Extract data using regex
            pattern = re.compile(
                r"(\d+\.\d+)\s+(\d+\.\d+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+|NaN)?\s+([-\d.eE+]+|NaN)?\s+([-\d.eE+]+|NaN)?\s+(\d+\.\d+)\s+(\d+\.\d+)"
            )
            match = pattern.match(line)

            if match:
                mjd_start = float(match.group(1))
                year_frac_start = float(match.group(2))
                c20 = np.double(match.group(3))
                c20_mean_diff = np.double(match.group(4))
                c20_sigma = np.double(match.group(5))
                c30 = match.group(6)
                c30_mean_diff = match.group(7)
                c30_sigma = match.group(8)
                mjd_end = float(match.group(9))
                year_frac_end = float(match.group(10))

                # Only add C30 if it exists (not)
                if c30.lower() != "nan":
                    c30 = np.double(c30)
                    c30_mean_diff = np.double(c30_mean_diff)
                    c30_sigma = np.double(c30_sigma)
                else:
                    c30 = None
                    c30_mean_diff = None
                    c30_sigma = None

                # Use mjd as key but in yyyy-mm-dd format
                mjd_key = julian.from_jd(mjd_start, fmt="mjd").date().strftime("%Y-%m")
                data_dict[mjd_key] = {
                    "year_frac_start": year_frac_start,
                    "mjd_start": mjd_start,
                    "c20": c20,
                    "c20_mean_diff": c20_mean_diff,
                    "c20_sigma": c20_sigma,
                    "c30": c30,
                    "c30_mean_diff": c30_mean_diff,
                    "c30_sigma": c30_sigma,
                    "mjd_end": mjd_end,
                    "year_frac_end": year_frac_end,
                }
    # Print a sample of the data to check if it's parsed correctly
    # for key in sorted(data_dict.keys())[:5]:  # print first 5 entries
    #     print(f"{key}: {data_dict[key]}")
    return data_dict


# def clm2cs_new(data):
#     """This is an other implementation of clm2cs which uses the clm2sc and then converts using
#     sc2cs functions

#     Args:
#         data (_type_): _description_
#     """
#     # read the data from clm to sc format
#     sc_mat, devsc_mat = clm2sc(data)

#     # number of files stacked
#     num_files = np.shape(sc_mat)[0]
#     r, c = np.shape(sc_mat)[1:]

#     # cs will be a square matrix
#     cs_mat = np.zeros((num_files, r, r))
#     devcs_mat = np.zeros((num_files, r, r))

#     for ith_file in range(num_files):
#         cs_mat[ith_file, :, :] = sc2cs.sc2cs(sc_mat[ith_file, :, :])
#         devcs_mat[ith_file, :, :] = sc2cs.sc2cs(devsc_mat[ith_file, :, :])


#     return cs_mat, devcs_mat


def parse_jpl_file(file_path: str) -> tuple[dict, dict]:
    """Reads the spherical harmonic data provided by JPL.

    Args:
        file_path (str): Absolute path to the file.

    Returns:
        tuple: A tuple containing:
            - jpl_header (dict): Parsed header information.
            - jpl_data (dict): Extracted spherical harmonic coefficients data.
    """
    # ensure that the file path is valid then proceed

    # check if the file is ziped or not

    if file_path[-3:] == ".gz":
        # open the file and read the lines
        with gzip.open(file_path, "r") as file:
            # read the file line wise -> obtain a list of bytes
            info_lines = file.readlines()

            for i in range(len(info_lines)):
                # find the end of header sentence in the text file
                if str(info_lines[i]) == str(
                    b"# End of YAML header\n",
                ):
                    end_of_header_idx = i
                    break

        # everything after the header is the numerical data
        header_info = info_lines[:end_of_header_idx]

        # parse the header strings to extract relavant metadata info
        jpl_header = parse_jpl_header(header_info)

        # parse the header strings to extract relavant metadata info
        jpl_data = extract_SH_data(file_path, source="itsg")

    return jpl_header, jpl_data


def parse_jpl_header(header_info: list) -> dict:
    """Parse JPL GRACE file header and return a dictionary of metadata.

    Args:
        header_info (list): List of byte strings read from the JPL file header.

    Returns:
        dict: Dictionary containing parsed header fields including title, institution,
            product_version, processing_level, normalization, permanent_tide_flag,
            degree, order, earth_gravity_param, and mean_equator_radius.
    """
    # parse the header info passed by the reader in as list of bytes
    # create a dictionary with key = important params from header file

    header: dict[str, Any] = {}

    # important info from header file
    # Dimension - Degree and Order
    # normalization info
    # permanent_tide_flag
    # earth_gravity_param - G*M_earth with units and value
    # mean_equator_radius - units and value
    # processing_level
    # product_version
    # conventions
    # institution
    # title
    # time_coverage_start
    # time_coverage_end
    # unused_days

    normal_keys = [
        "title",
        "institution",
        "product_version",
        "processing_level",
        "normalization",
        "permanent_tide_flag",
    ]
    dimension_keys = ["degree", "order"]
    date_time_keys = ["time_coverage_start", "time_coverage_end", "unused_days"]
    physical_constant_keys = ["earth_gravity_param", "mean_equator_radius"]

    for key in normal_keys:
        key_index_in_header = find_word(header_info, key)
        # print(f"{key} - header line = {key_index_in_header +1} value= {' '.join(parse_lines(header_info[key_index_in_header])[3:])[: -3]}")
        header[key] = " ".join(parse_lines(header_info[key_index_in_header])[3:])[:-3]

    for key in dimension_keys:
        key_index_in_header = find_word(header_info, key)
        val = int(
            " ".join(
                parse_lines(header_info[key_index_in_header], parse_fmt=r"\s+")[3:]
            )[:-3]
        )
        # print(f"{key} - {val}")
        header[key] = val

    for key in date_time_keys:
        # TODO: Look back and find what you meant....
        key_index_in_header = find_word(header_info, key)
        # find a way to make it date time object so it can be used later
        pass

    for key in physical_constant_keys:
        key_index_in_header = find_word(header_info, key)

        const_units = " ".join(parse_lines(header_info[key_index_in_header + 2])[3:])[
            :-3
        ]
        const_value = float(
            " ".join(parse_lines(header_info[key_index_in_header + 3])[3:])[:-3]
        )
        const_dict = {"units": const_units, "value": const_value}
        # returning a dict with value and corresponding units
        header[key] = const_dict

    return header


def parse_lines(line: str | bytes, parse_fmt: str = r"\s+") -> list:
    """Split a line string using the given regex pattern.

    Args:
        line: The input line to parse (string or bytes).
        parse_fmt (str, optional): Regex pattern used to split the line. Defaults to whitespace.

    Returns:
        list: List of substrings after splitting.
    """
    #  parses the liness and reutrns an array
    # '\s+' returns array with no whitespace

    parsed_array = re.split(parse_fmt, str(line))

    return parsed_array


def find_word(info_lines: list, search_key: str) -> int:
    """Find the index of the first line containing the given keyword.

    Args:
        info_lines (list): List of lines (strings or bytes) to search through.
        search_key (str): The keyword to locate within the lines.

    Returns:
        int: Index of the first line that contains the search key.
    """
    # finding the target word in the read lines

    for i in range(len(info_lines)):
        parsed_array = parse_lines(info_lines[i], parse_fmt=r"\s+")
        if search_key in parsed_array:
            search_idx = i
            break

    return search_idx


def parse_csr_file(file_path: str) -> tuple[list, dict]:
    """Read a CSR GRACE .gz file and return the header and spherical harmonic data.

    Args:
        file_path (str): Path to the CSR GRACE .gz file.

    Returns:
        tuple: A tuple containing:
            - csr_header (list): Raw header lines from the file.
            - csr_data (dict): Parsed spherical harmonic coefficients and metadata.
    """
    # ensure that the file path is valid then proceed

    # check if the file is ziped or not
    if file_path[-3:] == ".gz":
        # open the file and read the lines
        with gzip.open(file_path, "r") as file:
            # read the file line wise -> obtain a list of bytes
            info_lines = file.readlines()

            for i in range(len(info_lines)):
                # find the index of line which indicates end of header info
                if str(info_lines[i]) == str(
                    b"# End of YAML header\n",
                ):
                    end_of_header_idx = i
                    break

        header_info = info_lines[:end_of_header_idx]

        # parse the header strings to extract relavant metadata info
        csr_header = header_info  # parse_jpl_header(header_info)

        # parse the data strings to extract numeric data in suitable matrix fmt
        csr_data = extract_SH_data(file_path, source="csr")

        # Organize the data into either matrix, dataframe or dictionary format

    return csr_header, csr_data


def parse_csr_header() -> None:
    """Parse CSR GRACE file header (not yet implemented).

    Raises:
        NotImplementedError: This function is not yet implemented; use parse_jpl_header as reference.
    """
    # similar to JPL one

    raise NotImplementedError(
        "Similar to `parse_jpl_header`... not yet implemented seperately."
    )


def parse_itsg_file(file_path: str) -> tuple[list, dict]:
    """Read an ITSG GRACE .gfc file and return the header and spherical harmonic data.

    Args:
        file_path (str): Path to the ITSG GRACE .gfc file.

    Returns:
        tuple: A tuple containing:
            - istg_header (list): Raw header lines from the file.
            - itsg_data (dict): Parsed spherical harmonic coefficients and metadata.
    """
    # ensure that the file path is valid then proceed

    # check if the file is ziped or not

    # open the file and read the lines
    if file_path[-4:] == ".gfc":
        with open(file_path, "r") as file:
            # read the file line wise -> obtain a list of bytes
            info_lines = file.readlines()

            for i in range(len(info_lines)):
                if str(info_lines[i]) == str(
                    "end_of_head ==================================================================================\n",
                ):
                    end_of_header_idx = i
                    break

        istg_header = info_lines[:end_of_header_idx]

        # parse the header strings to extract relavant metadata info
        itsg_data = extract_SH_data(file_path, source="itsg")

    return istg_header, itsg_data


def parse_itsg_header(header_info: list) -> tuple[dict, str]:
    """Parse ITSG GRACE file header and return a metadata dictionary and model date string.

    Args:
        header_info (list): List of header lines read from the ITSG .gfc file.

    Returns:
        tuple: A tuple containing:
            - header_dict (dict): Dictionary of parsed header metadata.
            - date_str (str): Date string extracted from the model name (YYYY-MM format).
    """
    normal_keys = [
        "modelname",
        "product_type",
        "norm",
        "tide_system",
        "errors",
        "earth_gravity_constant",
        "radius",
        "max_degree",
    ]

    # physical_constant_keys = ['earth_gravity_constant', 'radius', ]

    for key in normal_keys:
        key_index_in_header = find_word(header_info, key)
        # print(f"{key} - header line = {key_index_in_header} value= {parse_lines(header_info[key_index_in_header])[1]}")

    model_name_idx = find_word(header_info, "modelname")
    date_str = parse_lines(header_info[model_name_idx])[1][-7:]

    header_dict: dict = {}
    """
    for key in physical_constant_keys:
        key_index_in_header = find_word(header_info, key)
        
        const_long_name = parse_lines(header_info[key_index_in_header])[0]
        const_value = float(parse_lines(header_info[key_index_in_header])[1])
        const_dict = {'long_name': const_long_name, 'value': const_value}
        print(const_dict)

    """
    return header_dict, date_str


def parse_tn13_header(header_info: list) -> tuple[str, str]:
    """Parse TN-13 replacement coefficient file header and return title and last reported date.

    Args:
        header_info (list): List of header lines read from the TN-13 file.

    Returns:
        tuple: A tuple containing:
            - title (str): Title string extracted from the header.
            - last_reported_date (str): Last reported data point date as a string.
    """
    # IMP Info
    # - Title
    # - Last reported data point
    # Special Notes
    # - 1, 2, 3, 4, 5

    # finding the index of important sub-headers like Title and Notes
    for i in range(len(header_info)):
        if "TITLE" in header_info[i]:
            title_idx = i
            break

    for i in range(len(header_info)):
        if "SPECIAL NOTES" in header_info[i]:
            break

    # The tile is
    title = (
        " ".join(re.split(r"\s+", header_info[title_idx + 1])[1:-1])
        + " ".join(re.split(r"\s+", header_info[title_idx + 2])[1:-1])
        + " ".join(re.split(r"\s+", header_info[title_idx + 3])[1:-1])
    )

    # TODO: later convert the str object to a date-time object
    last_reported_date = (re.split(r"\s+", header_info[title_idx + 3])[-2])[:-1]

    return title, last_reported_date


def parse_tn14_header() -> None:
    """Parse TN-14 replacement coefficient file header (not yet implemented).

    Raises:
        NotImplementedError: This function is not yet implemented.
    """
    # Key info
    # - Title
    # - Version
    # - Date Span
    # - Notes:

    # Constants
    # - Mean C20
    # - Mean C30
    # - GM
    # R

    pass


def find_date_in_replacemnt_file(
    replacemnt_mat: np.ndarray, file_type: str, epoch_begin: date | str, epoch_end: date | None = None
) -> list:
    """Find row indices in a replacement coefficient matrix matching the given epoch dates.

    Args:
        replacemnt_mat (numpy.ndarray): Matrix of replacement coefficients with date columns.
        file_type (str): Type of replacement file; either 'tn-13' or 'tn-14'.
        epoch_begin: Start date of the GRACE data epoch (date object or YYYY-MM string for ITSG).
        epoch_end (optional): End date of the GRACE data epoch. Required for JPL/CSR sources. Defaults to None.

    Returns:
        list: List of row indices in the replacement matrix that match the given epoch.

    Raises:
        ValueError: If file_type is not 'tn-13' or 'tn-14'.
    """
    # epoch_begin and epoch_end -> date from the grace data file
    # begin_date and end_data -> date from the replacement file (tn-13 or tn-14)

    rows, cols = replacemnt_mat.shape

    if file_type == "tn-13":
        time_buffer_itsg = timedelta(days=23)
        date_idxs = set()
        # think of a rather efficient searching scheme
        for i in range(rows):
            begin_date = datetime.strptime(
                str(int(replacemnt_mat[i][-2])), "%Y%m%d"
            ).date()
            end_date = datetime.strptime(
                str(int(replacemnt_mat[i][-1])), "%Y%m%d"
            ).date()

            if epoch_end:
                # for jpl and csr
                if begin_date == epoch_begin and end_date == epoch_end:
                    date_idxs.add(i)
                    print(
                        f"epoch-begin: {epoch_begin}, epoch-end: {epoch_end}, start: {begin_date}, end: {end_date}"
                    )
            else:
                # for itsg
                # begin_date = f"{begin_date.year}-{str(begin_date.month).zfill(2)}"
                if isinstance(epoch_begin, str):
                    epoch_begin = datetime.strptime(epoch_begin, "%Y-%m").date()

                if (
                    begin_date - time_buffer_itsg
                    <= epoch_begin
                    <= begin_date + time_buffer_itsg
                ):
                    date_idxs.add(i)
                    print(
                        f"start: {begin_date - time_buffer_itsg}, epoch-begin: {epoch_begin}, UB:{begin_date + time_buffer_itsg}"
                    )

            # Add bit more error handling statments
            # rest is fine -> if inputs's right - output is right

    elif file_type == "tn-14":
        # there will be only one row per month -> for sake of consistency using set
        # print("TN-14 Replacement file")
        date_idxs = set()
        # think of a rather efficient searching scheme
        time_buffer = timedelta(days=5)
        time_buffer_itsg = timedelta(days=23)
        for i in range(rows):
            begin_date = julian.from_jd(replacemnt_mat[i][-2], fmt="mjd").date()
            end_date = julian.from_jd(replacemnt_mat[i][-1], fmt="mjd").date()

            if epoch_end:
                if (
                    begin_date >= epoch_begin - time_buffer  # type: ignore[operator]
                    and end_date <= epoch_end + time_buffer
                ):
                    date_idxs.add(i)
                    print(
                        f"start: {begin_date}, epoch-begin: {epoch_begin}, LB:{epoch_begin - time_buffer}, UB: {epoch_end + time_buffer}, end: {end_date}, epoch-end: {epoch_end}"  # type: ignore[operator]
                    )
            else:
                # for itsg
                if isinstance(epoch_begin, str):
                    epoch_begin = datetime.strptime(epoch_begin, "%Y-%m").date()
                if (
                    begin_date - time_buffer_itsg
                    <= epoch_begin
                    <= begin_date + time_buffer_itsg
                ):
                    date_idxs.add(i)
                    print(
                        f"start: {begin_date - time_buffer_itsg}, epoch-begin: {epoch_begin}, UB:{begin_date + time_buffer_itsg}"
                    )

                    # Add bit more error handling statments
            # rest is fine -> if inputs's right - output is right

    else:
        raise ValueError(
            "Technical Note-13 (tn-13) and Technical Note 14 (tn-14) supported..."
        )

    return list(date_idxs)


def extract_C10_11_replcmnt_coeff(data_tn13: np.ndarray, source: str, epoch_begin: date | str, epoch_end: date | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Extract degree-1 replacement coefficients C10 and C11 from TN-13 data for the given epoch.

    Args:
        data_tn13 (numpy.ndarray): TN-13 replacement coefficient matrix.
        source (str): Data centre identifier ('jpl', 'csr', or 'itsg').
        epoch_begin: Start date of the GRACE epoch (date object or YYYY-MM string for ITSG).
        epoch_end (optional): End date of the GRACE epoch. Required for JPL/CSR. Defaults to None.

    Returns:
        tuple: A tuple containing:
            - C10 (numpy.ndarray): Replacement C10 coefficient row.
            - C11 (numpy.ndarray): Replacement C11 coefficient row.

    Raises:
        ValueError: If source is not 'jpl', 'csr', or 'itsg'.
    """
    # match the date
    file_type = "tn-13"

    if epoch_end is not None:
        end_epoch = epoch_end
    else:
        end_epoch = None

    if source == "jpl" or source == "csr":
        # find the necessary indxes
        replcmnt_idxs = find_date_in_replacemnt_file(
            data_tn13, file_type, epoch_begin, end_epoch
        )
        # extract the coeff from tn13 for required dates

        C10 = data_tn13[replcmnt_idxs[0], :-2]
        # extract the coeff from tn13 for required dates

        C11 = data_tn13[replcmnt_idxs[1], :-2]

    elif source == "itsg":
        replcmnt_idxs = find_date_in_replacemnt_file(
            data_tn13,
            file_type,
            f"{epoch_begin.year}-{str(epoch_begin.month).zfill(2)}",  # type: ignore[union-attr]
            end_epoch,
        )
        # extract the coeff from tn13 for required dates

        C10 = data_tn13[replcmnt_idxs[0], :-2]
        # extract the coeff from tn13 for required dates

        C11 = data_tn13[replcmnt_idxs[1], :-2]

    else:
        raise ValueError(
            "Invalid Source. The sources recoginized are CSR, ITSG and JPL"
        )

    return C10, C11


def extract_C20_replcmnt_coeff(data_tn14: np.ndarray, source: str, epoch_begin: date | str, epoch_end: date | None = None) -> np.ndarray:
    """Extract the C20 replacement coefficient from TN-14 data for the given epoch.

    Args:
        data_tn14 (numpy.ndarray): TN-14 replacement coefficient matrix.
        source (str): Data centre identifier ('jpl', 'csr', or 'itsg').
        epoch_begin: Start date of the GRACE epoch (date object or YYYY-MM string for ITSG).
        epoch_end (optional): End date of the GRACE epoch. Required for JPL/CSR. Defaults to None.

    Returns:
        numpy.ndarray: Array of shape (6,) containing [l, m, Clm, Slm, Clm_sdev, Slm_sdev] for C20.
    """
    # For JPL
    # generating a CLM array for C20 and C30
    # NOTE: Zonal coeff. does not have Slm - its taken as 0
    if source == "jpl" or source == "csr":
        replcmnt_idxs = find_date_in_replacemnt_file(
            data_tn14, "tn-14", epoch_begin, epoch_end
        )

        C20 = np.array(
            [
                2,
                0,
                data_tn14[replcmnt_idxs[0], 0:2][0],
                data_tn14[replcmnt_idxs[0], 0:2][1],
                0,
                0,
            ]
        )

    elif source == "itsg":
        replcmnt_idxs = find_date_in_replacemnt_file(
            data_tn14, "tn-14", epoch_begin, epoch_end=None
        )

        C20 = np.array(
            [
                2,
                0,
                data_tn14[replcmnt_idxs[0], 0:2][0],
                data_tn14[replcmnt_idxs[0], 0:2][1],
                0,
                0,
            ]
        )

    return C20


def extract_C30_replcmnt_coeff(data_tn14: np.ndarray, source: str, epoch_begin: date | str, epoch_end: date | None = None) -> np.ndarray:
    """Extract the C30 replacement coefficient from TN-14 data for the given epoch.

    Args:
        data_tn14 (numpy.ndarray): TN-14 replacement coefficient matrix.
        source (str): Data centre identifier ('jpl', 'csr', or 'itsg').
        epoch_begin: Start date of the GRACE epoch (date object or YYYY-MM string for ITSG).
        epoch_end (optional): End date of the GRACE epoch. Required for JPL/CSR. Defaults to None.

    Returns:
        numpy.ndarray: Array of shape (6,) containing [l, m, Clm, Slm, Clm_sdev, Slm_sdev] for C30,
            with NaN values replaced by zero.
    """
    if source == "jpl" or source == "csr":
        replcmnt_idxs = find_date_in_replacemnt_file(
            data_tn14, "tn-14", epoch_begin, epoch_end
        )

        # think about handling nan values while replacing and its impact
        # handle the nan issue

        C30 = np.array(
            [
                3,
                0,
                data_tn14[replcmnt_idxs[0], 2:4][0],
                data_tn14[replcmnt_idxs[0], 2:4][1],
                0,
                0,
            ]
        )

        # replace nan values with zeros
        C30[np.isnan(C30)] = 0

    elif source == "itsg":
        replcmnt_idxs = find_date_in_replacemnt_file(
            data_tn14, "tn-14", epoch_begin, epoch_end=None
        )

        C30 = np.array(
            [
                3,
                0,
                data_tn14[replcmnt_idxs[0], 2:4][0],
                data_tn14[replcmnt_idxs[0], 2:4][1],
                0,
                0,
            ]
        )

        # replace nan values with zeros
        C30[np.isnan(C30)] = 0

    return C30


def sub2ind(array_shape: tuple, rows: int | np.ndarray, cols: int | np.ndarray) -> int | np.ndarray:
    """Convert row and column subscripts to linear indices.

    Args:
        array_shape (tuple): Shape of the array as a tuple (num_rows, num_cols).
        rows (int or array-like): Row indices.
        cols (int or array-like): Column indices.

    Returns:
        int or array-like: Linear indices corresponding to the row and column subscripts.
    """
    # rows, list need to be linear array
    return rows * array_shape[1] + cols


# def klm2sc_new(data_mat, lmax: int):
#     sc_mat = np.zeros((lmax+1, 2*lmax + 2))
#     dev_sc_mat = np.zeros((lmax+1, 2*lmax + 2))
#     clm = data_mat[:, 2]
#     slm = data_mat[:, 3]
#     clm_std_dev = data_mat[:, 4]
#     slm_std_dev = data_mat[:, 5]

#     # first place the slm and then clm
#     index2 =0
#     for index1 in range(0,lmax+1,1):
#         sc_mat[index1:, lmax-index1] = slm[(index2):(index2 + lmax-index1+1)]
#         sc_mat[index1:, index1+lmax] = clm[(index2):(index2 + lmax-index1+1)]

#         dev_sc_mat[index1:, lmax-index1] = slm_std_dev[(index2):(index2 + lmax-index1+1)]
#         dev_sc_mat[index1:, index1+lmax] = clm_std_dev[(index2):(index2 + lmax-index1+1)]

#         index2 = index2 + lmax-index1+1

#     sc_mat=np.delete(sc_mat,lmax,axis=1)
#     dev_sc_mat=np.delete(dev_sc_mat,lmax,axis=1)

#     return sc_mat, dev_sc_mat
