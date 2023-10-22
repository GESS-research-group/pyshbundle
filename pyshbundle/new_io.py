# new implementationhek for reading files - combining some asspects from existing reader_replacer codes
# Author: Abhishek Mhamane, MS-Research Geoinformatics, IIT Kanpur


from copy import deepcopy
import julian
from tqdm import tqdm, trange
from pyshbundle import sc2cs, clm2sc
from datetime import datetime, timedelta

import gzip
import os
import numpy as np
import re


def clm2cs_new(data):
    """This is an other implementation of clm2cs which uses the clm2sc and then converts using
    sc2cs functions
    
    Args:
        data (_type_): _description_
    
    Returns:
        numpy.ndarray: Spherical harmonic coefficients in |C\S| format
        numpy.ndarray: Standard deviations associated with SH data arranged in |C\S| format
    """
    # read the data from clm to sc format
    sc_mat, devsc_mat = clm2sc.clm2sc(data)

    # number of files stacked
    num_files = np.shape(sc_mat)[0]
    r, c = np.shape(sc_mat)[1:]

    # cs will be a square matrix
    cs_mat = np.zeros((num_files, r, r))
    devcs_mat = np.zeros((num_files, r, r))

    for ith_file in range(num_files):
        cs_mat[ith_file, :, :] = sc2cs.sc2cs(sc_mat[ith_file, :, :])
        devcs_mat[ith_file, :, :] = sc2cs.sc2cs(devsc_mat[ith_file, :, :])

    return cs_mat, devcs_mat


def read_jpl(file_path: str):
    """Reads the spherical harmonic data provided by JPL

    Args:
        file_path (str): Absolute path to the file
    
    Returns:
        dict: Header info in a structured dictionary
        np.ndarray: SH coefficient data in form of a numpy n-dim array. Format is CLM.
        set: Epoch start and end times as per the GRACE datafile. (eg. 20020519.0)

    """
    # ensure that the file path is valid then proceed

    source = 'JPL'

    # check if the file is ziped or not

    if file_path[-3:] == '.gz':
        # open the file and read the lines
        with gzip.open(file_path, 'r') as file:

            # read the file line wise -> obtain a list of bytes
            info_lines = file.readlines()
            num_lines = len(info_lines)

            for i in range(len(info_lines)):
                # find the end of header sentence in the text file
                if str(info_lines[i]) == str(b'# End of YAML header\n',):
                    end_of_header_idx = i
                    break

        # everything after the header is the numerical data
        header_info = info_lines[:end_of_header_idx]
        jpl_data = info_lines[end_of_header_idx+1:]

        # parse the header strings to extract relavant metadata info
        jpl_header = parse_jpl_header(header_info)

        # call the parse data function and create a dictionary or matrix
        clm_mat, start_date, end_date = parse_jpl_data(jpl_data)

        # build a dictionary

        # Organize the data into either matrix, dataframe or dictionary format

    return jpl_header, clm_mat, start_date, end_date


def parse_jpl_header(header_info: list):
    """Parses the header info and returns important metadata as a dictionary 

    Args:
        header_info (list): list of string representing each of the lines 
    
    Returns:
        dict: Important/Relevant header info as a structured dictionary
    """

        # parse the header info passed by the reader in as list of bytes
    # create a dictionary with key = important params from header file

    header = {}

    # important info from header file

    normal_keys = ['title', 'institution', 'product_version',
                   'processing_level', 'normalization', 'permanent_tide_flag',
                   ]
    dimension_keys = ['degree', 'order']
    date_time_keys = ['time_coverage_start',
                      'time_coverage_end', 'unused_days']
    physical_constant_keys = ['earth_gravity_param', 'mean_equator_radius']

    for key in normal_keys:
        key_index_in_header = find_word(header_info, key)
        # print(f"{key} - header line = {key_index_in_header +1} value= {' '.join(parse_lines(header_info[key_index_in_header])[3:])[: -3]}")
        header[key] = ' '.join(parse_lines(
            header_info[key_index_in_header])[3:])[: -3]

    for key in dimension_keys:
        key_index_in_header = find_word(header_info, key)
        val = int(" ".join(parse_lines(
            header_info[key_index_in_header], parse_fmt='\s+')[3:])[: -3])
        # print(f"{key} - {val}")
        header[key] = val

    for key in date_time_keys:
        # TODO: Look back and find what you meant....
        key_index_in_header = find_word(header_info, key)
        # find a way to make it date time object so it can be used later
        pass

    for key in physical_constant_keys:
        key_index_in_header = find_word(header_info, key)

        const_long_name = ' '.join(parse_lines(
            header_info[key_index_in_header + 1])[3:])[: -3]
        const_units = ' '.join(parse_lines(
            header_info[key_index_in_header + 2])[3:])[: -3]
        const_value = float(' '.join(parse_lines(
            header_info[key_index_in_header + 3])[3:])[: -3])
        const_dict = {'units': const_units, 'value': const_value}
        # returning a dict with value and corresponding units
        header[key] = const_dict

    return header


def parse_jpl_data(jpl_data: list):
    """parses the numerical data from the list of strings read from the file

    Args:
        jpl_data (list): linewise data strings read from the file

    Returns:
        np.ndarray: A (n x 6) 2-d matrix representing the [degrer_l, order_m, C_lm, S_lm, sig_C_lm, sig_S_lm]
        set: Epoch begin and stop dates stored as float (eg. 20020519.0)
    """

    column_headers = ['record_key', 'degree_index', 'order_index', 'clm', 'slm',
                      'clm_std_dev', 'slm_std_dev', 'epoch_begin_time',
                      'epoch_stop_time', 'solution_flags', 'solution_comment']

    degree_index = []
    order_index = []
    clm = []
    slm = []
    clm_std_dev = []
    slm_std_dev = []

    # date format used is YYYYMMDD
    # using set to avoid duplicates
    epoch_begin_time = set()
    epoch_stop_time = set()

    for line in (jpl_data):
        parsed_array = parse_lines(line, parse_fmt='\s+')

        degree_index.append(int(parsed_array[1]))
        order_index.append(int(parsed_array[2]))
        clm.append(float(parsed_array[3]))
        slm.append(float(parsed_array[4]))
        clm_std_dev.append(float(parsed_array[5]))
        slm_std_dev.append(float(parsed_array[6]))
        # NOTE: YYYYMMDD strings typecasted to float values
        epoch_begin_time.add(float(parsed_array[7]))
        epoch_stop_time.add(float(parsed_array[8]))

    clm_matrix = np.array([degree_index, order_index, clm, slm,
                           clm_std_dev, slm_std_dev])

    # output will be a [x, 6] rectangular matrix
    # Taking transpose to achieve required shape

    return clm_matrix.T, epoch_begin_time, epoch_stop_time


def parse_lines(line, parse_fmt='\s+'):
    # HELPER FUNCTIONS

    #  parses the liness and reutrns an array
    # '\s+' returns array with no whitespace

    parsed_array = re.split(parse_fmt, str(line))

    return parsed_array


def find_word(info_lines: list, search_key: str):
    # HELPER FUNCTIONS
    # finding the target word in the read lines

    for i in range(len(info_lines)):
        parsed_array = parse_lines(info_lines[i], parse_fmt='\s+')
        if search_key in parsed_array:
            search_idx = i
            break

    return search_idx


def read_csr(file_path: str):
    """Reads the spherical harmonic data provided by CSR

    Args:
        file_path (str): Absolute path to the GRACE data file from CSR
    
    Returns:
        dict: Header info in a structured dictionary
        np.ndarray: SH coefficient data in form of a numpy n-dim array. Format is KLM
        set: Epoch begin and stop dates stored as float (eg. 20020519.0)
    """
    # ensure that the file path is valid then proceed

    source = 'CSR'

    # check if the file is ziped or not
    if file_path[-3:] == '.gz':
        # open the file and read the lines
        with gzip.open(file_path, 'r') as file:

            # read the file line wise -> obtain a list of bytes
            info_lines = file.readlines()
            num_lines = len(info_lines)

            for i in range(len(info_lines)):
                # find the index of line which indicates end of header info
                if str(info_lines[i]) == str(b'# End of YAML header\n',):
                    end_of_header_idx = i
                    break

        header_info = info_lines[:end_of_header_idx]
        # everything below header info is the sought after numeric data
        csr_data = info_lines[end_of_header_idx+1:]

        # parse the header strings to extract relavant metadata info
        csr_header = parse_jpl_header(header_info)

        # parse the data strings to extract numeric data in suitable matrix fmt
        klm_mat, start_time, stop_time = parse_csr_data(csr_data)

        # Organize the data into either matrix, dataframe or dictionary format

    return csr_header, klm_mat, start_time, stop_time


def parse_csr_header(header_info):
    """Parses the CSR header info and returns important metadata as a dictionary 

    Args:
        header_info (list): list of string representing each of the lines 
    
    Returns:
        dict: Header info in a structured dictionary
    """

    # similar to JPL one

    # parse the header info passed by the reader in as list of bytes
    # create a dictionary with key = important params from header file

    header = {}

    # important info from header file

    normal_keys = ['title', 'institution', 'product_version',
                   'processing_level', 'normalization', 'permanent_tide_flag',
                   ]
    
    dimension_keys = ['degree', 'order']
    
    date_time_keys = ['time_coverage_start',
                      'time_coverage_end', 'unused_days']
    
    physical_constant_keys = ['earth_gravity_param', 'mean_equator_radius']

    for key in normal_keys:
        key_index_in_header = find_word(header_info, key)
        # print(f"{key} - header line = {key_index_in_header +1} value= {' '.join(parse_lines(header_info[key_index_in_header])[3:])[: -3]}")
        header[key] = ' '.join(parse_lines(
            header_info[key_index_in_header])[3:])[: -3]

    for key in dimension_keys:
        key_index_in_header = find_word(header_info, key)
        val = int(" ".join(parse_lines(
            header_info[key_index_in_header], parse_fmt='\s+')[3:])[: -3])
        # print(f"{key} - {val}")
        header[key] = val

    for key in date_time_keys:
        # TODO: Look back and find what you meant....
        key_index_in_header = find_word(header_info, key)
        # find a way to make it date time object so it can be used later
        pass

    for key in physical_constant_keys:
        key_index_in_header = find_word(header_info, key)

        const_long_name = ' '.join(parse_lines(
            header_info[key_index_in_header + 1])[3:])[: -3]
        const_units = ' '.join(parse_lines(
            header_info[key_index_in_header + 2])[3:])[: -3]
        const_value = float(' '.join(parse_lines(
            header_info[key_index_in_header + 3])[3:])[: -3])
        const_dict = {'units': const_units, 'value': const_value}
        # returning a dict with value and corresponding units
        header[key] = const_dict

    return header


def parse_csr_data(csr_data: list):
    """parses the numerical data from the list of strings read from the file

    Args:
        csr_data (list): linewise data strings read from the file

    Returns:
        np.ndarray: A (n x 6) 2-d matrix representing the [degrer_l, order_m, C_lm, S_lm, sig_C_lm, sig_S_lm]
        set: Epoch begin and stop dates stored as float (eg. 20020519.0)
    """

    file_column_headers = ['record_key', 'degree_index', 'order_index', 'clm', 'slm',
                      'clm_std_dev', 'slm_std_dev', 'epoch_begin_time',
                      'epoch_stop_time', 'solution_flags', 'solution_comment']

    degree_index = []
    order_index = []
    clm = []
    slm = []
    clm_std_dev = []
    slm_std_dev = []
    epoch_begin_time = set()
    epoch_stop_time = set()

    for line in (csr_data):
        parsed_array = parse_lines(line, parse_fmt='\s+')

        degree_index.append(int(parsed_array[1]))
        order_index.append(int(parsed_array[2]))
        clm.append(float(parsed_array[3]))
        slm.append(float(parsed_array[4]))
        clm_std_dev.append(float(parsed_array[5]))
        slm_std_dev.append(float(parsed_array[6]))
        # begin and end time being a set object
        epoch_begin_time.add(float(parsed_array[7]))
        epoch_stop_time.add(float(parsed_array[8]))

    klm_matrix = np.array([degree_index, order_index, clm, slm,
                           clm_std_dev, slm_std_dev])

    # output will be a [..., 6] rectangular matrix

    return klm_matrix.T, epoch_begin_time, epoch_stop_time


def read_itsg(file_path: str):
    """Reads the spherical harmonic data provided by ITSG

    Args:
        file_path (str): Absolute path to the file
    
    Returns:
        dict: Header info in a structured dictionary
        np.ndarray: SH coefficient data in form of a numpy n-dim array. Format is KLM
        str: date of the data acquisition
    """

    # ensure that the file path is valid then proceed

    source = 'ITSG'

    # check if the file is ziped or not

    # open the file and read the lines
    if file_path[-4:] == '.gfc':
        with open(file_path, 'r') as file:

            # read the file line wise -> obtain a list of bytes
            info_lines = file.readlines()
            num_lines = len(info_lines)

            for i in range(len(info_lines)):
                if str(info_lines[i]) == str('end_of_head ==================================================================================\n',):
                    end_of_header_idx = i
                    break

        header_info = info_lines[:end_of_header_idx]
        itsg_data = info_lines[end_of_header_idx+1:]

        # call the parse header info functon and create a metadata dictionary
        header_dict, date_str = parse_itsg_header(header_info)

        # call the parse data function and create a dictionary or matrix
        itsg_clm_mat = parse_itsg_data(itsg_data)

        # Organize the data into either matrix, dataframe or dictionary format

    return header_dict, itsg_clm_mat, date_str


def parse_itsg_header(header_info: list):
    """Parses the ITSG header info and returns important metadata as a dictionary

    Args:
        header_info (list): list of string representing each of the lines

    Returns:
        dict: Important/Relevant header info as a structured dictionary
        str: date of the ITSG file straight out of the file name
    """

    normal_keys = ['modelname', 'product_type',
                   'norm', 'tide_system', 'errors', 'earth_gravity_constant', 'radius',
                   'max_degree']

    # physical_constant_keys = ['earth_gravity_constant', 'radius', ]

    for key in normal_keys:

        key_index_in_header = find_word(header_info, key)
        #print(f"{key} - header line = {key_index_in_header} value= {parse_lines(header_info[key_index_in_header])[1]}")

    model_name_idx = find_word(header_info, 'modelname')
    date_str = parse_lines(header_info[model_name_idx])[1][-7:]

    title = header_info[0].split('\n')[0]
    cite_info = ''
    for line in header_info[3:7]:
        cite_info += (line.split('\n')[0] + ' ')
    
    model_name = " ".join(str(header_info[10]).split()).split(" ")[-1]
    product_type = " ".join(str(header_info[11]).split()).split(" ")[-1]
    earth_gravity_constant = " ".join(str(header_info[12]).split()).split(" ")[-1]
    radius = " ".join(str(header_info[13]).split()).split(" ")[-1]
    max_degree = " ".join(str(header_info[14]).split()).split(" ")[-1]
    norm = " ".join(str(header_info[15]).split()).split(" ")[-1]
    tide_system = " ".join(str(header_info[16]).split()).split(" ")[-1]
    errors = " ".join(str(header_info[17]).split()).split(" ")[-1]
    itsg_header_dict = {'modelname': model_name, 'product': product_type,
                    'EarthGravityConst': earth_gravity_constant, 'radius': radius,
                    'MaxDegree': max_degree, 'norm': norm, 'tide_system': tide_system, 'error': errors, 'units': "SI",
                      'citeInfo': cite_info}

    return itsg_header_dict, date_str


def parse_itsg_data(itsg_data: list):
    """parses the numerical data from the list of strings read from the file

    Args:
        itsg_data (list): linewise data strings read from the file

    Returns:
        np.ndarray: A (n x 6) 2-d matrix representing the [degrer_l, order_m, C_lm, S_lm, sig_C_lm, sig_S_lm]
    """

    data_column_headers = ['key', 'L', 'M', 'C', 'S', 'sigma_C', 'sigma_S']

    order_index = []
    degree_index = []
    clm = []
    slm = []
    clm_std_dev = []
    slm_std_dev = []

    for line in (itsg_data):
        parsed_array = parse_lines(line, parse_fmt='\s+')

        degree_index.append(int(parsed_array[1]))
        order_index.append(int(parsed_array[2]))
        clm.append(float(parsed_array[3]))
        slm.append(float(parsed_array[4]))
        clm_std_dev.append(float(parsed_array[5]))
        slm_std_dev.append(float(parsed_array[6]))

    clm_matrix = np.array([degree_index, order_index, clm, slm,
                           clm_std_dev, slm_std_dev])

    # output will be a [..., 6] rectangular matrix

    return clm_matrix.T


def read_tn13(file_path):
    """_summary_

    Args:
        file_path (_type_): _description_

    Returns:
        _type_: _description_
    """

    # check if file exists at given path

    # check the file format - txt of zipped

    if file_path[-4:] == '.txt':
        with open(file_path, 'r') as file:

            # read the file line wise -> obtain a list of bytes
            info_lines = file.readlines()
            num_lines = len(info_lines)

            # find the end of header
            for i in range(len(info_lines)):
                if info_lines[i] == 'end of header ===============================================================================\n':
                    end_of_header_idx = i
                    break

        header_info = info_lines[:end_of_header_idx]
        tn13_data = info_lines[end_of_header_idx+1:]

        # call the parse header info functon and create a metadata dictionary

        # call the parse data function and create a dictionary or matrix
        tn13_clm_mat = parse_tn13_data(tn13_data)

        # Organize the data into either matrix, dataframe or dictionary format

    return tn13_clm_mat


def parse_tn13_header(header_info):
    """_summary_

    Args:
        header_info (_type_): _description_

    Returns:
        _type_: _description_
    """

    # IMP Info
    # - Title
    # - Last reported data point
    # Special Notes
    # - 1, 2, 3, 4, 5

    # finding the index of important sub-headers like Title and Notes
    for i in range(len(header_info)):
        if 'TITLE' in header_info[i]:
            title_idx = i
            break

    for i in range(len(header_info)):
        if 'SPECIAL NOTES' in header_info[i]:
            notes_idx = i
            break

    # The tile is
    title = ' '.join(re.split("\s+", header_info[title_idx + 1])[1:-1]) + ' '.join(re.split(
        "\s+", header_info[title_idx+2])[1:-1]) + ' '.join(re.split("\s+", header_info[title_idx+3])[1:-1])

    # TODO: later convert the str object to a date-time object
    last_reported_date = (re.split("\s+", header_info[title_idx+3])[-2])[:-1]

   #special_notes = []
    simplified_header = {'title': title, "last_updated": last_reported_date}

    # add parsing for special notes later

    return simplified_header



def parse_tn13_data(tn13_data):
    """_summary_

    Args:
        tn13_data (_type_): _description_

    Returns:
        _type_: _description_
    """

    data_column_headers = ['record_key', 'degree_index', 'order_index', 'clm', 'slm',
                           'clm_std_dev', 'slm_std_dev', 'epoch_begin_time', 'epoch_stop_time',
                           'solution_flag']

    # initiate empty lists for each of the column
    order_index = []
    degree_index = []
    clm = []
    slm = []
    clm_std_dev = []
    slm_std_dev = []
    epoch_begin_time = []
    epoch_stop_time = []

    for line in (tn13_data):
        parsed_array = parse_lines(line, parse_fmt='\s+')

        degree_index.append(int(parsed_array[1]))
        order_index.append(int(parsed_array[2]))
        clm.append(float(parsed_array[3]))
        slm.append(float(parsed_array[4]))
        clm_std_dev.append(float(parsed_array[5]))
        slm_std_dev.append(float(parsed_array[6]))
        # this will be important and will be helpful while searching
        # Converting to date-time object is avoided -> numpy matrix needs homogeneous data
        # datetime.strptime(parsed_array[7], '%Y/%m/%d').date() can be used later-on in later stages

        epoch_begin_time.append(float(parsed_array[7]))
        epoch_stop_time.append(float(parsed_array[8]))

    clm_matrix = np.array([degree_index, order_index, clm, slm,
                           clm_std_dev, slm_std_dev, epoch_begin_time, epoch_stop_time])

    # output will be a [..., 6] rectangular matrix

    return clm_matrix.T


def read_tn14(file_path):
    """_summary_

    Args:
        file_path (_type_): _description_

    Returns:
        _type_: _description_
    """

    # check the file extension - '.txt' or some zipped format

    # read the data into a list
    if file_path[-4:] == '.txt':
        with open(file_path, 'r') as file:

            # read the file line wise -> obtain a list of bytes
            info_lines = file.readlines()

            # find the end of header
            for i in range(len(info_lines)):
                if info_lines[i] == 'Product:\n':
                    end_of_header_idx = i
                    break

        header_info = info_lines[:end_of_header_idx]
        tn14_data = info_lines[end_of_header_idx+1:]

        # parse the data
        tn14_raplacement_mat = parse_tn14_data(tn14_data)

    return tn14_raplacement_mat



def parse_tn14_header(header_info):

    # Key info
    # - Title
    title = header_info[0].split(':')[1].split('\n')[0]
    # - Version
    version = header_info[1].split(':')[1].split('\n')[0]
    # - Date Span
    created_date = header_info[2].split(':')[1].split('\n')[0]
    date_span = header_info[3].split(':')[1].split('\n')[0]
    # - Notes:
    description = ''
    for line in header_info[7:13]:
        description += line.lstrip().split("\n")[0] + ' '

    # Constants
    # - Mean C20
    # - Mean C30
    # - GM
    # R
    header_dict = {'title': title, 'version': version, 'created_date': created_date,
                   'date_span': date_span, 'description': description}
    return header_dict


def parse_tn14_data(tn14_data):
    """_summary_

    Args:
        tn14_data (_type_): _description_

    Returns:
        _type_: _description_
    """

    data_column_headers = ['begin_MJD', 'begin_frac_date', 'C20', 'C20_sub_mean_C20',
                           'sigma_C20', 'C30', 'C30_sub_mean_C30', 'sigma_C30',
                           'end_MJD', 'end_frac_date']

    # initialise placeholder lists
    # actually either the MJD or fraction date can be used
    # Going ahead with MJD only for sake of consistency
    begin_MJD, end_MJD = [], []
    C20, C20_sub_mean_C20, sigma_C20 = [], [], []
    C30, C30_sub_mean_C30, sigma_C30 = [], [], []

    for data_line in tn14_data:
        parsed_array = parse_lines(data_line, parse_fmt='\s+')
        begin_MJD.append(float(parsed_array[0]))

        C20.append(float(parsed_array[2]))
        C20_sub_mean_C20.append(float(parsed_array[3]) * 1e-10)
        sigma_C20.append(float(parsed_array[4]) * 1e-10)

        # some of the values are NaN look for it
        C30.append(float(parsed_array[5]))
        C30_sub_mean_C30.append(float(parsed_array[6]) * 1e-10)
        sigma_C30.append(float(parsed_array[7]) * 1e-10)

        end_MJD.append(float(parsed_array[8]))

    # reorganizing the matrix so its similar to TN-13 matrix
    replacement_mat = np.array(
        [C20, sigma_C20, C30, sigma_C30, C20_sub_mean_C20, C30_sub_mean_C30, begin_MJD, end_MJD])

    return replacement_mat.T


def find_date_in_replacemnt_file(replacemnt_mat,file_type: str, source: str, epoch_begin):
    """_summary_

    Args:
        replacemnt_mat (_type_): _description_
        file_type (str): _description_
        source (str): _description_
        epoch_begin (_type_): _description_

    Raises:
        ValueError: _description_

    Returns:
        _type_: _description_
    """

    # epoch_begin and epoch_end -> date from the grace data file
    # begin_date and end_data -> date from the replacement file (tn-13 or tn-14)

    
    rows, cols = replacemnt_mat.shape
    if file_type == 'tn-13':
        time_buffer_itsg = timedelta(days=23)
        time_buffer_csr = timedelta(days=5)
        date_idxs = set()
        # think of a rather efficient searching scheme
        for i in range(rows):
            begin_date = datetime.strptime(
                str(int(replacemnt_mat[i][-2])), '%Y%m%d').date()
            end_date = datetime.strptime(
                str(int(replacemnt_mat[i][-1])), '%Y%m%d').date()

            if source == 'jpl':
                # for jpl and csr
                if begin_date == epoch_begin:
                    date_idxs.add(i)
                    #print(
                    #    f"epoch-begin: {epoch_begin}, start: {begin_date}, end: {end_date}")
            
            elif source == 'itsg':
                # for itsg
                #begin_date = f"{begin_date.year}-{str(begin_date.month).zfill(2)}"
                
                if type(epoch_begin) == str:
                    epoch_begin = datetime.strptime(
                        epoch_begin, "%Y-%m").date()
                
                if begin_date - time_buffer_itsg <= epoch_begin <= begin_date + time_buffer_itsg:
                    date_idxs.add(i)
                    print(f"Data Date - {epoch_begin}, replacemnt date(tn-13) = {begin_date}")
            
            elif source == 'csr':
                if begin_date-time_buffer_csr <= epoch_begin <= begin_date+time_buffer_csr:
                    date_idxs.add(i)
                    print(
                        f"GRACE Data File - epoch-begin: {epoch_begin}, epoch-end:; TN13-Products start: {begin_date}, end: {end_date}")

                pass
            # Add bit more error handling statments
            # rest is fine -> if inputs's right - output is right

    elif file_type == "tn-14":
        time_buffer_itsg = timedelta(days=23)
        time_buffer_csr = timedelta(days=5)
        # there will be only one row per month -> for sake of consistency using set
        # print("TN-14 Replacement file")
        date_idxs = set()
        # think of a rather efficient searching scheme
        for i in range(rows):
            begin_date = julian.from_jd(
                replacemnt_mat[i][-2], fmt='mjd').date()
            end_date = julian.from_jd(replacemnt_mat[i][-1], fmt='mjd').date()

            if source == 'jpl':
                if begin_date == epoch_begin:
                    date_idxs.add(i)
                    print(
                        f"start: {begin_date}, epoch-begin: {epoch_begin}, LB:{epoch_begin }, UB: , end: {end_date}, epoch-end:")
            
            elif source == 'itsg':
                # for itsg
                if type(epoch_begin) == str:
                    epoch_begin = datetime.strptime(
                        epoch_begin, "%Y-%m").date()
                    
                if begin_date - time_buffer_itsg <= epoch_begin <= begin_date + time_buffer_itsg:
                    date_idxs.add(i)
                    print(f"Data Date - {epoch_begin}, replacemnt date (tn-14)_ = {begin_date}")

            elif source == 'csr':
                if begin_date-time_buffer_csr <= epoch_begin <= begin_date+time_buffer_csr:
                    date_idxs.add(i)
                    print(
                        f"start: {begin_date}, epoch-begin: {epoch_begin}, LB:{epoch_begin - time_buffer_csr}, UB: , end: {end_date}, epoch-end:")
            
                    # Add bit more error handling statments
            # rest is fine -> if inputs's right - output is right

    else:
        raise ValueError(
            "Technical Note-13 (tn-13) and Technical Note 14 (tn-14) supported...")

    return list(date_idxs)




def extract_C10_11_replcmnt_coeff(data_tn13, source, epoch_begin, epoch_end=None):
    """_summary_

    Args:
        data_tn13 (_type_): _description_
        source (_type_): _description_
        epoch_begin (_type_): _description_
        epoch_end (_type_, optional): _description_. Defaults to None.

    Raises:
        ValueError: _description_

    Returns:
        _type_: _description_
    """

    # match the date
    file_type = 'tn-13'

    if epoch_end is not None:
        end_epoch = epoch_end
    else:
        end_epoch = None

    if source == 'jpl':

        # find the necessary indxes
        replcmnt_idxs = find_date_in_replacemnt_file(
            data_tn13, file_type='tn-13', source='jpl', epoch_begin=epoch_begin)
        # extract the coeff from tn13 for required dates

        C10 = data_tn13[replcmnt_idxs[0], :-2]
        # extract the coeff from tn13 for required dates

        C11 = data_tn13[replcmnt_idxs[1], :-2]

    elif source == 'itsg':

        replcmnt_idxs = find_date_in_replacemnt_file(
            data_tn13, file_type='tn-13', source='itsg', epoch_begin=f"{epoch_begin.year}-{str(epoch_begin.month).zfill(2)}")
        # extract the coeff from tn13 for required dates

        C10 = data_tn13[replcmnt_idxs[0], :-2]
        # extract the coeff from tn13 for required dates

        C11 = data_tn13[replcmnt_idxs[1], :-2]
    
    elif source == 'csr':

        # find the necessary indxes
        replcmnt_idxs = find_date_in_replacemnt_file(
            data_tn13, file_type='tn-13', source='csr', epoch_begin=epoch_begin)
        # extract the coeff from tn13 for required dates

        C10 = data_tn13[replcmnt_idxs[0], :-2]
        # extract the coeff from tn13 for required dates

        C11 = data_tn13[replcmnt_idxs[1], :-2]


    else:
        raise ValueError(
            "Invalid Source. The sources recoginized are CSR, ITSG and JPL")

    return C10, C11


def extract_C20_replcmnt_coeff(data_tn14, source, epoch_begin, epoch_end=None):
    """_summary_

    Args:
        data_tn14 (_type_): _description_
        source (_type_): _description_
        epoch_begin (_type_): _description_
        epoch_end (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """
    # For JPL
    # generating a CLM array for C20 and C30
    # NOTE: Zonal coeff. does not have Slm - its taken as 0
    if source == 'jpl':

        replcmnt_idxs = find_date_in_replacemnt_file(
            data_tn14, file_type='tn-14', source='jpl', epoch_begin=epoch_begin)

        C20 = np.array([2, 0, data_tn14[replcmnt_idxs[0], 0:2]
                       [0], data_tn14[replcmnt_idxs[0], 0:2][1], 0, 0])

    elif source == 'itsg':
        replcmnt_idxs = find_date_in_replacemnt_file(
            data_tn14, file_type='tn-14', source='itsg', epoch_begin=epoch_begin)

        C20 = np.array([2, 0, data_tn14[replcmnt_idxs[0], 0:2]
                       [0], data_tn14[replcmnt_idxs[0], 0:2][1], 0, 0])
    
    elif source == 'csr':
        replcmnt_idxs = find_date_in_replacemnt_file(
            data_tn14, file_type='tn-14', source='csr', epoch_begin=epoch_begin)

        C20 = np.array([2, 0, data_tn14[replcmnt_idxs[0], 0:2]
                       [0], data_tn14[replcmnt_idxs[0], 0:2][1], 0, 0])

    return C20


def extract_C30_replcmnt_coeff(data_tn14, source, epoch_begin, epoch_end=None):
    """_summary_

    Args:
        data_tn14 (_type_): _description_
        source (_type_): _description_
        epoch_begin (_type_): _description_
        epoch_end (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """
    if source == 'jpl':

        replcmnt_idxs = find_date_in_replacemnt_file(
            data_tn14, file_type='tn-14', source='jpl', epoch_begin=epoch_begin)

        # think about handling nan values while replacing and its impact
        # handle the nan issue

        # C30 data in TN-14 files for initial few years does not require any replacemnt indicated by NaN value
        if np.isnan(data_tn14[replcmnt_idxs[0], 2:4][0]) and np.isnan(data_tn14[replcmnt_idxs[0], 2:4][1]):
            return None
        else:
            C30 = np.array([3, 0, data_tn14[replcmnt_idxs[0], 2:4]
                       [0], data_tn14[replcmnt_idxs[0], 2:4][1], 0, 0])       

    elif source == 'itsg':
        replcmnt_idxs = find_date_in_replacemnt_file(
            data_tn14, file_type='tn-14', source='itsg', epoch_begin=epoch_begin)

        if np.isnan(data_tn14[replcmnt_idxs[0], 2:4][0]) and np.isnan(data_tn14[replcmnt_idxs[0], 2:4][1]):
            return None
        else:
            C30 = np.array([3, 0, data_tn14[replcmnt_idxs[0], 2:4][0], data_tn14[replcmnt_idxs[0], 2:4][1], 0, 0])       
    
    elif source == 'csr':
        replcmnt_idxs = find_date_in_replacemnt_file(
            data_tn14, file_type='tn-14', source='csr', epoch_begin=epoch_begin)

        # think about handling nan values while replacing and its impact
        # handle the nan issue

        if np.isnan(data_tn14[replcmnt_idxs[0], 2:4][0]) and np.isnan(data_tn14[replcmnt_idxs[0], 2:4][1]):
            return None
        else:
            C30 = np.array([3, 0, data_tn14[replcmnt_idxs[0], 2:4]
                       [0], data_tn14[replcmnt_idxs[0], 2:4][1], 0, 0])   

    return C30


def replace_zonal_coeff(data_mat, source, lmax, data_tn13, data_tn14, epoch_begin: float, epoch_end: float):
    """_summary_

    Args:
        data_mat (_type_): _description_
        source (_type_): _description_
        lmax (_type_): _description_
        data_tn13 (_type_): _description_
        data_tn14 (_type_): _description_
        epoch_begin (float): _description_
        epoch_end (float): _description_

    Returns:
        _type_: _description_
    """

    data_mat_copy = deepcopy(data_mat)

    if source == 'jpl':
        assert epoch_end is not None, "epoch_end argument cannot be None"
        # convert the float YYYYMMDD into datetime.date object
        epoch_begin = datetime.strptime(str(int(epoch_begin)), '%Y%m%d').date()
        epoch_end = datetime.strptime(str(int(epoch_end)), '%Y%m%d').date()

        # Extract the C10, C11, C20 and C30 from TN-13 and TN-14
        C10, C11 = extract_C10_11_replcmnt_coeff(
            data_tn13, 'jpl', epoch_begin, epoch_end)
        C20 = extract_C20_replcmnt_coeff(
            data_tn14, source, epoch_begin, epoch_end)
        C30 = extract_C30_replcmnt_coeff(
            data_tn14, source, epoch_begin, epoch_end)

        # For easy replacement purpose
        # [l, m, clm, slm, clm_dev, slm_dev]
        C00 = np.array([0, 0, 0, 0, 0, 0])

        # C30 is  at index - 3 in original matrix
        if C30 is not None:
            data_mat_copy[3, :] = C30

        # C20 is at index - 0 in original matrix
        data_mat_copy[0, :] = C20

        # stack the matrix row-wise
        data_mat_copy = np.row_stack([C11, data_mat_copy])
        data_mat_copy = np.row_stack([C10, data_mat_copy])
        data_mat_copy = np.row_stack([C00, data_mat_copy])

    elif source == 'csr':
        epoch_begin = datetime.strptime(str(int(epoch_begin)), '%Y%m%d').date()
        epoch_end = datetime.strptime(str(int(epoch_end)), '%Y%m%d').date()

        C10, C11 = extract_C10_11_replcmnt_coeff(
            data_tn13, 'csr', epoch_begin, epoch_end)

        C20 = extract_C20_replcmnt_coeff(
            data_tn14, 'csr', epoch_begin, epoch_end)
        C30 = extract_C30_replcmnt_coeff(
            data_tn14, 'csr', epoch_begin, epoch_end)

        # C10 is at index - 1
        # C20 is at index - 2
        # C30 is at index - 3
        # C11 is at index - lmax + 1
        data_mat_copy[lmax+1, :] = C11
        if C30 is not None:
            data_mat_copy[3, :] = C30        
        data_mat_copy[2, :] = C20
        data_mat_copy[1, :] = C10

    elif source == 'itsg':
        # the CSR dates are strings to begin with
        begin_date = datetime.strptime((epoch_begin), '%Y-%m').date()

        C10, C11 = extract_C10_11_replcmnt_coeff(
            data_tn13, 'itsg', epoch_begin=begin_date, epoch_end=None)
        
        print(C10, C11)

        C20 = extract_C20_replcmnt_coeff(
            data_tn14, 'itsg', epoch_begin=begin_date, epoch_end=None)
        C30 = extract_C30_replcmnt_coeff(
            data_tn14, 'itsg', epoch_begin=begin_date, epoch_end=None)

        # For easy replacement purpose
        # C10 is at index 1
        # C11 is at index 2
        # C20 is at index 3
        # C30 is at index 6
        if C30 is not None:
            data_mat_copy[6, :] = C30
        data_mat_copy[3, :] = C20
        data_mat_copy[2, :] = C11
        data_mat_copy[1, :] = C10

    return data_mat_copy


def klm2sc_new(data_mat, lmax: int):
    sc_mat = np.zeros((lmax+1, 2*lmax + 2))
    dev_sc_mat = np.zeros((lmax+1, 2*lmax + 2))
    clm = data_mat[:, 2]
    slm = data_mat[:, 3]
    clm_std_dev = data_mat[:, 4]
    slm_std_dev = data_mat[:, 5]
    
    # first place the slm and then clm
    index2 =0
    for index1 in range(0,lmax+1,1):
        sc_mat[index1:, lmax-index1] = slm[(index2):(index2 + lmax-index1+1)]
        sc_mat[index1:, index1+lmax] = clm[(index2):(index2 + lmax-index1+1)]

        dev_sc_mat[index1:, lmax-index1] = slm_std_dev[(index2):(index2 + lmax-index1+1)]
        dev_sc_mat[index1:, index1+lmax] = clm_std_dev[(index2):(index2 + lmax-index1+1)]
        
        index2 = index2 + lmax-index1+1

    sc_mat=np.delete(sc_mat,lmax,axis=1)
    dev_sc_mat=np.delete(dev_sc_mat,lmax,axis=1)

    return sc_mat, dev_sc_mat



def sub2ind(array_shape, rows, cols):
    # rows, list need to be linear array
    return rows*array_shape[1] + cols


def cklm2sc_new(clm_mat, lmax: int):
    """Transforms the spherical harmonics coefficients data in clm or klm format into a /S|C\ matrix

        clm data - [l, m, c_lm, s_lm]

    Args:
        clm_mat (np.ndarray): _description_
        lmax (int): maximum degree of spherical harmonic expansion

    Returns:
        _type_: _description_
    """

    # initialise an empty sc matrix
    sc_mat = np.zeros([lmax+1, 2*lmax + 1])
    dev_sc_mat = np.zeros([lmax+1, 2*lmax + 1])

    # navigating the maze of indices

    # Use logical indices

    # sc mat requires padding - Taken care of by the earlier initialisation
    #
    # filling the value at appropriate locaation is the key
    #
    # Approach-1
    # run through rows(degree) and fill the cols(order) respectively

    # Approach -2
    # create a row_func s.t. [....., C22, C21, C20, S21, S22, .......]
    # then stack the rows

    # First flatten the SC matrix - column wise aka Fortran style
    # get the flattented idx to be raplaced using sub2ind
    # replace the indices at those locations using
    # unflatten the matrix

    shape_sc = sc_mat.shape

    # following the approach similar to Octave implementation
    # using matrix operations to improve the time efficiency as compared to looping
    idx_s = sub2ind(sc_mat.shape, clm_mat[:, 0].astype(
        'i'), (lmax - clm_mat[:, 1]).astype('i')).astype('i')
    idx_c = sub2ind(sc_mat.shape, clm_mat[:, 0].astype(
        'i'), (lmax + clm_mat[:, 1]).astype('i')).astype('i')

    flat_sc = sc_mat.flatten("F")
    # Attention first place the slm coeff. or else it will relace zonal clm coeff.
    flat_sc[idx_s] = clm_mat[:, 3]
    flat_sc[idx_c] = clm_mat[:, 2]

    flat_sc2 = dev_sc_mat.flatten("F")
    flat_sc2[idx_s] = clm_mat[:, 5]
    flat_sc2[idx_c] = clm_mat[:, 4]

    dev_scmat = flat_sc2.reshape(shape_sc)

    scmat = flat_sc.reshape(shape_sc)

    # with one flag include for

    return scmat, dev_scmat

def check_format(scmat):
    
    pass