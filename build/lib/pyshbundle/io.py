# new implementationhek for reading files - combining some asspects from existing reader_replacer codes
# Author: Abhishek Mhamane, MS-Research Geoinformatics, IIT Kanpur
# 2024-06-10, cleaned, updated: Vivek Kumar Yadav, IISc Bengaluru


from tqdm import tqdm, trange
from . import sc2cs
from datetime import datetime, timedelta
import julian
import gzip
import os
import numpy as np
import re

def clm2cs_new(data):
    """This is an other implementation of clm2cs which uses the clm2sc and then converts using
    sc2cs functions

    Args:
        data (_type_): _description_
    """
    # read the data from clm to sc format
    sc_mat, devsc_mat = clm2sc(data)

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




def parse_jpl_file(file_path: str):
    """Reads the spherical harmonic data provided by JPL

    Args:
        file_path (str): Absolute path to the file
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

        # parse the header strings to extract relavant metadata info
        jpl_header = parse_jpl_header(header_info)

        # parse the header strings to extract relavant metadata info
        jpl_data = extract_SH_data(file_path, source='itsg')

    return jpl_header, jpl_data

def parse_jpl_header(header_info):

    # parse the header info passed by the reader in as list of bytes
    # create a dictionary with key = important params from header file

    header = {}
    
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
    
    normal_keys = ['title', 'institution', 'product_version',
                     'processing_level', 'normalization', 'permanent_tide_flag',
                    ]
    dimension_keys = ['degree', 'order']
    date_time_keys = ['time_coverage_start', 'time_coverage_end', 'unused_days']
    physical_constant_keys = ['earth_gravity_param', 'mean_equator_radius']
    
    for key in normal_keys:
        key_index_in_header = find_word(header_info, key)
        # print(f"{key} - header line = {key_index_in_header +1} value= {' '.join(parse_lines(header_info[key_index_in_header])[3:])[: -3]}")
        header[key] = ' '.join(parse_lines(header_info[key_index_in_header])[3:])[: -3]

    for key in dimension_keys:
        key_index_in_header = find_word(header_info, key)
        val = int(" ".join(parse_lines(header_info[key_index_in_header], parse_fmt='\s+')[3:])[: -3])
        # print(f"{key} - {val}")
        header[key] = val
    
    for key in date_time_keys:
        # TODO: Look back and find what you meant.... 
        key_index_in_header = find_word(header_info, key)
        # find a way to make it date time object so it can be used later
        pass
    
    for key in physical_constant_keys:
        key_index_in_header = find_word(header_info, key)
        
        const_long_name = ' '.join(parse_lines(header_info[key_index_in_header + 1])[3:])[: -3]
        const_units = ' '.join(parse_lines(header_info[key_index_in_header + 2])[3:])[: -3]
        const_value = float(' '.join(parse_lines(header_info[key_index_in_header + 3])[3:])[: -3])
        const_dict = {'units': const_units, 'value': const_value}
        # returning a dict with value and corresponding units
        header[key] = const_dict

    return header


def parse_lines(line, parse_fmt='\s+'):
    #  parses the liness and reutrns an array
    # '\s+' returns array with no whitespace

    parsed_array = re.split(parse_fmt, str(line))

    return parsed_array

def find_word(info_lines, search_key):
    # finding the target word in the read lines

    for i in range(len(info_lines)):
        parsed_array = parse_lines(info_lines[i], parse_fmt='\s+')
        if search_key in parsed_array:
            search_idx = i
            break
    
    return search_idx


def parse_csr_file(file_path: str):
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

        # parse the header strings to extract relavant metadata info
        csr_header = header_info # parse_jpl_header(header_info)

        # parse the data strings to extract numeric data in suitable matrix fmt
        csr_data = extract_SH_data(file_path, source='csr')

        # Organize the data into either matrix, dataframe or dictionary format      

    return csr_header, csr_data


def parse_csr_header():

    # similar to JPL one
    
    raise NotImplementedError("Similar to `parse_jpl_header`... not yet implemented seperately.")


def parse_itsg_file(file_path):

    # ensure that the file path is valid then proceed
    
    source = 'CSR'
    
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

        istg_header = info_lines[:end_of_header_idx]

        # parse the header strings to extract relavant metadata info
        itsg_data = extract_SH_data(file_path, source='itsg')

    return istg_header, itsg_data

def parse_itsg_header(header_info: list):

    normal_keys = ['modelname', 'product_type',
                     'norm', 'tide_system', 'errors', 'earth_gravity_constant', 'radius',
                     'max_degree'
                    ]
    
    # physical_constant_keys = ['earth_gravity_constant', 'radius', ]
    
    for key in normal_keys:
        
        key_index_in_header = find_word(header_info, key)
        #print(f"{key} - header line = {key_index_in_header} value= {parse_lines(header_info[key_index_in_header])[1]}")

    model_name_idx = find_word(header_info, 'modelname')
    date_str = parse_lines(header_info[model_name_idx])[1][-7:]

    header_dict = {}
    '''
    for key in physical_constant_keys:
        key_index_in_header = find_word(header_info, key)
        
        const_long_name = parse_lines(header_info[key_index_in_header])[0]
        const_value = float(parse_lines(header_info[key_index_in_header])[1])
        const_dict = {'long_name': const_long_name, 'value': const_value}
        print(const_dict)

    '''
    return header_dict, date_str


def parse_tn13_header(header_info):

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
    title = ' '.join(re.split("\s+", header_info[title_idx +1])[1:-1]) + ' '.join (re.split("\s+", header_info[title_idx+2])[1:-1]) + ' '.join(re.split("\s+", header_info[title_idx+3])[1:-1])
    
    # TODO: later convert the str object to a date-time object
    last_reported_date = (re.split("\s+", header_info[title_idx+3])[-2])[:-1]

    special_notes = []

    # add parsing for special notes later
    
    return title, last_reported_date

def parse_tn14_header():

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


def find_date_in_replacemnt_file(replacemnt_mat, file_type: str, epoch_begin, epoch_end=None):

    # epoch_begin and epoch_end -> date from the grace data file
    # begin_date and end_data -> date from the replacement file (tn-13 or tn-14)

    rows, cols = replacemnt_mat.shape

    if file_type == 'tn-13':
        time_buffer_itsg = timedelta(days=23)
        date_idxs = set()
        # think of a rather efficient searching scheme
        for i in range(rows):
            begin_date = datetime.strptime(str(int(replacemnt_mat[i][-2])), '%Y%m%d').date()
            end_date = datetime.strptime(str(int(replacemnt_mat[i][-1])), '%Y%m%d').date()

            if epoch_end:
                # for jpl and csr
                if begin_date == epoch_begin and end_date == epoch_end:
                    date_idxs.add(i)
                    print(f"epoch-begin: {epoch_begin}, epoch-end: {epoch_end}, start: {begin_date}, end: {end_date}")
            else:
                # for itsg
                #begin_date = f"{begin_date.year}-{str(begin_date.month).zfill(2)}"
                if type(epoch_begin) == str:
                    epoch_begin = datetime.strptime(epoch_begin, "%Y-%m").date()

                if begin_date - time_buffer_itsg <= epoch_begin <= begin_date + time_buffer_itsg:
                    date_idxs.add(i)
                    print(f"start: {begin_date - time_buffer_itsg}, epoch-begin: {epoch_begin}, UB:{begin_date + time_buffer_itsg}")

            
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
            begin_date = julian.from_jd(replacemnt_mat[i][-2], fmt='mjd').date()
            end_date = julian.from_jd(replacemnt_mat[i][-1], fmt='mjd').date()


            if epoch_end:
                if begin_date >= epoch_begin - time_buffer and end_date <= epoch_end + time_buffer:
                    date_idxs.add(i)
                    print(f"start: {begin_date}, epoch-begin: {epoch_begin}, LB:{epoch_begin - time_buffer}, UB: {epoch_end + time_buffer}, end: {end_date}, epoch-end: {epoch_end}")
            else:
                # for itsg
                if type(epoch_begin) == str:
                    epoch_begin = datetime.strptime(epoch_begin, "%Y-%m").date()
                if begin_date - time_buffer_itsg <= epoch_begin <= begin_date + time_buffer_itsg:
                    date_idxs.add(i)
                    print(f"start: {begin_date - time_buffer_itsg}, epoch-begin: {epoch_begin}, UB:{begin_date + time_buffer_itsg}")

                        # Add bit more error handling statments 
            # rest is fine -> if inputs's right - output is right
            
    else:
        raise ValueError("Technical Note-13 (tn-13) and Technical Note 14 (tn-14) supported...")
    
        
    return list(date_idxs)
    

def extract_C10_11_replcmnt_coeff(data_tn13, source, epoch_begin, epoch_end=None):

    # match the date
    file_type = 'tn-13'

    if epoch_end is not None:
        end_epoch = epoch_end
    else:
        end_epoch = None
    
    if source == 'jpl' or source == 'csr':

        # find the necessary indxes
        replcmnt_idxs = find_date_in_replacemnt_file(data_tn13, file_type, epoch_begin, end_epoch)
        # extract the coeff from tn13 for required dates

        C10 = data_tn13[replcmnt_idxs[0], :-2]
        # extract the coeff from tn13 for required dates

        C11 = data_tn13[replcmnt_idxs[1], :-2]

    elif source == 'itsg':
        
        replcmnt_idxs = find_date_in_replacemnt_file(data_tn13, file_type, f"{epoch_begin.year}-{str(epoch_begin.month).zfill(2)}", end_epoch)
        # extract the coeff from tn13 for required dates

        C10 = data_tn13[replcmnt_idxs[0], :-2]
        # extract the coeff from tn13 for required dates

        C11 = data_tn13[replcmnt_idxs[1], :-2]

    else:
        raise ValueError("Invalid Source. The sources recoginized are CSR, ITSG and JPL")

    return C10, C11


def extract_C20_replcmnt_coeff(data_tn14, source, epoch_begin, epoch_end=None):
    # For JPL 
    # generating a CLM array for C20 and C30
    # NOTE: Zonal coeff. does not have Slm - its taken as 0
    if source == 'jpl' or source == 'csr':
        
        replcmnt_idxs  = find_date_in_replacemnt_file(data_tn14, 'tn-14', epoch_begin, epoch_end)

        C20 = np.array([2, 0, data_tn14[replcmnt_idxs[0], 0:2][0], data_tn14[replcmnt_idxs[0], 0:2][1], 0, 0])
    
    elif source == 'itsg':
        replcmnt_idxs  = find_date_in_replacemnt_file(data_tn14, 'tn-14', epoch_begin, epoch_end=None)

        C20 = np.array([2, 0, data_tn14[replcmnt_idxs[0], 0:2][0], data_tn14[replcmnt_idxs[0], 0:2][1], 0, 0])


    return C20

def extract_C30_replcmnt_coeff(data_tn14, source, epoch_begin, epoch_end=None):
    if source == 'jpl' or source == 'csr':
            
        replcmnt_idxs  = find_date_in_replacemnt_file(data_tn14, 'tn-14', epoch_begin, epoch_end)

        # think about handling nan values while replacing and its impact
        # handle the nan issue

        C30 = np.array([3, 0, data_tn14[replcmnt_idxs[0], 2:4][0], data_tn14[replcmnt_idxs[0], 2:4][1], 0, 0])

        # replace nan values with zeros
        C30[np.isnan(C30)] = 0

    elif source == 'itsg':
        replcmnt_idxs  = find_date_in_replacemnt_file(data_tn14, 'tn-14', epoch_begin, epoch_end=None)

        C30 = np.array([3, 0, data_tn14[replcmnt_idxs[0], 2:4][0], data_tn14[replcmnt_idxs[0], 2:4][1], 0, 0])

        # replace nan values with zeros
        C30[np.isnan(C30)] = 0
    
    return C30

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
    idx_s = sub2ind(sc_mat.shape, clm_mat[:, 0].astype('i'), (lmax - clm_mat[:, 1]).astype('i')).astype('i')
    idx_c = sub2ind(sc_mat.shape, clm_mat[:, 0].astype('i'), (lmax + clm_mat[:, 1]).astype('i')).astype('i')

    
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

def extract_SH_data(file_path, source):
    # Initialize an empty dictionary to store the coefficients and dates
    data = {
        'coefficients': {},
        'time_coverage_start': None,
        'time_coverage_end': None
    }


    # Regular expression pattern to match the lines with coefficients
    coeff_pattern_csr = re.compile(r'^GRCOF2\s+(\d+)\s+(\d+)\s+([-+]?\d*\.\d+E[-+]?\d+)\s+([-+]?\d*\.\d+E[-+]?\d+)\s+([-+]?\d*\.\d+E[-+]?\d+)\s+([-+]?\d*\.\d+E[-+]?\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+')
    coeff_pattern_jpl = re.compile(r'^GRCOF2\s+(\d+)\s+(\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)')
    coeff_pattern_itsg = re.compile(r'^gfc\s+(\d+)\s+(\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)$')


    if source=='jpl': coeff_pattern=coeff_pattern_jpl
    elif source=='csr': coeff_pattern=coeff_pattern_csr
    elif source=='itsg': coeff_pattern=coeff_pattern_itsg


    # Regular expression patterns to match the time coverage start and end lines
    start_pattern = re.compile(r'time_coverage_start\s*:\s*([\d\-T:.]+)')
    end_pattern = re.compile(r'time_coverage_end\s*:\s*([\d\-T:.]+)')
    timeindex_itsg = re.compile(r'^modelname\s+(.+)$')


    # Open and read the gzipped file to extract the time coverage start and end dates
    if source=='itsg':
        with open(file_path, 'rt') as file:
            for line in file:
                # Strip any leading/trailing whitespace characters
                line = line.strip()

                # Search for time coverage start
                start_match = timeindex_itsg.search(line)
                if start_match:
                    data['time_coverage_start'] = start_match.group(1)
                
                # Break the loop if both dates are found
                if data['time_coverage_start']:
                    break
            # File is automatically closed here due to the 'with' statement
        with open(file_path, 'rt') as file:
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
                    clm = np.longdouble(coeff_match.group(3))
                    slm = np.longdouble(coeff_match.group(4))
                    clm_sdev = np.longdouble(coeff_match.group(5))
                    slm_sdev = np.longdouble(coeff_match.group(6))
                    
                    # Store the coefficients in the dictionary
                    data['coefficients'][(degree, order)] = {'Clm': clm, 'Slm': slm,
                                                            'Clm_sdev': clm_sdev, 'Slm_sdev': slm_sdev}



    elif source=='csr' or source=='jpl':
        with gzip.open(file_path, 'rt') as file:   # gzip.open
            for line in file:
                # Strip any leading/trailing whitespace characters
                line = line.strip()

                # Search for time coverage start
                start_match = start_pattern.search(line)
                if start_match:
                    data['time_coverage_start'] = start_match.group(1)
                
                # Search for time coverage end
                end_match = end_pattern.search(line)
                if end_match:
                    data['time_coverage_end'] = end_match.group(1)
                
                # Break the loop if both dates are found
                if data['time_coverage_start'] and data['time_coverage_end']:
                    break
            # File is automatically closed here due to the 'with' statement


        # Open and read the gzipped file again to extract the coefficients
        with gzip.open(file_path, 'rt') as file:
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
                    clm = np.longdouble(coeff_match.group(3))
                    slm = np.longdouble(coeff_match.group(4))
                    clm_sdev = np.longdouble(coeff_match.group(5))
                    slm_sdev = np.longdouble(coeff_match.group(6))
                    
                    # Store the coefficients in the dictionary
                    data['coefficients'][(degree, order)] = {'Clm': clm, 'Slm': slm,
                                                            'Clm_sdev': clm_sdev, 'Slm_sdev': slm_sdev}
    return data


def extract_deg1_coeff_tn13(file_path):
    data_dict = {}
    
    with open(file_path, 'rt') as file:
        lines = file.readlines()
        for line in lines:

            # Extract data using regex
            pattern = re.compile(r'^GRCOF2\s+(\d+)\s+(\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)')
            match = pattern.match(line)
            
            if match:
                degree = int(match.group(1))
                order = int(match.group(2))
                Clm = float(match.group(3))
                Slm = float(match.group(4))
                Clm_sdev = np.longdouble(match.group(5))
                Slm_sdev = np.longdouble(match.group(6))
                epoch_begin = match.group(7)
                epoch_end = match.group(8)

                # Use epoch start as key but in yyyy-mm-dd format
                epoch_key=datetime.strptime(epoch_begin, '%Y%m%d.%H%M%S').strftime('%Y-%m')
                data_dict[epoch_key, degree, order] = {
                    'degree': degree,
                    'order': order,
                    'Clm': Clm,
                    'Slm': Slm,
                    'Clm_sdev': Clm_sdev,
                    'Slm_sdev': Slm_sdev,
                    'epoch_begin': epoch_begin,
                    'epoch_end': epoch_end,
                }
    # Print a sample of the data to check if it's parsed correctly
    # for key in sorted(data_dict.keys())[:5]:  # print first 5 entries
    #     print(f"{key}: {data_dict[key]}")
    return data_dict

def extract_deg2_3_coeff_tn14(file_path):
    data_dict = {}
    
    with open(file_path, 'rt') as file:
        lines = file.readlines()
        for line in lines:

            # Extract data using regex
            pattern = re.compile(
                r'(\d+\.\d+)\s+(\d+\.\d+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+|NaN)?\s+([-\d.eE+]+|NaN)?\s+([-\d.eE+]+|NaN)?\s+(\d+\.\d+)\s+(\d+\.\d+)')
            match = pattern.match(line)
            
            if match:
                mjd_start = float(match.group(1))
                year_frac_start = float(match.group(2))
                c20 = np.longdouble(match.group(3))
                c20_mean_diff = np.longdouble(match.group(4))
                c20_sigma = np.longdouble(match.group(5))
                c30 = match.group(6)
                c30_mean_diff = match.group(7)
                c30_sigma = match.group(8)
                mjd_end = float(match.group(9))
                year_frac_end = float(match.group(10))

                # Only add C30 if it exists (not)
                if c30.lower() != 'nan':
                    c30 = np.longdouble(c30)
                    c30_mean_diff = np.longdouble(c30_mean_diff)
                    c30_sigma = np.longdouble(c30_sigma)
                else:
                    c30 = None
                    c30_mean_diff = None
                    c30_sigma = None

                # Use mjd as key but in yyyy-mm-dd format
                mjd_key = julian.from_jd(mjd_start, fmt='mjd').date().strftime('%Y-%m')
                data_dict[mjd_key] = {
                    'year_frac_start': year_frac_start,
                    'mjd_start': mjd_start,
                    'c20': c20,
                    'c20_mean_diff': c20_mean_diff,
                    'c20_sigma': c20_sigma,
                    'c30': c30,
                    'c30_mean_diff': c30_mean_diff,
                    'c30_sigma': c30_sigma,
                    'mjd_end': mjd_end,
                    'year_frac_end': year_frac_end
                }
    # Print a sample of the data to check if it's parsed correctly
    # for key in sorted(data_dict.keys())[:5]:  # print first 5 entries
    #     print(f"{key}: {data_dict[key]}")
    return data_dict