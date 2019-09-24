# -*- coding: utf-8 -*-.
"""Provides default routines for integrating AMISR instruments into pysat

"""

from __future__ import absolute_import

import numpy as np
import xarray as xr

import pysat
import visuamisr as visr

def amisr_rules():
    """ General acknoledgement statement for AMISR data

    Returns
    -------
    ackn : string
        String with general acknowledgement for all AMISR data

    """
    contact_location = "https://amisr.com/amisr/requests/"
    ackn = "".join(["For all data requests or questions please contact either",
                    " the PI using the contact info available at: ",
                    contact_location])

    return ackn

def load(fnames, tag=None, sat_id=None,
         xarray_coords=['ave_times', 'beamcodes', 'range_gate', 'uncor_gate'],
         xarray_attrs=['site_latitude', 'site_longitude', 'site_altitude',
                       'site_name', 'site_code']):
    """ Loads data from the AMISR data files into xarray

    Parameters
    ----------
    fnames : array-like
        iterable of filename strings, full path, to data files to be loaded.
        This input is nominally provided by pysat itself.
    tag : string ('')
        tag name used to identify particular data set to be loaded.
        This input is nominally provided by pysat itself. While
        tag defaults to None here, pysat provides '' as the default
        tag unless specified by user at Instrument instantiation.
    sat_id : string ('')
        Satellite ID used to identify particular data set to be loaded.
        This input is nominally provided by pysat itself.
    xarray_coords : list
        List of keywords to use as coordinates, with time as the first
        key in the list (default=['ave_times', 'beamcodes', 'range_gate',
        'uncor_gate'])
    xarray_attrs : list
        List of keywords to include as attributes (default=['site_latitude',
        'site_longitude', 'site_altitude', 'site_name', 'site_code'])

    Returns
    -------
    data : xr.DataSet
        An xarray DataSet holding the data from the HDF5 file
    meta : pysat.Meta
        Meta data from the HDF5 file, as well as default values from pysat

    """
    time_key = xarray_coords[0]

    # Open only the first file, output is a dict
    file_data = visr.analyze.read_data(fnames[0])

    # Define the dimensions of each coordiante
    file_data['range_gate'] = np.arange(0, file_data['range'].shape[1], 1)
    file_data['uncor_gate'] = np.arange(0,
                                        file_data['altitude_uncor'].shape[1], 1)
    coords = {coord: {'dims': coord, 'data': file_data[coord]}
              for coord in xarray_coords}

    # Define the attributes for this datas set
    attrs = {attr: file_data[attr] for attr in xarray_attrs}

    # Transpose the magnetic coordinates
    file_data['babs'] = file_data['babs'].transpose()
    file_data['kvec'] = file_data['kvec'].transpose()

    # Extend the keys that have extra dimensions
    extended_keys = ['times', 'kvec']
    for dkey in extended_keys:
        for dd in range(file_data[dkey].shape[-1]):
            ekey = "{:s}_{:d}".format(dkey, dd)
            file_data[ekey] = file_data[dkey].take(dd, axis=-1)

    # Extend the list of keys that are not data keys
    extended_keys.extend(xarray_attrs)
    extended_keys.extend(xarray_coords)

    # Extract the data variables
    data_vars = dict()
    for dkey in file_data.keys():
        if dkey not in extended_keys:
            fdata = file_data[dkey]

            # Determine the dimensions of each data type
            dims = list(xarray_coords)
            for i, j in enumerate(fdata.shape):
                if file_data[dims[i]].shape[0] != j:
                    while len(dims) > i and file_data[dims[i]].shape[0] != j:
                        dims.pop(i)

            while(len(dims) > len(fdata.shape) or
                  file_data[dims[-1]].shape[0] != fdata.shape[-1]):
                dims.pop()

            data_vars[dkey] = {'dims': dims, 'data': fdata}

    # Convert from dictionary to xarray
    xdata = xr.Dataset.from_dict({'coords': coords, 'attrs': attrs,
                                  'dims': xarray_coords,
                                  'data_vars': data_vars})

    # Rename the time coordinate
    if time_key != 'time':
        xdata.rename({time_key: 'time'}, inplace=True)

    # Assign the meta data
    meta = pysat.Meta()
    meta.info = {'acknowledgements': amisr_rules()}
    for dkey in data_vars.keys():
        meta[dkey] = get_metadata(dkey)

    for dkey in xarray_coords:
        meta[dkey] = get_metadata(dkey)

    for dkey in xarray_attrs:
        meta[dkey] = get_metadata(dkey)

    return xdata, meta

def get_metadata(data_key):
    """ Get metadata for a specified data key

    Parameters
    ----------
    data_key : string
        Data key corresponding to known meta data information

    Returns
    -------
    meta_dict : dict
        Dictionary with desired meta data

    """

    long_name = {"az": "Azimuth",
                 "el": "Elevation",
                 "density": "Corrected Neutral Density",
                 "edensity": "Corrected Electron Density",
                 "Te": "Electron Temperature",
                 "Ti": "Ion Temperature",
                 "vel": "Ion Velocity",
                 "eTe": "Electorn Temperature Error",
                 "eTi": "Ion Temperature Error",
                 "evel": "Ion Velocity Error",
                 "range": "Radar range gate",
                 "altitude": "Altitude",
                 "density_uncor": "Uncorrected Neutral Density",
                 "edensity_uncor": "Uncorrected Electron Density",
                 "altitude_uncor": "Uncorrected Altitude",
                 "latitude": "Latitude",
                 "longitude": "Longitude",
                 "babs": "Magnetic Field Magnitude",
                 "times_0": "Signal Start Time",
                 "times_1": "Signal End Time",
                 "kvec_0": "Magnetic field vector (1/3)",
                 "kvec_1": "Magnetic field vector (2/3)",
                 "kvec_2": "Magnetic field vector (3/3)",
                 "site_latitude": "Radar Latitude",
                 "site_longitude": "Radar Longitude",
                 "site_altitude": "Radar Altitude",
                 "site_name": "Radar Name",
                 "site_code": "Radar Code",
                 "ave_time": "Average Time",
                 "beamcodes": "Beam Code",
                 "range_gate": "Range Gate",
                 "uncor_gate": "Uncorrected Range Gate"}

    units = {"az": "degrees", "el": "degrees", "density": "N/A",
             "edensity": "N/A", "Te": "K", "Ti": "K", "vel": "m/s",
             "eTe": "K",  "eTi": "K", "evel": "m/s", "range": "",
             "altitude": "m",  "density_uncor": "N/A",
             "edensity_uncor": "N/A", "altitude_uncor": "m",
             "latitude": "degrees", "longitude": "degrees", "babs": "nT",
             "times_0": "", "times_1": "", "kvec_0": "", "kvec_1": "",
             "kvec_2": "", "site_latitude": "degrees",
             "site_longitude": "degrees", "site_altitude": "m",
             "site_name": "", "site_code": "", "ave_time": "", "beamcodes": "",
             "range_gate": "", "uncor_gate": ""}

    scale = {"az": "linear", "el": "linear", "density": "log",
             "edensity": "log", "Te": "linear", "Ti": "linear", "vel": "linear",
             "eTe": "linear",  "eTi": "linear", "evel": "linear",
             "range": "linear", "altitude": "linear",  "density_uncor": "log",
             "edensity_uncor": "log", "altitude_uncor": "linear",
             "latitude": "linear", "longitude": "linear", "babs": "linear",
             "times_0": "linear", "times_1": "linear", "kvec_0": "linear",
             "kvec_1": "linear", "kvec_2": "linear", "site_latitude": "linear",
             "site_longitude": "linear", "site_altitude": "linear",
             "site_name": "", "site_code": "", "ave_time": "linear",
             "beamcodes": "linear", "range_gate": "linear",
             "uncor_gate": "linear"}

    # Set the default output
    meta_dict = {'units': '', 'long_name': '', 'desc': '', 'label': '',
                 'scale': 'linear', 'notes': '', 'value_min': np.nan,
                 'value_max': np.nan, 'fill': np.nan}

    # Assign the particular data
    if data_key in units.keys():
        meta_dict = {'units': units[data_key], 'long_name': long_name[data_key],
                     'desc': long_name[data_key], 'label': long_name[data_key],
                     'scale': scale[data_key]}

    return meta_dict

def list_remote_files():
    """ List remote files on AMISR database
    """

    return list()

def download(date_array, tag='', sat_id='', data_path=None, user=None,
             password=None):
    """ Downloads data from the AMISR data base

    Parameters
    ----------
    date_array : array-like
        list of datetimes to download data for. The sequence of dates need not
        be contiguous.
    tag : string ('')
        Tag identifier used for particular dataset. This input is provided by
        pysat.
    sat_id : string  ('')
        Satellite ID string identifier used for particular dataset. This input
        is provided by pysat.
    data_path : string (None)
        Path to directory to download data to.
    user : string (None)
        User string input used for download. Provided by user and passed via
        pysat. If an account
        is required for dowloads this routine here must error if user not
        supplied.
    password : string (None)
        Password for data download.

    Notes
    -----
    The user's names should be provided in field user. Ruby Payne-Scott should
    be entered as Ruby+Payne-Scott

    The password field should be the user's email address. These parameters
    are passed to AMISR when downloading.

    The affiliation field is set to pysat to enable tracking of pysat
    downloads.

    """

    return

