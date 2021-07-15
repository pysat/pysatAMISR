# -*- coding: utf-8 -*-.
"""Provides default routines for integrating AMISR instruments into pysat

"""

import numpy as np
import warnings
import xarray as xr

import pysat
import visuamisr as visr

def amisr_rules():
    """ General acknowledgement statement for AMISR data

    Returns
    -------
    ackn : str
        String with general acknowledgement for all AMISR data

    """
    contact_location = "https://amisr.com/amisr/requests/"
    ackn = "".join(["For all data requests or questions please contact either",
                    " the PI using the contact info available at: ",
                    contact_location])

    return ackn

def load(fnames, tag=None, inst_id=None,
         xarray_coords=['ave_times', 'beamcodes', 'range_gate', 'uncor_gate'],
         xarray_attrs=['site_latitude', 'site_longitude', 'site_altitude',
                       'site_name', 'site_code']):
    """ Loads data from the AMISR data files into xarray

    Parameters
    ----------
    fnames : array-like
        iterable of filename strings, full path, to data files to be loaded.
        This input is nominally provided by pysat itself.
    tag : str or NoneType
        tag name used to identify particular data set to be loaded.
        This input is nominally provided by pysat itself. (default=None)
        tag unless specified by user at Instrument instantiation.
    inst_id : str or NoneType
        Instrument ID used to identify particular data set to be loaded.
        (default=None)
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
    for dkey in data_vars.keys():
        meta[dkey] = get_metadata(dkey, meta.labels)

    for dkey in xarray_coords:
        meta[dkey] = get_metadata(dkey, meta.labels)

    for dkey in xarray_attrs:
        meta[dkey] = get_metadata(dkey, meta.labels)

    return xdata, meta

def get_metadata(data_key, mlabels):
    """ Get metadata for a specified data key

    Parameters
    ----------
    data_key : str
        Data key corresponding to known meta data information
    mlabels : pysat.MetaLabels
        MetaLabels class object with Instrument meta data labels

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

    # TODO HERE: fill these in based on data
    min_val = {data_key: np.nan}
    max_val = {data_key: np.nan}
    fill_val = {data_key: np.nan}

    # Assign the particular data
    if data_key in units.keys():
        meta_dict = {mlabels.units: units[data_key],
                     mlabels.name: long_name[data_key],
                     mlabels.desc: long_name[data_key],
                     mlabels.notes: '',
                     mlabels.min_val: min_val[data_key],
                     mlabels.max_val: max_val[data_key],
                     mlabels.fill_val: fill_val[data_key]}

    return meta_dict

def list_remote_files():
    """ List remote files on AMISR database
    """

    return list()

def download(date_array, tag='', inst_id='', data_path=None):
    """ Downloads data from the AMISR data base

    Parameters
    ----------
    date_array : array-like
        list of datetimes to download data for. The sequence of dates need not
        be contiguous.
    tag : str
        Tag identifier used for particular dataset. This input is provided by
        pysat. (default='')
    inst_id : str
        Satellite ID str identifier used for particular dataset. This input
        is provided by pysat. (default='')
    data_path : str
        Path to directory to download data to. (default=None)

    Warnings
    --------
    Scripted downloads not currently supported by SRI

    """

    # Link address:
    # https://data.amisr.com/database/dbase_site_media/PFISR/Experiments/20190202.004/DataFiles/20190202.004_ac_1min-fitcal.h5
    # https://data.amisr.com/database/dbase_site_media/PFISR/Experiments/20190202.002/DataFiles/20190202.002_lp_5min-fitcal.h5

    remoteaccess = {'method': 'http', 'host': 'data.amisr.com',
                    'path': 'database/dbase_site_media',}

    warnings.warn('script downloads are currently not supported by SRI')

    return


