# -*- coding: utf-8 -*-.
"""Supports the Incoherent Scatter Radar at Poker Flat

The Poker Flat Incoherent Scatter Radar (PFISR) observes electrion density
and temperature, ion temperature, ion-neutral collision frequency,
O+ composition and ion velocity through various experiments.

Downloads data from the SRI Madrigal Database.

Parameters
----------
platform : string
    'isr'
name : string
    'pf'
tag : string
    'lp', 'ac'

Example
-------
    import pysat
    import datetime as dt
    dmsp = pysat.Instrument('isr', 'pf', 'lp', clean_level='clean')
    dmsp.download(dt.datetime(2010, 2, 19), dt.datetime(2010, 2, 24))
    dmsp.load(2010,52)

Note
----
    Please provide name and email when downloading data with this routine.

"""

from __future__ import print_function
from __future__ import absolute_import
import datetime as dt
import functools
import logging
import numpy as np

import pysat

from . import methods

logger = logging.getLogger('pysatAMISER_logger')

platform = 'isr'
name = 'pf'
tags = {'lp': 'Long Pulse', 'ac': 'Alternating Code',
        'vvelsLat': 'Vector velocities'}
sat_ids = {'': list(tags.keys())}
test_dates = {'': {'lp': dt.datetime(2014, 2, 19),
                   'ac': dt.datetime(2014, 2, 19)}}

# support list files routine
# use the default CDAWeb method
pf_fname = '{year:4d}{month:02d}{day:02d}.{version:03d}_'
supported_tags = {ss: {'lp': pf_fname + "lp_5min-cal.h5",
                       'ac': pf_fname + "ac_5min-cal.h5",
                       'vvelsLat': pf_fname + "??_5min-cal-vvelsLat-300sec.hf"}
                  for ss in sat_ids.keys()}
list_files = functools.partial(pysat.instruments.methods.nasa_cdaweb.list_files,
                               supported_tags=supported_tags)

# support listing files currently available on remote server
list_remote_files = functools.partial(methods.list_remote_files,
                                      supported_tags=supported_tags)

# let pysat know that data is spread across more than one file
# multi_file_day=True

# Set to False to specify using xarray (not using pandas)
# Set to True if data will be returned via a pandas DataFrame
pandas_format = False

# support load routine
load = functools.partial(methods.load)

# Support download routine
download = functools.partial(methods.download)

# ISRs will sometimes include multiple days within a file labeled with a single
# date. Filter out this extra data using the pysat nanokernel processing queue.
# To ensure this function is always applied first, we set the filter
# function as the default function for (PF).
# Default function is run first by the nanokernel on every load call.
default = pysat.instruments.methods.madrigal.filter_data_single_date


def init(self):
    """Initializes the Instrument object with values specific to PF ISR

    Runs once upon instantiation.

    Parameters
    ----------
    self : pysat.Instrument
        This object

    Returns
    --------
    Void : (NoneType)
        Object modified in place.


    """

    print(methods.amisr_rules())

    return


def clean(self):
    """Routine to return PF ISR data cleaned to the specified level

    Returns
    --------
    Void : (NoneType)
        data in inst is modified in-place.

    Notes
    --------
    Supports 'clean', 'dusty', 'dirty'
    'Clean' is over 120 km for drifts
    'Dusty' and 'Dirty' default to 'Clean'
    'None' None

    Routine is called by pysat, and not by the end user directly.

    """

    min_alt = 120000.0 # 120 km in meters

    # Default to selecting all of the finite (not Inf or NaN) heights
    if self.clean_level in ['clean', 'dusty', 'dirty']:
        if self.clean_level in ['clean', 'dusty']:
            print('WARNING: this level 2 data has no quality flags')

        ida_cor = (self.data['altitude'].values < min_alt)
        ida_uncor = (self.data['altitude_uncor'].values < min_alt)

        # downselect altitude-based data using the cleaning conditions above
        for dkey in self.data.data_vars.keys():
            if self.data[dkey].values.shape == ida_cor.shape:
                self.data[dkey].values[ida_cor] = np.nan
            elif self.data[dkey].values.shape == ida_uncor.shape:
                self.data[dkey].values[ida_uncor] = np.nan
    else:
        # downselection not necessary, provide warning
        if np.any(data['altitude'].values < min_alt):
            logger.WARNING("densities below 120 km require reprocessing")

    return
