# -*- coding: utf-8 -*-.
"""Supports the Incoherent Scatter Radar at Poker Flat

The Poker Flat Incoherent Scatter Radar (PFISR) observes electrion density
and temperature, ion temperature, ion-neutral collision frequency,
O+ composition and ion velocity through various experiments.

Downloads data from the SRI Madrigal Database.

Properties
----------
platform : string
    'isr'
name : string
    'pf'
tag : string
    'lp', 'ac'

Example
-------
::

    import pysat
    import pysatAMISR as pyisr
    pfisr = pysat.Instrument(inst_module=pyisr.instruments.isr_pf, tag='ac')
    pfisr.load(2019, 33)

"""

import datetime as dt
import functools
import numpy as np

import pysat
from pysatMadrigal.instruments.methods import general as mad_gen

from pysatAMISR.methods import general

logger = pysat.logger

# ----------------------------------------------------------------------------
# Instrument attributes

platform = 'isr'
name = 'pf'
tags = {'lp': 'Long Pulse', 'ac': 'Alternating Code',
        'vvelsLat': 'Vector velocities'}
inst_ids = {'': list(tags.keys())}

pandas_format = False

# ----------------------------------------------------------------------------
# Instrument test attributes

_test_dates = {inst_id: {tag: dt.datetime(2014, 2, 19)
                         for tag in inst_ids[inst_id]}
               for inst_id in inst_ids.keys()}
_test_download_travis = {inst_id: {tag: False for tag in inst_ids[inst_id]}
                         for inst_id in inst_ids.keys()}

# ----------------------------------------------------------------------------
# Instrument methods

# ISRs will sometimes include multiple days within a file labeled with a single
# date. Filter out this extra data using the pysat nanokernel processing queue.
# To ensure this function is always applied first, we set the filter
# function as the default function for (PF).
# Default function is run first by the nanokernel on every load call.
preprocess = mad_gen.filter_data_single_date


def init(self):
    """Initializes the Instrument object with values specific to PF ISR
    """

    self.acknowledgements = general.amisr_rules()
    self.references = "ADD REFERENCES"
    logger.info(self.acknowledgements)

    return


def clean(self):
    """Routine to return PF ISR data cleaned to the specified level

    Note
    ----
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
            logger.warning('this level 2 data has no quality flags')

        ida_cor = (self.data['altitude'].values < min_alt)
        ida_uncor = (self.data['altitude_uncor'].values < min_alt)

        # Downselect altitude-based data using the cleaning conditions above
        for dkey in self.data.data_vars.keys():
            if self.data[dkey].values.shape == ida_cor.shape:
                self.data[dkey].values[ida_cor] = np.nan
            elif self.data[dkey].values.shape == ida_uncor.shape:
                self.data[dkey].values[ida_uncor] = np.nan
    else:
        # Downselection not necessary, provide warning
        if np.any(data['altitude'].values < min_alt):
            logger.warning("densities below 120 km require reprocessing")

    return


# ----------------------------------------------------------------------------
# Instrument functions

# Support list files routine; use the default CDAWeb method
pf_fname = '{year:4d}{month:02d}{day:02d}.{version:03d}_'
supported_tags = {ss: {'lp': pf_fname + "lp_5min-cal.h5",
                       'ac': pf_fname + "ac_5min-cal.h5",
                       'vvelsLat': pf_fname + "??_5min-cal-vvelsLat-300sec.hf"}
                  for ss in inst_ids.keys()}
list_files = functools.partial(pysat.instruments.methods.general.list_files,
                               supported_tags=supported_tags)

# Support listing files currently available on remote server
list_remote_files = functools.partial(general.list_remote_files,
                                      supported_tags=supported_tags)

# Support load routine
load = functools.partial(general.load)

# Support download routine
download = functools.partial(general.download)

