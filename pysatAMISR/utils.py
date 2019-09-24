# -*- coding: utf-8 -*-.
"""Provides useful routines for AMISR instruments

"""

from __future__ import print_function
from __future__ import absolute_import

import pysat
from pysat import utils as pyutils


def calc_geoloc(inst):
    """ Calculate the instrument measurement location in geodetic coordinates

    Parameters
    ----------
    inst : (pysat.Instrument)
        AMISR pysat Instrument object

    Notes
    -----
    Adds 'gdlat' and 'gdlon' to the Instrument.  Intended to be added as a
    set in the custom processing when loading the AMISR data.

    """

    # Get the gedetic latitude, longitude, and radial height for the corrected
    # and uncorrected data
    for skey in ['', '_uncor']:
        akey = "altitude{:s}".format(skey)
        gout = pyutils.coords.local_horizontal_to_global_geo(
            inst.data['az'], inst.data['el'], inst.data[akey]/1000.0,
            inst.data.attrs['site_latitude'], inst.data.attrs['site_longitude'],
            inst.data.attrs['site_altitude']/1000.0, geodetic=True)

        lat_key = "gdlat{:s}".format(skey)
        lon_key = "gdlon{:s}".format(skey)
        rad_key = "gdrad{:s}".format(skey)

        # Assigning as data, to ensure that the number of coordinates match
        # the number of data dimensions
        inst.data = inst.data.assign(lat_key=gout[0], lon_key=gout[1],
                                     rad_key=gout[2])
        inst.data.rename({"lat_key": lat_key, "lon_key": lon_key,
                          "rad_key": rad_key}, inplace=True)

        # Add metadata for the new data values
        sname = "Corrected" if len(skey) == 0 else "Uncorrected"
        lat_name = "{:s} Latitude".format(sname)
        lon_name = "{:s} Longitude".format(sname)
        rad_name = "{:s} Radial Height".format(sname)
        inst.meta[lat_key] = {inst.meta.units_label: 'degrees',
                              inst.meta.name_label: lat_name,
                              inst.meta.desc_label: lat_name,
                              inst.meta.plot_label: lat_name,
                              inst.meta.axis_label: lat_name,
                              inst.meta.scale_label: 'linear',
                              inst.meta.min_label: -90.0,
                              inst.meta.max_label: 90.0,
                              inst.meta.fill_label: np.nan}
        inst.meta[lon_key] = {inst.meta.units_label: 'degrees',
                              inst.meta.name_label: lon_name,
                              inst.meta.desc_label: lon_name,
                              inst.meta.plot_label: lon_name,
                              inst.meta.axis_label: lon_name,
                              inst.meta.scale_label: 'linear',
                              inst.meta.fill_label: np.nan}
        inst.meta[rad_key] = {inst.meta.units_label: 'km',
                              inst.meta.name_label: rad_name,
                              inst.meta.desc_label: rad_name,
                              inst.meta.plot_label: rad_name,
                              inst.meta.axis_label: rad_name,
                              inst.meta.scale_label: 'linear',
                              inst.meta.fill_label: np.nan}

    return
