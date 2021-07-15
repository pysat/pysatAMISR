# -*- coding: utf-8 -*-.
"""Module for loading AMISR data into pysat
"""
import os as _os

from pysatAMISR import instruments, utils  # noqua F401

# Set the package version
with open(_os.path.join(_os.path.abspath(_os.path.dirname(__file__)),
                        "version.txt"), "r") as fin:
    __version__ = fin.read().strip()

# Clean up
del fin
