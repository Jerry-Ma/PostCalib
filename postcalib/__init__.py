# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
PostCalib is a tool that takes calibrated ODI data from QuickReduce, and
delivers science-ready mosaics.
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *  # noqa: F403
# ----------------------------------------------------------------------------

if not _ASTROPY_SETUP_:  # noqa: F405
    # For egg_info test builds to pass, put package imports here.

    from . import qa
    from . import main_cli
