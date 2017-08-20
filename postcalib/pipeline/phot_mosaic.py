#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2017-08-20 06:03
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
phot_mosaic.py

This script contains utility routines that helps create the mosaic.
"""


import os
import re
import sys
import glob
import numpy as np
from astropy.io import fits
from postcalib.apus.common import get_log_func


def get_header(*args, **kwargs):
    if not args:
        args = sys.argv[1:]
    image_file, jobdir, header_file = args
    log = get_log_func(default_level='debug', **kwargs)
    log("link header for {}".format(image_file))
    parsed_filename = re.match(kwargs['reg_inputs'],
                               os.path.basename(image_file)).groupdict()
    hdr = glob.glob(os.path.join(
        jobdir, kwargs['phot_hdr_glob'].format(**parsed_filename)))
    if len(hdr) == 0:
        log("no header found for {}".format(image_file))
        if os.path.exists(header_file):
            log("purge {}".format(header_file))
            os.remove(header_file)
    else:
        linkname = os.path.basename(hdr[0])
        if not os.path.exists(header_file):
            log("link {} to {}".format(linkname, header_file))
            os.symlink(linkname, header_file)
        else:
            log("link exists {}".format(header_file))


def apply_weight(*args, **kwargs):
    if not args:
        args = sys.argv[1:]
    image_file, weight_file, masked_file = args
    log = get_log_func(default_level='debug', **kwargs)
    log("apply weight to {}".format(image_file))
    hlin = fits.open(image_file)
    hlwht = fits.open(weight_file)
    for i in range(len(hlin)):
        if isinstance(hlin[i], (fits.PrimaryHDU, fits.ImageHDU)):
            badmask = hlwht[i].data == 0
            hlin[i].data[badmask] = np.nan
    hlin.writeto(masked_file, overwrite=True)
