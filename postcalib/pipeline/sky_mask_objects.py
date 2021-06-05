#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-08-03 15:00
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
sky_mask_objects.py

The script is called to mask out objects on a science exposure. The mask
is typically obtained via running of SExtractor and set the CHECKIMAGE
type to SEGMENTATION.

The program takes an science exposure and its segmentation map, combines
them and creates an sky image with all the detected objects masked as
NAN.

A dilation filter is run on the mask prior to the masking so that
the edge pixels of objects with elevated values are also masked. However,
the size of the dilation filter should not be too large especially
when the number of the total frames is low, otherwise there may be
holes in the combined image due to lack of sufficient coverage.

Inputs
------
image: fits image
    The science image containing objects to be removed
segmentation: fits image
    The segmentation check-image of the input image, typically created
    from running SExtractor. The pixels with objects on have value > 1
    and pixels of the sky have value = 0

Outputs
-------
sky_image: fits image
    The image with all objects masked as NAN
"""

from __future__ import (absolute_import, division, print_function)
import os
import re
import sys
import glob
import warnings

import pyregion
import numpy as np
from scipy.ndimage.morphology import binary_dilation
from astropy.io import fits

from postcalib.wiyn import get_layout
from postcalib import qa
from postcalib.apus.common import get_log_func


def main(*args, **kwargs):
    if not args:
        args = sys.argv[1:]
        image_file, segment_file, out_file = args
    else:
        (image_file, segment_file), out_file = args
    log = get_log_func(default_level='debug', **kwargs)

    with fits.open(image_file, memmap=True) as image:

        # look for region masks
        skymask_dir = kwargs['skymask_dir']
        obsid = image[0].header['OBSID']
        regions = glob.glob(os.path.join(skymask_dir, 'skymask.reg'))
        regions.extend(
                glob.glob(
                    os.path.join(skymask_dir, '*{}*.reg'.format(obsid))))
        # looks for additional mask stats with skymask_{key}.reg
        _regions = glob.glob(os.path.join(skymask_dir, "skymask_*.reg"))
        # print(regions)
        # print(_regions)
        add_regions = []
        for reg in sorted(_regions, key=lambda x: len(x)):
            regkey = re.match(
                    "skymask_(.+)\.reg", os.path.basename(reg)).groups()[0]
            if regkey in os.path.basename(image_file):
                add_regions.append(reg)
        if len(add_regions) > 1:
            log('warning', "more than one skymask_*.reg found for {}, use"
                "the last one {}".format(image_file, add_regions[-1]))
        elif len(add_regions) == 1:
            regions.append(add_regions[-1])
        else:
            pass
        log("use DS9 region masks\n{}".format('\n'.join(regions)))

        with fits.open(segment_file, memmap=True) as segment:
            layout = get_layout(image)
            for ext, hdu in layout.enumerate(image):
                otaxy = layout.get_ext_ota(ext)
                log("working on OTA {0}".format(otaxy))
                apply_segment_mask(
                        hdu, segment[ext].data,
                        **kwargs)
                # apply region mask
                apply_region_mask(hdu, regions)
                # image[ext].data = data
        log("save to sky image {}".format(out_file))
        image.writeto(out_file, overwrite=True)
        pr = qa.create_preview(
                hdulist=image, filename=out_file, delete_data=True)
        pr.save()


def apply_segment_mask(hdu, segdata, **kwargs):
    # log = get_log_func(default_level='debug', **kwargs)
    # dilation
    mask = binary_dilation(segdata, iterations=10)
    hdu.data[mask] = np.nan
    data = hdu.data
    with warnings.catch_warnings():
        warnings.filterwarnings(
                'ignore', ".+", RuntimeWarning)
        # also mask out the low value edges
        bkg = np.nanmedian(data) * 3 - np.nanmean(data) * 2
        std = np.nanstd(data)
        hdu.data[np.isnan(data) | (data < bkg - 10 * std)] = np.nan


def apply_region_mask(hdu, regions, **kwargs):
    log = get_log_func(default_level='debug', **kwargs)
    for region in regions:
        log("apply mask region {0}".format(region))
        # create a fresh wcs from keys
        # tmp_hdu = fits.ImageHDU(data=hdu.data)
        # # add basic wcs
        # for key in ['EQUINOX', 'CTYPE1', 'CTYPE2', 'CRPIX1', 'CRPIX2',
        #             'CRVAL1', 'CRVAL2', 'CUNIT1', 'CUNIT2', 'CD1_1', 'CD2_1',
        #             'CD1_2', 'CD2_2']:
        #     tmp_hdu.header[key] = hdu.header[key]
        try:
            mask = pyregion.open(region).get_mask(hdu=hdu)
            hdu.data[mask] = np.nan
        except ValueError:
            log("unable to apply region mask {},"
                " please check the format of the file".format(region))


if __name__ == "__main__":
    main()
