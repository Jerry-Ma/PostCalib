#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2017-08-19 22:08
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
sky_subtract.py

This script subtracts a template from an image. The template can
be subtracted in difference modes:

    fringe:
        Fringe mask is created via a binary filtering process,
        and the mask is used for measuring the amplitude used
        for scaling the template before the subtraction.
        In this mode, the region affected by pupil ghost is
        excluded from calculating the scaling factor, enabling
        a separation of the fringe and pupil ghost.
        A second round of combine-subtract could remove the
        residual pupil ghost if necessary.
    pupil:
        Pupil mask is created via a region file, and the
        mask is used for measuring the amplitude of the pupil
        ghost, which is used for the consequent subtraction.
        In this mode the subtraction only affects the central
        OTAs where pupil ghost is present.

Inputs
------
image: fits file
    The image to be subtracted.
template: fits file
    The template to used.

Outputs
-------
out_image: fits file
    The image with the features in the template removed.
"""


from __future__ import (absolute_import, division, print_function)
import os
import sys
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import itertools
import numpy as np
import pyregion

from postcalib.wiyn import get_layout
from postcalib.apus.common import get_log_func
from postcalib import qa


def main(*args, **kwargs):
    log = get_log_func(default_level='debug', **kwargs)
    if not args:
        args = sys.argv[1:]
        image, template, out_file = args
    else:
        (image, template), out_file = args
    log("subtract {} from {}".format(template, image))
    hdulist = subtract_fringe(image, template, kwargs)

    log("write to {}".format(out_file))
    hdulist.writeto(out_file, overwrite=True)
    pr = qa.create_preview(
            hdulist=hdulist, filename=out_file, delete_data=True)
    pr.save()


def subtract_fringe(image, template, kwargs):
    log = get_log_func(default_level='debug', **kwargs)

    fringe = fits.open(template, memmap=True)
    hdulist = fits.open(image, memmap=True)
    layout = get_layout(hdulist)
    for ext, hdu in layout.enumerate(hdulist):
        ota = layout.get_ext_ota(ext)
        log("work on ext {} OTA {}".format(ext, ota))
        # get pupilmask
        regmask_file = os.path.join(
                os.path.dirname(__file__),
                'pupilmask', 'pg_large{}.reg'.format(ota))
        if os.path.exists(regmask_file):
            log("read pupil region mask {}".format(regmask_file))
            hdu = fits.ImageHDU(data=hdu.data)
            pmask = pyregion.open(regmask_file).get_mask(hdu=hdu)
            pmask = ~pmask  # outside is True
        else:
            pmask = np.ones_like(hdu.data, dtype=bool)
        # get fringe mask
        fmask = mask_fringe(fringe[ext].data, layout, thresh=0.8)
        fmask[~pmask] = np.nan  # inside pupil is nan
        hdulist[ext].data = de_fringe(
                hdu.data, fringe[ext].data, fmask, ota, log)
    return hdulist


def mask_fringe(data, wl, thresh=0.8):
    # bin each cell by nbin
    mask = np.empty_like(data) * np.nan
    nbin = 2
    bw = wl.CW / nbin
    bh = wl.CH / nbin
    for cj, ci in itertools.product(range(wl.NCX), range(wl.NCY)):
        (cl, cr), (cb, ct) = wl.get_cell_rect(0, 0, cj, ci)
        for bj, bi in itertools.product(range(nbin), repeat=2):
            bl = cl + bj * bw
            br = bl + bw
            bb = cb + bi * bh
            bt = bb + bh
            # bx = 0.5 * (bl + br)
            # by = 0.5 * (bb + bt)
            bb, bt, bl, br = map(int, (bb, bt, bl, br))
            grid = data[bb:bt, bl:br]
            mask[bb:bt, bl:br][grid > np.nanpercentile(
                grid, thresh * 100)] = 1
            mask[bb:bt, bl:br][grid < np.nanpercentile(
                grid, 100 - thresh * 100)] = 0
    return mask


def de_fringe(data, fringe, fringemask, ota, log):
    # skip all nan extension
    if np.all(np.isnan(data)):
        log("skip all NAN OTA {}".format(ota))
        return data
    # measure amplitude by sampling from the mask
    sample_size = int(1e5)
    hi_index = np.where(fringemask == 1)
    hi_samp_index = np.random.choice(len(hi_index[0]), size=sample_size)
    hi_samp = (hi_index[0][hi_samp_index],
               hi_index[1][hi_samp_index])
    lo_index = np.where(fringemask == 0)
    lo_samp_index = np.random.choice(len(lo_index[0]), size=sample_size)
    lo_samp = (lo_index[0][lo_samp_index],
               lo_index[1][lo_samp_index])
    # with open("fringesample{}.reg".format(otaxy), 'w') as fo:
    #     fo.write("image\n")
    #     for i in range(len(hi_samp[0])):
    #         fo.write("line({},{},{},{})\n".format(
    #             hi_samp[1][i], hi_samp[0][i], lo_samp[1][i], lo_samp[0][i]))
    scale = (data[hi_samp[0], hi_samp[1]] - data[lo_samp[0], lo_samp[1]]) / (
             fringe[hi_samp[0], hi_samp[1]] - fringe[lo_samp[0], lo_samp[1]])
    # sigmaclip
    scale = scale[(~np.isnan(scale)) & (~np.isinf(scale))]
    mean, med, stdev = sigma_clipped_stats(scale, sigma=2, maxiters=10)

    # import matplotlib.pyplot as plt
    # plt.hist(scale, np.linspace(
    #     med - stdev * 5, med + stdev * 5, 100))
    # plt.axvline(med)
    # # plt.imshow(fringemask)
    # plt.show()

    # print("scaling stats:", mean, med, stdev)
    # _, med, std = sigma_clipped_stats(scale, sigma=3.0)
    log("scaling factor: {0} +/- {1}".format(med, stdev))
    data -= fringe * med
    return data
