#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2017-08-19 17:56
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
sky_combine.py

This script combines a number of sky images to a fringe/pupil ghost
template.

The script makes use of the C extension written by R. Kotulla in his
QuickReduce package for the sigma-clipped image combine.

The script performs the operation in parallel, each process for one
OTA.

Inputs
------
images: fits files
    The sky images to be combined. The objects in these images should
    be masked.

Outputs
-------
template image: fits file
    The resultant image after combining the inputs
"""

from __future__ import (absolute_import, division, print_function)
import os
import sys

import numpy as np
from astropy.io import fits
# from functools import partial
import itertools

# from multiprocessing import cpu_count, Pool

from scipy.ndimage import uniform_filter  # , gaussian_filter, median_filter
# from scipy.stats import sigmaclip
import pyregion
from scipy import interpolate

from postcalib.utils import mp_traceback
from postcalib.wiyn import get_layout
from postcalib.apus.common import get_log_func
from postcalib import qa
from postcalib.qr.podi_cython import sigma_clip_median


def main(*args, **kwargs):
    log = get_log_func(default_level='debug', **kwargs)
    if not args:
        args = sys.argv[1:]
        images = args[:-1]
        out_file = args[-1]
    else:
        images, out_file = args
    hdulists = [fits.open(image, memmap=True) for image in images]
    layouts = [get_layout(h) for h in hdulists]
    # get 5odi if any, to enable combining podi with 5odi
    for hdulist, layout in zip(hdulists, layouts):
        if layout.instru == '5odi':
            break  # 5odi
    else:
        hdulist = hdulists[0]
        layout = layouts[0]  # podi

    otas = layout.ota_order

    # memory_limit = 4  # G
    # _cpu_count = int(memory_limit / (len(images) * 0.1))
    # _cpu_count = min(cpu_count(), _cpu_count)
    # if _cpu_count == 0:
    #     _cpu_count = 1
    # log("using {} CPUs".format(_cpu_count))
    # pool = Pool(_cpu_count)
    # data_dict = dict(pool.map_async(
    #         partial(mp_worker,
    #                 images=images, layout=layout, kwargs=kwargs),
    #         otas).get(9999999))
    data_dict = []
    for ota in otas:
        data_dict.append(mp_worker(
            ota, images=images, layout=layout, kwargs=kwargs))
    data_dict = dict(data_dict)
    for ota in otas:
        ext = layout.get_ota_ext(ota)
        log("write to ext {} OTA {}".format(ext, ota))
        hdulist[ext].data = data_dict[ota]
    hdulist.writeto(out_file, overwrite=True)
    pr = qa.create_preview(
            hdulist=hdulist, filename=out_file, delete_data=True)
    pr.save()


def smooth(*args, **kwargs):
    log = get_log_func(default_level='debug', **kwargs)
    if not args:
        args = sys.argv[1:]
    in_file, out_file = args
    log("smoothing {}".format(in_file))
    hdulist = fits.open(in_file, memmap=True)
    layout = get_layout(hdulist)
    otas = layout.ota_order

    # pool = Pool(cpu_count())
    # data_dict = dict(pool.map_async(
    #         partial(smooth_tile,
    #                 image=in_file, layout=layout,
    #                 width=8, kwargs=kwargs),
    #         otas).get(9999999))
    data_dict = []
    for ota in otas:
        data_dict.append(smooth_tile(
            ota, image=in_file, layout=layout, width=8, kwargs=kwargs))
    data_dict = dict(data_dict)

    for ota in otas:
        ext = layout.get_ota_ext(ota)
        log("write to ext {} OTA {}".format(ext, ota))
        hdulist[ext].data = data_dict[ota]
    hdulist.writeto(out_file, overwrite=True)
    pr = qa.create_preview(
            hdulist=hdulist, filename=out_file, delete_data=True)
    pr.save()


@mp_traceback
def smooth_tile(ota, image, layout, width, kwargs):
    log = get_log_func(default_level='debug', **kwargs)
    ext = layout.get_ota_ext(ota)
    log("work on ext {} OTA {}".format(ext, ota))
    hdulist = fits.open(image, memmap=True)
    data = hdulist[ext].data
    # min_count = max(int(width * 0.1), 3)
    # 2sigma clip
    for cj, ci in itertools.product(
            range(layout.NCX), range(layout.NCY)):
        # print "smoothing OTA {0} cell {1}{2}".format(otaxy, cj, ci)
        (l, r), (b, t) = layout.get_cell_rect(0, 0, cj, ci)
        l, r, b, t = map(int, (l, r, b, t))
        cell = data[b:t, l:r]
        if np.all(np.isnan(cell)):
            log("cell {0}{1} skip due to all NAN".format(cj, ci))
            continue
        # cell = sigma_clip(cell, sigma=2,
        #                   cenfunc=bn.nanmedian, stdfunc=bn.nanstd)
        data[b:t, l:r] = cell
    # box_2D_kernel = Box2DKernel(width)
    # _data = convolve(data, box_2D_kernel)
    # _data = convolve_fft(fix_data(data), box_2D_kernel, allow_huge=True,
    #                      fill_value=np.nanmedian(data))
    _data = uniform_filter(fix_data(data, layout, log),
                           width, mode='nearest')

    # _data = bn.move_mean(data, width, min_count=min_count, axis=1)
    # _data = bn.move_mean(_data, width, min_count=min_count, axis=0)
    for cj, ci in itertools.product(range(layout.NCX), range(layout.NCY)):
        # print "smoothing OTA {0} cell {1}{2}".format(otaxy, cj, ci)
        (l, r), (b, t) = layout.get_cell_rect(0, 0, cj, ci)
        l, r, b, t = map(int, (l, r, b, t))
        if np.all(np.isnan(data[b:t, l:r])):
            log("cell {0}{1} skip due to all NAN".format(cj, ci))
            continue
        else:
            data[b:t, l:r] = _data[b:t, l:r]
    # cell = data[b:t, l:r]
    # data[b:t, l:r] = cell
    # return data
    return ota, data


@mp_traceback
def mp_worker(ota, images, layout, kwargs):
    log = get_log_func(default_level='debug', **kwargs)
    ext = layout.get_ota_ext(ota)
    log("work on ext {} OTA {}".format(ext, ota))

    hdulists = [fits.open(image, memmap=True) for image in images]
    # get bkg mask
    regmask_file = os.path.join(
            os.path.dirname(__file__),
            'pupilmask', 'pg_large{}.reg'.format(ota))
    if os.path.exists(regmask_file):
        log("read pupil region mask {}".format(regmask_file))
        hdu = fits.ImageHDU(data=hdulists[0][ext].data)
        mask = pyregion.open(regmask_file).get_mask(hdu=hdu)
        mask = ~mask
    else:
        mask = None

    # scale
    data = []
    for hdulist in hdulists:
        data.append(scale_to_bkg(hdulist[ext].data, mask))
    data = np.dstack(data).astype('d')
    # free some memory
    for hdulist in hdulists:
        del hdulist[ext].data
    combined = np.empty(
            (data.shape[0] * data.shape[1]),
            dtype=data.dtype)
    sigma_clip_median(data.reshape(-1, len(hdulists)), combined)

    return ota, combined.reshape(data.shape[:2])


def scale_to_bkg(data, bkgmask):
    if bkgmask is not None:
        _data = data[bkgmask]
    else:
        _data = data
    mode = np.nanmedian(_data) * 3. - np.nanmean(_data) * 2
    data /= mode
    return data


def fix_data(data, wl, log):
    """generate a LF interpolation of data for nan data"""
    # bin each cell by nbin
    nbin = 16
    bw = wl.CW / nbin
    bh = wl.CH / nbin
    # bs = np.hypot(bw, bh) * 0.5
    samp_v = np.empty((wl.NCY * nbin, wl.NCX * nbin))
    samp_i = []
    samp_j = []
    for cj, ci in itertools.product(range(wl.NCX), range(wl.NCY)):
        (cl, cr), (cb, ct) = wl.get_cell_rect(0, 0, cj, ci)
        for bj, bi in itertools.product(range(nbin), repeat=2):
            bl = cl + bj * bw
            br = bl + bw
            bb = cb + bi * bh
            bt = bb + bh
            bx = 0.5 * (bl + br)
            by = 0.5 * (bb + bt)
            bb, bt, bl, br = map(int, (bb, bt, bl, br))
            v = np.nanmedian(data[bb:bt, bl:br])
            samp_v[ci * nbin + bi, cj * nbin + bj] = v
            if cj == ci and bj == bi:
                samp_j.append(bx)
                samp_i.append(by)
    # pad value to handle edges
    samp_i.insert(0, 2 * samp_i[0] - samp_i[1] - 50)
    samp_i.append(2 * samp_i[-1] - samp_i[-2] + 50)
    samp_j.insert(0, 2 * samp_j[0] - samp_j[1] - 50)
    samp_j.append(2 * samp_j[-1] - samp_j[-2] + 50)
    # sigma clip the samp_v for the padding value
    # _samp_v, _, _ = sigmaclip(samp_v[~np.isnan(samp_v)], 2, 2)
    _samp_v = samp_v
    padval = np.nanmedian(_samp_v)
    if np.isnan(padval):
        log("warning", "not able to get an estimate of the padval")
    samp_v[np.isnan(samp_v)] = padval
    _samp_v = np.ones((samp_v.shape[0] + 2, samp_v.shape[1] + 2)) * padval
    _samp_v[1:-1, 1:-1] = samp_v
    samp_v = _samp_v
    # print(samp_i, samp_j)

    samp_i, samp_j = np.meshgrid(np.array(samp_i),
                                 np.array(samp_j),
                                 indexing='ij')
    # smooth the anchor array before spline
    # box_2D_kernel = Box2DKernel(3)
    # samp_v = convolve(samp_v, box_2D_kernel, boundary='extend')
    # samp_v = samp_v - np.nanmedian(samp_v)
    m = ~np.isnan(samp_v)
    kx = ky = 5
    log(
        "size of sampling array: {0} {1}".format(samp_v[m].size, samp_v.size))
    if len(samp_v[m]) < (kx + 1) * (ky + 1):
        # raise RuntimeError("unable to create LF map")
        log("unable to create LF map, use median")
        _lf = np.ones_like(data) * np.median(data)
    else:
        # do spline interpolation
        spline = interpolate.bisplrep(
                samp_i[m], samp_j[m], samp_v[m], kx=kx, ky=ky)
        ii, jj = np.mgrid[0:data.shape[0], 0:data.shape[1]]
        _lf = interpolate.bisplev(ii[:, 0], jj[0, :], spline).reshape(
                data.shape)
        # mask bad cell off on lf
        # lf = np.ones_like(_lf) * np.nan
        # for cj, ci in itertools.product(range(wl.NCX), range(wl.NCY)):
        #     (cl, cr), (cb, ct) = wl.get_cell_rect(0, 0, cj, ci)
        #     if np.any(~np.isnan(data[cb:ct, cl:cr])):
        #         lf[cb:ct, cl:cr] = _lf[cb:ct, cl:cr]
    edgemask = np.zeros_like(data, dtype=bool)
    # edgemask[50:-50, 50:-50] = True
    edgemask[:, :] = True
    _lf[~np.isnan(data) & edgemask] = data[~np.isnan(data) & edgemask]
    return _lf
