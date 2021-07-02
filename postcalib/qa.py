#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2017-08-14 14:54
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
qa.py
"""


from .wiyn import get_layout
import os
import logging
import warnings
from astropy.visualization import ZScaleInterval  # , PercentileInterval
from astropy.visualization.mpl_normalize import ImageNormalize
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clip  # , mad_std
import itertools

import matplotlib
matplotlib.use("agg")
# import matplotlib.image as mimg
import matplotlib.pyplot as plt  # noqa: E402


class Preview(object):
    """
    A thin wrapper for previewing a fits image
    """
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)

    def save(self, name=None):

        if name is None:
            savename = self.default_savename
        elif os.path.isdir(name):
            savename = os.path.join(name, os.path.basename(
                self.default_savename))
        else:
            savename = name
        self.fig.savefig(savename, pad_inches=0.0, bbox_inches='tight')
        self.logger.info('save preview to {0}'.format(savename))

    def show(self):
        # plt.show()
        raise NotImplementedError()


def create_preview(hdulist=None, binning=8, filename=None, delete_data=True):
    logger = logging.getLogger("qa.preview")

    if hdulist is None and filename is not None:
        hdulist = fits.open(filename, memmap=True)
        close_after = True
    elif hdulist is not None:
        if filename is None:
            filename = hdulist.filename()
        close_after = False
    else:
        raise ValueError("there has to be at least one of hdulist or filename")
    filenamebase = os.path.basename(filename).rsplit(".fz", 1)[0].rsplit(
            ".fits", 1)[0]
    dirname = os.path.dirname(filename)

    logger.info("create preview for {}".format(filename))
    layout = get_layout(hdulist, binning=binning)
    logger.info("instrument {}, number of OTAs {}, preview binning {}".format(
        layout.instru, layout.n_ota, binning))

    # binning mosaic
    (_, size_x), (_, size_y) = layout.get_ota_rect(
            layout.n_ota_x - 1, layout.n_ota_y - 1)
    size_x, size_y = map(int, (size_x, size_y))
    preview_data = np.empty((size_y, size_x), dtype='f') * np.NAN
    binned_data = []
    for ext, hdu in layout.enumerate(hdulist):
        bin_shape = list(
                map(int, (layout.OH, binning, layout.OW, binning)))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            binned = np.nanmean(
                    np.nanmean(
                        np.reshape(hdu.data, bin_shape), axis=-1), axis=1)
        # binned[0:100, 0:100] = np.nan
        interval = ZScaleInterval()
        # interval = PercentileInterval(99)
        try:
            vmin, vmax = interval.get_limits(binned)
        except IndexError:
            vmin, vmax = np.min(binned), np.max(binned)
        norm = ImageNormalize(vmin=vmin, vmax=vmax)
        _ox, _oy = layout.get_ext_ota(ext, return_tuple=True)
        ox = _ox - layout.ota_range_x[0]  # offset to corner ota
        oy = _oy - layout.ota_range_y[0]
        (l, r), (b, t) = layout.get_ota_rect(ox, oy)
        l, r, b, t = map(int, (l, r, b, t))
        norm_binned = np.full_like(binned, np.nan)
        norm_binned[~np.isnan(binned)] = norm(binned[~np.isnan(binned)])
        preview_data[b:t, l:r] = norm_binned
        binned_data.append(binned)
        if hdulist._file.memmap and delete_data:
            del hdu.data  # possibly free some memory
    binned_data = np.dstack(binned_data)

    guide_otas = find_guide_otas(binned_data, layout, thresh=10, logger=logger)

    # preview figure
    fig = plt.figure(figsize=(2 * layout.n_ota_x, 2 * layout.n_ota_y))
    ax = fig.add_subplot(1, 1, 1)
    ax.imshow(preview_data, vmin=0, vmax=1, origin='lower', aspect=1)
    tick_x, tick_y = layout.get_ota_bins()
    tick_x = [(x[0] + x[1]) * 0.5 for x in tick_x][:layout.n_ota_x]
    tick_y = [(y[0] + y[1]) * 0.5 for y in tick_y][:layout.n_ota_y]
    ax.yaxis.set_ticks(tick_y)
    ax.yaxis.set_ticklabels(['Y{0}'.format(j) for j
                             in range(*layout.ota_range_y)])
    ax.xaxis.set_ticks(tick_x)
    ax.xaxis.set_ticklabels(['X{0}'.format(j) for j
                             in range(*layout.ota_range_x)])
    ax.set_title(filenamebase)

    default_savename = os.path.join(dirname, filenamebase + '.png')
    pr = Preview(
            fig=fig, ax=ax, logger=logger, default_savename=default_savename,
            binning=binning, binned_data=binned_data, guide_otas=guide_otas,
            preview_data=preview_data,
            layout=layout)
    if close_after:
        hdulist.close()
    return pr


def find_guide_otas(data, layout, thresh=10, logger=None):
    if layout.instru == 'podi':
        return []
    centers = []
    corners = []
    cs = 32 / layout.binning
    for cj, ci in itertools.product(range(layout.NCX), range(layout.NCY)):
        (cl, cr), (cb, ct) = layout.get_cell_rect(0, 0, cj, ci)
        centers.append((cl + layout.CW * 0.35, cr - layout.CW * 0.35,
                        cb + layout.CH * 0.35, ct - layout.CH * 0.35))
        corners.extend([(cl, cl + cs, cb, cb + cs),
                        (cr - cs, cr, cb, cb + cs),
                        (cl, cl + cs, ct - cs, ct),
                        (cr - cs, cr, ct - cs, ct),
                        ])
    mask = np.ones_like(data[:, :, 0]) * np.nan
    for box in centers:
        l, r, b, t = map(int, box)
        mask[b:t, l:r] = 1
    for box in corners:
        l, r, b, t = map(int, box)
        mask[b:t, l:r] = 0
    otas = []
    offsets = []
    for i in range(layout.n_ota):
        _data = data[:, :, i]
        center_data = sigma_clip(
                _data[(mask == 1) & (~np.isnan(_data))],
                # stdfunc=mad_std,
                sigma=3)
        corner_data = sigma_clip(
                _data[(mask == 0) & (~np.isnan(_data))],
                # stdfunc=mad_std,
                sigma=3)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            offset = (np.median(corner_data) - np.median(center_data)) / \
                np.std(center_data)
        offsets.append(offset)
        otas.append(layout.get_ext_ota(i + 1))
        # print("OTA {} has corner excess {}".format(
        #     layout.get_ext_ota(i + 1), offset))
    guide_otas = [o for i, o in enumerate(otas) if offsets[i] > thresh]
    if len(guide_otas) == 0:
        if logger is None:
            def log(mesg):
                print(mesg)
        else:
            log = logger.warning
        guide_otas = sorted(otas, key=lambda o: offsets[otas.index(o)])[-1:]
        log("unable to locate guide ota,"
            " use the most likely one with thresh {}".format(
                offsets[otas.index(guide_otas[0])]))
        # log("thresholds:\n{}".format("\n".join(
        #     '{}: {}'.format(otaxy, thresh) for otaxy, thresh in zip(
        #         layout.ota_order, offsets))))
    return guide_otas


# if __name__ == "__main__":
#     import sys
#     # from astropy.io import fits
#     with fits.open(sys.argv[1]) as hdulist:
#         create_preview(hdulist)
