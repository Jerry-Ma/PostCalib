#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2017-08-20 01:17
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
phot_calib.py

This script handles the photometry calibration of ODI images.
It follows the algorithm detailed in Doc 01, namely, the
least chisq based self-calibration.

The routine takes an collection of matched catalogs with
sources have both the instrument magnitude and the reference
truth magnitude, tries to determine the zero points of each
OTA by minimizing the differences.

The zero points of each OTA consists of two components, one is
for the global transparency, which is uniform across the entire
image, the other is for per OTA variation, which differs OTA
by OTA, but the intra-OTA differences are the same for different
exposures. In addition, the color term is also a free parameter.

Inputs
------
matched catalogs: list of ASCII catalogs
    The entries in this catalog should have both the measured and
    the truth magnitude from reference catalog (SDSS, by default).
    Moreover, it should contain the extension number of each source
    By default, the comparison is done using the value of MAG_AUTO.

Outputs
-------
header: ASCII header
    The header file contains FLXSCALE and COLOR, to be used for
    SWarp.
"""


from __future__ import (absolute_import, division, print_function)
import os
import re
import sys
import glob

import numpy as np
import lmfit
import time
from scipy.stats import sigmaclip
from astropy.time import Time
from astropy.io import fits
from astropy.table import Table, Column, vstack
# from functools import partial
import itertools
import subprocess

# from multiprocessing import cpu_count, Pool

# from scipy.ndimage import uniform_filter  # , gaussian_filter, median_filter
import pyregion
# from scipy import interpolate

# from postcalib.utils import mp_traceback
from postcalib.wiyn import get_layout
from postcalib.apus.common import get_log_func
# from postcalib import qa

from cycler import cycler
import matplotlib
matplotlib.use("agg")
# import matplotlib.image as mimg
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import cm  # noqa: E402
from mpl_toolkits.axes_grid.inset_locator import inset_axes  # noqa: E402
import matplotlib.colors as mc  # noqa: E402
from matplotlib import rc  # noqa: E402


def main(*args, **kwargs):
    log = get_log_func(default_level='debug', **kwargs)
    if not args:
        args = sys.argv[1:]
        cat_files = args[0:-1:2]
        image_files = args[1:-1:2]
        master_file = args[-1]
    else:
        in_files, master_file = args
        cat_files, image_files = zip(*in_files)
    log("self calib on\n{}".format('\n'.join(cat_files)))
    # print(kwargs['reg_phot'], master_file)
    parsed_outfile = re.match(kwargs['reg_phot'],
                              os.path.basename(master_file)).groupdict()
    photgrp = parsed_outfile['ppflag']
    band = parsed_outfile['band']
    log("{} cats in phot group {} band {}".format(
        len(cat_files), photgrp, band))
    model_flags = kwargs['phot_model_flags']
    hdr_suffix = kwargs['phot_hdr_suffix']
    log("phot model flags {}, output to .{}".format(model_flags, hdr_suffix))

    with fits.open(image_files[0], memmap=True) as hdulist:
        layout = get_layout(hdulist)
    ck = get_context(photgrp, band, layout)
    master = get_calibrated_master(cat_files, ck, model_flags, kwargs)
    # generate fluxscale header
    catind = np.unique(master['catind'])
    for i in catind:
        subbulk = master[master['catind'] == i]
        catfile = subbulk['catfile'][0]
        log("processing catalog {}".format(catfile))
        airmass = subbulk[ck['airmass']][0]
        # figure out hdrfile name from catfile through indexing
        hdrfile = look_up_images(image_files, catfile).rsplit(
                '.fits', 1)[0] + '.' + hdr_suffix
        log("processing swarp header {0}".format(hdrfile))
        valid_s = []
        headers = []
        for ext in range(1, ck['n_ota_x'] * ck['n_ota_y'] + 1):
            ota_xy = layout.get_ext_ota(ext)
            extentry = subbulk[subbulk[ck['otaxy']] == ota_xy]
            if len(extentry) == 0:
                k = b = 99
                s = '<replace>'
            else:
                extentry = extentry[0]
                k = extentry['catkins']
                b = extentry['catbota'] + extentry['catbair'] * airmass + \
                    extentry['catbcat']
                s = str(10 ** ((25.0 - b) / 2.5))
                valid_s.append(float(s))
            headers.append("{0:8s}={1:>22s}{2:s}\n".format(
                "MYCOLOR", str(k), " / color term {0}-{1}"
                .format(*ck['cband'])))
            headers.append("{0:8s}={1:>22s}{2:s}\n".format(
                "MYZEROP", str(b), " / self-calibrated zero point"))
            headers.append("{0:8s}={1:>22s}{2:s}\n".format(
                "FLXSCALE", str(s), " / sc flux scale to zp=25"))
            headers.append("END     \n")
        mean_s = np.median(valid_s)
        with open(hdrfile, 'w') as fo:
            for hdr in headers:
                if '<replace>' in hdr:
                    fo.write("{0:8s}={1:>22s}{2:s}\n".format(
                        "FLXSCALE", str(mean_s),
                        " / sc flux scale to zp=25 (median)"))
                else:
                    fo.write(hdr)
    # create QA plots
    qa_master(master, ck, master_file.rsplit(".cat", 1)[0] + ".png", kwargs)

    log("master calib file saved: {0}".format(master_file))
    master.write(master_file, format='ascii.commented_header')


def match_refcat(*args, **kwargs):
    log = get_log_func(default_level='debug', **kwargs)
    if not args:
        args = sys.argv[1:]
        cat_file, ref_file, image_file, out_file = args
    else:
        (cat_file, ref_file, image_file), out_file = args
    log("match {} to {}".format(cat_file, ref_file))

    # figure out band
    parsed_filename = re.match(
            kwargs['reg_inputs'], os.path.basename(cat_file)
            ).groupdict()
    # print(parsed_filename)
    band = parsed_filename['band']
    log("use band {}".format(band))
    if band == 'u':
        log('no ps1 data for u band')
        refkey = 'sdss'
    else:
        # figure out coverage from recat
        refcat = Table.read(ref_file, format='ascii.commented_header')
        samp_lim = (17, 18)
        samp_counts = {}
        for key in ['sdss', 'ps1']:
            samp_counts[key] = len(
                refcat[
                    (refcat['{}_{}'.format(band, key)] > samp_lim[0]) &
                    (refcat['{}_{}'.format(band, key)] < samp_lim[1])
                    ])
        if samp_counts['sdss'] > 1.1 * samp_counts['ps1']:
            log("use sdss ({sdss}) over ps1 ({ps1})".format(**samp_counts))
            refkey = 'sdss'
        elif samp_counts['ps1'] > 1.1 * samp_counts['sdss']:
            log("use ps1 ({ps1}) over sdss ({sdss})".format(**samp_counts))
            refkey = 'ps1'
        else:
            log("use preferred sdss ({sdss})".format(**samp_counts))
            refkey = 'sdss'

    icmd1 = 'select "{band}_{refkey} < 90. && {band}_{refkey} > 0' \
            ' && err_{band}_{refkey} <= 0.10857"' \
        .format(band=band, refkey=refkey)
    icmd2 = 'select "MAG_AUTO < 90. && MAGERR_AUTO >= 0.001 && FLAGS == 0"'
    stilts_match = """{stilts_cmd}
                   tmatch2
                   matcher=sky
                   in1={ref_file}
                   ifmt1=ascii
                   icmd1={icmd1}
                   in2={cat_file}
                   ifmt2=ascii
                   icmd2={icmd2}
                   values1=ra_{refkey} dec_{refkey}
                   values2=ALPHA_J2000 DELTA_J2000
                   params=1.2
                   join=1and2
                   out={out_file}
                   ofmt=ascii""".format(
                        stilts_cmd=kwargs['stilts_cmd'], **locals())
    subprocess.check_call(map(str.strip, stilts_match.split("\n")))

    tbl = Table.read(out_file, format='ascii.commented_header')
    # get ready for self-calibration
    log("{} {} stars in total".format(len(tbl), refkey))

    copycol = [
            ('ra_{refkey}', 'REF_RA'), ('dec_{refkey}', 'REF_DEC'),
            ('u_{refkey}', 'REF_MAG_U'), ('err_u_{refkey}', 'REF_ERR_U'),
            ('g_{refkey}', 'REF_MAG_G'), ('err_g_{refkey}', 'REF_ERR_G'),
            ('r_{refkey}', 'REF_MAG_R'), ('err_r_{refkey}', 'REF_ERR_R'),
            ('i_{refkey}', 'REF_MAG_I'), ('err_i_{refkey}', 'REF_ERR_I'),
            ('z_{refkey}', 'REF_MAG_Z'), ('err_z_{refkey}', 'REF_ERR_Z'),
            ('ALPHA_J2000', 'ODI_RA'), ('DELTA_J2000', 'ODI_DEC'),
            ('MAG_AUTO', 'ODI_MAG_AUTO'), ('MAGERR_AUTO', 'ODI_ERR_AUTO'),
            ('XWIN_IMAGE', 'ODI_X'), ('YWIN_IMAGE', 'ODI_Y'),
            ]
    for oc, nc in copycol:
        oc = oc.format(refkey=refkey)
        if oc in tbl.colnames:
            tbl[nc] = tbl[oc]
        else:
            log('warning', 'column {} does not exist for refcat'.format(oc))
    # get layout
    with fits.open(image_file, memmap=True) as hdulist:
        layout = get_layout(hdulist)
        # header columns
        for key in ['AIRMASS', 'EXPMEAS']:
            col = Column([hdulist[0].header[key], ] * len(tbl), name=key)
            tbl.add_column(col)
        # get mjd time
        obstime = Time(
                hdulist[0].header['DATE-MID'], format='isot', scale='utc')
        col_time = Column([obstime.mjd, ] * len(tbl), name='MJD')
        tbl.add_column(col_time)
    # odixy column
    ota_xy = [layout.get_ext_ota(ext) for ext in tbl['EXT_NUMBER']]
    col_ota = Column(ota_xy, name='ODI_OTA')
    tbl.add_column(col_ota)

    log("save to matched cat {}".format(out_file))
    tbl.write(out_file, format='ascii.commented_header')


def cleanup(*args, **kwargs):
    """
    Cleanup the input catalog so that sources fall near the edges are
    rejected. The remaining sources are like to have good photometry.
    """
    log = get_log_func(default_level='debug', **kwargs)
    if not args:
        args = sys.argv[1:]
        cat_file, image_file, out_file = args
    else:
        (cat_file, image_file), out_file = args
    log("cleanup {}".format(cat_file))
    # out_reg = os.path.splitext(out_file)[0].'.reg'

    cell_edge_pad = bpm_edge_pad = 25
    ota_edge_pad = 50
    cat = Table.read(cat_file, format='ascii.sextractor')

    # figure out layout
    # mask out region
    with fits.open(image_file, memmap=True) as hdulist:
        layout = get_layout(hdulist)
        log("layout {}".format(layout.instru))
        skymask_dir = kwargs['skymask_dir']
        obsid = hdulist[0].header['OBSID']
        regions = glob.glob(os.path.join(skymask_dir, 'skymask.reg'))
        regions.extend(
                glob.glob(
                    os.path.join(skymask_dir, '*{}*.reg'.format(obsid))))
        log("use DS9 region masks\n{}".format('\n'.join(regions)))
        regmask = get_regmask(
                cat['XWIN_IMAGE'], cat['YWIN_IMAGE'], cat['EXT_NUMBER'],
                hdulist, regions, **kwargs)
    # remove edge source around cell
    edgemask = get_edgemask(
            cat['XWIN_IMAGE'], cat['YWIN_IMAGE'], cell_edge_pad, layout)
    # remove edge source around bpm
    bpmask = get_bpmask(
            cat['XWIN_IMAGE'], cat['YWIN_IMAGE'], cat['EXT_NUMBER'],
            bpm_edge_pad, layout)
    # remove edge source around OTA
    otaedgemask = get_ota_edgemask(
            cat['XWIN_IMAGE'], cat['YWIN_IMAGE'], ota_edge_pad, layout)
    cat = cat[edgemask & bpmask & otaedgemask & regmask]
    log("save to {}".format(out_file))
    cat.write(out_file, format='ascii.commented_header')


def get_regmask(xs, ys, exts, hdulist, regions, **kwargs):
    log = get_log_func(default_level='debug', **kwargs)
    layout = get_layout(hdulist)
    filters = {ota: [] for ota in layout.ota_order}
    for ext, hdu in layout.enumerate(hdulist):
        ota = layout.get_ext_ota(ext)
        # create a fresh wcs from keys
        tmp_hdu = fits.ImageHDU(data=hdu.data)
        # add basic wcs
        for key in ['EQUINOX', 'CTYPE1', 'CTYPE2', 'CRPIX1', 'CRPIX2',
                    'CRVAL1', 'CRVAL2', 'CUNIT1', 'CUNIT2', 'CD1_1', 'CD2_1',
                    'CD1_2', 'CD2_2']:
            tmp_hdu.header[key] = hdu.header[key]
        for region in regions:
            try:
                filter_ = pyregion.open(region).as_imagecoord(
                        tmp_hdu.header).get_filter()
                filters[ota].append(filter_)
            except ValueError as e:
                if "need more than 0 values to unpack" not in e.message:
                    log("unable to apply region mask {} due to '{}'"
                        " please check the format of the file".format(
                            region, e))
    # loop over catalog
    mask = np.ones_like(xs, dtype=bool)
    for i in range(len(xs)):
        otaxy = layout.get_ext_ota(exts[i])
        filter_list = filters[otaxy]
        for filter_ in filter_list:
            if filter_.inside1(xs[i], ys[i]):
                mask[i] = False
                break
    return mask


def get_edgemask(xs, ys, e, wl):
    mask = np.zeros_like(xs, dtype=bool)
    for cj, ci in itertools.product(range(wl.NCX), range(wl.NCY)):
        (cl, cr), (cb, ct) = wl.get_cell_rect(0, 0, cj, ci)
        mask = mask | ((xs > cl + e) & (xs < cr - e) &
                       (ys > cb + e) & (ys < ct - e))
    return mask


def get_ota_edgemask(xs, ys, e, wl):
    mask = (xs > e) & (xs < wl.OW - e) & (ys > e) & (ys < wl.OH)
    return mask


def get_bpmask(xs, ys, exts, pad, layout):
    bpmdir = os.path.join(os.path.dirname(__file__), 'bpm')
    if layout.instru == '5odi':
        bpmdir = os.path.join(bpmdir, 'odi_5x6')
    elif layout.instru == 'podi':
        bpmdir = os.path.join(bpmdir, 'podi')
    else:
        raise ValueError('ODI instru {0} not recognized'
                         .format(layout.instru))
    mask = np.ones_like(xs, dtype=bool)
    bpms = {}
    for otaxy in layout.ota_order:
        bpm = []
        bpm_file = os.path.join(bpmdir, 'bpm_xy{0}.reg'.format(otaxy))
        with open(bpm_file, 'r') as fo:
            for ln in fo.readlines():
                rect = re.match(r'box\(([0-9+-., ]+)\)', ln.strip())
                if rect is not None:
                    rect = list(map(float, rect.group(1).split(',')))
                    # print "box from bpm: {0}".format(rect)
                    bpm.append((
                        rect[0] - rect[2] * 0.5,
                        rect[0] + rect[2] * 0.5,
                        rect[1] - rect[3] * 0.5,
                        rect[1] + rect[3] * 0.5))
                else:
                    continue
        bpms[otaxy] = bpm
    for i in range(len(xs)):
        otaxy = layout.get_ext_ota(exts[i])
        bpm = bpms[otaxy]
        for l, r, b, t in bpm:
            if xs[i] > l - pad and xs[i] < r + pad and \
               ys[i] > b - pad and ys[i] < t + pad:
                mask[i] = False
    return mask


def get_context(photgrp, band, layout):
    cband = {'u': ('u', 'g'),
             'g': ('g', 'r'),
             'r': ('g', 'r'),
             'i': ('r', 'i'),
             'z': ('i', 'z'),
             }
    ck = {
        'smag': 'REF_MAG_{0}'.format(band.upper()),
        'semag': 'REF_ERR_{0}'.format(band.upper()),
        'smaglims': dict(z=(0, 19), u=(0, 19), i=(0, 20)).get(band, (0, 20)),
        'mag': 'ODI_MAG_AUTO',
        'emag': 'ODI_ERR_AUTO',
        # 'mag': 'MAG_APER_3',
        # 'emag': 'MAGERR_APER_3',
        'otaxy': 'ODI_OTA',
        'photgrp': photgrp,
        'instru': layout.instru,
        'band': band,
        'cband': cband[band],
        'cmag1': 'REF_MAG_{0}'.format(cband[band][0].upper()),
        'cmag2': 'REF_MAG_{0}'.format(cband[band][1].upper()),
        'clip0': 5.0,
        'clipsc': 3.0,
        'airmass': 'AIRMASS',
        'kins': dict(g=0.14, z=-0.1342).get(band, 0.),
        'n_ota_x': layout.n_ota_x,
        'n_ota_y': layout.n_ota_y,
        'ota_range_x': layout.ota_range_x,
        'ota_range_y': layout.ota_range_y,
        }
    return ck


def get_calibrated_master(in_cats, ck, model_flags, kwargs):
    log = get_log_func(default_level='debug', **kwargs)
    log("run self calibrate with photgrp {0}".format(ck['photgrp']))
    master = get_master_calib(in_cats, ck, kwargs)
    log("self calibration with flags: {0}".format(model_flags))
    start_time = time.time()
    sc_ret = self_calibrate(master, ck, model_flags=model_flags)
    log("self calibration finished after {0}s"
        .format(time.time() - start_time))
    # save the data to master_calib
    outpvals, outpuncs, clipmask, rstd, rcstd = sc_ret
    master['catkins'] = outpvals['kins']
    master['catkinsunc'] = outpuncs['kins']
    master['catkair'] = outpvals['kair']
    master['catkairunc'] = outpuncs['kair']
    master['catbota'] = np.array([outpvals['b{0:.0f}'.format(o)]
                                  for o in master[ck['otaxy']]])
    master['catbotaunc'] = np.array([outpuncs['b{0:.0f}'.format(o)]
                                     for o in master[ck['otaxy']]])
    master['catbair'] = outpvals['bair']
    master['catbairunc'] = outpuncs['bair']
    master['catbcat'] = np.array([outpvals['bcat{0:.0f}'.format(o)]
                                  for o in master['catind']])
    master['catbcatunc'] = np.array([outpuncs['bcat{0:.0f}'.format(o)]
                                     for o in master['catind']])
    master['catclipflag'][clipmask] = 1
    master['catres'] = rstd
    master['catclipres'] = rcstd
    return master


def get_master_calib(tablelist, ck, kwargs):
    '''get calibration table files from input file list'''
    log = get_log_func(default_level='debug', **kwargs)
    # load table and compile to a master calib table
    log('create master')
    master = []
    for i, fname in enumerate(tablelist):
        try:
            tbl = Table.read(fname, format='ascii.commented_header')
        except IOError:
            log("table {0} does not exist".format(fname))
            continue
        # global zero-th order zp by a 5 sigma clip
        zp = tbl[ck['smag']] - tbl[ck['mag']]
        clipped, lolim, uplim = sigmaclip(zp, ck['clip0'], ck['clip0'])
        zp0 = np.median(clipped)
        log('zero-th order zp clipping ({0:.1f} sigma):'
            ' {1:+.3f} {2:.3f} {3:+.3f}'
            .format(ck['clip0'], lolim - zp0, zp0, uplim - zp0))
        smag_lo, smag_up = ck['smaglims']
        qualitymask = (smag_lo < tbl[ck['smag']]) & \
                      (tbl[ck['smag']] < smag_up) & \
                      (zp > lolim) & (zp < uplim)  # & \
        tbl = tbl[qualitymask]
        n = len(tbl)
        if n < 10:
            log("table {0} contains < 10 objects".format(fname))
            continue
        log('read in {0} ({1})'.format(fname, n))
        # append the table level info
        col_ind = Column([i] * n, name='catind')
        col_cfn = Column([os.path.basename(fname)] * n, name='catfile')
        col_zp0 = Column([zp0] * n, name='catzp0')
        for col in [col_ind, col_cfn, col_zp0]:
            tbl.add_column(col)
        # place holder for the self calibration result
        # (kins + kair * X) * (color) + bota + bair * X + bcat
        for name in [
                'catkins', 'catkinsunc', 'catkair', 'catkairunc',
                'catbota', 'catbotaunc', 'catbair', 'catbairunc',
                'catbcat', 'catbcatunc', 'catclipflag',
                'catres', 'catclipres'
                ]:
            col = Column([0.0] * n, name=name)
            tbl.add_column(col)
        # col_kins = Column([0.0] * n, name='catkins')
        # col_kair = Column([0.0] * n, name='catkair')
        # col_bota = Column([0.0] * n, name='catbota')
        # col_bair = Column([0.0] * n, name='catbair')
        # col_bcat = Column([0.0] * n, name='catbcat')
        # for col in [col_ind, col_cfn, col_zp0,
        #             col_kins, col_kair, col_bota, col_bair, col_bcat]:
        #     tbl.add_column(col)
        master.append(tbl)
    master = vstack(master, join_type='outer')
    return master


def sc_func(params, color, zp, zperr, otaxy, catind, airmass):
    # (kins + kair * X) * (color) + bota + bair * X + bcat
    bota = np.array([params['b{0:.0f}'.format(o)].value for o in otaxy])
    bcat = np.array([params['bcat{0:.0f}'.format(o)].value for o in catind])
    model = (params['kins'].value + params['kair'].value * airmass) * color \
        + bota + params['bair'].value * airmass + bcat
    return (model - zp) / zperr


def self_calibrate(bulk, ck, model_flags=('color', 'ota', 'cat', 'airmass')):
    '''
    For each exposure/ota, fit an offset so that the over all
    dispersion is minimized
    '''
    color = bulk[ck['cmag1']] - bulk[ck['cmag2']]
    zp1 = bulk[ck['smag']] - bulk[ck['mag']]
    zp1err = np.hypot(bulk[ck['semag']], bulk[ck['emag']])
    params = lmfit.Parameters()
    params.add('kins', value=ck['kins'], vary='color' in model_flags)
    params.add('kair', value=0., vary=False)
    params.add('bair', value=0, vary='airmass' in model_flags)
    otaxy = bulk[ck['otaxy']]
    catind = bulk['catind']
    airmass = bulk[ck['airmass']]
    for xy in np.unique(otaxy):
        vary = False if xy == 33 else 'ota' in model_flags
        params.add('b{0:.0f}'.format(xy), value=np.mean(zp1), vary=vary)
    for ind in np.unique(catind):
        vary = 'cat' in model_flags
        params.add('bcat{0:.0f}'.format(ind), 0.0, vary=vary)
    outparams = lmfit.minimize(
            sc_func, params, args=(color, zp1, zp1err,
                                   otaxy, catind, airmass)
            ).params
    # do a sigma clipping on the de-trended data
    residue = sc_func(
            outparams, color, zp1, zp1err, otaxy, catind, airmass) * zp1err
    clipped, lo, up = sigmaclip(residue, ck['clipsc'], ck['clipsc'])
    clipmask = (residue >= lo) & (residue <= up)
    # do the fitting again
    outparams = lmfit.minimize(
            sc_func, outparams,
            args=(color[clipmask], zp1[clipmask], zp1err[clipmask],
                  otaxy[clipmask], catind[clipmask],
                  airmass[clipmask])).params

    parvals = outparams.valuesdict()
    paruncs = {k: outparams[k].stderr for k in parvals.keys()}
    return parvals, paruncs, clipmask, np.std(residue), np.std(clipped)


def look_up_images(images, cat):
    for image in images:
        if os.path.basename(image).rsplit('.fits', 1)[0] in cat:
            return image
    else:
        raise RuntimeError("unable to find image for {}".format(cat))


def qa_master(master, ck, savename, kwargs):
    log = get_log_func(default_level='debug', **kwargs)
    # load table and compile to a master calib table
    log('create QA plots')
    set_color()
    cband = ck['cband']
    fig = plt.figure(figsize=(18, 10))
    ax = fig.add_subplot(1, 1, 1)

    ax.set_xlabel(r'REF {0} - {1} (mag)'.format(*cband))
    ax.set_ylabel(r'$m_{REF} - m_{ins} - k_{ins} \times (%s - %s) '
                  r'- b_{OTA} - b_{X}X - b_{cat}$ (mag)' % cband)

    plot_zp_color(ax, master, ck, kwargs)
    log("save plot to {}".format(savename))
    plt.savefig(savename, bbox_inches='tight')


def plot_zp_color(ax, bulk, ck, kwargs):
    log = get_log_func(default_level='debug', **kwargs)
    ax.set_xlim((-0.2, np.max(bulk['catind']) + 1))
    ax.set_ylim((-0.7, 1.))

    color = bulk[ck['cmag1']] - bulk[ck['cmag2']]
    zp1 = bulk[ck['smag']] - bulk[ck['mag']]
    zp1err = np.hypot(bulk[ck['semag']], bulk[ck['emag']])
    airmass = bulk[ck['airmass']]
    otaxy = bulk[ck['otaxy']]
    catind = bulk['catind']

    pltkw = dict(fmt='o', ms=3, capsize=0, mew=0)
    pltkw_ex = dict(fmt='D', ms=2, capsize=0, mew=1, fillstyle='none')
    # legs = [mympl.get_dummy_leg() for _ in range(len(np.unique(catind)))]
    legs = []

    bota = bulk['catbota']
    botaunc = bulk['catbotaunc']
    bcat = bulk['catbcat']
    bcatunc = bulk['catbcatunc']
    kins = bulk['catkins'][0]
    kair = bulk['catkair'][0]
    bair = bulk['catbair'][0]
    log("residue: {0} ({1})".format(bulk['catres'][0],
                                    bulk['catclipres'][0]))
    model = (kins + kair * airmass) * color \
        + bota + bair * airmass + bcat
    res = zp1 - model
    clipmask = bulk['catclipflag'].astype(bool)
    # res = res - np.median(res)
    for ii, i in enumerate(np.unique(catind)):
        im = (bulk['catind'] == i) & clipmask
        ex = (bulk['catind'] == i) & (~clipmask)
        # (kins + kair * X) * (color) + bota + bair * X
        if np.any(im):
            leg = ax.errorbar(color[im] + bulk['catind'][im], res[im],
                              yerr=zp1err[im],
                              **pltkw)
            for lc in leg[2]:
                ecolor = lc.get_color()[0]
                eecolor = change_hsv(ecolor, s=0.2, v=0.9)
                lc.set_color(eecolor)
            shortname = re.match(r'.+_20(?:\d\d)(.+)_odi_\w\..+',
                                 bulk['catfile'][im][0]).group(1).replace(
                                         '_', '-')
            legs.append((leg, shortname))
            if np.any(ex):
                legex = ax.errorbar(color[ex] + bulk['catind'][ex], res[ex],
                                    yerr=zp1err[ex], color=eecolor, **pltkw_ex)
                for lc in legex[2]:
                    lc.set_color(eecolor)
    # xg = np.linspace(-2, 2, 10)
    ax.axhline(0.05 * np.log10(np.e) * 2.5, ls='--', color='gray')
    ax.axhline(-0.05 * np.log10(np.e) * 2.5, ls='--', color='gray')
    ax.text(-0.1, 0.0, (r'$\pm 5\%$'),
            verticalalignment='center', fontsize=20)

    # inset plot for the best fit parameters
    psize = 200. / 72.27
    px = inset_axes(ax,
                    width=psize,
                    height=psize,
                    loc=1,
                    bbox_to_anchor=(0.05, 0.05, 0.87, 0.9),
                    bbox_transform=ax.transAxes,
                    borderpad=0,
                    )
    cpx = inset_axes(ax,
                     width=psize * 0.05,
                     height=psize,
                     loc=1,
                     bbox_to_anchor=(0.05, 0.05, 0.87, 0.9),
                     bbox_transform=ax.transAxes,
                     borderpad=0,
                     )
    # get pval bota array
    bota = np.ones((ck['n_ota_y'], ck['n_ota_x']), dtype='d') * np.nan
    botaunc = np.ones((ck['n_ota_y'], ck['n_ota_x']), dtype='d') * np.nan
    for y in range(*ck['ota_range_y']):
        for x in range(*ck['ota_range_x']):
            if np.any(otaxy == 10 * x + y) > 0:
                bota[y - ck['ota_range_y'][0], x - ck['ota_range_x'][0]] \
                    = bulk['catbota'][
                                otaxy == 10 * x + y][0]
                botaunc[y - ck['ota_range_y'][0], x - ck['ota_range_x'][0]] \
                    = bulk['catbotaunc'][
                                otaxy == 10 * x + y][0]
            else:
                log("no data found for OTA {0}{1}".format(x, y))

    bota = bota[::-1, :]
    pleg = px.imshow(bota, interpolation='nearest', aspect=1, cmap=cm.coolwarm)
    cb = plt.colorbar(pleg, cax=cpx)
    cb.set_label('$b_{OTA}$')
    px.yaxis.set_ticks(range(ck['n_ota_y']))
    px.yaxis.set_ticklabels(['Y{0}'.format(j) for j
                             in reversed(range(*ck['ota_range_y']))])
    px.xaxis.set_ticks(range(ck['n_ota_x']))
    px.xaxis.set_ticklabels(['X{0}'.format(j) for j
                             in range(*ck['ota_range_x'])])
    # plot airmass scatter
    xsize = 200. / 72.27
    xx = inset_axes(ax,
                    width=xsize * 1.5,
                    height=xsize,
                    loc=9,
                    bbox_to_anchor=(0.05, 0.05, 0.87, 0.9),
                    bbox_transform=ax.transAxes,
                    borderpad=0,
                    )
    u_data = []
    for ii, i in enumerate(np.unique(catind)):
        u_x = airmass[bulk['catind'] == i][0]
        u_b = bcat[bulk['catind'] == i][0]
        u_bunc = bcatunc[bulk['catind'] == i][0]
        color = legs[ii][0][0].get_color()
        pltkw_xx = dict(pltkw, ms=9 if ii == 0 else 6, color=color)
        xx.errorbar([u_x], [u_b], yerr=[u_bunc], **pltkw_xx)
        u_data.append((u_x, u_b))
    xx.plot(*zip(*u_data), linestyle='-', color='gray', zorder=0)
    xx.set_xlim((np.min(airmass) - 0.1, np.max(airmass) + 0.1))
    xx.set_xlabel("Airmass")
    xx.set_ylabel("$b_{cat}$")

    # show label text
    label = [
            r'$k_{ins}$=%.4f' % (kins),
            r'${\langle}b_{ota}{\rangle}$=%.4f' % (np.mean(bota)),
            r'${\langle}{\delta}b_{ota}{\rangle}$=%.4f' % (np.mean(botaunc)),
            r'${\Delta}b_{ota}$=%.4f' % (np.std(bota)),
            r'$b_{X}$=%.4f' % (bair),
            r'$n_{obj}$=%d' % (len(zp1)),
            r'${\Delta}res.$=%.5f' % (np.std(res)),
            r'$n_{obj,clip}$=%d' % (len(zp1[clipmask])),
            r'${\Delta}res._{clip}$=%.5f' % (np.std(res[clipmask]))
            ]
    ax.text(0.05, 0.92, '\n'.join(label),
            transform=ax.transAxes,
            verticalalignment='top')
    ax.legend(*zip(*legs), loc='lower center', ncol=4)
    return legs


def change_hsv(c, h=None, s=None, v=None, frac=False):
    '''Quickly change the color in hsv space'''
    if isinstance(c, str):
        rgb = np.array([[mc.hex2color(c), ]])
    else:  # rgb
        rgb = np.array([[c[:3], ]])
    hsv = mc.rgb_to_hsv(rgb)
    # print c
    # print rgb
    # print hsv
    for i, j in enumerate([h, s, v]):
        if j is not None:
            if frac:
                if j < 1:
                    hsv[0][0][i] = hsv[0][0][i] * j
                else:
                    hsv[0][0][i] = hsv[0][0][i] + \
                            (1 - hsv[0][0][i]) * (1 - j)
            else:
                hsv[0][0][i] = j
    return mc.hsv_to_rgb(hsv)[0][0]


def set_color():
    kelly = np.array([
        [255, 179, 0], [128, 62, 117], [255, 104, 0], [166, 189, 215],
        [193, 0, 32], [206, 162, 98], [129, 112, 102], [0, 125, 52],
        [246, 118, 142], [0, 83, 138], [255, 122, 92], [83, 55, 122],
        [255, 142, 0], [179, 40, 81], [244, 200, 0], [127, 24, 13],
        [147, 170, 0], [89, 51, 21], [241, 58, 19], [35, 44, 22]]) / 255.
    rc('axes', prop_cycle=cycler('color', kelly))
