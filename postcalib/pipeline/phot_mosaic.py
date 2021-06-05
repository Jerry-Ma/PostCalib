#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2017-08-20 06:03
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
phot_mosaic.py

This script contains utility routines that helps create the mosaic.
"""


from __future__ import print_function, division
import os
import re
import sys
import glob
import numpy as np
from astropy.io import fits
from postcalib.apus.common import get_log_func
import subprocess
from astropy.table import Table


from scipy.stats import sigmaclip
from cycler import cycler
import matplotlib
matplotlib.use("agg")
# import matplotlib.image as mimg
import matplotlib.pyplot as plt  # noqa: E402
# from matplotlib import cm  # noqa: E402
# from mpl_toolkits.axes_grid.inset_locator import inset_axes  # noqa: E402
import matplotlib.colors as mc  # noqa: E402
from matplotlib import rc  # noqa: E402


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
            if os.path.islink(header_file):
                os.unlink(header_file)
            log("link {} to {}".format(linkname, header_file))
            os.symlink(linkname, header_file)
        else:
            log("link exists {}".format(header_file))


def _apply_weight(*args, **kwargs):
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


def apply_weight(*args):
    if len(args) != 2:
        raise RuntimeError("not suppose to happen")
    image_file, weight_file = args
    hlin = fits.open(image_file)
    hlwht = fits.open(weight_file)
    for i in range(len(hlin)):
        if isinstance(hlin[i], (fits.PrimaryHDU, fits.ImageHDU)):
            badmask = hlwht[i].data == 0
            hlin[i].data[badmask] = np.nan
    hlin.writeto(image_file, overwrite=True)


def sex_to_ascii(*args):
    if len(args) > 1:
        raise RuntimeError("not suppose to happen")
    filename = args[0]
    tbl = Table.read(filename, format='ascii.sextractor')
    tbl.write(filename, format='ascii.commented_header', overwrite=True)


def header_to_config(*args, **kwargs):
    if len(args) > 1:
        raise RuntimeError("not suppose to happen")
    filename = args[0]
    hdr = fits.Header.fromfile(filename)
    config = {
            'CENTER_TYPE': 'MANUAL',
            'PIXELSCALE_TYPE': 'MANUAL',
            }
    config['CENTER'] = '{}, {}'.format(hdr['CRVAL1'], hdr['CRVAL2'])
    config['PIXEL_SCALE'] = np.sqrt(
            hdr['CD1_1'] ** 2. + hdr['CD1_2'] ** 2.) * 3600.
    config['IMAGE_SIZE'] = '{}, {}'.format(hdr['NAXIS1'], hdr['NAXIS2'])
    config.update(kwargs)
    with open(filename, 'w') as fo:
        for k, v in config.items():
            fo.write("{:20s}    {}\n".format(k, v))


def match_refcat(*args, **kwargs):
    log = get_log_func(default_level='debug', **kwargs)
    if not args:
        args = sys.argv[1:]
        cat_file, image_file, ref_file, out_file = args
    else:
        (cat_file, image_file, ref_file), out_file = args
    log("match {} to {}".format(cat_file, ref_file))

    # figure out band
    parsed_filename = re.match(
            kwargs['reg_mosaic'], os.path.basename(cat_file)
            ).groupdict()
    # print(parsed_filename)
    band = parsed_filename['band']
    log("use sdss band {}".format(band))

    icmd1 = 'select "{band}_sdss < 90. && {band}_sdss > 0' \
            ' && err_{band}_sdss <= 0.10857"' \
        .format(band=band)
    icmd2 = 'select "MAG_AUTO < 90. && MAGERR_AUTO >= 0.0005 && FLAGS == 0"'
    ocmd = ('addcol sep_gaia "skyDistanceDegrees('
            'ALPHA_J2000, DELTA_J2000, ra_gaia, dec_gaia) * 3600"')
    stilts_match = """{stilts_cmd}
                   tmatch2
                   matcher=sky
                   in1={ref_file}
                   ifmt1=ascii
                   icmd1={icmd1}
                   in2={cat_file}
                   ifmt2=ascii
                   icmd2={icmd2}
                   values1=ra_sdss dec_sdss
                   values2=ALPHA_J2000 DELTA_J2000
                   params=1.2
                   join=1and2
                   out={out_file}
                   ocmd={ocmd}
                   ofmt=ascii""".format(
                        stilts_cmd=kwargs['stilts_cmd'], **locals())
    subprocess.check_call(map(str.strip, stilts_match.split("\n")))

    tbl = Table.read(out_file, format='ascii.commented_header')
    log("save to matched cat {}".format(out_file))
    tbl.write(out_file, format='ascii.commented_header')
    # create plot for catalog and matched catalog
    cat_tbl = Table.read(cat_file, format='ascii.commented_header')
    savename = cat_file.rsplit(".cat", 1)[0] + ".png"
    qa_catalog(cat_tbl, tbl, image_file, savename, kwargs)


def qa_catalog(cat_tbl, ref_tbl, image_file, savename, kwargs):
    log = get_log_func(default_level='debug', **kwargs)
    log('create QA plots')
    set_color()
    fig = plt.figure(figsize=(8, 16))

    # depth histo, astro histo, zp offset
    ax = fig.add_subplot(3, 1, 1)
    plot_depth_histo(ax, cat_tbl)
    ax = fig.add_subplot(3, 1, 2)
    plot_astro_histo(ax, ref_tbl)
    ax = fig.add_subplot(3, 1, 3)
    plot_zp_color(ax, ref_tbl, image_file)

    log("save plot to {}".format(savename))
    # plt.show()
    plt.savefig(savename, bbox_inches='tight')


def plot_zp_color(ax, ref_tbl, image_file):
    band = re.match(r'.+odi_([ugriz])\.fits$',
                    os.path.basename(image_file)).groups()[0]
    smag = '{}_sdss'.format(band)
    semag = 'err_{}_sdss'.format(band)
    smaglims = dict(z=(17, 19), u=(17, 19)).get(band, (17, 20))
    mag = 'MAG_AUTO'
    emag = 'MAGERR_AUTO'
    sigma = 3.0

    # figure out color bands
    hdr = fits.open(image_file, memmap=True)[0].header
    k = hdr['MYCOLOR']
    c1, c2 = re.match(r"color term ([ugriz])-([ugriz])",
                      hdr.comments['MYCOLOR']).groups()
    c1, c2 = '{}_sdss'.format(c1), '{}_sdss'.format(c2)
    # good mask
    cat = ref_tbl[
            (ref_tbl[smag] > smaglims[0]) &
            (ref_tbl[smag] < smaglims[1])
            ]
    color = cat[c1] - cat[c2]
    zp = cat[smag] - cat[mag] - k * color
    zp0 = np.median(zp)
    good = (zp > zp0 - 0.5) & (zp < zp0 + 0.5)
    cat = cat[good]
    zp = zp[good]
    color = color[good]
    zperr = np.hypot(cat[semag], cat[emag])
    clipped, lo, up = sigmaclip(zp, sigma, sigma)
    im = (zp < up) & (zp > lo)
    ex = (zp > up) | (zp < lo)

    xlo = np.percentile(color, 5)
    xup = np.percentile(color, 95)
    ax.set_xlim((xlo - 0.2 * (xup - xlo), xup + 0.2 * (xup - xlo)))
    ax.set_ylim((-0.35, 0.55))
    pltkw = dict(fmt='o', ms=3, capsize=0, mew=0, color='red')
    pltkw_ex = dict(fmt='D', ms=2, capsize=0, mew=1, fillstyle='none')
    legs = []
    eecolor = 'red'
    if np.any(im):
        leg = ax.errorbar(color[im], zp[im] - zp0, yerr=zperr[im], **pltkw)
        for lc in leg[2]:
            ecolor = lc.get_color()[0]
            eecolor = change_hsv(ecolor, s=0.2, v=0.9)
            lc.set_color(eecolor)
            legs.append((leg, 'clipped'))
    if np.any(ex):
        legex = ax.errorbar(color[ex], zp[ex] - zp0, yerr=zperr[ex],
                            color=eecolor, **pltkw_ex)
        for lc in legex[2]:
            lc.set_color(eecolor)
    ax.axhline(0.05 * np.log10(np.e) * 2.5, ls='--', color='gray')
    ax.axhline(-0.05 * np.log10(np.e) * 2.5, ls='--', color='gray')
    ax.text(xlo, 0.15, r'$\pm 5\%$',
            verticalalignment='center', fontsize=20)
    label = [
            r'$n_{obj}=%d$' % (len(zp)),
            r'${\Delta}res.=%.5f$' % (np.std(zp)),
            r'$n_{obj,clip}=%d$' % (len(zp[im])),
            r'${\Delta}res._{clip}=%.5f$' % (np.std(zp[im]))
            ]
    ax.text(0.05, 0.92, '\n'.join(label),
            transform=ax.transAxes,
            verticalalignment='top')


def plot_astro_histo(ax, ref_tbl):
    cliplim = 0.2
    lims = (0, cliplim * 2)  # arcsec

    ax.set_xlim(lims)
    ax.set_xlabel(r'Separation (arcsec)')
    ax.set_ylabel(r'Count')
    bins = np.arange(lims[0], lims[1], 0.005)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    legs = []
    lbls = []
    for key, name in [('sep_gaia', 'Gaia'), ('Separation', 'SDSS')]:
        sep = ref_tbl[key]
        # sep = sep[~sep.mask]
        sep = sep[sep >= 0]
        sep_clipped = sep[sep < cliplim]
        sep_mean = np.median(sep_clipped)
        sep_std = np.std(sep_clipped)
        hist, _ = np.histogram(sep, bins=bins)
        leg = ax.errorbar(bin_centers, hist, drawstyle='steps-mid')
        lbl = "To {}: ${:.2f} \\pm {:.2f}\\ arcsec$".format(
                name, sep_mean, sep_std)
        legs.append(leg)
        lbls.append(lbl)
    ax.legend(legs, lbls, loc='upper right')


def plot_depth_histo(ax, cat_tbl):
    lims = (16, 28)
    ax.set_xlim(lims)
    ax.set_xlabel(r'MAG_AUTO')
    ax.set_ylabel(r'Count')

    bins = np.arange(lims[0], lims[1], 0.1)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    cat = cat_tbl[cat_tbl['MAG_AUTO'] > -90]
    snr = 2.5 * np.log10(np.e) / cat['MAGERR_AUTO']
    legs = []
    lbls = []
    for cut in [2, 3, 5, 10, 20]:
        subcat = cat[snr >= cut]
        mag = subcat['MAG_AUTO']
        hist, _ = np.histogram(mag, bins=bins)
        leg = ax.errorbar(bin_centers, hist, drawstyle='steps-mid')
        legs.append(leg)
        lbls.append("$SNR \geq {}\ ({:.2f}\%)$".format(
            cut,
            100 * len(subcat) / len(cat)))
    ax.legend(legs, lbls, loc='upper left')


def set_color():
    kelly = np.array([
        [255, 179, 0], [128, 62, 117], [255, 104, 0], [166, 189, 215],
        [193, 0, 32], [206, 162, 98], [129, 112, 102], [0, 125, 52],
        [246, 118, 142], [0, 83, 138], [255, 122, 92], [83, 55, 122],
        [255, 142, 0], [179, 40, 81], [244, 200, 0], [127, 24, 13],
        [147, 170, 0], [89, 51, 21], [241, 58, 19], [35, 44, 22]]) / 255.
    rc('axes', prop_cycle=cycler('color', kelly))


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
