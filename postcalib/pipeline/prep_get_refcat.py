#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2017-08-18 01:57
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
prep_get_refcat.py

This script is called to retrieve reference catalog for a collection
of images.

The program looks into the footprints of all the input images and
combines them. The reference catalog server is queried and the
catalog within the combined footprint of the input images is
downloaded.

For now, the script queries SDSS, Pan-Starrs, and Gaia.
The results of the three are merged into a master catalog, refcat.cat.

Along with the catalogs, a DS9 region file is also generated, which
contains boxes showing the (combined) footprints of the input images.

Note: this script will produce an error if the queried sky area is larger
than 5x5 deg2.

Inputs
------
images: fits images
    The input images for which the reference catalog is to be retrieved

Outputs
-------
recat_{object}.cat: ASCII table
    The output master reference catalog retrieved.
sdss_{object}.cat: ASCII table
    The SDSS table
pnst_{object}.cat: ASCII table
    The Pan-Starr table
gaia_{object}.cat: ASCII table
    The Gaia table
"""


from __future__ import (absolute_import, division, print_function)
import os
import re
import sys
# import glob
import itertools
from astropy.io import fits
from astropy.coordinates import SkyCoord
import subprocess
import astropy.units as u
import numpy as np
from astroquery.sdss import SDSS
from astroquery.gaia import Gaia
from astropy.table import Table, vstack, unique


from postcalib.wiyn import get_layout
from postcalib.apus.common import get_log_func


def main(*args, **kwargs):
    if not args:
        args = sys.argv[1:]
        in_files = args[:-1]
        out_file = args[-1]
    else:
        in_files, out_file = args
    log = get_log_func(default_level='debug', **kwargs)

    outdir = os.path.dirname(out_file)
    outbase = os.path.basename(out_file)
    out_reg = os.path.splitext(out_file)[0] + '.reg'
    sdss_filter = re.match(r'.+_odi_([ugriz])\.+', in_files[0]).group(1)

    # determine the footprint of the images
    box, ds9reg = combined_sky_footprint(in_files)
    w, e, s, n = box
    cra, cdec, width, height = to_ds9_box(box)
    area = width * height
    log("{} footprint ra {}-{} dec {}-{}".format(out_file, w, e, s, n))
    log("estimated area {:.2f} deg2".format(area))
    if area > 25:
        log('error', "the combined footprint is too large!")
        raise RuntimeError("please check the image pointings and make sure "
                           "they are close-by to each other on the sky")
    # write ds9 region
    log("save to DS9 region {}".format(out_reg))
    with open(out_reg, 'w') as fo:
        fo.write(ds9reg)

    stilts_in = []
    cats_in = []
    # query sdss and get reference catalog
    sdss_cat = query_sdss(
            filter=sdss_filter,
            min_ra=w, max_ra=e,
            min_dec=s, max_dec=n,
            **kwargs
            )
    if sdss_cat is not None and len(sdss_cat) > 0:
        out_sdss = os.path.join(outdir, outbase.replace("refcat", 'sdss'))
        log("save to SDSS catalog {}".format(out_sdss))
        sdss_cat.write(out_sdss, format='ascii.commented_header')
        stilts_in.append(
            "ifmt%N%=ascii in%N%={out_sdss} values%N%=ra~dec suffix%N%=_sdss ")
        cats_in.append((sdss_cat, '_sdss'))

    # query gaia
    gaia_cat = query_gaia(
            min_ra=w, max_ra=e,
            min_dec=s, max_dec=n,
            **kwargs
            )
    if gaia_cat is not None and len(gaia_cat) > 0:
        out_gaia = os.path.join(outdir, outbase.replace("refcat", 'gaia'))
        log("save to Gaia catalog {}".format(out_gaia))
        gaia_cat.write(out_gaia, format='ascii.commented_header')
        stilts_in.append(
            "ifmt%N%=ascii in%N%={out_gaia} values%N%=ra~dec suffix%N%=_gaia ")
        cats_in.append((gaia_cat, '_gaia'))

    # query Pan-STARRS DR1
    out_ps1 = os.path.join(outdir, outbase.replace("refcat", 'ps1'))
    ps1_cat = query_panstarrs(
            min_ra=w, max_ra=e,
            min_dec=s, max_dec=n,
            out_file=out_ps1,
            cra=cra, cdec=cdec,
            width=width, height=height,
            **kwargs
            )
    if ps1_cat is not None and len(ps1_cat) > 0:
        log("save to PS1 catalog {}".format(out_ps1))
        ps1_cat.write(out_ps1, format='ascii.commented_header')
        stilts_in.append(
            "ifmt%N%=ascii in%N%={out_ps1} values%N%=ra~dec suffix%N%=_ps1 ")
        cats_in.append((ps1_cat, '_ps1'))

    nin = len(stilts_in)
    if nin == 0:
        raise RuntimeError("no reference catalog can be found")
    elif nin == 1:
        # add suffix to the catalog columns and save
        out_cat, suffix = cats_in[0]
        for col in out_cat.colnames:
            out_cat.rename_column(col, col + suffix)
        out_cat.write(out_file, format='ascii.commented_header',
                      overwrite=True)
    else:
        # merge the catalogs
        stilts_merge = (
                "{stilts_cmd} tmatchn nin={nin} "
                "out={out_file} ofmt=ascii multimode=group "
                "matcher=sky params=1 fixcols=all " +
                ' '.join(
                    in_.replace("%N%", str(i + 1))
                    for i, in_ in enumerate(stilts_in)) + ' '
                ' '.join(
                    'join{}=always'.format(i + 1) for i in range(nin))
                ).format(
                    stilts_cmd=kwargs['stilts_cmd'], **locals())

        subprocess.check_call(
            re.split(r'%+', stilts_merge.replace(" ", '%').replace('~', " ")))
    refcat = Table.read(out_file, format='ascii.commented_header')
    log("save to refcat {}".format(out_file))
    log("{} stars in total".format(len(refcat)))


def stack_cats(in_cats, out_cat, **kwargs):
    """
    Returns a vertically stacked catalog from multiple input catalogs.
    Any duplicated entry is removed.

    Inputs
    ------
    in_cats: list of ASCII tables
        The input catalogs to be stacked.

    Outputs
    -------
    out_cat: ASCII table
        The stacked catalog
    """
    tbls = [Table.read(f, format='ascii.commented_header') for f in in_cats]
    tbl = vstack(tbls, join_type='exact')
    # use ra_sdss as key
    # tbl_w_key = tbl[~tbl['ra_sdss'].mask]
    # tbl_wo_key = tbl[tbl['ra_sdss'].mask]
    # run unique separately
    tbls = []
    keys = ['ra_sdss', 'ra_gaia', 'ra_ps1']

    def get_tbl(tbl, ikey):
        key = keys[ikey]
        tbl_w_key = tbl[~tbl[key].mask]
        tbl_wo_key = tbl[tbl[key].mask]
        tbls.append(unique(tbl_w_key, keys=key))
        if ikey == len(keys) - 1:
            pass
            # tbls.append(unique(tbl_wo_key, keys=key))
        else:
            get_tbl(tbl_wo_key, ikey + 1)

    get_tbl(tbl, 0)
    tbl = vstack(tbls, join_type='exact')
    log = get_log_func(default_level='debug', **kwargs)
    log("save to master refcat {}".format(out_cat))
    log("{} stars in total".format(len(tbl)))
    tbl.write(out_cat, format='ascii.commented_header')


def query_sdss(**kwargs):
    log = get_log_func(default_level='debug', **kwargs)
    sql_query = [
        "SELECT ra,dec,raErr,decErr,u,err_u,g,err_g,r,err_r,i,err_i,z,err_z",
        "FROM Star WHERE",
        "    ra BETWEEN {min_ra:f} and {max_ra:f}",
        "AND dec BETWEEN {min_dec:f} and {max_dec:f}",
        "AND ((flags_{filter:s} & 0x10000000) != 0)",     # detected in BINNED1
        # not EDGE, NOPROFILE, PEAKCENTER, NOTCHECKED, PSF_FLUX_INTERP,
        # SATURATED, or BAD_COUNTS_ERROR"
        "AND ((flags_{filter:s} & 0x8100000c00a4) = 0)",
        # not DEBLEND_NOPEAK or small PSF error
        "AND (((flags_{filter:s} & 0x400000000000) = 0) or "
        "(psfmagerr_{filter:s} <= 0.2))",
        # not INTERP_CENTER or not COSMIC_RAY
        "AND (((flags_{filter:s} & 0x100000000000) = 0) or "
        "(flags_{filter:s} & 0x1000) = 0)"]
    sql_query = '\n'.join(sql_query).format(**kwargs)
    log("query SDSS with sql string\n{}".format(sql_query))
    cat = SDSS.query_sql(sql_query)
    log("{} stars found in SDSS".format(len(cat)))
    return cat


def query_gaia(**kwargs):
    log = get_log_func(default_level='debug', **kwargs)
    sql_query = [
        "SELECT ra,dec,ra_error,dec_error,phot_g_mean_mag,source_id",
        "FROM gaiadr1.gaia_source WHERE",
        "    ra BETWEEN {min_ra:f} and {max_ra:f}",
        "AND dec BETWEEN {min_dec:f} and {max_dec:f}", ]
    sql_query = '\n'.join(sql_query).format(**kwargs)
    log("query Gaia with sql string\n{}".format(sql_query))
    try:
        job = Gaia.launch_job_async(sql_query)
        cat = job.get_results()
        log("{} stars found in Gaia".format(len(cat)))
        return cat
    except Exception as e:
        log("unable to query gaia catalog")
        return None


def query_panstarrs(**kwargs):
    log = get_log_func(default_level='debug', **kwargs)
    url = 'http://archive.stsci.edu/panstarrs/search.php'
    querycols = ('objName,objiD,'
                 'raMean,decMean,raMeanErr,decMeanErr,'
                 'nDetections,objInfoFlag,qualityFlag,'
                 'ng,nr,ni,nz,ny,'
                 'gMeanPSFMag,gMeanPSFMagErr,gFlags,'
                 'rMeanPSFMag,rMeanPSFMagErr,rFlags,'
                 'iMeanPSFMag,iMeanPSFMagErr,iFlags,'
                 'zMeanPSFMag,zMeanPSFMagErr,zFlags,'
                 'yMeanPSFMag,yMeanPSFMagErr,yFlags')
    payload = {
            'outputformat': 'VOTable',
            'equinox': 'J2000',
            'nDetections': '>4',
            'selectedColumnsCsv': querycols.lower(),
            'max_records': 50001,
            'max_rpp': 500,
            'action': 'Search',
            'skipformat': 'on',
            }
    # check the target size to compute request params
    ps1_size = 30. * 2 ** 0.5  # arcmin
    target_width = kwargs['width'] * 60
    target_height = kwargs['height'] * 60

    # get the query radius and number pointings
    def get_pts(width, min_, max_):
        n_pts = 1
        while True:
            w = width / n_pts
            if w <= ps1_size:
                borders = np.linspace(min_, max_, n_pts + 1)
                return (borders[:-1] + borders[1:]) * 0.5, w / 2. ** 0.5
            else:
                n_pts += 1

    ra_centers, ra_rad = get_pts(
            target_width, kwargs['min_ra'], kwargs['max_ra'])
    dec_centers, dec_rad = get_pts(
            target_height, kwargs['min_dec'], kwargs['max_dec'])
    # print(ra_centers, ra_rad, dec_centers, dec_rad)
    if ra_rad > dec_rad:
        query_rad = ra_rad
        query_size = query_rad * 2 ** 0.5
        dec_n_pts = int(((0.999 * target_height) // query_size) + 1)
        dec_pts_off = (
                dec_n_pts * query_size - target_height) * 0.5 / 60.
        dec_borders = np.linspace(
                kwargs['min_dec'] - dec_pts_off,
                kwargs['min_dec'] - dec_pts_off + dec_n_pts * query_size / 60.,
                dec_n_pts + 1)
        dec_centers = (dec_borders[:-1] + dec_borders[1:]) * 0.5
    else:
        query_rad = dec_rad
        cosdec = np.cos(np.deg2rad(kwargs['cdec']))
        query_size = query_rad * 2 ** 0.5
        ra_n_pts = int(((0.999 * target_width) // query_size) + 1)
        ra_pts_off = (
                ra_n_pts * query_size - target_height) * 0.5 / 60.
        ra_borders = np.linspace(
                kwargs['min_ra'] - ra_pts_off / cosdec,
                kwargs['min_ra'] -
                (ra_pts_off + ra_n_pts * query_size / 60.) / cosdec,
                ra_n_pts + 1)
        ra_centers = (ra_borders[:-1] + ra_borders[1:]) * 0.5
        ra_n_pts = int(((0.999 * target_width) // query_rad) + 1)
        ra_pts_off = (
                ra_n_pts * query_rad * 2 ** 0.5 - target_width) * 0.5 / 60.
        ra_borders = np.linspace(
                kwargs['min_ra'] - ra_pts_off,
                kwargs['max_ra'] + ra_pts_off, ra_n_pts)
        ra_centers = (ra_borders[:-1] + ra_borders[1:]) * 0.5
    # print(ra_centers, dec_centers, query_rad)
    # import pdb
    # pdb.set_trace()
    # do the request
    import requests
    from astropy.table import Table

    cats = []
    for i, (cra, cdec) in enumerate(
            itertools.product(ra_centers, dec_centers)):
        _payload = dict(payload, ra=cra, dec=cdec, radius=query_rad)
        # from astropy.table import Column
        log("query PAN-STARRS DR1 chunk {} with payload\n{}".format(
            i, _payload))
        response = requests.get(url=url, params=_payload)
        if response.status_code == requests.codes.ok:
            out_file = kwargs['out_file'] + 'chunk{}'.format(i)
            with open(out_file, 'w') as fo:
                fo.write(response.text)
            # import pdb
            # pdb.set_trace()
            # filter the catalog
            cat = Table.read(out_file, format='votable')
            log("{} PS1 sources before filtering".format(len(cat)))
            # good source flags
            goodflag = 0x00000004 | 0x00000010 | 0x00000020
            # cat['qualityFlag'].fill_value = 0xFFFFFFFF
            # qualityflag = cat['qualityFlag'].astype(int)
            qualityflag = cat['qualityflag']
            cat = cat[((qualityflag & goodflag) == goodflag) &
                      (cat['ramean'] >= kwargs['min_ra']) &
                      (cat['ramean'] <= kwargs['max_ra']) &
                      (cat['decmean'] >= kwargs['min_dec']) &
                      (cat['decmean'] <= kwargs['max_dec'])
                      ]
            log("{} PS1 sources after filtering".format(len(cat)))
            cats.append(cat)
    if len(cats) == 0:
        log("unable to found PS1 sources")
        return None
    cat = unique(vstack(cats, join_type='exact'), keys='objid')
    outcols = [
            None, None, 'ra', 'dec', 'raErr', 'decErr',
            None, None, None,
            None, None, None, None, None,
            'g', 'err_g', 'flag_g',
            'r', 'err_r', 'flag_r',
            'i', 'err_i', 'flag_i',
            'z', 'err_z', 'flag_z',
            'y', 'err_y', 'flag_y',
            ]
    for outcol, col in zip(outcols, querycols.split(',')):
        if outcol is not None:
            cat.rename_column(col.lower(), outcol)
    return cat


def combined_sky_footprint(images):
    bbox = None
    ds9reg = ["global color=yellow", ]
    for image in images:
        with fits.open(image, memmap=True) as hdulist:
            coord = SkyCoord(hdulist[0].header['RA'],
                             hdulist[0].header['DEC'],
                             unit=(u.hourangle, u.degree))
            ra = coord.ra.degree
            dec = coord.dec.degree
            layout = get_layout(hdulist)
            box = layout.get_sky_footprint(center=(ra, dec))
            icra, icdec, iwidth, iheight = to_ds9_box(box)
            ds9reg.append('fk5; box({0},{1},{2},{3}, 0)\n'.format(
                    icra, icdec, iwidth, iheight))
            ds9reg.append(
                'fk5; point({0},{1}) # color=red text={{{2}}}\n'.format(
                        icra, icdec, hdulist[0].header['OBSID'][-4:]))
            bbox = merge_bbox(bbox, box)
    cra, cdec, width, height = to_ds9_box(bbox)
    ds9reg.append('fk5; box({0},{1},{2},{3}, 0) # color=red'.format(
            cra, cdec, width, height))
    ds9reg = "\n".join(ds9reg)
    return bbox, ds9reg


def merge_bbox(box1, box2):
    if box1 is None:
        return box2
    elif box2 is None:
        return box1
    else:
        l1, r1, b1, t1 = box1
        l2, r2, b2, t2 = box2
        ramin, ramax = get_min_max([l1, l2, r1, r2])
        decmin, decmax = get_min_max([b1, b2, t1, t2])
        return ramin, ramax, decmin, decmax


def get_min_max(l):
    minl = np.inf
    maxl = -np.inf
    dist = 0
    for i in l:
        for j in l:
            # TODO get handle of degree wrapping
            if np.abs(i - j) > dist:
                minl = min(i, j)
                maxl = max(i, j)
                dist = np.abs(i - j)
    if dist > 180.:
        raise RuntimeError(
            "there is degree wrapping that that code is not able to handle "
            "for input sky coord values {}".format(l))
    return minl, maxl


def to_ds9_box(box):
    w, e, s, n = box
    cra = (w + e) * 0.5
    cdec = (s + n) * 0.5
    dra = (e - w) * np.cos(cdec * np.pi / 180.)
    ddec = (n - s)
    return cra, cdec, dra, ddec


if __name__ == "__main__":
    main()
