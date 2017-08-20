#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2017-08-10 09:43
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
core.py

The core functionalities
"""

from __future__ import (absolute_import, division, print_function)
import re
import os
import stat
import logging
from datetime import datetime
import shutil
from shutilwhich import which
from astropy.io.misc import yaml
from astropy.io import fits
from astropy.table import Table, Column, unique
from multiprocessing import Pool, cpu_count
from functools import partial
import StringIO


from . import qa
from . import pipeline
from .apus import core as apuscore
from .apus.common import touch_file


def setup_workdir(workdir=".", overwrite_dir=False, backup_config=True):
    """
    Setup execution environment and dump default configuration file

    Parameters
    ----------
    workdir: str
        The directory to store the created dirs and files.

    require_empty: bool
        when set True, raise RuntimeError if the workdir is not empty

    backup_config: bool
        when not require_empty, this controls whether the files are
        get backed up when an overwriting is about to occur.
    """
    logger = logging.getLogger(name="setup")

    logdir = 'logs'  # stores logs for the apus pipelines
    externdir = 'extern'  # stores external dependencies
    tmpdir = 'tmp'  # stores temporary files
    configfile = 'postcalib.yaml'  # main configuration file

    # validate the workdir
    workdir_content = os.listdir(workdir)
    if workdir_content:
        if not overwrite_dir:
            # skip setup when all files are present
            if all(d in workdir_content for d in [
                    logdir, externdir, tmpdir, configfile]):
                raise RuntimeError(
                    "it seems that the workdir has been setup already. "
                    "Re-run with -f to proceed anyway")
            else:
                raise RuntimeError(
                    "the workdir is not empty, re-run with -f to proceed "
                    "anyway")
        else:
            logger.warning(
                    "the workdir is not empty but a forced setup is requested")
    logger.info("setup workdir {}".format(workdir))

    # external dependencies
    externdir = os.path.join(workdir, externdir)
    which_path = os.path.abspath(externdir) + ':' + os.environ['PATH']
    if os.path.isdir(externdir):
        logger.warning("use existing extern dir {}".format(externdir))
    else:
        os.makedirs(externdir)
        logger.info("create extern dir {}".format(externdir))
    logger.info("check external dependencies")
    astromatic_prefix = []
    for name, cmds, datadir in [
            ("SExtractor", ("sex", 'ldactoasc'), "sextractor"),
            ("SCAMP", ("scamp", ), "scamp"),
            ("SWarp", ("swarp", ), "swarp")]:
        cmds = [which(cmd, path=which_path) for cmd in cmds]
        if any(c is None for c in cmds):
            raise RuntimeError("not able to locate {}"
                               .format(name))
        prefix = os.path.normpath(
                os.path.join(os.path.dirname(cmds[0]), '../'))
        datadir = os.path.join(prefix, "share", datadir)
        if not os.path.exists(datadir):
            raise RuntimeError(
                "not able to locate data files for {0}. It is likely that {0} "
                "is compiled within the source directory but without proper "
                "installation. To resolve the issue, either run `make install`"
                " in the source directory, or manually link the data "
                "directory to {1}.".format(name, datadir))
        logger.info("{0:10s} ... OK".format(name))
        astromatic_prefix.append(prefix)
    if len(set(astromatic_prefix)) > 1:
        raise RuntimeError(
            "it seems that the SExtractor, SCAMP and SWarp are not installed "
            "into the same prefix. PostCalib for now does not deal with this "
            "situation. Try re-configure SExtractor, SCAMP and SWarp with "
            "--prefix=<prefixpath>")
    astromatic_prefix = os.path.normpath(astromatic_prefix[0])
    logger.info("use shared astromatic prefix {}".format(
        astromatic_prefix))
    stilts_cmd = which("stilts", path=which_path)
    if stilts_cmd is None:
        logger.warning("not able to find stilts. Get from internet")
        # retrieve stilts
        from astropy.utils.data import download_file
        stilts_jar_tmp = download_file(
                "http://www.star.bris.ac.uk/%7Embt/stilts/stilts.jar",
                cache=True)
        stilts_jar = os.path.join(externdir, 'stilts.jar')
        shutil.copyfile(stilts_jar_tmp, stilts_jar)
        stilts_cmd = os.path.join(externdir, 'stilts')
        with open(stilts_cmd, 'w') as fo:
            fo.write("""#!/bin/sh
java -Xmx4000M -classpath "{0}:$CLASSPATH" uk.ac.starlink.ttools.Stilts "$@"
""".format(os.path.abspath(stilts_jar)))
        os.chmod(stilts_cmd, os.stat(stilts_cmd).st_mode | stat.S_IEXEC)
    logger.info("{0:10s} ... OK".format("stilts"))
    logger.info("use stilts {}".format(stilts_cmd))
    funpack_cmd = which("funpack", path=which_path)
    if funpack_cmd is None:
        logger.warning("not able to find funpack. Get from internet")
        # retrieve stilts
        from astropy.utils.data import download_file
        funpack_tmp = download_file(
                "http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/"
                "cfitsio_latest.tar.gz",
                cache=True)
        funpack_src = os.path.join(externdir, 'cfitsio')
        import tarfile
        with tarfile.open(funpack_tmp) as tar:
            tar.extractall(path=externdir)
        logger.warning("try compiling funpack")
        import subprocess
        try:
            for command in [
                    './configure', 'make', 'make fpack', 'make funpack']:
                subprocess.check_call(command.split(), cwd=funpack_src)
        except subprocess.CalledProcessError:
            raise RuntimeError("unable to compile funpack")
        funpack_cmd = os.path.join(externdir, 'funpack')
        shutil.copy(
                os.path.join(externdir, 'cfitsio', 'funpack'),
                funpack_cmd)
    logger.info("{0:10s} ... OK".format("funpack"))
    logger.info("use funpack {}".format(funpack_cmd))

    # create logging directory
    logdir = os.path.join(workdir, logdir)
    if os.path.isdir(logdir):
        logger.warning("use existing log dir {}".format(logdir))
    else:
        os.makedirs(logdir)
        logger.info("create log dir {}".format(logdir))

    # setup scratch space
    tmpdir = os.path.join(workdir, tmpdir)
    # def freespace_GiB(path):
    #     stat = os.statvfs(path)
    #     return stat.f_bfree * stat.f_frsize / 1024 ** 3

    # logger.info("{:.2f} GiB free in {}".format(
    #     freespace_GiB(scratch_dir), scratch_dir))
    if os.path.isdir(tmpdir):
        logger.warning("use existing tmp dir {}".format(tmpdir))
    else:
        os.makedirs(tmpdir)
        logger.info("create tmp dir {}".format(tmpdir))

    # dump default config
    time_fmt = "%b-%d-%Y_%H-%M-%S"
    base_fmt = ("{obsid}_{object}_odi_{band}")
    config = """## config file of PostCalib
## {time:s}

# setting
phot_model_flags: 'color,ota,cat'

# qa inputs
qa_headers: [
    'OBSID', 'CRVAL1', 'CRVAL2', 'OBJECT', 'EXPMEAS', 'AIRMASS',
    'SEEING', 'SKY_MEDI', 'SKY_STD',
    'FILTER', 'INSTRUME', 'MJD-OBS'
    ]

# naming
reg_ppa: '{reg_ppa:s}'  # regex to parse the default filenames from PPA
fmt_orig: 'orig_{base_fmt}.{{ext}}'
reg_orig: '{reg_orig:s}'  # regex to parse the filenames in job-in dir

reg_inputs: '{reg_inputs:s}'  # regex to parse the filenames in jobdir
fmt_inputs: 'orig_{base_fmt}.{{ext}}'  # format string for image in jobdir
sel_inputs: ['orig_*.fits', 'orig_*.fits.fz']

fmt_masked: 'masked_{base_fmt:s}.fits'  # format string of masked images
sel_masked: 'masked_*.fits'   # selection string of masked images

fmt_selected: '{{ppflag}}{{imflag}}_{base_fmt:s}.fits'
sel_fcomb: 'fcomb[0-9]*_masked_*.fits'
fmt_objcat: '{{ppflag}}objcat_{base_fmt:s}.cat'  # format of object catalog
fmt_objmask: '{{ppflag}}objmask_{base_fmt:s}.fits'  # format of object masks
sel_objmask: 'fcomb[0-9]*_objmask_*.fits'   # selection string of object masks
fmt_sky: '{{ppflag}}sky_{base_fmt:s}.fits'  # format string of sky images
fmt_fcomb: '{{ppflag}}combined.fits'
reg_fcomb: '{reg_fcomb:s}'  # regex to parse grouped master
fmt_fsmooth: '{{ppflag}}smoothed.fits'

sel_fsub: 'fsub[0-9]*_masked_*.fits'
fmt_fsub_fsmooth: 'fcomb{{grpid}}_smoothed.fits'  # fmt for subtract fcomb
fmt_fsub: 'fsub_{base_fmt:s}.fits'

sel_fsubed: 'fsub_*.fits'
sel_phot: 'phot[0-9]*_*_*.fits'
fmt_photcat: '{{ppflag}}{{imflag}}_{base_fmt:s}.cat'
fmt_photcat_cleaned: '{{ppflag}}{{imflag}}_{base_fmt:s}.cln.cat'
fmt_photcat_matched: '{{ppflag}}{{imflag}}_{base_fmt:s}.zp.cat'
fmt_phot: '{{ppflag}}odi_{{band}}.cat'
reg_phot: '{reg_phot:s}'  # regex to parse grouped phot master
phot_hdr_suffix: 'hdr_phot'
phot_hdr_glob: 'phot[0-9]_*{{obsid}}*odi_{{band}}.hdr_phot'

sel_mosaic: 'mosaic[0-9]*_*_*.fits'
fmt_mosaic_orig: 'swarp{{grpid}}_{{imflag}}_odi_{{band}}.fits'
fmt_mosaic: 'coadd{{grpid}}_{{imflag}}_odi_{{band}}.fits'
fmt_mosaic_wht: 'coadd{{grpid}}_{{imflag}}_odi_{{band}}.wht.fits'
reg_mosaic: {reg_mosaic:s}

reg_grp: '{reg_grp:s}'  # regex to parse grouped images
fmt_grp: '{{ppflag}}{{grpid}}_{{imflag}}_{base_fmt:s}.fits'


# environ
workdir: {workdir}
tmpdir: {tmpdir}
logdir: {logdir}
astromatic_prefix: {astromatic_prefix}
stilts_cmd: {stilts_cmd}
funpack_cmd: {funpack_cmd}
""".format(time=datetime.now().strftime(time_fmt),
           version="0.0",
           reg_ppa=r'(?P<obsid>20\d{6}T\d{6}\.\d)_(?P<object>.+?)'
                   r'_odi_(?P<band>[ugriz])\.(?P<procid>\d{4})\.'
                   r'(?P<ext>fits|fits\.fz)$',
           reg_orig=r'(?P<ppflag>[^_/]+_)?(?P<imflag>[^_/]+)_'
                    r'(?P<obsid>20\d{6}T\d{6}\.\d)_(?P<object>.+?)'
                    r'_odi_(?P<band>[ugriz])'
                    r'\.(?P<ext>fits|fits\.fz)$',
           reg_inputs=r'(?P<ppflag>[^_/]+_)?(?P<imflag>[^_/]+)_'
                    r'(?P<obsid>20\d{6}T\d{6}\.\d)_(?P<object>.+?)'
                    r'_odi_(?P<band>[ugriz])'
                    r'\.(?P<ext>[^/]+)$',
           reg_fcomb=r'(?P<ppflag>[^_/]+_)(?P<imflag>[^/]+)\.fits',
           reg_phot=r'(?P<ppflag>[^_/]+)_odi_(?P<band>[ugriz])'
                    r'\.(?P<ext>[^/]+)',
           reg_mosaic=r'(?P<ppflag>[a-z]+)(?P<grpid>\d+)_(?P<imflag>[^_/]+)'
                      r'_odi_(?P<band>[ugriz]).fits',
           reg_grp=r'(?P<ppflag>[a-z]+)(?P<grpid>\d+)_(?P<imflag>[^_/]+)_'
                    r'(?P<obsid>20\d{6}T\d{6}\.\d)_(?P<object>.+?)'
                    r'_odi_(?P<band>[ugriz])'
                    r'\.(?P<ext>[^/]+)$',
           # reg_inputs=r'(?P<ppflag>[^_/]+_)?(?P<imflag>[^_/]+)_'
           #            r'(?P<obsid>20\d{6}T\d{6}\.\d)_(?P<object>.+?)'
           #            r'_odi_(?P<band>[ugriz])'
           #            r'_(?P<featgrp>\d+)_(?P<photgrp>\d+)_(?P<mscgrp>\d+)'
           #            r'\.(?P<ext>fits|fits\.fz)$',
           **locals())
    configfile = os.path.join(workdir, configfile)
    if os.path.exists(configfile):
        if backup_config:
            timestamp = datetime.fromtimestamp(
                os.path.getmtime(configfile)).strftime(time_fmt)
            bakfile = os.path.join(
                workdir,
                "{1}_{0}{2}".format(timestamp, *os.path.splitext(
                    os.path.basename(configfile)))
                )
            logger.warning("backup existing config file to {}".format(
                bakfile))
            os.rename(configfile, bakfile)
        else:
            logger.warning(
                    "overwrite existing config file {}".format(configfile))
    with open(configfile, 'w') as fo:
        fo.write(config)


def _qa_worker(image, config):
    """Create a QA summary for given input image"""
    headers = config['qa_headers']
    values = []
    with fits.open(image, memmap=True) as hdulist:
        for key in headers:
            try:
                val = hdulist[0].header[key]
            except KeyError:
                val = hdulist[1].header[key]
            except KeyError:
                val = ""
            if key == "OBJECT":
                # escape name
                val = re.sub(r'\s+', '_', val.strip())
            values.append(val)
        # QA image
        preview = qa.create_preview(hdulist=hdulist)
        preview.save()
        # mask guide ota
        values.insert(0, ','.join(map(str, preview.guide_otas)))
        # append filename
        values.append(image)
    keys = ['mask_otas', ] + headers + ['filename', ]
    return keys, values


def init_job(config_file, jobkey, images,
             overwrite_dir=False, backup_jobfile=True, dirs_only=False):
    logger = logging.getLogger(jobkey + ".init")

    # load config file
    with open(config_file, 'r') as fo:
        logger.info("use config file {}".format(config_file))
        logger.info("\n\n{}\n".format(fo.read()))
        fo.seek(0)
        config = yaml.load(fo)

    if dirs_only:
        logger.warning("init job {} dirs only".format(jobkey))

    # create or get input dir
    workdir = config['workdir']
    jindir = os.path.join(workdir, jobkey + ".in")
    logger.info("initialize inputs of job {}".format(jobkey))
    if os.path.exists(jindir):
        if overwrite_dir:
            logger.warning("overwrite existing job inputs dir {}".format(
                jindir))
        else:
            raise RuntimeError(
                "job inputs dir {} exists. Re-run with -f to overwrite".format(
                    jindir))
    else:
        os.makedirs(jindir)
        logger.info("create job inputs dir {}".format(jindir))

    if not dirs_only:
        # purge any existing good links that matches the linked name regex
        reg_orig = re.compile(config['reg_orig'])
        for filename in os.listdir(jindir):
            filename = os.path.join(jindir, filename)
            if os.path.islink(filename) and re.match(
                    reg_orig, os.path.basename(filename)):
                logger.warning("x {}".format(filename))
                os.unlink(filename)
                # remove preview files if any
                previewname = filename.rstrip(".fz").rstrip(".fits") + '.png'
                if os.path.exists(previewname) and not os.path.islink(
                        previewname):
                    logger.warning("x {}".format(previewname))
                    os.remove(previewname)

        # link over to job-in dir
        reg_ppa = re.compile(config['reg_ppa'])
        linked = []
        for i, image in enumerate(images):
            imagebase = os.path.basename(image)
            parsed_filename = re.match(reg_ppa, imagebase)
            if parsed_filename is None:
                logger.info("skip unmatched file {}".format(image))
                continue
            parsed_filename = parsed_filename.groupdict()
            linkname = config['fmt_orig'].format(**parsed_filename)
            link = os.path.join(jindir, linkname)
            logger.info("{} -> {}".format(imagebase, link))
            os.symlink(os.path.abspath(image), link)
            linked.append(link)  # update to the link path
        images = linked
        logger.info("{} input images linked".format(len(images)))

        # run QA
        pool = Pool(cpu_count())
        qa_keys, qa_vals = zip(*pool.map_async(
                partial(_qa_worker,  config=config),
                images).get(9999999))
        # tabulate QA results
        qa_tbl = Table(
                rows=qa_vals, names=qa_keys[0],
                )
        # prepend grouping table
        qa_tbl.add_column(Column(range(len(qa_tbl)), name='numid'), index=0)
        for i, (group, group_key) in enumerate([
                ('fcomb', "FILTER"), ('fsub', 'FILTER'),
                ('phot', 'FILTER'),
                ('mosaic', 'OBJECT')]):
            colname = group + "_group"
            keys = qa_tbl[group_key]
            unique_keys = list(unique(qa_tbl[[group_key]])[group_key])
            col = Column([unique_keys.index(k) for k in keys])
            qa_tbl.add_column(col, index=i + 1, name=colname)  # after numid

        logger.info("summary of input images\n{}".format(
            "\n".join(qa_tbl.pformat(max_lines=-1, max_width=-1))))
        # save the table
        tblname = os.path.join(jindir, '{}.txt.in'.format(jobkey))

        def quote(x):
            x = str(x)
            if ' ' in x or x == "":
                return "\"{}\"".format(x)
            else:
                return x
        _fo = StringIO.StringIO()
        qa_tbl.write(_fo, format='ascii.fixed_width', delimiter=" ",
                     formats={c: quote for c in qa_tbl.colnames})
        with open(tblname, 'w') as fo:
            fo.write('#' + _fo.getvalue()[1:])
        # copy table to work_dir for the next step
        jobfile = os.path.join(workdir, '{}.txt'.format(jobkey))

        if os.path.exists(jobfile):
            if backup_jobfile:
                time_fmt = "%b-%d-%Y_%H-%M-%S"
                timestamp = datetime.fromtimestamp(
                    os.path.getmtime(jobfile)).strftime(time_fmt)
                bakfile = os.path.join(
                    workdir,
                    "{1}_{0}{2}".format(timestamp, *os.path.splitext(
                        os.path.basename(jobfile)))
                    )
                logger.warning("backup existing job file to {}".format(
                    bakfile))
                os.rename(jobfile, bakfile)
            else:
                logger.warning(
                        "overwrite existing job file {}".format(jobfile))

        shutil.copy(tblname, jobfile)

    # create skymask dir
    skymask_dir = os.path.join(workdir, jobkey + ".skymask")
    if os.path.exists(skymask_dir):
        logger.warning("overwrite existing skymask dir {}".format(skymask_dir))
    else:
        os.makedirs(skymask_dir)
        logger.info("create sky mask dir {}".format(skymask_dir))
    with open(os.path.join(skymask_dir, 'README'), 'w') as fo:
        fo.write("""# Auto-generated README file, do not edit.

User can create and add DS9 region files to mask out
very bright stars that are hard to get rid of automatically.

The naming of the created DS9 region files could be
    1.  "skymask.reg"
        This file serves as the "global" mask, i.e., entries of this file
        are used to mask all images.

    2. *{obsid}*.reg
        This file is used when creating sky image from the image with "OBSID"
        keyword = obsid. The OBSID is a UTC time-like string that identifies
        uniquely the image, e.g. "20140925T201215.6".
        Example of valid region filenames:
            20140925T201215.6.reg
            20140925T201215.6_bright.reg
            gsc_20140925T201215.6.reg

The content of the region file should typically contain entries that are
specified in sky coordinates (RA, Dec). The file is parsed and processed with
python package "pyregion", which have good support of a variety of region
shapes. Please check https://pyregion.readthedocs.io/en/latest/ for more
details.""")
    touch_file(os.path.join(skymask_dir, 'skymask.reg'))

    # create job config
    # jobconfig =


def run_pipeline(config_file, jobfile, apus_args=None):
    jobkey = os.path.splitext(os.path.basename(jobfile))[0]
    logger = logging.getLogger(jobkey)
    logger.info("run pipeline {}".format(jobkey))

    # load config file
    with open(config_file, 'r') as fo:
        logger.info("use config file {}".format(config_file))
        # logger.info("\n\n{}\n".format(fo.read()))
        fo.seek(0)
        config = yaml.load(fo)

    # create jobdir
    workdir = config['workdir']
    jobdir = os.path.join(workdir, jobkey)
    if os.path.exists(jobdir):
        logger.warning("use existing jobdir {}".format(jobdir))
    else:
        os.makedirs(jobdir)
        logger.info("create job inputs dir {}".format(jobdir))

    config['jobdir'] = jobdir
    config['jobfile'] = jobfile
    config['jobkey'] = jobkey
    config['skymask_dir'] = jobdir + ".skymask"

    logger.info("arguments passed to Apus {}".format(apus_args))
    logger.propagate = False
    # call apus
    apuscore.bootstrap(
            dict(
                jobkey=jobkey,
                jobdir=jobdir,
                logdir=config['logdir'],
                env_overrides={
                    'path_prefix': config['astromatic_prefix'],
                    'tmpdir': config['tmpdir']},
                tlist=pipeline.get_tlist(config),
                ),
            apus_args
            )
    # # create or get input dir
    # workdir = config['workdir']
    # jobdir = os.path.join(workdir, jobkey)
    # if os.path.exists(jobdir):
    #     logger.warning("use existing jobdir {}".format(jobdir))
    # else:
    #     os.makedirs(jobdir)
    #     logger.info("create job inputs dir {}".format(jobdir))
