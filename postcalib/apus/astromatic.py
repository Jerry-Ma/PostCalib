#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2017-08-16 21:59
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
astromatic.py
"""


from __future__ import (absolute_import, division, print_function,
                        )
import os
import re
import logging


def get_scamp_checkplot_name(dirname, prefix=None):
    """
    Get string specifying paths of a full set of SCAMP diagnostic plots,
    to be used within SCAMP configuration files.

    Parameters
    ----------
    dirname: str
        The directory to save the plots

    prefix: str or None, optional
        The prefix to be add to the filenames of the plots

    Returns
    -------
    names: str
        Paths of diagnostic plots joined with ',' to be used in SCAMP
        configuration files.
    """
    scamp_checkplot = ['fgroups', 'distort', 'astr_interror2d',
                       'astr_interror1d', 'astr_referror2d',
                       'astr_referror1d', 'astr_chi2', 'psphot_error']
    if prefix is None:
        prefix = ""
    else:
        prefix = prefix + "_"
    return ','.join([os.path.join(dirname, '{0}{1}'.format(prefix, i))
                     for i in scamp_checkplot])


def parse_astromatic_conf(*conf_files):
    """
    Parse astromatic-style configuration files into a dict

    Parameters
    ----------
    positional args: str
        The list of configuration files to be read and parsed.
        Duplicated entries are overwritten following the order supplied.

    Returns
    -------
    params: dict
        Dict that contains the configuration (key, value) pairs from the
        given configuration files.
    """
    params = {}
    for fname in conf_files:
        with open(fname, 'r') as fo:
            for ln in fo.readlines():
                ln = ln.strip()
                if ln == '' or ln.startswith('#'):
                    continue
                else:
                    key, rest = map(str.strip, ln.split(None, 1))
                    value = list(map(str.strip, rest.split('#')))[0]
                    params[key] = value
    return params


def dump_astromatic_conf(infile, outfile, overwrite=False, **kwargs):
    """
    Genearte a modified version the input configuration file according
    to input kwargs

    Parameters
    ----------
    infile: str
        The template configuration file acts as the base
    outfile: str
        Output filename of the modified configuration file
    keyword args: str, optional
        The (key, value) pairs specifying how the input configuration
        file is to be modified.

    Returns
    -------
    outfile: str
        The dumped new configuration file
    """

    logger = logging.getLogger(__name__)
    if os.path.isfile(outfile) and not overwrite:
        raise RuntimeError("file exist:{0}".format(outfile))
    with open(outfile, 'w') as fo:
        for oln in infile.readlines():
            ln = oln.strip()
            if len(ln) == 0 or ln.startswith("#"):
                fo.write(oln)
                continue
            else:
                try:
                    keyval = list(map(str.strip, ln.split('#', 1)))[0]
                except ValueError:
                    keyval = ln
                try:
                    key, val = map(str.strip, keyval.split(None, 1))
                except ValueError:
                    key, val = keyval, None
                newval = kwargs.get(key, None)
                if newval is None:
                    fo.write(oln)
                    continue
                else:
                    if val is None:
                        fo.write(oln.replace(key, '{0}  {1}'.format(
                            key, newval)))
                    else:
                        # only replace the value part
                        iv = oln.index(key) + len(key)
                        if '#' in oln:
                            jv = oln.index('#')
                        fo.write(oln[:iv] + oln[iv:jv].replace(val, newval) +
                                 oln[jv:])
    logger.info("+> {0:s}".format(outfile))
    return outfile


def dump_sex_param(infile, outfile, *args, **kwargs):
    """
    Genearte a modified version the input parameter file according
    to input args and kwargs

    Parameters
    ----------
    infile: str
        The template parameter file acts as the base
    outfile: str
        Output filename of the modified parameter file
    positional args: str, optional

    keyword args: str, optional
        The (key, value) pairs specifying how the input parameter
        file is to be modified.

    Returns
    -------
    outfile: str
        The dumped new parameter file
    """

    logger = logging.getLogger(__name__)
    err = 'Not a valid output para: {0:s}'
    if os.path.isfile(outfile) and not kwargs.get('overwrite', False):
        raise RuntimeError("file exist:{0}".format(outfile))
    content = infile.getvalue()
    # merge common with args: handle array-like keys
    re_key = re.compile('(\w+)\s*(\(\s*\d+\s*\))?')
    keys = []
    for arg in args:
        for key in arg:
            match = re.match(re_key, key.strip())
            if match is not None:
                _key = match.groups()[0]
                if _key not in content:
                    raise RuntimeError(err.format(key))
                if _key in keys:
                    keys.index(_key)
                    keys[keys.index(_key)] = key
                else:
                    keys.append(key)
            else:
                raise RuntimeError(err.format(key))
    logger.info('output params: {0}'.format(', '.join(keys)))
    with open(outfile, 'w') as fo:
        for key in keys:
            fo.write('{0:23s}  #\n'.format(key))
        fo.write('#' * 26 + '\n')
        fo.write(content)
    logger.info("+> {0:s}".format(outfile))
    return keys


class AmConfig(object):
    """
    Hold definitions and defaults for astromatic tools
    """

    sexparam_default = [
        # coord
        'ALPHA_J2000', 'DELTA_J2000', 'X_IMAGE', 'Y_IMAGE',
        'NUMBER', 'EXT_NUMBER',
        # phot
        'MAG_AUTO', 'MAGERR_AUTO', 'MAG_APER', 'MAGERR_APER',
        'FLUX_AUTO', 'FLUXERR_AUTO', 'FLUX_APER', 'FLUXERR_APER',
        'BACKGROUND', 'THRESHOLD',
        # scamp
        'XWIN_IMAGE', 'YWIN_IMAGE',
        'ERRAWIN_IMAGE', 'ERRBWIN_IMAGE', 'ERRTHETAWIN_IMAGE',
        'FLAGS', 'FLAGS_WEIGHT', 'FLAGS_WIN',  # 'IMAFLAGS_ISO',
        'FLUX_RADIUS',
        # ref key
        'X_WORLD', 'Y_WORLD',
        'ERRA_WORLD', 'ERRB_WORLD', 'ERRTHETA_WORLD',
        # PSF shape
        'FWHM_WORLD', 'FWHM_IMAGE',
        'A_IMAGE', 'B_IMAGE',
        'THETA_IMAGE', 'ELLIPTICITY',
        'CLASS_STAR'
        ]
    sex_default = {
            'FILTER_NAME': 'default.conv',
            'STARNNW_NAME': 'default.nnw',
            'WRITE_XML': 'N',
            'BACKPHOTO_TYPE': 'LOCAL',
            'PIXEL_SCALE': 0,
            'HEADER_SUFFIX': '.none',
            'GAIN_KEY': 'bug_of_sex_219',
            'NTHREADS': 0,
            }
    scamp_default = {
            'CHECKPLOT_RES': '1024',
            'SAVE_REFCATALOG': 'Y',
            'WRITE_XML': 'N',
            }
    swarp_default = {
            'INTERPOLATE': 'N',
            'FSCALASTRO_TYPE': 'VARIABLE',
            'DELETE_TMPFILES': 'Y',
            'NOPENFILES_MAX': '1000000',
            'WRITE_XML': 'N',
            }
    tmpdir = '/tmp'
    path_prefix = '/usr'

    def __init__(self, **kwargs):
        """
        Populate properties, with optional overrides from kwargs
        """
        self.set_overrides(kwargs)

    def get(self, prop):
        return self.overrides.get(prop, getattr(self, prop))

    def set_overrides(self, overrides):
        self.overrides = overrides
        self.share_dir = os.path.join(self.get('path_prefix'), 'share')
        self.bin_dir = os.path.join(self.get('path_prefix'), 'bin')
        for i in ['bin', 'share']:
            for sname, lname in [('sex', 'sextractor'),
                                 ('scamp', 'scamp'),
                                 ('swarp', 'swarp')]:
                if sname == 'sex' and i == 'bin':
                    lname = 'sex'  # sextractor binary naming
                setattr(self, '{0}{1}'.format(sname, i),
                        os.path.join(self.get('{0}_dir'.format(i)), lname))


def normalize_params(params):
    """
    Return a dict with all values being str type, suitable for
    writing into astromatic configuration files.
    """
    outparams = {}
    for k, v in params.items():
        if isinstance(v, (list, tuple)):
            v = ','.join(map(str, v))
        else:
            v = str(v)
        outparams[k] = v
    return outparams
