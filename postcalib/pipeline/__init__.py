#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2017-08-17 15:35
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
__init__.py

PostCalib pipelines
"""

from __future__ import (absolute_import, division, print_function)
import os
import re
import sys
from . import prep_grouping
from . import prep_mask_otas
from . import prep_get_refcat
from . import sky_mask_objects
from . import sky_combine
from . import sky_subtract
from . import phot_calib
from . import phot_mosaic


def get_tlist(config):
    t00 = dict(
        name='select images',
        func=prep_mask_otas.select_images,
        pipe='transform',
        in_=os.path.abspath(config['jobfile']),
        extras=config['jobdir'],
        out=config['jobkey'] + '.checker',
        jobs_limit=1,
        kwargs={
            "reg_orig": config['reg_orig'],
            "fmt_inputs": config['fmt_inputs'],
            }
            )
    t01 = dict(
        name='mask OTAs',
        # func=pyfunc('prep_mask_otas.py', '{in} {out}'),
        func=prep_mask_otas.main,
        pipe='transform',
        in_=(config['sel_inputs'], config['reg_inputs']),
        extras=os.path.abspath(config['jobfile']),
        # add_inputs=[
        #     os.path.abspath(config['jobfile']),
        #     t00],
        out=fmtname(config['fmt_masked']),
        kwargs={
            "reg_inputs": config['reg_inputs'],
            "fmt_masked": config['fmt_masked'],
            "funpack_cmd": config['funpack_cmd'],
            },
        follows=t00,
        )
    # reference catalogs
    t10 = dict(
        name='get refcats',
        func=prep_get_refcat.main,
        pipe='collate',
        in_=(t01, config['reg_inputs']),
        out='refcat_{object[0]}.cat',
        kwargs={
            "stilts_cmd": config['stilts_cmd'],
            }
            )
    t15 = dict(
        name='merge refcats',
        func=prep_get_refcat.stack_cats,
        pipe='collate',
        in_=(t10, r'(?P<kind>refcat)_(?P<object>.+)\.cat'),
        out='{kind[0]}.cat',
            )
    # fringe/pupil removal
    t20 = dict(
        name='select fcomb',
        func=prep_grouping.main,
        pipe='transform',
        in_=os.path.abspath(config['jobfile']),
        extras=config['jobdir'],
        out=config['jobkey'] + '.fcomb_group',
        jobs_limit=1,
        kwargs={
            "grpkey": 'fcomb',
            "sel_inputs": config['sel_masked'],
            "reg_inputs": config['reg_inputs'],
            "fmt_selected": config['fmt_selected'],
            },
        follows=t01,
            )
    t21 = dict(
        name='get objmask',
        func='sex',
        pipe='transform',
        in_=(config['sel_fcomb'], config['reg_inputs']),
        out_keys=['CATALOG_NAME', 'CHECKIMAGE_NAME'],
        out=[fmtname(config['fmt_objcat']),
             fmtname(config['fmt_objmask'])],
        params={'CATALOG_TYPE': 'ASCII_HEAD',
                'DETECT_MINAREA': 3,
                'DETECT_THRESH': 3,
                'ANALYSIS_THRESH': 3,
                'DEBLEND_MINCONT': 0.05,
                'BACK_SIZE': 128,
                'CHECKIMAGE_TYPE': 'SEGMENTATION',
                },
        follows=t20,
            )
    t22 = dict(
        name='mask objects',
        func=sky_mask_objects.main,
        pipe='transform',
        in_=(config['sel_fcomb'], config['reg_inputs']),
        add_inputs=fmtname(config['fmt_objmask']),
        out=fmtname(config['fmt_sky']),
        follows=[t20, t21],
        kwargs={
            'skymask_dir': config['skymask_dir']
            }
            )
    t23 = dict(
        name='create ftemp',
        func=sky_combine.main,
        pipe='collate',
        in_=(t22, config['reg_inputs']),
        out=fmtname(config['fmt_fcomb']),
        jobs_limit=1,
            )
    t24 = dict(
        name='smooth ftemp',
        func=sky_combine.smooth,
        pipe='transform',
        in_=(t23, config['reg_fcomb']),
        out=fmtname(config['fmt_fsmooth']),
        jobs_limit=1,
            )
    # fringe/pupil sub
    t30 = dict(
        name='select fsub',
        func=prep_grouping.main,
        pipe='transform',
        in_=os.path.abspath(config['jobfile']),
        extras=config['jobdir'],
        out=config['jobkey'] + '.fsub_group',
        jobs_limit=1,
        kwargs={
            "grpkey": 'fsub',
            "sel_inputs": config['sel_masked'],
            "reg_inputs": config['reg_inputs'],
            "fmt_selected": config['fmt_selected'],
            },
        follows=t01,
            )
    t31 = dict(
        name='subtract ftemp',
        func=sky_subtract.main,
        pipe='transform',
        in_=(config['sel_fsub'], config['reg_grp']),
        add_inputs=fmtname(config['fmt_fsub_fsmooth']),
        out=fmtname(config['fmt_fsub']),
        follows=[t30, t24],
            )
    # phot calib
    t40 = dict(
        name='select phot',
        func=prep_grouping.main,
        pipe='transform',
        in_=os.path.abspath(config['jobfile']),
        extras=config['jobdir'],
        out=config['jobkey'] + '.phot_group',
        jobs_limit=1,
        kwargs={
            "grpkey": 'phot',
            "sel_inputs": config['sel_fsubed'],
            "reg_inputs": config['reg_inputs'],
            "fmt_selected": config['fmt_selected'],
            "fallbacks": [('masked', config['sel_masked']), ]
            },
        follows=[t01, t31],
            )
    t41 = dict(
        name='get photcat',
        func='sex',
        pipe='transform',
        in_=(config['sel_phot'], config['reg_inputs']),
        out=fmtname(config['fmt_photcat']),
        params={'CATALOG_TYPE': 'ASCII_HEAD',
                'DETECT_MINAREA': 3,
                'DETECT_THRESH': 3,
                'ANALYSIS_THRESH': 3,
                'DEBLEND_MINCONT': 0.05,
                'PHOT_APERTURES': [18, 27, 36, 45, 54, 72, 90, 109],
                'BACK_SIZE': 128,
                # 'VERBOSE_TYPE': 'QUIET',
                },
        outparams=['MAG_APER(8)', 'MAGERR_APER(8)', 'FLUX_MAX',
                   'AWIN_IMAGE', 'BWIN_IMAGE', 'ELONGATION'],
        follows=t40,
            )
    t42 = dict(
        name='cleanup photcat',
        func=phot_calib.cleanup,
        pipe='transform',
        in_=(t41, config['reg_inputs']),
        add_inputs="{basename[0]}.fits",
        out=fmtname(config['fmt_photcat_cleaned']),
        kwargs={
            'skymask_dir': config['skymask_dir']
            }
            )
    t43 = dict(
        name='match refcat',
        func=phot_calib.match_refcat,
        pipe='transform',
        in_=(t41, config['reg_inputs']),
        replace_inputs=[t42['out'], t10['out'], '{basename[0]}.fits'],
        out=fmtname(config['fmt_photcat_matched']),
        follows=t10,
        kwargs={
            'reg_inputs': config['reg_inputs'],
            'stilts_cmd': config['stilts_cmd'],
            }
            )
    t44 = dict(
        name='get flxscale',
        func=phot_calib.main,
        pipe='collate',
        in_=(t41, config['reg_inputs']),
        replace_inputs=[t43['out'], '{basename[0]}.fits'],
        out=fmtname(config['fmt_phot']),
        kwargs={
            'reg_phot': config['reg_phot'],
            'phot_model_flags': config['phot_model_flags'],
            'phot_hdr_suffix': config['phot_hdr_suffix'],
            },
        follows=[t43, t40]
            )
    t50 = dict(
        name='select mosaic',
        func=prep_grouping.main,
        pipe='transform',
        in_=os.path.abspath(config['jobfile']),
        extras=config['jobdir'],
        out=config['jobkey'] + '.mosaic_group',
        jobs_limit=1,
        kwargs={
            "grpkey": 'mosaic',
            "sel_inputs": config['sel_fsubed'],
            "reg_inputs": config['reg_inputs'],
            "fmt_selected": config['fmt_selected'],
            "fallbacks": [('masked', config['sel_masked']), ],
            },
        follows=[t01, t31],
            )
    t51 = dict(
        name='get mschdr',
        func=phot_mosaic.get_header,
        pipe='transform',
        in_=(config['sel_mosaic'], config['reg_inputs']),
        extras=config['jobdir'],
        out="{{basename[0]}}.{}".format(config['phot_hdr_suffix']),
        kwargs={
            'reg_inputs': config['reg_inputs'],
            'phot_hdr_glob': config['phot_hdr_glob'],
            },
        follows=[t50, t44]
            )
    t52 = dict(
        name='create mosaic',
        func='swarp',
        pipe='collate',
        in_=(t51, config['reg_grp']),
        add_inputs='{basename[0]}.fits',
        in_keys=[('dummy', 'in')],
        out=[fmtname(config['fmt_mosaic_orig']),
             fmtname(config['fmt_mosaic_wht'])],
        params={
            'HEADER_SUFFIX': '.{}'.format(config['phot_hdr_suffix']),
            'PIXELSCALE_TYPE': 'MANUAL',
            'PIXEL_SCALE': 0.20,
            'DELETE_TMPFILES': 'Y',
            'FSCALE_DEFAULT': '0/0',
            'BACK_SIZE': 64,
                },
        follows=t50,
        jobs_limit=1,
            )
    t53 = dict(
        name='apply whtmap',
        func=phot_mosaic.apply_weight,
        pipe='transform',
        in_=(t52, config['reg_mosaic']),
        out=fmtname(config['fmt_mosaic']),
            )

    tlist = [
            t00, t01,   # mask
            t10, t15,   # refcat
            t20, t21, t22, t23, t24,          # comb
            t30, t31,       # sub
            t40, t41, t42, t43, t44,      # phot
            t50, t51, t52, t53
            ]
    return tlist


def pyfunc(script, arg):
    return "{} -u {} {}".format(
            sys.executable,
            os.path.join(os.path.dirname(__file__), script),
            arg)


def fmtname(string):
    return re.sub(r'{(.+?)}', r'{\1[0]}', string)
