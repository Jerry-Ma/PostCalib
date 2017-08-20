#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2017-08-19 19:54
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
setup_package.py
"""

from __future__ import absolute_import

import os
from distutils.extension import Extension

ROOT = os.path.relpath(os.path.dirname(__file__))


def get_extensions():
    sources = ["podi_cython.pyx", "sigma_clip_mean.c", "sigma_clip_median.c"]
    include_dirs = ['numpy', os.path.dirname(__file__)]

    exts = [
        Extension(name='postcalib.qr.podi_cython',
                  sources=[os.path.join(ROOT, source) for source in sources],
                  include_dirs=include_dirs,
                  libraries=['gsl', 'gslcblas',  "m"]
                  )
    ]

    if os.environ.get('READTHEDOCS', False):
        return []
    else:
        return exts


def requires_2to3():
    return False
