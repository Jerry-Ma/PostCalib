#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2017-08-17 12:08
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
common.py
"""


from __future__ import (absolute_import, division, print_function)
from functools import wraps
import os


def get_log_func(default_level='debug', **kwargs):
    """return a convenient logger function"""
    logger = kwargs.pop('logger', None)
    logger_mutex = kwargs.pop('logger_mutex', None)

    def logfunc(*args):
        if len(args) == 1:
            level, mesg = default_level, args[0]
        else:
            level, mesg = args
        if logger_mutex is None or logger is None:
            print("[{0}] {1}".format(level, mesg))
        else:
            with logger_mutex:
                getattr(logger, level)(mesg)
    return logfunc


def touch_file(out_file):
    with open(out_file, 'a'):
        os.utime(out_file, None)


def create_symbolic_link(in_file, out_file, **kwargs):
    """create relative symbolic link

    :in_file: filename of original file
    :out_file: filename of symlink to create

    """
    # print(in_file)
    # print(out_file)
    log = get_log_func(**kwargs)
    in_file = os.path.abspath(in_file)
    out_file = os.path.abspath(out_file)
    in_dir = os.path.dirname(in_file)
    out_dir = os.path.dirname(out_file)
    if kwargs.get('relpath', True):
        src_path = os.path.relpath(in_file, out_dir)
    else:
        src_path = os.path.abspath(in_file)
    log('info', 'link {0} -> {1}'.format(src_path, out_file))
    if in_dir == out_dir or in_file == out_file:
        raise RuntimeError('target and source are from the same directory {0}'
                           .format(in_dir))
    if os.path.lexists(out_file):
        if not os.path.islink(out_file):
            raise RuntimeError("%s exists and is not a link" % out_file)
        try:
            os.unlink(out_file)
        except Exception as e:
            log('info', 'cannot unlink {0}: {1}'.format(out_file, e))
    os.symlink(src_path, out_file)
