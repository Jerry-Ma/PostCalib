#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2017-08-16 21:19
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
utils.py
"""


from __future__ import (absolute_import, division, print_function,
                        )
import re
from functools import wraps    # enable pickling of decorator


def ensure_list(value, tuple_ok=False):
    """
    Wrap the value into a list if not.

    Parameters
    ----------
    value: list, tuple, dict, or str
        The value to be wrapped
    tuple_ok: bool
        If True, tuple is treated as list and returned. Otherwise
        wrapped as a single value.

    Returns
    -------
    result: list
        The value if it is list or [value ] if not.
    """

    if tuple_ok:
        listclass = (list, tuple)
        elemclass = (str, dict, )
    else:
        listclass = list
        elemclass = (str, tuple, dict)
    if isinstance(value, elemclass) or callable(value):
        value = [value, ]
    elif value is None:
        value = []
    elif isinstance(value, listclass):
        if isinstance(value, tuple):
            value = list(value)
    else:
        raise RuntimeError("not able to ensure list type for {0}"
                           .format(value))
    return value


def unwrap_if_len_one(container):
    """
    Unwrap the container if there is only one element.

    Parameters
    ----------
    container: list or tuple
        The container to be unwrapped.

    Returns
    -------
    container: list or element-like
        The container or the first element of the container if that is
        the only one.
    """
    return container if len(container) > 1 else container[0]


def ensure_args_as_list(*iargs, **ikwargs):
    """
    Return a decorator that could Wrap string argument as list

    Parameters
    ----------
    positional args: int
        the indices of the arguments to be wrapped
    keyword args:
        keyword args that get passed to ensure_list function

    Returns
    -------
    wrapper: decorator
        Decorator used to wrap functions
    """
    def wrapper(func):
        @wraps(func)
        def wrapped_func(*args, **kwargs):
            newargs = list(args)
            for i in iargs:
                newargs[i] = ensure_list(args[i], **ikwargs)
            return func(*newargs, **kwargs)
        return wrapped_func
    return wrapper


def alert(string):
    """
    Highlight string in terminal
    """
    return '\033[91m{0}\033[0m'.format(string)


def de_alert(string):
    """
    Remove any escape sequence in string
    """
    return re.sub(r'\033\[\d+m', '', string)


class Namespace(object):
    """
    Hold attributes for named access
    """
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __repr__(self):
        keys = sorted(self.__dict__)
        items = ("{}={!r}".format(k, self.__dict__[k]) for k in keys)
        return "{}({})".format(type(self).__name__, ", ".join(items))

    def __eq__(self, other):
        return self.__dict__ == other.__dict__


# # monkey patch the multiprocess Pool
# import multiprocessing  # noqa: E402
# # We must import this explicitly, it is not imported by the top-level
# # multiprocessing module.
# import multiprocessing.pool  # noqa: E402


# class NoDaemonProcess(multiprocessing.Process):
#     # make 'daemon' attribute always return False
#     def _get_daemon(self):
#         return False

#     def _set_daemon(self, value):
#         pass

#     daemon = property(_get_daemon, _set_daemon)


# # We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# # because the latter is only a wrapper function, not a proper class.
# class MyPool(multiprocessing.pool.Pool):
#     Process = NoDaemonProcess


# multiprocessing.Pool = MyPool
