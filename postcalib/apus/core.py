#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-05-19 23:23
# Python Version :  2.7.10
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
core.py

"""

from __future__ import (absolute_import, division, print_function)

import os
import re
import sys
import glob
import time
import logging
import pickle
from io import StringIO
import subprocess
from copy import copy
from functools import wraps    # enable pickling of decorator
from datetime import timedelta
from functools import reduce
# from collections import Iterable
# from tempfile import NamedTemporaryFile
from . import utils
from .utils import ensure_list, unwrap_if_len_one, ensure_args_as_list
from . import astromatic as am
from . import common

import ruffus
import ruffus.cmdline as cmdline
from ruffus.proxy_logger import make_shared_logger_and_proxy
# from ruffus import pipeline_printout_graph
# from ruffus import mkdir
from ruffus import Pipeline
from ruffus import formatter
# from ruffus import add_inputs
from ruffus import output_from
from ruffus.ruffus_exceptions import error_ambiguous_task


class ApusConfig(object):
    """
    Hold configurations for running the pipeline
    """

    dirs = [
        ('jobdir', None), ('confdir', '{jobdir:s}'), ('diagdir', '{jobdir:s}'),
        ('task_io_default_dir', '{jobdir:s}'),
        ('logdir', ''),
        ]
    specs = [('jobkey', None), ('inputs', []), ('tlist', [])]
    runtime = [
        ('env_overrides', {}),
        ('logger', None), ('logger_mutex', None),
        ('log_file', '{jobkey:s}.log'), ('history_file', '{jobkey:s}.ruffus')
        ]

    def __init__(self, config=None, **kwargs):
        # mark the begin time
        self.timestamp = time.time()
        missingkeyerror = "missing key in apus config: {0}"
        if config is None:
            config = utils.Namespace(**kwargs)
        for key, defval in self.dirs + self.specs + self.runtime:
            try:
                val = getattr(config, key)
            except AttributeError:
                if defval is None and key not in dict(self.runtime).keys():
                    raise RuntimeError(missingkeyerror.format(key))
                elif isinstance(defval, str):
                    val = defval.format(**self.__dict__)
                else:
                    val = defval
                setattr(config, key, val)
            setattr(self, key, val)

    def get_dirs(self):
        return [d for d in set(getattr(self, k)
                for k in dict(self.dirs).keys()) if d != '']

    def get_task_names(self):
        return [t['name'] for t in self.tlist]


def configure(config, args):
    """
    Setup runtime from config module/dict and command line args

    Parameters
    ----------
    config: dict or Namespace
        Hold configurations used to initialize ApusConfig object
    args: list
        list of arguments to be passed to Ruffus.cmdline module

    Returns
    -------
    apusconf: ApusConfig
        Hold configurations of the Apus
    option: Namespace
        Hold parsed command line arguments
    """
    if isinstance(config, dict):
        apusconf = ApusConfig(**config)
    else:
        apusconf = ApusConfig(config=config)

    parser = cmdline.get_argparse(description="""
+- Astronomy Pipeline Using ruffuS, specifically tweaked for PostCalib -+
""", version=ruffus.__version__, prog='postcalib run ... -a ')
    parser.add_argument(
            '-r', '--redo-all', action='store_true',
            help='force redo all tasks')
    parser.add_argument(
            '-l', '--list-tasks', action='store_true',
            help='list the task names and exit')

    parser.set_defaults(
            verbose=['0', ],
            log_file=os.path.join(apusconf.logdir, apusconf.log_file),
            history_file=os.path.join(apusconf.logdir, apusconf.history_file)
            )
    option = parser.parse_args(args)
    # handle logger
    logger, logger_mutex = make_shared_logger_and_proxy(
            logger_factory, apusconf.jobkey, [option.log_file, option.verbose])
    apusconf.logger = logger
    apusconf.logger_mutex = logger_mutex
    return apusconf, option


def build_pipeline(config):
    """
    Assemble the pipeline

    Parameters
    ----------
    config: ApusConfig object
        Hold the configurations

    Returns
    -------
    pipe: Pipeline object
        The pipeline object with tasks setup
    """
    pipe = Pipeline(name=config.jobkey)
    # mkdirs
    # t00 = {'name': 'create jobdir'}
    # dirs = config.get_dirs()
    # pipe.mkdir(dirs, name=t00['name'])
    for task in config.tlist:
        create_ruffus_task(pipe, config, task)
    return pipe


def bootstrap(config=None, option=None):
    """entry point; parse command line argument, create pipeline object,
    and run it
    """
    print("+- Apus powered by Ruffus ver {0} -+".format(ruffus.__version__))
    if config is None:
        config = sys.modules['__main__']
    if option is None:
        option = sys.argv[1:]
    config, option = configure(config, option)
    if option.list_tasks:
        tlist = config.get_task_names()
        if not tlist:
            print("no tasks found")
        else:
            for t in tlist:
                print("{0}".format(t))
        sys.exit(0)
    # set up astromatic config
    config.am = am.AmConfig(**config.env_overrides)
    build_pipeline(config)
    # handle redo-all
    if option.redo_all:
        task_list = ruffus.pipeline_get_task_names()
        option.forced_tasks.extend(task_list)
    if len(option.forced_tasks) > 0:
        for t in option.forced_tasks:
            config.logger.info("forced redo: {0}".format(utils.alert(t)))
    cmdline.run(option, checksum_level=1)


def get_task_func(func):
    """
    Return pipeline function to be executed for a given task func
    """
    errmsg = 'not a meaningful task func value: {0}'.format(func)
    if isinstance(func, str):
        if func.lower() in ['sex', 'scamp', 'swarp']:
            return astromatic_task
        else:
            command = func.split()
            if '{in}' in command or '{out}' in command:
                return subprocess_task
            else:
                raise RuntimeError(errmsg)
    elif callable(func):
        return callable_task
    else:
        raise RuntimeError(errmsg)


def aggregate_task_inputs(inputs):
    """
    Sort-out the inputs.

    Parameters
    ----------
    inputs: list of (list|tuple|str|callable)
        The inputs specified in task inputs. Each element represent
        one set of input. list: a list of filenames;
        tuple: (filename(s), formatter) pair; str: treated as a list
        of one; callable: generator that yield filenames, maximum of
        one generator input is allowed

    Returns
    -------
    formatter_inputs: list of tuple
        Inputs that have formatter along with
    simple_inputs: list of str
        Inputs without formatter
    generator_inputs: list with length <= 1
        Input specified as a generator if any
    """
    tuple_num_elem = 2  # allow pair-like tuple
    formatter_inputs = []
    simple_inputs = []
    generator_inputs = []
    for in_ in inputs:
        if callable(in_):
            generator_inputs.append(in_)
        elif isinstance(in_, tuple) and len(in_) == tuple_num_elem:
            if callable(in_[0]):
                raise RuntimeError('generator input should not have formatter'
                                   ' {0}'.format(in_))
            formatter_inputs.append(in_)
        elif isinstance(in_, list):
            simple_inputs.extend(in_)
        elif isinstance(in_, (str, dict, )):
            simple_inputs.append(in_)
        else:
            raise RuntimeError('invalid input {0}'.format(in_))
    if len(generator_inputs) > 1:
        raise RuntimeError('found multiple generator inputs {0}'.format(
            generator_inputs))
    return formatter_inputs, simple_inputs, generator_inputs


def normalize_taskname(name):
    """
    Returns escaped task name.
    """
    return name.replace(' ', '_').lower()


def create_ruffus_task(pipe, config, task, **kwargs):
    """
    Create Ruffus task from the task dictionary and add to pipe

    Task keys: name, func, pipe, [in_, out]
    Optional task keys: add_inputs, replace_inputs,
                        in_keys, out_keys,
                        ruffus_kwargs

    Parameters
    ----------
    pipe: Ruffus Pipeline object
        The pipeline the created task will be added to
    config: ApusConfig object
        The configuration object
    task: dict
        The dict defines the task
    keyword args: optional
        Additional config entries to be updated to the
        configuration object
    """

    # validate the task dict first
    missingkeyerror = "missing key in task dict: {0}"
    for key in ['name', 'func', 'pipe', ['in_', 'out']]:
        if all(k not in task.keys() for k in ensure_list(key)):
            raise RuntimeError(missingkeyerror.format(key))
    config = copy(config)
    config.__dict__.update(**kwargs)

    # process task dict
    task_name_list = config.get_task_names()

    pipe_func = getattr(pipe, task['pipe'])
    task_name = task['name']
    task_func = get_task_func(task['func'])

    task_args = []
    task_kwargs = {'name': task_name, 'task_func': task_func}
    # handle follows
    task_follows = [t['name'] if isinstance(t, dict) else t
                    for t in ensure_list(task.get('follows', []))]
    # additional logic for astromatic tasks
    if task_func.__name__ == 'astromatic_task':
        task['params'] = am.normalize_params(task.get('params', {}))
        # generate configuration file if not supplied
        if 'conf' not in ensure_list(task.get('in_keys', None)) or \
                task.get('auto_conf', True):
            pre_task = {
                'name': 'autoconf {0}'.format(task_name),
                'func': dump_config_files,
                'pipe': 'originate',
                'out': os.path.join(config.confdir, 'conf.{0}'.format(
                        normalize_taskname(task_name))),
                'extras': os.path.join(
                    config.confdir,
                    'conf.{0}.checker'.format(normalize_taskname(task_name))),
                'params': task.pop('params'),  # remove params from this task
                'outparams': task.get('outparams', []),
                'verbose': False,
                'check_if_uptodate': check_config_uptodate,
                'diagdir': config.diagdir,
                'prog': get_am_prog(task['func']),
                'jobs_limit': 1
                }
            create_ruffus_task(
                pipe, config, pre_task, task_io_default_dir='')
            task_follows.append(pre_task['name'])
            # connect pre_task to this
            # conf_inputs = task.get('extras', [])
            # conf_inputs.append(os.path.abspath(pre_task['out']))
            # task_inkeys.append('conf')
            conf_inputs = ensure_list(task.get('add_inputs', None))
            conf_inputs.insert(0, os.path.abspath(pre_task['out']))
            task_inkeys = copy(task.get('in_keys', ['in', ]))
            _inkeys = ['in', 'in+']
            for i, key in enumerate(task_inkeys):
                if key in _inkeys:
                    task_inkeys[i] = (key, 'conf')
                    break
                elif isinstance(key, tuple) and \
                        any(j in _inkeys for j in key):
                    _keys = list(key)
                    _keys.insert(1 - len(conf_inputs), 'conf')
                    task_inkeys[i] = tuple(_keys)
                    break
                elif isinstance(key, list) and \
                        any(j in _inkeys for j in key):
                    _keys = tuple(key + ['conf', ])
                    task_inkeys[i] = _keys
                    break
            task['add_inputs'] = conf_inputs
            task['in_keys'] = task_inkeys

    # handle input
    def handle_input(in_, in_key, formatter_key):
        formatter_inputs, simple_inputs, generator_inputs = \
            aggregate_task_inputs(
                ensure_list(task.get(in_, []), tuple_ok=False))
        if len(generator_inputs) > 0:  # generator_inputs goes to unnamed arg
            task_args.extend(generator_inputs)
        # simple_inputs get common general formatter
        if len(simple_inputs) > 0:
            formatter_inputs.append((simple_inputs, r'.+'))
        # handle formatter_inputs
        task_inputs = []
        task_formatters = []
        for in_, reg in formatter_inputs:
            in_ = [i['name'] if isinstance(i, dict) else i
                   for i in ensure_list(in_)]
            temp_in = []
            temp_reg = reg
            for i in in_:
                if i in task_name_list:
                    temp_in.append(output_from(i))
                    continue
                elif not os.path.isabs(i):  # prepend default io dir
                    i = os.path.join(config.task_io_default_dir, i)
                    temp_reg = r'(?:[^/]*/)*' + reg
                else:
                    pass
                if re.search(r'[?*,\[\]{}]', i) is not None:
                    config.logger.info('{0:^27s} (glob): {1}'.format(
                        task_name, i))
                    temp_in.append(i)
                else:
                    config.logger.info('{0:^27s} (file): {1}'.format(
                        task_name, i))
                    temp_in.append(i)
            task_inputs.append(temp_in)  # list of list
            task_formatters.append(temp_reg)  # list of regex
        if len(task_inputs) > 0:
            task_inputs = reduce(lambda a, b: a + b, task_inputs)  # flatten
            if len(task_inputs) > 0:
                task_kwargs[in_key] = unwrap_if_len_one(task_inputs)
            if task['pipe'] != 'merge':  # require formatter for non-merge pipe
                task_kwargs[formatter_key] = formatter(*task_formatters)
    handle_input('in_', 'input', 'filter')
    if 'in2' in task.keys():
        handle_input('in2', 'input2', 'filter2')

    def resolve_task_name(s):
        if isinstance(s, dict):
            s = s['name']
        if s in task_name_list:
            try:  # have to replace task name with task
                s, = pipe.lookup_task_from_name(s, "__main__")
            except (ValueError, error_ambiguous_task):
                pass
        return s
    # handle additional inputs and replace_inputs
    for inkey in ['add_inputs', 'replace_inputs']:
        task_inkey = []
        for in_ in ensure_list(task.get(inkey, None)):
            in_ = resolve_task_name(in_)
            if isinstance(in_, str) and not os.path.isabs(in_):
                in_ = os.path.join(config.task_io_default_dir, in_)
            task_inkey.append(in_)
        if len(task_inkey) > 0:
            task_kwargs[inkey] = tuple(task_inkey)  # ruffus req.
    # handle outputs
    task_output = []
    for out in ensure_list(task.get('out', None)):
        if not os.path.isabs(out):
            out = os.path.join(config.task_io_default_dir, out)
        task_output.append(out)
    if len(task_output) > 0:
        task_kwargs['output'] = unwrap_if_len_one(task_output)
    else:  # flag file for checkpointing
        task_kwargs['output'] = os.path.join(
            config.task_io_default_dir,
            normalize_taskname(task_name) + '.success')
    # handle context as extra
    task_extras = []
    for extra in ensure_list(task.get('extras', None)):
        if isinstance(extra, str) and not os.path.isabs(extra):
            task_extras.append(os.path.join(config.task_io_default_dir, extra))
        else:
            task_extras.append(extra)
    # context is parameters that are passed to the task function
    context_exclude_task_keys = [
            'in_', 'out', 'add_inputs', 'allow_slice']
    context_key_defaults = {
            'verbose': True,
            'kwargs': {},
            }
    context = {k: v for k, v in list(task.items()) +
               list(context_key_defaults.items())
               if k not in context_exclude_task_keys}
    for key, defval in context_key_defaults.items():
        context[key] = task.get(key, defval)
    if 'follows' in context.keys():
        # for cleaner debug info
        context['follows'] = unwrap_if_len_one(task_follows)
    task_context = {
            'task': context,
            'logger': config.logger,
            'logger_mutex': config.logger_mutex,
            'am': config.am
            }
    task_extras.append(task_context)
    task_kwargs['extras'] = task_extras
    # create ruffus task
    ruffus_task = pipe_func(*task_args, **task_kwargs)
    if len(task_follows) > 0:
        ruffus_task.follows(*task_follows)
    # handle job_limit
    jobs_limit = task.get('jobs_limit', None)
    if jobs_limit is not None:
        ruffus_task.jobs_limit(jobs_limit)
    # add finish signal
    ruffus_task.posttask(task_finish_signal(task_name, config))
    # handle forced run
    if task.get('check_if_uptodate', None) is not None:
        ruffus_task.check_if_uptodate(task['check_if_uptodate'])
    return ruffus_task


def logger_factory(logger_name, args):
    """
    Provide logging modules, adapted from Ruffus, with some small tweaks
    """
    log_file_name, verbose = args
    new_logger = logging.getLogger(logger_name)

    class debug_filter(logging.Filter):
        """ignore INFO messages"""
        def filter(self, record):
            return logging.INFO != record.levelno

    class NullHandler(logging.Handler):
        """for when there is no logging"""
        def emit(self, record):
            pass

    new_logger.setLevel(logging.DEBUG)
    has_handler = False

    # log to file if that is specified
    if log_file_name:
        handler = logging.FileHandler(log_file_name, delay=False)

        class stipped_down_formatter(logging.Formatter):
            def format(self, record):
                prefix = ""
                if not hasattr(self, "first_used"):
                    self.first_used = True
                    prefix = "\n" + self.formatTime(record, "%Y-%m-%d")
                    prefix += " %(name)s\n" % record.__dict__
                self._fmt = " %(asctime)s - %(levelname)-7s - %(message)s"
                old_msg = record.msg
                record.msg = utils.de_alert(record.msg)
                out = prefix + logging.Formatter.format(self, record)
                record.msg = old_msg
                return out
        handler.setFormatter(stipped_down_formatter(
            "%(asctime)s - %(name)s - %(levelname)6s - %(message)s", "%H:%M:%S"
            ))
        handler.setLevel(logging.DEBUG)
        new_logger.addHandler(handler)
        has_handler = True

    # log to stderr if verbose
    if verbose:
        stderrhandler = logging.StreamHandler(sys.stderr)
        stderrhandler.setFormatter(logging.Formatter("[%(name)s] %(message)s"))
        stderrhandler.setLevel(logging.DEBUG)
        # if log_file_name:
        #     stderrhandler.addFilter(debug_filter())
        new_logger.addHandler(stderrhandler)
        has_handler = True

    # no logging
    if not has_handler:
        new_logger.addHandler(NullHandler())

    return new_logger


def task_finish_signal(task_name, config):
    """
    Print task execution time upon finishing
    """
    def ret_task_finish_signal():
        # refresh flowchart
        # subprocess.call(['python', sys.argv[0],
        #                  '-T', config.tlist[-1]['name'],
        #                  '--flowchart', 'test.png'])
        with config.logger_mutex:
            elapsed_time = timedelta(
                    seconds=time.time() - config.timestamp)
            config.logger.debug('task {0:^45s} finished @ {1:s}'
                                .format(utils.alert(task_name),
                                        utils.alert(elapsed_time))
                                )
    return ret_task_finish_signal


def get_flag_file(out_files):
    suffix = '.success'
    if len(out_files) >= 1 and out_files[-1][-len(suffix):] == suffix:
        flag = out_files[0]
        out_files = out_files[:-1]
    else:
        flag = None
    return out_files, flag


def get_am_prog(func):
    program = [i for i in ['sex', 'scamp', 'swarp']
               if func.lower().startswith(i)][0]
    return program


def to_callable_task_args(conv_func):
    def wrapper(func):
        @wraps(func)
        def wrapped_func(in_files, out_files, *extras):
            if len(extras) == 0:  # from originate, rename the variables
                in_files, out_files, extras = [], in_files, out_files
            in_files += extras[:-1]
            out_files, flag_file = get_flag_file(out_files)
            # flatten any third level list
            for i, in_ in enumerate(in_files):
                if isinstance(in_, tuple):
                    _in_ = []
                    for item in in_:
                        if isinstance(item, list):
                            _in_.extend(item)
                        else:
                            _in_.append(item)
                    in_files[i] = tuple(_in_)
            # print(in_files, out_files)
            try:
                overlap = list(set(in_files).intersection(set(out_files)))
                if len(overlap) > 0:
                    raise RuntimeError(
                            'danger: output {0} has same filename as inputs'
                            .format(unwrap_if_len_one(overlap)))
            except TypeError:
                pass
            context = copy(extras[-1])
            context['flag_file'] = flag_file
            if conv_func is not None:
                # copy task as well
                context['task'] = dict(context['task'], func=conv_func(
                        in_files, out_files, context))
            return func(in_files, out_files, context)
        return wrapped_func
    return wrapper


def _subprocess_callable(in_files, out_files, context):
    # split tuples from string
    flattened_files = []
    for in_ in in_files:
        if isinstance(in_, tuple):
            flattened_files.extend(list(in_))
        else:
            flattened_files.append(in_)
    in_files = flattened_files
    command = context['task']['func'].split()
    for key, value in zip(['{in}', '{out}'], [in_files, out_files]):
        if key in command:
            i = command.index(key)
            command[i:i + 1] = value
    return documented_subprocess_call(command, flag_file=context['flag_file'])


@ensure_args_as_list(0, 1, tuple_ok=True)
@to_callable_task_args(_subprocess_callable)
def subprocess_task(in_files, out_files, context):
    """run command as subprocess by replacing func key in task context
    with callable"""
    return callable_task(in_files, out_files, context)


def _astromatic_callable(in_files, out_files, context):
    """return the command for executing astromatic task"""
    task = context['task']
    amconf = context['am']
    # get program
    prog = get_am_prog(task['func'])
    # split up inputs types
    rectified_inputs = get_astromatic_inputs(in_files, task['in_keys'])
    command = [amconf.get('{0}bin'.format(prog)), ]
    params = {}
    for key, val in rectified_inputs:
        # here val is a list of values
        if key == 'in':
            command.extend(val)
        elif key == 'conf':
            if len(val) < 1:
                raise RuntimeError('no configuration file specified for {0}'
                                   .format(task['func']))
            else:
                command.extend(['-c', val[0]])
                params.update(am.parse_astromatic_conf(*val[1:]))
        elif key == 'dummy':
            pass
        else:  # values should be concat by comma
            params[key] = ','.join(val)
    params.update(task.get('params', {}))
    for key, val in params.items():
        command.extend(['-{0}'.format(key), "{0}".format(val)])
    # handle outkeys
    default_outkeys = {
            'sex': ['CATALOG_NAME', ],
            'scamp': [],
            'swarp': ['IMAGEOUT_NAME', 'WEIGHTOUT_NAME']
            }
    out_keys = ensure_list(task.get('out_keys', default_outkeys[prog]))
    for i, key in enumerate(out_keys):
        command.extend(['-{0}'.format(key), out_files[i]])
    return documented_subprocess_call(command, flag_file=context['flag_file'])


@ensure_args_as_list(0, 1)
@to_callable_task_args(_astromatic_callable)
def astromatic_task(in_files, out_files, context):
    """create astromatic command to run by subprocess"""
    return callable_task(in_files, out_files, context)


def get_astromatic_inputs(inputs, in_keys):
    """identify the inputs by mapping to the in_keys"""

    # collate, no add_input: [(in1, in2, ...), ex1, ex2]
    # key: [key1, ekey1, ekey2]

    # collate, add_input: [((in1, add1, add2), (in2, add1, add2)), ex1, ex2]
    # key: [(key1, key2 key3), ekey1, ekey2]

    # transform, no add_input: [ in, ex1, ex2]
    # key: [key1, ekey1, ekey2]

    # transform, add_input: [(in, add1), ex1, ex2]
    # key: [(key1, key2 key3), ekey1, ekey2]

    # ret:
    #   [(key1, [in1, in2, .. ]), (key2, [add1_1, add1_2 ..])
    # print(inputs)
    # print(in_keys)
    ret_keys = []
    ret_vals = []
    # expand in+
    if 'in+' in in_keys:
        nin = len(inputs) - len(in_keys)
        iin = in_keys.index('in+')
        # print("expand in+ to")
        in_keys[iin:iin + 1] = ['in', ] * (nin + 1)
    # first level zip for exkeys
    for key, val in zip(in_keys, inputs):
        if isinstance(key, tuple):  # deal with tuple keys:
            if all(isinstance(v, tuple) for v in val):
                val = list(zip(*val))  # tuple of tuple, need to be zipped
            if len(key) == len(val):
                for k, vv in zip(key, val):
                    if isinstance(vv, tuple):
                        _v = tuple(set(vv))
                        vv = _v[0] if len(_v) == 1 else vv
                    ret_keys.append(k)
                    ret_vals.append(vv)
            else:
                raise RuntimeError(
                    "mismatched number of"
                    " keys ({0}) and inputs {1}".format(len(key), len(val)))
        else:
            val = unwrap_if_len_one(val)
            if isinstance(val, str) and re.search('[*?]', val) is not None:
                val = tuple(glob.glob(val))
            ret_keys.append(key)
            ret_vals.append(val)
    # aggregate duplicated keys
    ag_keys = list(set(ret_keys))
    ag_vals = [[] for _ in range(len(ag_keys))]
    for i, key in enumerate(ret_keys):
        ag_vals[ag_keys.index(key)].append(ret_vals[i])
    for i, val in enumerate(ag_vals):
        # print(i, val)
        # validate list-type keys
        if any([isinstance(v, tuple) for v in val]):
            if len(val) > 1:
                raise RuntimeError("list-type key {0} should only appear once"
                                   " in in_keys".format(ag_keys[i]))
            else:
                ag_vals[i] = val[0]
    return zip(ag_keys, ag_vals)


@ensure_args_as_list(0, 1)
@to_callable_task_args(None)
def callable_task(in_files, out_files, context):
    log = common.get_log_func(**context)
    task = context['task']
    func = task['func']
    args = in_files + out_files
    verbose = task.get('verbose', True)
    if func.__doc__ is not None and func.__doc__.startswith('subprocess: '):
        caller_string = func.__doc__[len('subprocess: '):]
    else:
        caller_string = "{0}({1})".format(
                func.__module__ + "." + func.__name__, ', '.join(
                    map(str, args)))
    if verbose:
        log('debug', caller_string)
    kwargs = dict(
            task['kwargs'],
            task=task,
            am=context['am'],
            logger=context['logger'],
            logger_mutex=context['logger_mutex'])
    output = func(*args, **kwargs)
    if task.get('after_func', None) is not None:
        task['after_func'](*out_files)
    flag_file = context['flag_file']
    if flag_file is not None:
        common.touch_file(flag_file)
    if output or verbose:
        output = "finished silently" if not output else 'finished'
        log('debug', output)


def dump_config_files(conf_file, checker_file, **kwargs):
    """create configuration file"""

    logger, logger_mutex = kwargs['logger'], kwargs['logger_mutex']
    task = kwargs['task']
    amconf = kwargs['am']
    prog = task['prog']
    am_bin = amconf.get('{0}bin'.format(prog))
    am_share = amconf.get('{0}share'.format(prog))
    # handle parameters
    conf_params = dict(amconf.get('{0}_default'.format(prog)),
                       **task.get('params', {}))
    if prog == 'sex':
        if 'PARAMETERS_NAME' not in task.get('in_keys', []):
            params_file = conf_params.get('PARAMETERS_NAME', None)
            if params_file is None or not os.path.isfile(params_file):
                params_file_keys = task.get('outparams', [])
                params_file = conf_file + '.sexparam'
                with logger_mutex:
                    logger.info('params file: {0}'.format(params_file))
                fo = StringIO(
                        subprocess.check_output([am_bin, '-dp']).decode())
                am.dump_sex_param(
                        fo, params_file, amconf.get('sexparam_default'),
                        params_file_keys, overwrite=True)
                # with logger_mutex:
                #     logger.info("keys: {0}".format(', '.join(keys)))
                fo.close()
                conf_params['PARAMETERS_NAME'] = params_file
    elif prog == 'scamp':
        conf_params['CHECKPLOT_NAME'] =\
            am.get_scamp_checkplot_name(
                    os.path.abspath(task['diagdir']),
                    prefix=normalize_taskname(task['name']))
    elif prog == 'swarp':
        if conf_params.get('RESAMPLE_DIR', None) is None:
            conf_params['RESAMPLE_DIR'] = amconf.get('tmpdir')
    # convert to strings
    # use absolute file paths whenever possible
    for k, v in conf_params.items():
        if isinstance(v, list):
            v = ', '.join(map(str, v))
        else:
            v = str(v)
        conf_params[k] = v
        if re.search('[/.]', v) is not None:  # only replace file-like vals
            if os.path.isfile(v):
                conf_params[k] = os.path.abspath(v)
            elif os.path.isfile(os.path.join(am_share, v)):
                conf_params[k] = os.path.abspath(os.path.join(am_share, v))
            else:
                conf_params[k] = v
    # create conf file
    fo = StringIO(subprocess.check_output([am_bin, '-dd']).decode())
    am.dump_astromatic_conf(fo, conf_file, overwrite=True, **conf_params)
    with logger_mutex:
        logger.info('conf file: {0}'.format(conf_file))
        for k, v in conf_params.items():
            if os.path.isfile(v):
                v = os.path.relpath(v)
            logger.info('{0:>20s}: {1:s}'.format(k, v))
    fo.close()
    # write checker file
    with open(checker_file, 'wb') as fo:
        pickle.dump(task.get('params', {}), fo)
        pickle.dump(task.get('outparams', {}), fo)


def check_config_uptodate(*args, **kwargs):
    conf_file, checker_file, context = args[-3:]
    for f in [conf_file, checker_file]:
        if not os.path.isfile(f):
            return True, "missing file {0}".format(f)
    task = context['task']
    with open(checker_file, 'rb') as fo:
        try:
            old_params = pickle.load(fo)
            old_outparams = pickle.load(fo)
        except EOFError:
            return True, "corrupted checker"
    new_params = task.get('params', {})
    new_outparams = task.get('outparams', [])
    if len(new_params) != len(old_params) or \
            len(new_outparams) != len(old_outparams):
        return True, 'params/outparams changed its size'
    for key, val in new_params.items():
        if key in old_params.keys() and old_params[key] == val:
            continue
        else:
            return True, 'params dict changed its content'
    else:
        overlap = list(set(new_outparams).intersection(old_outparams))
        if len(overlap) == len(new_outparams):
            # print out params
            return False, "no change of params/outparams"
        else:
            return True, 'outparams list changed its content'


def documented_subprocess_call(command, flag_file=None):
    def call(*args, **kwargs):
        # handle scamp refcatalog suffix
        # if '-ASTREFCAT_NAME' in command:
        #     ikey = command.index('-ASTREFCAT_NAME') + 1
        #     refcatkey = command[ikey]
        #     refcatfiles = glob.glob(
        #             "{0}_?{1}".format(*os.path.splitext(refcatkey)))
        #     if len(refcatfiles) == 0:
        #         refcatfiles = glob.glob(refcatkey)
        #     if len(refcatfiles) > 0:
        #         command[ikey] = ','.join(refcatfiles)
        log = common.get_log_func(**kwargs)

        # with NamedTemporaryFile() as fo:
        #     subprocess.check_call(command, stdout=fo,
        #                           stderr=subprocess.STDOUT)
        #     fo.seek(0)
        #     output = fo.read()
        # output = subprocess.check_output(command)
        proc = subprocess.Popen(command,
                                stdout=subprocess.PIPE,
                                bufsize=1,
                                # stderr=subprocess.PIPE
                                )
        has_output = False
        for ln in iter(proc.stdout.readline, b''):
            if ln:
                has_output = True
            log('debug', ln.strip('\n'))
        if proc.poll() is not None and proc.returncode != 0:
            err_msg = "subprocess failed with code {0}".format(proc.returncode)
            # an error happened!
            # err_msg = "%s\nsubprocess failed with code: %s" % (
            #         proc.stderr.read(),
            #         proc.returncode)
            raise RuntimeError(err_msg)
        if flag_file is not None:
            common.touch_file(flag_file)
        return has_output
        # return output
    # print(' '.join(command))
    # quote items with string
    call.__doc__ = 'subprocess: ' + ' '.join(
            ['"{0}"'.format(c) if ' ' in c else c for c in command])
    return call
