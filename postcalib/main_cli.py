#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2017-08-10 09:23
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
main_cli.py

The command-line interface
"""


def main(args=None):

    import os
    # import sys
    import argparse
    import inspect
    import logging.config
    from .utils import PathType, RecursiveHelpAction  # , PlainHelpAction
    from . import core

    parser = argparse.ArgumentParser(
        add_help=False,
        description="""
A pipeline to post-process and stack QR calibrated WIYN ODI images"""
        )
    parser.add_argument(
            '-h', '--help', action=RecursiveHelpAction,
            help='show this (and more) help message and exit')
    parser.add_argument(
            "-l", "--log", dest="loglevel", default='DEBUG', nargs=None,
            choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
            help="set the logging level")

    subparsers = parser.add_subparsers(
            title="actions",
            help="available actions")

    # create the parser for the "setup" command
    parser_setup = subparsers.add_parser("setup", help="setup workdir")
    parser_setup.add_argument(
            "dir", type=PathType(exists=True, type="dir"),
            metavar="DEST_DIR",
            nargs='?', default=".",
            help="the dir to setup as workdir, default is the current dir",
            )
    parser_setup.add_argument(
            "-f", "--force", action="store_true",
            help="force the setup even if TARDIR is not empty",
            )
    parser_setup.add_argument(
            "-o", "--overwrite", action="store_true",
            help="overwrite any existing file without bakcup in case "
                 "a forced setup is requested"
            )

    def f_setup(option):
        core.setup_workdir(workdir=os.path.abspath(option.dir),
                           overwrite_dir=option.force,
                           backup_config=not option.overwrite)
    parser_setup.set_defaults(func=f_setup)

    # create the parser for the "init" command
    parser_init = subparsers.add_parser(
            "init", help="collect input images")
    jobkey_arg = parser_init.add_argument(
            "jobkey",
            metavar="JOBKEY", nargs='*',
            help="the string identifier of the job to be initialized. "
                 "It will also be the basename of the jobdir to be created. "
                 "This can be omitted if an input list is used, in which "
                 "case the jobkey is the file basename"
            )
    parser_init.add_argument(
            "-f", "--force", action="store_true",
            help="force the init even if the jobdir exists",
            )
    parser_init.add_argument(
            "-o", "--overwrite", action="store_true",
            help="overwrite  existing job file without bakcup in case "
                 "a forced init is requested"
            )

    config_file_arg = parser_init.add_argument(
            "-c", "--config-file", type=PathType(exists=True, type='file'),
            metavar="CONFIG_FILE",
            nargs=None,
            help="the config file to use. If omitted, look into the current "
                 "directory for one",
            )
    inputs_group = parser_init.add_mutually_exclusive_group(required=True)
    inputs_group.add_argument(
            "-i", "--input-images", type=PathType(exists=True, type='file'),
            metavar="IMAGE",
            nargs="*",
            help="the input images",
            )
    input_list_arg = inputs_group.add_argument(
            "-l", "--input-list", type=PathType(exists=True, type='file'),
            metavar="IMAGE_LIST",
            nargs=None,
            help="the input images specified via a text file, one per line",
            )
    inputs_group.add_argument(
            "-d", "--dirs-only", action='store_true',
            help="if specified, only the job directories are created",
            )

    def f_init(option):
        # parse the input file list if specified
        if option.input_images is not None:
            images = option.input_images
        elif option.input_list is not None:
            images = []
            with open(option.input_list, 'r') as fo:
                for ln in fo:
                    image = ln.strip()
                    if image.startswith("#") or image == "":
                        continue
                    if os.path.exists(image):
                        images.append(image)
                    else:
                        raise argparse.ArgumentError(
                                input_list_arg,
                                "invalid image path in image list: {}".format(
                                    image))
        elif option.dirs_only:
            images = []
        else:
            raise  # not possible
        # look for the config file to use
        if option.config_file is None:
            config_file = "postcalib.yaml"
            if not os.path.exists(config_file):
                raise argparse.ArgumentError(
                        config_file_arg,
                        "no valid config file found. Either specify one via "
                        " -c or run the command in a dir that has been setup "
                        "as workdir")
        else:
            config_file = option.config_file
        # infer the jobkey if possible
        if not option.jobkey:
            if option.input_list is None:
                raise argparse.ArgumentError(
                        jobkey_arg,
                        "no jobkey specified or can be inferred")
            else:
                jobkey = os.path.basename(option.input_list).split('.')[0]
        else:
            jobkey = option.jobkey[0].strip(os.path.sep)
        core.init_job(os.path.abspath(config_file), jobkey, images,
                      overwrite_dir=option.force,
                      backup_jobfile=not option.overwrite,
                      dirs_only=option.dirs_only)
    parser_init.set_defaults(func=f_init)

    # create the parser for the "run" command
    parser_run = subparsers.add_parser("run", help="run the pipeline")
    parser_run.add_argument(
            "jobfile",
            type=PathType(exists=True, type='file'),
            metavar="INPUT_TABLE", nargs=None,
            help="the ASCII table generated from the init step containing "
                 "info of the inputs for the stacking pipeline",
            )
    parser_run.add_argument(
            '-a', '--apus-args', nargs=argparse.REMAINDER, default=[],
            help="the arguments that get passed to apus")
    config_file_arg = parser_run.add_argument(
            "-c", "--config-file", type=PathType(exists=True, type='file'),
            metavar="CONFIG_FILE",
            nargs=None,
            help="the config file to use. If omitted, look into the current "
                 "directory for one",
            )

    def f_run(option):
        # look for the config file to use
        if option.config_file is None:
            config_file = "postcalib.yaml"
            if not os.path.exists(config_file):
                raise argparse.ArgumentError(
                        config_file_arg,
                        "no valid config file found. Either specify one via "
                        " -c or run the command in a dir that has been setup "
                        "as workdir")
        else:
            config_file = option.config_file
        core.run_pipeline(
                os.path.abspath(config_file),
                option.jobfile, apus_args=option.apus_args)
    parser_run.set_defaults(func=f_run)

    # create an example command to print out example workflow
    parser_example = subparsers.add_parser(
            "example",
            add_help=False,
            usage=argparse.SUPPRESS,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=" " """
  $ mkdir workdir                 # create an empty dir
  $ cd workdir                    # get into it
  $ postcalib setup               # setup the dir as postcalib workdir

  $ postcalib init testjob -i <paths to QR calibrated images>
                                  # init a job named testjob with given
                                  # QR calibrated images. This will generate
                                  # a testjob.in dir and a testjob.txt table

  $ vi testjob.txt                # edit job table to config how the
                                  # pipeline is run

  $ postcalib run testjob.txt ...
                                  # run the pipeline. The running of the
                                  # pipeline is powered by Apus,
                                  # run `postcalib run --apus-help` to
                                  # print a full help"""
                                  )  # noqa: W291
    parser_example.add_argument(
            '-h', '--help', action=argparse._HelpAction,
            help=argparse.SUPPRESS)

    def f_example(option):
        print(parser_example.description)
    parser_example.set_defaults(func=f_example)

    option = parser.parse_args(args)
    # logger level
    logging.config.dictConfig({
        'version': 1,
        'formatters': {
            'standard': {
                'format': '%(asctime)s [%(levelname)s] %(name)s: %(message)s'
            },
            'short': {
                'format': '[%(levelname)s] %(name)s: %(message)s'
            },
        },
        'handlers': {
            'default': {
                'class': 'logging.StreamHandler',
                'formatter': 'short',
            },
        },
        'loggers': {
            '': {
                'handlers': ['default'],
                'level': option.loglevel,
                'propagate': False
            },
        }
    })
    # execute the action
    if hasattr(option, 'func'):
        try:
            option.func(option)
        except RuntimeError as e:
            logger = inspect.trace()[-1][0].f_locals.get("logger", None)
            if logger is not None:
                logger.exception(e)
            else:
                print(e)
    else:
        parser.print_help()
