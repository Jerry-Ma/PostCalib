Command Line Interface
======================

PostCalib CLI
-------------

Full output of PostCalib help document.


.. code-block:: bash

    $ postcalib -h

    usage: postcalib [-h] [-l {DEBUG,INFO,WARNING,ERROR}]
                     {setup,init,run,example} ...

    A pipeline to post-process and stack QR calibrated WIYN ODI images

    optional arguments:
      -h, --help            show this (and more) help message and exit
      -l {DEBUG,INFO,WARNING,ERROR}, --log {DEBUG,INFO,WARNING,ERROR}
                            set the logging level

    actions:
      {setup,init,run,example}
                            available actions
        setup               setup workdir
        init                collect input images
        run                 run the pipeline

    setup:
      usage: postcalib setup [-h] [-f] [-o] [DEST_DIR]

      positional arguments:
        DEST_DIR         the dir to setup as workdir, default is the current dir

      optional arguments:
        -h, --help       show this help message and exit
        -f, --force      force the setup even if TARDIR is not empty
        -o, --overwrite  overwrite any existing file without bakcup in case a forced
                         setup is requested

    init:
      usage: postcalib init [-h] [-f] [-o] [-c CONFIG_FILE]
                            (-i [IMAGE [IMAGE ...]] | -l IMAGE_LIST | -d)
                            [JOBKEY [JOBKEY ...]]

      positional arguments:
        JOBKEY                the string identifier of the job to be initialized. It
                              will also be the basename of the jobdir to be created.
                              This can be omitted if an input list is used, in which
                              case the jobkey is the file basename

      optional arguments:
        -h, --help            show this help message and exit
        -f, --force           force the init even if the jobdir exists
        -o, --overwrite       overwrite existing job file without bakcup in case a
                              forced init is requested
        -c CONFIG_FILE, --config-file CONFIG_FILE
                              the config file to use. If omitted, look into the
                              current directory for one
        -i [IMAGE [IMAGE ...]], --input-images [IMAGE [IMAGE ...]]
                              the input images
        -l IMAGE_LIST, --input-list IMAGE_LIST
                              the input images specified via a text file, one per
                              line
        -d, --dirs-only       if specified, only the job directories are created

    run:
      usage: postcalib run [-h] [-a ...] [-c CONFIG_FILE] INPUT_TABLE

      positional arguments:
        INPUT_TABLE           the ASCII table generated from the init step
                              containing info of the inputs for the stacking
                              pipeline

      optional arguments:
        -h, --help            show this help message and exit
        -a ..., --apus-args ...
                              the arguments that get passed to apus
        -c CONFIG_FILE, --config-file CONFIG_FILE
                              the config file to use. If omitted, look into the
                              current directory for one

    example:

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
                                        # print a full help


Apus/Ruffus CLI
---------------


Full output of Apus command line options. These options can be
supplied when running `postcalib run ..` with `-a ...`.
Anything after `-a` will be passed to Apus.

.. code-block:: bash

    $ postcalib run fls201409.txt -a -h

    [INFO] fls201409: arguments passed to Apus ['-h']
    +- Apus powered by Ruffus ver 2.6.3 -+
    usage: postcalib run ... -a  [-h] [--verbose [VERBOSE]] [--version] [-L FILE]
                                 [-T JOBNAME] [-j N] [--use_threads] [-n]
                                 [--touch_files_only] [--recreate_database]
                                 [--checksum_file_name FILE] [--flowchart FILE]
                                 [--key_legend_in_graph]
                                 [--draw_graph_horizontally]
                                 [--flowchart_format FORMAT]
                                 [--forced_tasks JOBNAME] [-r] [-l]

    +- Astronomy Pipeline Using ruffuS, specifically tweaked for PostCalib -+

    optional arguments:
      -h, --help            show this help message and exit
      -r, --redo-all        force redo all tasks
      -l, --list-tasks      list the task names and exit

    Common options:
      --verbose [VERBOSE], -v [VERBOSE]
                            Print more verbose messages for each additional
                            verbose level.
      --version             show program's version number and exit
      -L FILE, --log_file FILE
                            Name and path of log file

    pipeline arguments:
      -T JOBNAME, --target_tasks JOBNAME
                            Target task(s) of pipeline.
      -j N, --jobs N        Allow N jobs (commands) to run simultaneously.
      --use_threads         Use multiple threads rather than processes. Needs
                            --jobs N with N > 1
      -n, --just_print      Don't actually run any commands; just print the
                            pipeline.
      --touch_files_only    Don't actually run the pipeline; just 'touch' the
                            output for each task to make them appear up to date.
      --recreate_database   Don't actually run the pipeline; just recreate the
                            checksum database.
      --checksum_file_name FILE
                            Path of the checksum file.
      --flowchart FILE      Don't run any commands; just print pipeline as a
                            flowchart.
      --key_legend_in_graph
                            Print out legend and key for dependency graph.
      --draw_graph_horizontally
                            Draw horizontal dependency graph.
      --flowchart_format FORMAT
                            format of dependency graph file. Can be 'pdf', 'svg',
                            'svgz' (Structured Vector Graphics), 'pdf', 'png'
                            'jpg' (bitmap graphics) etc
      --forced_tasks JOBNAME
                            Task(s) which will be included even if they are up to
                            date.
