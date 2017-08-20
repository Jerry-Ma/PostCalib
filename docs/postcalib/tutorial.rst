
Tutorial
========

This tutorial aims at give the user a step by step instruction
of how to use PostCalib.

Overview
--------
The set-up and running of PostCalib is control via various
command lines.

Below is (part of) the listing of the help message from the built-in
main help:

.. code-block:: bash

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


In the most general case, one would follow the order of

1. setup workdir

2. collect input images

3. run the pipeline

Workflow
--------

.. _setup-details:

1. Setup
^^^^^^^^

.. code-block:: bash

        $ postcalib setup -h

        setup:
      usage: postcalib setup [-h] [-f] [-o] [DEST_DIR]

      positional arguments:
        DEST_DIR         the dir to setup as workdir, default is the current dir

      optional arguments:
        -h, --help       show this help message and exit
        -f, --force      force the setup even if TARDIR is not empty
        -o, --overwrite  overwrite any existing file without bakcup in case a forced
                         setup is requested

This step is to setup runtime environments/configurations for PostCalib
in an empty directory. If the `DEST_DIR` is not empty, one need to supply
`-f` switch to force the setup.

The setup will first check the external, non-python command line tools
that PostCalib depends on. To be specific, the external dependencies
that are checked during setup are

    1. SExtractor, SCAMP, and SWarp
    2. stilts
    3. fpack/funpack

For the SExtractor, SCAMP and SWarp suite, the program further requires
that they are properly installed in to a common prefix path. This can
be done by supplying `--prefix=<path>` to `./configure` when compiling
the codes, and issue `make && make install`.

For the latter two items, if they are not present in the `PATH` environment
variable, PostCalib will download the package from online and automatically
store them in the `<workdir>/extern` folder. This folder is by default
in the PATH of PostCalib so this downloading will occur only once.

In addition to the directory structures and the extern dependencies,
the setup step also creates a runtime configuration file named `postcalib.yaml`
in the work directory. This file by default contains entries that
defines variables that affect how the pipeline is run.

For most cases, the default values are sufficient, and users do not need to
modify the file. If running setup in a directory that already contains
such a file, it will be by default backed up. Supplying `-o` to the
command line arguments could force the file being overwritten.


2. Init
^^^^^^^^

.. code-block:: bash

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


This step is to collect input files to a single directory and make ready
for the execution of pipeline in the next step.

`JOBKEY` is a required argument for `init` and it is used as the unique
identifier of the job to be run. The collected input images are stored
in a subdirectory named `{JOBKEY}.in' in the `{workdir}`, and an empty
directory named `{JOBKEY}.skymask` is also created along side, which is
to be used for supplying user defined mask regions.

The inputs of the job can either be specified via `-i` followed by
a list of images in the command line, or via `-l` followed by a
text file containing the valid paths pointing to input image files.
If `-d` is specified instead, no inputs are read and only the files
and directories are re-generated.

The inputs file specified are symlinked to the `{JOBKEY}.in` directory
and preview png images are created for visual check. Finally, a
table of the collected images are created, one image per entry, with
the columns typically containing:

.. code-block:: bash

 numid fcomb_group fsub_group phot_group mosaic_group mask_otas
 OBSID CRVAL1 CRVAL2 OBJECT  EXPMEAS  AIRMASS  SEEING
 SKY_MEDI   SKY_STD   FILTER   INSTRUME       MJD-OBS
 filenam


The columns ending with `_group` are used as input flags for
controlling the execution of the pipeline. Additional columns
could be added to the table provided that it is stored as a
header keywords of the input images.

We refer to this table as the **job table**


3. Run
^^^^^^

.. code-block:: bash

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

Once the `init` step is finished, we are ready to run the pipeline.
As mentioned in the previous subsection, the execution of the pipeline
is configured via the job table. Here we show an example of the job table,
omitting the header keywords part for cleaner look:


.. code-block:: bash

    # numid   fcomb_group   fsub_group   phot_group   mosaic_group   mask_otas               OBSID
          0             0            0            0              0      56,11   20170401T022619.1 
          1             0            0            0              0      56,53   20170401T022619.2 
          2             0            0            0              0      56,53   20170401T022619.3 
          3             0            0            0              0      56,53   20170401T022619.4 
          4             0            0            0              0      56,53   20170401T022619.5 
          5             0            0            0              0      56,21   20170401T022619.6 
          6             0            0            0              0      56,11   20170401T022619.7 
          7             0            0            0              1      56,53   20170401T033917.1 
          8             0            0            0              1      56,52   20170401T033917.2 
          9             0            0            0              1      56,11   20170401T033917.3 
         10             0            0            0              1      56,21   20170401T033917.4 
         11             0            0            0              1      56,11   20170401T033917.5 
         12             0            0            0              1      56,53   20170401T043947.1 
         13             0            0            0              1      56,53   20170401T043947.2 

Let's start with `mask_otas` column. This column is used to mask-out certain
OTAs from the input image or even the entire one (i.e., excluding the file).
By default, one thing the `init` does is to examine automatically all the
OTAs and identify which OTA is used as the guiding OTA. In addition
to that, upon visually examine the inputs, we decided to also mask out the
OTA number 56, which has more or less a skewed QE at the corner of the
cells.

The `mask_otas` support wild card. For example, if `mask_otas` == "1*"
then mask_otas  will be expends to [10, 11, 12, 13, ... 19]. A value of
empty string represent no OTA is to be masked. If all OTAs of an image
is masked, the image is considered to be remoed.

For the columns ends with `_group`, they are used as the group id for
collating the inputs in performing various operations:

    1. `fcomb_group`:
        This group is used to indicate how the images are grouped for
        creating the fringe/pupil ghost template.
        In the example above, all the values are 0 therefore only
        one template will be created (identified as `fcomb0`) using
        all the images listed as inputs.
        A value of -1 indicate that the image does not go into this
        the combination process.
    2. `fsub_group`:
        Similar to `fcomb_group`, but for specifing which image to
        subtract which template. The values in this column should
        always be present in the `fcomb_group` column.
    3. `phot_group`:
        This group is used to select images that performs the
        Mininizedr-chisq self-calibration together.
    4. `mosaic-Group` is use to specify what images are to be
           combined to craete the final stack(s).

Thanks to the help of Apus/Ruffus, which provides the check-pointing
for the pipeline execution. This means that only the files
that are affected by the changed values of the table are get re-computed
and updated, which allows the user spend more time in checking and
optimizing there selection of images.

Once we finished editing the job tablea, it is time to issue the
command `$ postcalib run <JOB_TASBLE>` to execute the pipeline.

The pipeline can be controled by a more-or-less independent, but rich
command line arguments, provided by `Ruffus`, the pipeline backend.

Here are some exampled usage:

.. code-block:: bash

    # pring out the names of all the available tasks built-in to the pipeline
    $ postcalib 201703z.txt -a -l
        +- Apus powered by Ruffus ver 2.6.3 -+
        select images
        mask OTAs
        get refcats
        merge refcats
        select fcomb
        get objmask
        mask objects
        create ftemp
        smooth ftemp
        select fsub
        subtract ftemp
        select phot
        get photcat
        cleanup photcat
        match refcat
        get flxscale
        select mosaic
        get mschdr
        create mosaic
        apply whtmap


    # create a DAG showing the topology of the tasks:
    $ postcalib 201703z.txt -a --flowchart flowchart.png

    # run the pipeline until task "get flxscale", and forced the
    # task "select phot" to re-run
    $ postcalib 201703z.txt -a --T "get flxscale" --forced_tasks 'select phot'

Further Readings
----------------

Example
^^^^^^^

For a walk-through of a real world example, check:

.. toctree::
    :maxdepth: 1

    example.rst


Pipeline
^^^^^^^^

For the details of the pipeline steps, refer to:

.. toctree::
    :maxdepth: 1

    pipeline.rst



