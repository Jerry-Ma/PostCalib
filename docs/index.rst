.. doctest-skip-all


===========
PostCalib
===========


This is the documentation for the **PostCalib** package, a pipeline for
creating science-ready mosaics from WIYN ODI images.

Code and issue tracker are on `GitHub <https://github.com/Jerry-Ma/PostCalib>`_.


Introduction
==============


PostCalib is a pipeline for creating science-ready mosaics from a collection of
WIYN ODI images.

This tool works as the next step (post-processing) of the WIYN ODI data
reduction pipeline `QuickReduce <https://github.com/WIYN-ODI/QuickReduce>`_ .
PostCalib typically takes flat-fielded and WCS-corrected images produced
from QR, performs various collection-based operations such as
fringe/pupil ghost removal, photometry calibration, and image stacking.

PostCalib is written on top of `Apus <https://github.com/Jerry-Ma/apus>`__
and facilitates a full suite of command line switches/arguments to help
the users control how the pipeline is run.


Installation
==============


The latest/development version of PostCalib can either be download directory
from Github or pip installed.


Using pip
----------


.. code-block:: bash

   $ pip install https://github.com/Jerry-Ma/PostCalib/archive/master.zip



Building from source
---------------------



.. code-block:: bash

    $ # If you have a github account:
    $ git clone git@github.com:Jerry-Ma/PostCalib.git
    $ # If you do not:
    $ git clone https://github.com/Jerry-Ma/PostCalib.git
    $ cd PostCalib
    $ python setup.py install



Requirements
----------------


PostCalib works with Python 2.7 and 3.4 later

The following packages are required for PostCalib installation & use:

* `Cython`
* `numpy`
* `scipy`
* `matplotlib`
* `cycler`
* `astropy`
* `PyYAML`
* `astroquery`
* `pyregion`
* `shutilwhich`
* `ruffus`
* `lmfit`

The following non-python softwares are required as well:

* SExtractor, SCAMP, and SWarp
* stilts
* fpack/funpack

See :ref:`setup-details` for more details.


Get Started
===============


Tutorial
----------

The primary way of running PostCalib is through command line.
The recommended way of getting started is to follow the tutorial:

.. toctree::
    :maxdepth: 1

    postcalib/tutorial.rst


Example
-----------

For a real-world example, check:

.. toctree::
    :maxdepth: 1

    postcalib/example.rst


Pipeline
-------------
For more on the details of the pipeline steps, check the pipeline document:

.. toctree::
    :maxdepth: 2

    postcalib/pipeline.rst


Reference/API
===============


For a comprehensive listing of the command line options, see:

.. toctree::
    :maxdepth: 1

    postcalib/commandline.rst

For the documentation of PostCalib, check:

.. toctree::
    :maxdepth: 2

    postcalib/index.rst


License
========


PostCalib is licensed under a 3-clause BSD style license - see :doc:`license`.
