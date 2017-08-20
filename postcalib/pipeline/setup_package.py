# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import


def get_package_data():
    return {"postcalib.pipeline": [
        'bpm/odi_5x6/*.reg',
        'bpm/podi/*.reg',
        'bpm/pupilmask/*.reg'
        ]}
