#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2017-08-19 15:47
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
prep_grouping.py

This file contains routines the handles the grouping of
input images from the job table.

With the help of these routines, the users are able to
define the inputs of the various steps of calibrations.

For example, fcomb_group determines what images goes into
sky image construction and subtraction of additive features
such as fringe and pupil ghost.

The value in the groups are explained as follows:

    fcomb_group: non-negative int or -1
        The value serves as the id of the group. All images
        that share the same group id are gathered together to
        create a combined composite template
        (see Doc02 for detailed discussions).
        By default, the images are grouped based on filter, i.e,
        images with the same filter have the same fcomb_group.
        A special value -1 indicates that the entry does not
        go into the consequent processes.
        If all the images have fcomb_group=-1 then no template
        will be created.

    fsub_group: non-negative int or -1
        The value indicate which composite template will be used
        for subtracting the fringe/pupilghost.
        The values here should be a valid group id of fcomb_group.
        Otherwise an error will occur.
        A special value -1 indicates that no subtraction will be
        done for that entry.

    phot_group: non-negative int or -1
        The photometry calibration group. By default, this is
        set to be grouped with respect to the filter. A phot
        group of value -1 indicate that this image does not
        go into the self-calibration fitting, hence no FLXSCALE
        header will be created. Such image will not be able to
        get into the mosaic step.

    mosaic_group: non-negative int or -1
        The value indicate how the images in the job file are
        grouped to create mosaic. If an image have phot_group = -1,
        it will force the mosaic_group = -1.
        If mosaic_group is set to -1, the image is excluded from
        creating mosaics.

Inputs
------
jobfile: ASCII table
    Contains the relevant columns described above.
jobdir: path
    The path to the jobdir to locate the files

Outputs
-------
checkfile: ASCII table
    Save content as jobfile, serves as the checker for diff-ed updating.
"""


import os
import re
import sys
import glob
from astropy.table import Table
# from postcalib.wiyn import get_layout
from postcalib.apus.common import get_log_func, touch_file
# from postcalib import qa


def main(*args, **kwargs):
    if not args:
        args = sys.argv[1:]
    log = get_log_func(default_level='debug', **kwargs)
    jobfile, jobdir, checkfile = args
    grpkey = kwargs['grpkey']
    grpcol = grpkey + '_group'

    job_table = Table.read(jobfile, format='ascii.commented_header')
    job_table.pprint(max_width=None, max_lines=-1)

    # cases:
    # no checker_table
    #   purge all links first {prefix}_*_{obsid}_*
    #   entries in job_table, add if >= 0
    # has checker_table
    #   entry is in checker, and group the same, skip
    #   otherwise purge {prefix}_*_{obsid}_*
    #   add if >= 0

    if not os.path.exists(checkfile):
        checker_table = []
    else:
        checker_table = Table.read(
                checkfile, format='ascii.commented_header')

    purge_glob = os.path.join(jobdir, grpkey + "[0-9]*_{obsid}_*")
    # glob the images
    add_count = skip_count = rm_count = 0
    images = glob.glob(os.path.join(jobdir, kwargs['sel_inputs']))
    # process the images to merge with fallbacks
    fb_images = []
    for entry in job_table:
        image = get_image(images, entry['OBSID'])
        if image is None:
            log("no images available for {}".format(entry['OBSID']))
            if 'fallbacks' in kwargs:
                for fkey, fallback in kwargs['fallbacks']:
                    fb_image = get_image(
                            glob.glob(os.path.join(jobdir, fallback)),
                            entry['OBSID'])
                    if fb_image is None:
                        continue
                    else:
                        log("fallback to {}".format(fb_image))
                        fb_images.append(fb_image)
                        break
        else:
            pass
    # if len(images) == 0:
    #     log("no images available for selection, try using fallbacks")
    #     if 'fallbacks' in kwargs:
    #         for fkey, fallback in kwargs['fallbacks']:
    #             fb_images = glob.glob(os.path.join(jobdir, fallback))
    #             if len(fb_images) == 0:
    #                 continue
    #             else:
    #                 log("fallback to {}".format(fkey))
    #                 images = fb_images
    #                 break
    rm_grp = set()
    add_grp = set()
    skipped_grp = {}
    for inname in images + fb_images:
        # parse image name
        parsed_filename = re.match(
                kwargs['reg_inputs'], os.path.basename(inname)
                ).groupdict()
        entry = get_entry(job_table, OBSID=parsed_filename['obsid'])
        ppflag = "{}{}_".format(grpkey, entry[grpcol])
        outname = os.path.join(jobdir, kwargs['fmt_selected'].format(
                **dict(parsed_filename, ppflag=ppflag)
                ))
        # purge all but outname
        skip_count_flag = True
        rm_count_flag = True
        for f in glob.glob(purge_glob.format(**parsed_filename)):
            bf = os.path.basename(f)
            if entry[grpcol] >= 0 and os.path.exists(outname) and \
                    bf.startswith(ppflag):
                if entry[grpcol] not in skipped_grp.keys():
                    skipped_grp[entry[grpcol]] = inname
                log("existing file {} for link {}".format(f, outname))
                if skip_count_flag:
                    skip_count += 1
                    skip_count_flag = False
                continue
            else:
                # query the checkfile for the old grpid
                old_entry = [e for e in checker_table
                             if e['OBSID'] == parsed_filename['obsid']]
                if not old_entry:
                    log("warning", "something wrong with the checker file")
                rm_grp.add(old_entry[0][grpcol])
                log("purge obsolete file {}".format(f))
                os.remove(f)
                if rm_count_flag:
                    rm_count += 1
                    rm_count_flag = False
                # also purge no id product
                for f in glob.glob(os.path.join(
                        jobdir, grpkey + "_{obsid}_*".format(
                            **parsed_filename))):
                    log("purge obsolete file {}".format(f))
                    os.remove(f)
        # create link of >= 0
        # print(inname, outname)
        if not os.path.exists(outname) and entry[grpcol] >= 0:
            add_grp.add(entry[grpcol])
            log("link for {} to {}".format(grpkey, outname))
            os.symlink(os.path.basename(inname), outname)
            add_count += 1
    # for grp that has rm but no add, touch one file to trigger update
    touch_grp = rm_grp - add_grp
    print(rm_grp, add_grp, touch_grp)
    for grpid in touch_grp:
        if grpid in skipped_grp.keys():
            touch_file(skipped_grp[grpid])
            log("touch group {} of {} for triggering update".format(
                grpid, skipped_grp[grpid]))

    log("{} entries + {}, - {}, skipped {}".format(
        grpkey, add_count, rm_count, skip_count))
    job_table.write(checkfile, format='ascii.commented_header')


def get_entry(tbl, **kwargs):
    k, v = list(kwargs.items())[0]
    if tbl is None:
        return None
    else:
        entry = [e for e in tbl if e[k] == v]
        if len(entry) > 0:
            return entry[0]
        else:
            return None


def get_image(images, obsid):
    image = [i for i in images if obsid in os.path.basename(i)]
    if len(image) > 0:
        return image[0]
    else:
        return None
