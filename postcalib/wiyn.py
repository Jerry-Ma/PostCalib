#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2017-08-14 15:40
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
wiyn.py
"""


from astropy.io import fits
import numpy as np


class ODILayout(object):

    def __init__(self, binning=1.0):

        self.NCX = 8
        self.NCY = 8
        self.NOX = 8
        self.NOY = 8

        self.CW = 480 / binning   # cell width
        self.CH = 494 / binning   # cell height
        self.CGW = 28 / binning   # cell gap width
        self.CGH = 11 / binning   # cell gap height

        # self.OW = self.CW * self.NCX + self.CGW * (self.NCX - 1)  # ota wd
        # self.OH = self.CH * self.NCY + self.CGH * (self.NCY - 1)  # ota ht
        self.OW = 4096 / binning  # ota width
        self.OH = 4096 / binning  # ota height
        self.OG = 200 / binning   # ota gap (width and height)
        self.ps = 0.11 * binning  # pixel scale
        self.binning = binning

    def get_xy_from_oxy(self, ox, oy, x, y):
        '''Return global x and y with given ota id and ota x and y'''
        gx = ox * (self.OW + self.OG) + x
        gy = oy * (self.OH + self.OG) + y
        return gx, gy

    def get_ota_rect(self, ox, oy):
        '''Return rect (in global x and y: l, r, b, t) of ota with given id'''
        left, bottom = self.get_xy_from_oxy(ox, oy, 0, 0)
        right, top = self.get_xy_from_oxy(ox, oy, self.OW, self.OH)
        return (left, right), (bottom, top)

    def get_cell_rect(self, ox, oy, cx, cy):
        '''Return rect of (in global x and y: l, r, b, t) of cell with
        given ota id and cell id'''
        left, bottom = self.get_xy_from_oxy(ox, oy, cx * (self.CW + self.CGW),
                                            cy * (self.CH + self.CGH))
        right = left + self.CW
        top = bottom + self.CH
        return (left, right), (bottom, top)

    def get_ota_bins(self):
        '''Return two list of tuples, for x and y direction, respectively.
        The each tuple in each list is the bound left and right global
        coordinates for that ota'''
        xbin = []
        ybin = []
        for o in range(max(self.NOX, self.NOY)):
            rect = self.get_ota_rect(o, o)
            xbin.append(rect[0])
            ybin.append(rect[1])
        return xbin[:self.NOX], ybin[:self.NOY]

    def get_cell_bins(self):
        '''Return two list of tuples, for x and y direction, respectively.
        The each tuple in each list is the bound left and right global
        coordinates for that cell'''
        xbin = []
        ybin = []
        for o in range(max(self.NOX, self.NOY)):
            _xbin = []
            _ybin = []
            for c in range(max(self.NCX, self.NCY)):
                rect = self.get_cell_rect(o, o, c, c)
                _xbin.append(rect[0])
                _ybin.append(rect[1])
            xbin += _xbin[:self.NCX]
            ybin += _ybin[:self.NCY]
        return xbin[:self.NOX], ybin[:self.NOY]


class ODILayoutInstruMixin(object):

    def get_data_extent(self):
        n_ota_x = self.ota_x_range[1] - self.ota_x_range[0]
        n_ota_y = self.ota_y_range[1] - self.ota_y_range[0]
        return {
                'n_ota_x': n_ota_x,
                'n_ota_y': n_ota_y,
                'ota_range_x': self.ota_range_x,
                'ota_range_y': self.ota_range_y,
                }

    def get_ext_ota(self, ext, return_tuple=False):
        '''return the ota xy of ext'''
        ota = self.ota_order[ext - 1]
        if return_tuple:
            return ota // 10, ota % 10
        else:
            return ota

    def get_ota_ext(self, ota):
        '''return the extension of ota'''
        return self.ota_order.index(int(ota)) + 1

    def get_sky_footprint_bbox(self, center=None):
        '''return range of ra and dec for image'''
        if center is None:
            ra = dec = 0
        else:
            ra, dec = center
        cosdec = np.cos(dec * np.pi / 180.)
        w, e, s, n = self.bbox
        return ra + w / cosdec, ra + e / cosdec, dec + s, dec + n

    def get_sky_footprint(self, center=None):
        return self.get_sky_footprint_bbox(center=center)

    def get_sky_footprint_cbox(self, center=None):
        w, e, s, n = self.get_sky_footprint_bbox(center=center)
        cra = (w + e) * 0.5
        cdec = (s + n) * 0.5
        dra = (e - w) * np.cos(cdec * np.pi / 180.)
        ddec = (n - s)
        return cra, cdec, dra, ddec

    def enumerate(self, hdulist):
        ext = 1
        while ext < self.n_ota + 1:
            yield ext, hdulist[ext]
            ext += 1


class ODILayout5ODI(ODILayout, ODILayoutInstruMixin):

    ota_order = (33, 34, 43, 44, 32,
                 23, 24, 42, 35, 53,
                 45, 54, 22, 25, 52,
                 55, 31, 13, 41, 14,
                 36, 46, 21, 12, 15,
                 51, 26, 56, 11, 16)
    ota_range_x = (1, 6)
    ota_range_y = (1, 7)
    n_ota_x = 5
    n_ota_y = 6
    n_ota = 30
    instru = '5odi'
    bbox = (-27. / 60, 19. / 60, -27. / 60, 28. / 60)

    def __init__(self, binning=1.0):
        super(ODILayout5ODI, self).__init__(binning=binning)


class ODILayoutPODI(ODILayout, ODILayoutInstruMixin):

    ota_order = (33, 34, 44,
                 43, 42, 32,
                 22, 23, 24)
    ota_range_x = (2, 5)
    ota_range_y = (2, 5)
    n_ota_x = 3
    n_ota_y = 3
    n_ota = 9
    instru = 'podi'
    bbox = (-18. / 60, 9.5 / 60, -18. / 60, 9.5 / 60)

    def __init__(self, binning=1.0):
        super(ODILayoutPODI, self).__init__(binning=binning)


def get_layout(image=None, binning=1.0, instru=None):
    if image is not None:
        if isinstance(image, fits.HDUList):
            instru = image[0].header['INSTRUME']
        else:
            with fits.open(image, memmap=True) as hdulist:
                instru = hdulist[0].header['INSTRUME']
    if instru == '5odi':
        return ODILayout5ODI(binning=binning)
    elif instru == 'podi':
        return ODILayoutPODI(binning=binning)
    else:
        raise ValueError("invalid instrument {}".format(instru))
