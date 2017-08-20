#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2017-08-17 21:12
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
test_prep_mask_otas.py
"""


def test_parse_mask_otas():
    from ..prep_mask_otas import parse_mask_otas
    assert parse_mask_otas("*") == [
            i + j for i in map(str, range(10)) for j in map(str, range(10))]
