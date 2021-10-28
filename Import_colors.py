#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 15:33:04 2021

@author: rnoah
"""


def import_rgb_colors(file_dir):
    color_scheme = {}
    name_scheme = []
    divider = 255 # From Illustrator -> reformat as 0 <= val <= 1
    with open(file_dir, 'r') as f:
        for index, line in enumerate(f):
            if index > 0:
                seg = int(line.split()[0])
                name = line.split()[1]
                r = int(line.split()[2])/divider
                g = int(line.split()[3])/divider
                b = int(line.split()[4])/divider
        
                color_scheme[seg] = [r, g, b]
                name_scheme.append(name)
    return color_scheme, name_scheme