#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 11:42:08 2020

@author: rnoah
"""
import numpy as np
import pywt


def denoising(y_BLC, ppnm, iterations=3, thr=0.1, level_factor=4, wave = 'db4', thr_mode = 'soft'):
    """
    Data-Input: Baseline-corrected y-data, and points per nm (ppnm)
    User defined parameters:
        iterations: Levels of decomposition
        thr = multiplier of denoising threshold - is multiplied with ppnm.
            The applied threshold is scaled on each level by level_factor.
        level_factor = scaling of denoising threshold per level
        wave = type of wavelet
        thr_mode = choose soft or hard threshold
        for more details check pywt documentation
    """
    coeffs = pywt.wavedec(y_BLC, wave, level = iterations)
    level_thr = thr*ppnm
    for i in range(iterations):
        coeffs[i+1] = pywt.threshold(coeffs[i+1], np.std(coeffs[i+1])*level_thr, thr_mode)
        level_thr = level_thr * level_factor
    return pywt.waverec(coeffs, wave)