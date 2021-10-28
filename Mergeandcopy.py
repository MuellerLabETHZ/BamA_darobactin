#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 15:28:44 2020

@author: rnoah
"""
import os
from shutil import copy


def mergeandcopy(FD_Curves, FD_Curves_prev, prev_files_dir, directory, saveselected):
    if saveselected == True:
        srcdir = directory
        destdir = os.path.join(directory, 'Selected_Curves')
        os.makedirs(destdir, exist_ok=True)

        # Copy selected curves from previous set to destination dir
        prev_filenames = os.listdir(prev_files_dir)
        for filename in prev_filenames:
            filedir = os.path.join(prev_files_dir, filename)
            copy(filedir, destdir)

    # Copy selected curves of new set to destination dir
    # Rename files if same name exists in prev set (adapt filename and entry in FD_Curves)
    #1. Make a list of prev names (only selected curves)
    prev_names = []
    for i in range(len(FD_Curves_prev)):
        if FD_Curves_prev[i]['Selected'] == True:
            prev_names.append(FD_Curves_prev[i]['name'])
    #2. Check if name is in prev name list -> add number if present
    for j in range(len(FD_Curves)):
        if FD_Curves[j]['name'] in prev_names:
            filenamecounter = 2
            init_name = FD_Curves[j]['name']
            while FD_Curves[j]['name'] in prev_names:
                new_name = init_name+'_'+str(filenamecounter)
                FD_Curves[j]['name'] = new_name
                filenamecounter += 1
            #3. Copy file with new name if selected
            if saveselected == True:
                if FD_Curves[j]['Selected'] == True:
                    # 3.1. rename original file
                    init_path = os.path.join(srcdir, init_name + '.txt')
                    new_path = os.path.join(srcdir, new_name + '.txt')
                    os.rename(init_path, new_path)
                    # 3.2. copy renamed file to destdir
                    copy(new_path, destdir)
                    # 3.3. rename original files back to inital names
                    os.rename(new_path, init_path)
        # If filename does not exist yet, just copy to destdir
        else:
            if saveselected == True:
                if FD_Curves[j]['Selected'] == True:
                    filename = FD_Curves[j]['name'] + '.txt'
                    filedir = os.path.join(srcdir, filename)
                    copy(filedir, destdir)
    # Combine FD_curves_prev and FD_Curves (with updated filenames if necessary)
    FD_Curves = FD_Curves_prev + FD_Curves
    return FD_Curves