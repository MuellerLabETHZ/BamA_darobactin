#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np

class FD(object):
    """This module creates a class for FD_Curves.
The input is the abs_path of a file (together with filename).
Important: a global variable named directory,
which contains the absolute path to the file
needs to be initialized before calling FDcurve.
The get-function: returns np-array of Vertical Deflection (VD) or
Tip-Sample_Separation from any segment of any curve:
FDcurve(Filename).get(segmentNr,'VD'/'TSS'),
alternatively segNr + segName can be used."""
    def __init__(self, drct):
        self.filepath = drct
        self.filename = drct.split('/')[-1].rstrip('.txt')
        self.segName = []
        self.segNum = {}
        self.segInfo = {}
        self.singlesegment = True
        self.data_length = 0
        
        InHeader = True
        with open(self.filepath, 'r') as f:
            for line in f:
                ### Here Values are extracted only from the first Header of the file f.
                
                # Extraction of fancyNames to self.fn
                if line.startswith('# fancyNames:') and InHeader == True:
                    self.fn = line.replace('# fancyNames: "','').replace('"\n','').split('" "')
                                           
                # Go through multiple lines and fill dict self.segInfo
                elif line.startswith('# force-settings.segment.') and InHeader == True:
                    #Splitted contains line, starting from seg-number, splitted by first '.'
                    self.singlesegment = False
                    splitted = line.rstrip().split('.',2)[2].split('.',1)
                    
                    # Extracts 'style' of segment (extend/pause/retract)
                    if splitted[-1].startswith('style:'):
                        self.segInfo[splitted[0] + 'style'] = splitted[-1].split()[-1]
                    # Extracts duration of segment
                    elif splitted[-1].startswith('duration:'):
                        self.segInfo[splitted[0] + 'duration'] = float(splitted[-1].split()[-1])
                    # Extracts number of points of segment
                    elif splitted[-1].startswith('num-points:'):
                        self.segInfo[splitted[0] + 'num-points'] = int(splitted[-1].split()[-1])
                    # Extracts z-start of segment in nm
                    elif splitted[-1].startswith('z-start:'):
                        self.segInfo[splitted[0] + 'z-start'] = float(splitted[-1].split()[-1])*10**9
                    # Extracts z-end of segment in nm
                    elif splitted[-1].startswith('z-end:'):
                        self.segInfo[splitted[0] + 'z-end'] = float(splitted[-1].split()[-1])*10**9
                
                # Extracts segmentIndex of current segment
                elif line.startswith('# segmentIndex:'):
                    segNr = line.strip()[-1]
                # Extracts segmentName of current segment and together with segNr key to empty list in segNum si created.
                elif line.startswith('# segment: '):
                    self.segName.append(segNr + line.split()[-1])
                    self.segNum[self.segName[-1]] = []
                
                # Extracts numerical data of each segment and stores it in segNum with key segName
                elif not line.startswith('#') and not line.startswith('\n'):
                    InHeader = False
                    self.data_length += 1
                    self.segNum[self.segName[-1]].append([float(i) for i in line.strip().split()])
                    
                # For unfolding curves, where only the retraction segment was stored
                if self.singlesegment == True:
                    if line.startswith('# force-settings.relative-z-start:'):
                        self.zstart = float(line.split(' ')[-1])*10**9
                    elif line.startswith('# force-settings.relative-z-end:'):
                        self.zend = float(line.split(' ')[-1])*10**9
                
        # Finds position of VD-column
        self.VD_pos = self.fn.index('Vertical Deflection')
        # Finds position of TSS-column
        if 'Vertical Tip Position' in self.fn:
            self.TSS_pos = self.fn.index('Vertical Tip Position')
        elif 'Tip-Sample Separation' in self.fn:
            self.TSS_pos = self.fn.index('Tip-Sample Separation')
        
    # This method returns VD or TSS of any segment
    def get(self,segname,fancyname):
        seg = []
        if fancyname == 'VD':
            pos = self.VD_pos
            Scaling = 10**12 # Convert N to pN
        elif fancyname == 'TSS':
            pos = self.TSS_pos
            Scaling = 10**9 # Convert m to nm
        if str(segname) not in self.segNum:
            for segment in self.segNum:
                if segment[0] == str(segname):
                    segname = segment
        for subseg in self.segNum[segname]:
            seg.append(subseg[pos]*Scaling)
        seg = np.asarray(seg)
        return seg

    # Method to return additional parameters and info extracted from Header.
    def info(self, segment, param):
        if param == 'ppnm':
            if self.singlesegment == False:
                numpoints = self.segInfo[str(segment) + 'num-points']
                z_dist = self.segInfo[str(segment) + 'z-start'] - self.segInfo[str(segment) + 'z-end']
                ppnm = numpoints / z_dist
                self.segInfo[str(segment) + 'ppnm'] = ppnm
                key = str(segment) + str(param)
                return self.segInfo[key]
            else:
                numpoints = self.data_length
                zdist = self.zstart-self.zend
                ppnm = numpoints / zdist
                return ppnm

