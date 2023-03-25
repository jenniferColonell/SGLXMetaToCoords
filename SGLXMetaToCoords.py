# -*- coding: utf-8 -*-
"""
Requires python 3

The main() function at the bottom of this file can run from an
interpreter, or, the helper functions can be imported into a
new module or Jupyter notebook.

Standalone program to generate a coordinate file from a SpikeGLX
metadata file for 3A, NP1.0, or NP2 single shank or multishank probes.

The output is set by the outType parameter:
                  0 for text coordinate file; 
                  1 for Kilosort or Kilosort2 channel map file;
                  2 for strings to paste into JRClust .prm file 
                


@author: Jennifer Colonell, Janelia Research Campus

"""

import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from pathlib import Path
from tkinter import Tk
from tkinter import filedialog
import shutil


# Parse ini file returning a dictionary whose keys are the metadata
# left-hand-side-tags, and values are string versions of the right-hand-side
# metadata values. We remove any leading '~' characters in the tags to match
# the MATLAB version of readMeta.
#
# The string values are converted to numbers using the "int" and "float"
# fucntions. Note that python 3 has no size limit for integers.
#
def readMeta(metaPath):
    metaDict = {}
    if metaPath.exists():
        # print("meta file present")
        with metaPath.open() as f:
            mdatList = f.read().splitlines()
            # convert the list entries into key value pairs
            for m in mdatList:
                csList = m.split(sep='=')
                if csList[0][0] == '~':
                    currKey = csList[0][1:len(csList[0])]
                else:
                    currKey = csList[0]
                metaDict.update({currKey: csList[1]})
    else:
        print("no meta file")
        
    return(metaDict)
    
   
# Return counts of each imec channel type that composes the timepoints
# stored in the binary files.
#
def ChannelCountsIM(meta):
    chanCountList = meta['snsApLfSy'].split(sep=',')
    AP = int(chanCountList[0])
    LF = int(chanCountList[1])
    SY = int(chanCountList[2])
    
    return(AP, LF, SY)


# =========================================================
# Return geometry paramters for supported probe types
# These are used to calculate positions from metadata
# that includes only ~snsShankMap
#
def getGeomParams(meta):
# many part numbers have the same geometry parameters ;
# define those sets in lists
# [nShank, shankWidth, shankPitch, even_xOff, odd_xOff, horizPitch, vertPitch, rowsPerShank, elecPerShank]
# offset and pitch values in um
    np1_stag_70um       = [1,  70,   0,   27,   11,  32,  20,  480,  960]
    nhp_lin_70um        = [1,  70,   0,   27,   27,  32,  20,  480,  960]
    nhp_stag_125um_med  = [1, 125,   0,   27,   11,  87,  20, 1368, 2496]
    nhp_stag_125um_long = [1, 125,   0,   27,   11,  87,  20, 2208, 4416]
    nhp_lin_125um_med   = [1, 125,   0,   11,   11, 103,  20, 1368, 2496]
    nhp_lin_125um_long  = [1, 125,   0,   11,   11, 103,  20, 2208, 4416]
    uhd_8col_1bank      = [1,  70,   0,   14,   14,   6,   6,   48,  384]
    uhd_8col_16bank     = [1,  70,   0,   14,   14,   6,   6,  768, 6144]
    np2_ss              = [1,  70,   0,   27,   27,  32,  15,  640, 1280]
    np2_4s              = [4,  70, 250,   27,   27,  32,  15,  640, 1280]
    NP1120              = [1,  70,   0, 6.75, 6.75, 4.5, 4.5,  192,  384]
    NP1121              = [1,  70,   0, 6.25, 6.25,   3,   3,  384,  384]
    NP1122              = [1,  70,   0, 6.75, 6.75, 4.5, 4.5,   24,  384]
    NP1123              = [1,  70,   0,10.25,10.25,   3,   3,   32,  384]
    NP1300              = [1,  70,   0,   11,   11,  48,  20,  480,  960]
    NP1200              = [1,  70,   0,   27,   11,  32,  20,   64,  128]
    NXT3000             = [1,  70,   0,   53,   53,   0,  15,  128,  128]
    
    M = dict([
    ('3A',np1_stag_70um),
    ('PRB_1_4_0480_1',np1_stag_70um),
    ('PRB_1_4_0480_1_C',np1_stag_70um),
    ('NP1010',np1_stag_70um), 
    ('NP1011',np1_stag_70um),
    ('NP1012',np1_stag_70um),
    ('NP1013',np1_stag_70um),
    
    ('NP1015',nhp_lin_70um),
    ('NP1015',nhp_lin_70um),
    ('NP1016',nhp_lin_70um),
    ('NP1017',nhp_lin_70um),
       
    ('NP1020',nhp_stag_125um_med),
    ('NP1021',nhp_stag_125um_med),
    ('NP1030',nhp_stag_125um_long),
    ('NP1031',nhp_stag_125um_long),
    
    ('NP1022',nhp_lin_125um_med),
    ('NP1032',nhp_lin_125um_long),
    
    ('NP1100',uhd_8col_1bank),
    ('NP1110',uhd_8col_16bank),
    
    ('PRB2_1_4_0480_1',np2_ss),
    ('PRB2_1_2_0640_0',np2_ss),
    ('NP2000',np2_ss),
    ('NP2003',np2_ss),
    ('NP2004',np2_ss),
    
    ('PRB2_4_2_0640_0',np2_4s),
    ('NP2010',np2_4s),
    ('NP2013',np2_4s),
    ('NP2014',np2_4s),
    
    ('NP1120',NP1120),
    ('NP1121',NP1121),
    ('NP1122',NP1122),
    ('NP1123',NP1123),
    ('NP1300',NP1300),
    
    ('NP1200',NP1200),
    ('NXT3000',NXT3000)
    ])
    
    # get probe part number; if absent, this is a 3A
    if 'imDatPrb_pn' in meta:   
        pn = meta['imDatPrb_pn']
    else:
        pn = '3A';


    if pn in M:
        geomList = M[pn]
    else:
        print('unsupported probe part number\n')
        geomList = []
    
    return geomList


# =========================================================
# Parse snsGeomMap for XY coordinates
#
def geomMapToGeom(meta):

    # read in the shank map   
    geomMap = meta['snsGeomMap'].split(sep=')')
    
    # there is an entry in the map for each saved channel
    # number of entries in map -- subtract 1 for header, one for trailing ')'
    nEntry = len(geomMap) - 2

    shankInd = np.zeros((nEntry,))
    xCoord = np.zeros((nEntry,))
    yCoord = np.zeros((nEntry,))
    connected = np.zeros((nEntry,))
    
    for i in range(nEntry):
        # get parameter list from this entry, skipping first header entry
        currEntry = geomMap[i+1]
        currEntry = currEntry[1:len(currEntry)]
        currList = currEntry.split(sep=':')
        shankInd[i] = int(currList[0])
        xCoord[i] = float(currList[1])
        yCoord[i] = float(currList[2])
        connected[i] = int(currList[3])

    # parse header for number of shanks
    currList = geomMap[0].split(',')
    nShank = int(currList[1]);
    shankWidth = float(currList[2]);
    shankPitch = float(currList[3]);
    
    return nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected


# =========================================================
# Build snsGeomMap from xy coordinates
#
def snsGeom(meta, shankInd, xCoord, yCoord, use):
    # header
    pn = meta['imDatPrb_pn']
    geomList = getGeomParams(meta)
    # geomList =
    # [nShank, shankWidth, shankPitch, even_xOff, odd_xOff, horizPitch, vertPitch, rowsPerShank, elecPerShank]
    
    snsGeomStr = '~snsGeomMap=(' + pn
    snsGeomStr = snsGeomStr + ',{:d},{:g},{:g})'.format(geomList[0], geomList[2], geomList[1])
    nEntry = shankInd.shape[0]
    
    for i in range(0,nEntry):
        snsGeomStr = snsGeomStr + '({:g}:{:g}:{:g}:{:g})'.format(shankInd[i], xCoord[i], yCoord[i], use[i])
    
    return snsGeomStr


# =========================================================
# Get XY coordinates from snsShankMap plus hard coded geom values
#
def shankMapToGeom(meta):
    # get number of saved AP channels (some early metadata files have a
    # SYNC entry in the snsChanMap)
    AP, LF, SY = ChannelCountsIM(meta)

    shankMap = meta['snsShankMap'].split(sep=')')
    
    shankInd = np.zeros((AP,))
    colInd = np.zeros((AP,))
    rowInd = np.zeros((AP,))
    connected = np.zeros((AP,))
    xCoord = np.zeros((AP,))
    yCoord = np.zeros((AP,))
    
    for i in range(AP):
        # get parameter list from this entry, skipping first header entry
        currEntry = shankMap[i+1]
        currEntry = currEntry[1:len(currEntry)]
        currList = currEntry.split(sep=':')
        shankInd[i] = int(currList[0])
        colInd[i] = int(currList[1])
        rowInd[i] = int(currList[2])
        connected[i] = int(currList[3])
   
    geomList = getGeomParams(meta);
    # geomList = 
    # [nShank, shankWidth, shankPitch, even_xOff, odd_xOff, horizPitch, vertPitch, rowsPerShank, elecPerShank]
    
    oddRows = np.bool_(rowInd%2);
    evenRows = ~oddRows;
    xCoord = colInd*float(geomList[5]);
    xCoord[evenRows] = xCoord[evenRows] + geomList[3]
    xCoord[oddRows] = xCoord[oddRows] + geomList[4]
    yCoord = rowInd*float(geomList[6])
    
    nShank = geomList[0]
    shankWidth = geomList[1]
    shankPitch = geomList[2]
    
    return nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected


# =========================================================
# Plot x z positions of all electrodes and saved channels
#
def plotSaved(xCoord, yCoord, shankInd, meta):
    
    geomList = getGeomParams(meta)
    # geomList = 
    # [nShank, shankWidth, shankPitch, even_xOff, odd_xOff, horizPitch, vertPitch, rowsPerShank, elecPerShank]
    
    # calculate positions on one shank
    nCol = geomList[8]/geomList[7]
    rowInd = np.arange(geomList[8])
    rowInd = np.floor(rowInd/nCol)
    oddRows = np.bool_(rowInd%2);
    evenRows = ~oddRows;
    
    colInd = np.arange(geomList[8])
    colInd = (colInd % nCol)
    
    xall = colInd*geomList[5]
    xall[evenRows] = xall[evenRows] + geomList[3]
    xall[oddRows] = xall[oddRows] + geomList[4]
    
    yall = rowInd*geomList[6];
    
    
    fig = plt.figure(figsize=(2,12))
    
    shankSep = geomList[2]
    
    # loop over shanks
    for sI in range(geomList[0]):
        
        # plot all positions   
        marker_style = dict(c='w', edgecolor = 'k', linestyle='None', marker='s', s=5)                 
        plt.scatter(shankSep*sI + xall, yall, **marker_style)
        
        # plot selected positions
        currInd = np.argwhere(shankInd == sI)
        marker_style = dict(c='b', edgecolor = 'g', linestyle='None', marker='s', s=15) 
        plt.scatter(shankSep*sI + xCoord[currInd], yCoord[currInd], **marker_style)
   
   # after looping over all shanks, show plot 
    plt.show()

    return


def CoordsToText(meta, chans, xCoord, yCoord, connected, shankInd, shankSep, baseName, savePath, buildPath ):

    if buildPath:
        newName = baseName + '_siteCoords.txt'
        saveFullPath = Path(savePath / newName)
    else:
        saveFullPath = savePath

    # Note that the channel index written is the index of that channel in the saved file

    with open(saveFullPath, 'w') as outFile:
        for i in range(0,chans.size):
            currX = shankInd[i]*shankSep + xCoord[i]
            currLine = '{:d}\t{:g}\t{:g}\t{:g}\n'.format(i, currX, yCoord[i], shankInd[i])
            outFile.write(currLine)


def CoordsToNPY(meta, chans, xCoord, yCoord, connected, shankInd, shankSep, baseName, savePath, buildPath ):

    if buildPath:
        newName = baseName + '_siteCoords.npy'
        saveFullPath = Path(savePath / newName)
    else:
        saveFullPath = savePath

    # write an npy file of nChanx2 
    nchan = xCoord.shape[0]
    geom = np.zeros((nchan,2))
    geom[:,0] = xCoord + shankInd*shankSep
    geom[:,1] = yCoord
    
    np.save(saveFullPath, geom)

            
def CoordsToJRCString(meta, chans, xCoord, yCoord, connected, shankInd, shankSep, baseName, savePath, buildPath ):
    
    if buildPath:
        newName = baseName +'_forJRCprm.txt'
        saveFullPath = Path(savePath / newName)
    else:
        saveFullPath = savePath
    
    # siteMap, equivalent of chanMap in KS, is the order of channels in the saved file
    # rather than original channel indicies.
    nChan = chans.size
    siteMap = np.arange(0,nChan, dtype = 'int')
    siteMap = siteMap + 1   # convert to 1-based for MATLAB
    
    shankInd = shankInd + 1 # conver to 1-based for MATLAB
 
    shankStr = 'shankMap = ['
    coordStr = 'siteLoc = ['
    siteMapStr = 'siteMap = ['
    
    xCoord = shankInd*shankSep + xCoord
    
    for i in range(0,chans.size-1):
        shankStr = shankStr + '{:g},'.format(shankInd[i])   # convert to 1-based for MATLAB
        coordStr = coordStr + '{:g},{:g};'.format(xCoord[i], yCoord[i])
        siteMapStr = siteMapStr + '{:d},'.format(siteMap[i])    # convert to 1-based for MATLAB
    
    # final entries
    shankStr = shankStr + '{:g}];\n'.format(shankInd[nChan-1])
    coordStr = coordStr + '{:g},{:g}];\n'.format(xCoord[nChan-1], yCoord[nChan-1])
    siteMapStr = siteMapStr + '{:d}];\n'.format(siteMap[nChan-1])
    
    with open(saveFullPath, 'w') as outFile:
        outFile.write(shankStr) 
        outFile.write(coordStr)
        outFile.write(siteMapStr)


def CoordsToKSChanMap(meta, chans, xCoord, yCoord, connected, shankInd, shankSep, baseName, savePath, buildPath ): 
    
    if buildPath:
        newName = baseName +'_kilosortChanMap.mat'
        saveFullPath = Path(savePath / newName)
    else:
        saveFullPath = savePath
    
    nChan = chans.size
    # channel map is the order of channels in the file, rather than the 
    # original indicies of the channels
    chanMap0ind = np.arange(0,nChan,dtype='float64')
    chanMap0ind = chanMap0ind.reshape((nChan,1))
    chanMap = chanMap0ind + 1
    
    connected = (connected==1)
    connected = connected.reshape((nChan,1))
    
    xCoord = shankInd*shankSep + xCoord
    xCoord = xCoord.reshape((nChan,1))
    yCoord = yCoord.reshape((nChan,1))
    
    kcoords = shankInd + 1
    kcoords = kcoords.reshape((nChan,1))
    kcoords = kcoords.astype('float64')
    
    name = baseName
    
    mdict = {
            'chanMap':chanMap,
            'chanMap0ind':chanMap0ind,
            'connected':connected,
            'name':name,
            'xcoords':xCoord,
            'ycoords':yCoord,
            'kcoords':kcoords,
            }
    scipy.io.savemat(saveFullPath, mdict)

    
def CoordsToGeomMap(meta, chans, xCoord, yCoord, connected, shankInd, shankSep, baseName, savePath, buildPath ):

    if buildPath:
        newName = baseName + '_plusGeom.meta'
        saveFullPath = Path(savePath / newName)
    else:
        print('Can only make new metadata in same directory as original.')
        return
    
    origPath = Path(savePath / (baseName + '.meta'))
    shutil.copy(origPath, saveFullPath)
    
    snsGeomStr = snsGeom(meta, shankInd, xCoord, yCoord, connected)
    
    with open(saveFullPath, 'a') as outFile:
        outFile.write(snsGeomStr)


# Given a path to a SpikeGLX metadata file, write out coordinates
# in formats for analysis software to consume
# Input params:
#   metaFullPath: full path, including the file name
#   outType:  format for the output
#   badChan:  channels other than reference channels to exclude
#   destFullPath: 
    
def MetaToCoords(metaFullPath, outType, badChan= np.zeros((0), dtype = 'int'), destFullPath = '', showPlot=False):
       
    # Read in metadata; returns a dictionary with string for values
    meta = readMeta(metaFullPath)
    
    # Get coordinates for saved channels from snsGeomMap, if present,
    # otherwise from snsShankMap
    if 'snsGeomMap' in meta:
        [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = geomMapToGeom(meta);
    else:   
        [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = shankMapToGeom(meta);
    
    if showPlot:
        plotSaved(xCoord, yCoord, shankInd, meta)
        
    [AP,LF,SY] = ChannelCountsIM(meta)    
    chans = np.arange(AP)
    
    # Channels identified as noisy by kilosort helper indexed
    # according to position in the file
    # since these can include the SYNC channel, remove any from
    # list that are outside the range of AP channels
    badChan = badChan[badChan < AP]
    connected[badChan] = 0    

    baseName = metaFullPath.stem
    
    if outType >= 0:
        if len(destFullPath) == 0:
            savePath = metaFullPath.parent
            buildPath = True
        else:
            buildPath = False
            savePath = destFullPath
        outputSwitch = {
                0: CoordsToText,
                1: CoordsToKSChanMap,
                2: CoordsToJRCString,
                3: CoordsToGeomMap,
                4: CoordsToNPY
        }
        
        writeFunc = outputSwitch.get(outType)
        writeFunc(meta, chans, xCoord, yCoord, connected, shankInd, shankPitch, baseName, savePath, buildPath )
    
    return xCoord, yCoord, shankInd
    
# Sample calling program to get a metadata file from the user,
# output a file set by outType
#   0 = tab delimited text file of coordinates in um, index, x, y, shank index
#   1 = KS2 chan map .mat file
#   2 = strings of channel map, shank index, and x,y pairs for JRClust
#   3 = create a new metadata file with snsGeomMap appended (for converting old metadata to new)
#   4 = npy file with (nchan,2) matrix of xy coordinates, e.g. for YASS and related sorters
# file is saved to the same path as the metadata file.
#
def main():
    
    outType = 1
    
    # Get file from user
    root = Tk()         # create the Tkinter widget
    root.withdraw()     # hide the Tkinter root window

    # Windows specific; forces the window to appear in front
    root.attributes("-topmost", True)

    metaFullPath = Path(filedialog.askopenfilename(title="Select meta file"))
    root.destroy()      # destroy the Tkinter widget
   
    MetaToCoords( metaFullPath=metaFullPath, outType=outType, showPlot=True  )

        

if __name__ == "__main__":
    main()