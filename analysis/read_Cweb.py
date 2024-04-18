#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 10:12:17 2020

@author: aknebe
"""

import numpy as np
import sys
import struct

#==============================================================================
# routine to read a Cweb ASCII file
#==============================================================================
def read_CwebASCII(filename):
    
    print('o read_CwebASCII():')
    print('   reading',filename)
    Cweb = np.loadtxt(filename)
    print('   done')
    
    return Cweb
    
#==============================================================================
# routine to read a Cweb binary file
#   check
#         https://docs.python.org/3/library/struct.html#format-characters
#   for format options
#==============================================================================
def read_Cweb(filename, PWEB=False):
    print('o read_Cweb():')
    print('   reading',filename)
    with open(filename, "rb") as file:
        
        # read all the header information
        #---------------------------------
        fileContent = file.read(4)                        # int32_t
        one = struct.unpack('i', fileContent)[0]
        print('       one     =',one)
        
        fileContent = file.read(8)                        # uint64_t
        Nnodes = struct.unpack('L', fileContent)[0]
        print('       Nnodes  =',Nnodes)

        fileContent = file.read(8)                        # uint64_t
        L = struct.unpack('L', fileContent)[0]
        print('       L       =',L)
        
        fileContent = file.read(4)                        # float
        BoxSize = struct.unpack('f', fileContent)[0]
        print('       BoxSize =',BoxSize)
        
        # read the actual Cweb[] matrix
        #-------------------------------
        Cweb = np.fromfile(file, dtype=np.float32)
        if (PWEB==False):
            Cweb = np.reshape(Cweb,(Nnodes,22))           # 22 values per node
        else:
            Cweb = np.reshape(Cweb,(Nnodes,34))           # 34 values per node
          
        # how to access ==> Cweb[inode,ivalue]
        
    print('   done')
    return Cweb


#==============================================================================
# Weiguang's routine to read a Cweb binary file
#==============================================================================
def readCweb(fileall, endian=None, UonGrid=None, PWEB=None, quiet=None, selected=None):
    """
    readCweb(fileall,endian=None,UonGrid=None):
    read Cweb binary out puts
    fileall: path + filename
             last 4 charaters "Cweb" means read only one file
             else try to read multi-files by adding ".000x.Cweb" to fileall
    endian:  default '='
    UonGrid: default 34
    quiet  : default False
    selected:select Return data colomn.
        e.g. selected=[0,1,2] will Return x,y,z
    Return : data
    Structure of returned data:
    typedef struct {
        float x;0
        float y;1    // cell-centre [kpc/h]
        float z;2
        float dens;3 // overdensity respected to mean density
        float Vx;4
        float Vy;5   // velocity [km/sec]
        float Vz;6
        float Wx;7
        float Wy;8   // vorticity [km/sec]
        float Wz;9
        float lambda1;10
        float lambda2;11
        float lambda3;12
        float local_shear[3][3];13:22
    #ifdef PWEB
        float lambda1;22
        float lambda2;23
        float lambda3;24
        float local_shear[3][3];25:34
    #endif
    #ifdef UonGrid
            float u;34
    #endif
    } Cweb_t;
    """

    if PWEB is None:
        print('o readCweb():')
    else:
        print('o readCweb(PWEB=True):')
    print('   reading',fileall)


    if endian is None:
        endian = '='
    dims = 22
    if PWEB:
        dims += 12
    if UonGrid:
        dims += 1

    if fileall[-4:] == "Cweb":  # read only one file
        opf = open(fileall, 'rb')

        swap = struct.unpack(endian + 'i', opf.read(4))[0]
        Nnodes = struct.unpack(endian + 'q', opf.read(8))[0]
        L = struct.unpack(endian + 'q', opf.read(8))[0]
        BoxSize = struct.unpack(endian + 'f', opf.read(4))[0]
        if quiet:
            quiet = True
        else:
            print("       Nnodes  = ", Nnodes)
            print("       L       = ", L)
            print("       Boxsize = ", BoxSize)

        data = np.ndarray(shape=(Nnodes, dims), dtype="float32",
                          buffer=opf.read(Nnodes * dims * 4))
        opf.close()

    else:  # try to read multi-files
        exts = "0000"
        rnod, rfil = 0, 0
        ttnd = 10
        while (rnod < ttnd - 1):
            tmpe = exts + str(rfil)
            newf = fileall + "." + tmpe[-4:] + ".Cweb"
            try:
                opf = open(newf, 'rb')
            except IOError:
                print('       cannot open', newf)
                raise AttributeError('Please check and try again')

            swap = struct.unpack(endian + 'i', opf.read(4))[0]
            Nnodes = struct.unpack(endian + 'q', opf.read(8))[0]
            L = struct.unpack(endian + 'q', opf.read(8))[0]
            BoxSize = struct.unpack(endian + 'f', opf.read(4))[0]
            if quiet:
                quiet = True
            else:
                print("       Nnodes  = ", Nnodes)
                print("       L       = ", L)
                print("       Boxsize = ", BoxSize)

            if rnod == 0:
                ttnd = L**3
                data = np.zeros((ttnd, dims), dtype='float32')

            data[rnod:rnod + Nnodes] = \
                np.ndarray(shape=(Nnodes, dims), dtype="float32",
                           buffer=opf.read(Nnodes * dims * 4))
            rnod += Nnodes
            rfil += 1
            opf.close()

    del(swap)
    if selected is None:
        return(data)
    else:
        return(data[:, selected])
