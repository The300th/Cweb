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
def read_Cweb(filename):
    print('o read_Cweb():')
    print('   reading',filename)
    with open(filename, "rb") as file:
        
        # read all the header information
        #---------------------------------
        fileContent = file.read(4)                        # int32_t
        one = struct.unpack('i', fileContent)[0]
        print('       one     =',one)
        
        fileContent = file.read(4)                        # int32_t
        Pweb = struct.unpack('i', fileContent)[0]
        print('       Pweb    =',Pweb)
        
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
        if (Pweb==0):
            Cweb = np.reshape(Cweb,(Nnodes,22))           # 22 values per node
        else:
            Cweb = np.reshape(Cweb,(Nnodes,35))           # 35 values per node
          
        # how to access ==> Cweb[inode,ivalue]
        
    print('   done')
    return Cweb


#==============================================================================
# Weiguang's routine to read a Cweb binary file
#==============================================================================
def readCweb(fileall, endian=None, UonGrid=None, quiet=None, selected=None):
    """
    readCweb(fileall,endian=None,UonGrid=None):
    read Cweb binary out puts
    fileall: path + filename
             last 4 charaters "Cweb" means read only one file
             else try to read multi-files by adding ".000x.Cweb" to fileall
    endian:  default '='
    UonGrid: default 35
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
        float pot;23
        float lambda1;24
        float lambda2;25
        float lambda3;26
        float local_shear[3][3];27:35
    #endif
    #ifdef UonGrid
            float u;36
    #endif
    } Cweb_t;
    """

    print('   reading',fileall)


    if endian is None:
        endian = '='
    dims = 22
    if UonGrid:
        dims += 1

    if fileall[-4:] == "Cweb":  # read only one file
        opf = open(fileall, 'rb')

        swap = struct.unpack(endian + 'i', opf.read(4))[0]
        Pweb = struct.unpack(endian + 'i', opf.read(4))[0]
        Nnodes = struct.unpack(endian + 'q', opf.read(8))[0]
        L = struct.unpack(endian + 'q', opf.read(8))[0]
        BoxSize = struct.unpack(endian + 'f', opf.read(4))[0]
        
        if Pweb==1:
            dims += 13
            
        if quiet:
            quiet = True
        else:
            print("       Pweb    = ", Pweb)
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
            Pweb = struct.unpack(endian + 'i', opf.read(4))[0]
            Nnodes = struct.unpack(endian + 'q', opf.read(8))[0]
            L = struct.unpack(endian + 'q', opf.read(8))[0]
            BoxSize = struct.unpack(endian + 'f', opf.read(4))[0]
            
            if Pweb==1:
                dims += 13
                
            if quiet:
                quiet = True
            else:
                print("       Pweb    = ", Pweb)
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


#==============================================================================
# routine to read a mesh binary file (generated by write_griddata())
#   check
#         https://docs.python.org/3/library/struct.html#format-characters
#   for format options
#==============================================================================
def read_mesh(filename,Pweb=0):
    print('o read_mesh():')
    print('   reading',filename)
    with open(filename, "rb") as file:
        
        # read all the header information
        #---------------------------------
        fileContent = file.read(4)                        # int
        Lx = struct.unpack('i', fileContent)[0]
        print('       Lx      =',Lx)
        fileContent = file.read(4)                        # int
        Ly = struct.unpack('i', fileContent)[0]
        print('       Ly      =',Ly)
        fileContent = file.read(4)                        # int
        Lz = struct.unpack('i', fileContent)[0]
        print('       Lz      =',Lz)

        fileContent = file.read(8)                        # double
        zred = struct.unpack('d', fileContent)[0]
        print('       zred    =',zred)

        fileContent = file.read(8)                        # double
        boxsize = struct.unpack('d', fileContent)[0]
        print('       boxsize =',boxsize)

        fileContent = file.read(8)                        # double
        omega0 = struct.unpack('d', fileContent)[0]
        print('       omega0  =',omega0)

        fileContent = file.read(8)                        # double
        lambda0 = struct.unpack('d', fileContent)[0]
        print('       lambda0 =',lambda0)

        fileContent = file.read(8)                        # double
        pmass = struct.unpack('d', fileContent)[0]
        print('       pmass   =',pmass)

        fileContent = file.read(8)                        # long
        npart = struct.unpack('L', fileContent)[0]
        print('       npart   =',npart)
        
        Nnodes = Lx*Ly*Lz

        # read the actual node[] matrix
        #-------------------------------
        node = np.fromfile(file, dtype=np.float32)
        if (Pweb==0):
            node = np.reshape(node,(Nnodes,4))           # 4 values per node (dens, densVx, densVy, densVz)
        else:
            node = np.reshape(node,(Nnodes,5))           # 5 values per node (+pot)
          
        # how to access ==> node[inode,ivalue]
        
    print('   done')
    return node


