#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 15:08:48 2023

@author: qhong7
"""
import re
import numpy as np
#from matplotlib import pyplot as plt
import spglib

def line2vec(s):
    list_s = re.split('([- ])', s)
    #list_s = re.split('(\W)', s)
    counter = 0
    vec = []
    neg = 1.
    for c in list_s:
        if len(c) > 1:
            vec.append(neg * float(c))
            neg = 1.
            counter += 1
            if counter == 3:  return vec
        elif c == '-':
            neg = -1.
            
def pos_match(pos0,pos1,latt):
    for s in range(pos1.shape[0]):
        dist_min = np.linalg.norm(pos1[s] - pos0[s])
        p_min = pos1[s]
        if dist_min < 1.0:
            pass
        else:
            for i in range(-1,2):
                for j in range(-1,2):
                    for k in range(-1,2):
                        p = pos1[s] + np.matmul(np.asarray([i,j,k]), latt)
                        dist = np.linalg.norm(p - pos0[s])
                        if dist < dist_min:
                            p_min = p
                            dist_min = dist
        pos1[s] = p_min
    return pos1

def get_pos(strct_name):
    #flag_simple = False
    #try: 
    #    file = open('../diff/simple_OUTCAR_collect_'+strct_name, 'r')
    #    flag_simple = True
    #except: 
    #    file = open('../diff/OUTCAR_collect_'+strct_name, 'r')
    #    file_s = open('../diff/simple_OUTCAR_collect_'+strct_name, 'w')
    file = open('OUTCAR_collect'+strct_name, 'r')
    
    i = 0
    flag_lattice, flag_position = False, False
    counter = 0
    dist_his = []
    time = 0
    time_his = []
    position_his = []
    lattice_his = []
    natom = 0
    print('Reading OUTCAR......')
    while i==0 or Line:# and i<10000:
        i += 1
        Line = file.readline()
        
        if flag_lattice:
            lattice.append(line2vec(Line))
            remain_lattice -= 1
            if remain_lattice == 0: 
                flag_lattice = False
                lattice = np.asarray(lattice)
                lattice_his.append(lattice)
        
        if flag_position:
            if remain_position<natom+1: position.append(line2vec(Line))
            remain_position -= 1
            if remain_position == 0: 
                flag_position = False
                position = np.asarray(position)
                if counter == 0:
                    position0 = position
                else:
                    # match this position with last one
                    position = pos_match(position_last, position, lattice)
                    # calculate distance
                    dist_mat = position - position0
                    dist = []
                    for v in dist_mat: dist.append(np.linalg.norm(v))
                    dist_his.append(dist)
                    time += stepsize
                    time_his.append(time)
                counter += 1
                if counter%4000 == 0: print("%d MD steps read in" %counter)
                position_last = position
                position_his.append(position)

        if Line.startswith( '      direct lattice vectors'): 
            flag_lattice = True
            remain_lattice = 3
            lattice = []
        if Line.startswith( ' POSITION'):  
            flag_position = True
            remain_position = natom+1  # 120 80
            position = []
        if Line.startswith( '   POTIM  ='):                                     
            stepsize = float(Line.split(' ')[6])
        if Line.startswith( '   number of dos'):            
            natom = int(Line.split('NIONS =    ')[1])
        #else:
        #    continue
        #if not flag_simple: file_s.write(Line)
        
            
    file.close()
    #if not flag_simple: file_s.close()
    position_his = np.asarray(position_his)
    lattice_his = np.asarray(lattice_his)
    
    ##%%
    #print()
    #ll = position_his.shape[0]
    #position_tmp = position_his[ll//2:ll,:,:]
    #position_tmp = position_his
    pos = np.mean(position_his[position_his.shape[0]//2:,:,:],axis=0)
    #print(pos)
    #print(position_his[position_his.shape[0]//2:ll,:,:].shape)
    #lattice_his = lattice_his[lattice_his.shape[0]//2:,:,:]
    latt = np.mean(lattice_his[lattice_his.shape[0]//2:,:,:],axis=0)
    #cart = frac*latt; frac= cart * inv(latt)
    pos = np.matmul(pos,np.linalg.inv(latt))
    #print(pos)
    lattice = [tuple(lat) for lat in latt]
    l = len(pos)
    #pos = [tuple(pos[i]) for i in range(int(l*0.6),l)]
    pos = [tuple(pos[i]) for i in range(l)]
    
    return [position_his,lattice_his,l]
#%%
import time

def spg_anal(strct_name):
    #return 0

#strct_name = 'Nd_H_2300'
    cutoff = [2,4,8]
    #cutoff = [2]
    [position_his,lattice_his,natoms] = get_pos(strct_name)
    ll = position_his.shape[0]
    
    sym = dict()
    if strct_name == '':
        print(" The material is current OUTCAR_collect.")
    else:
        print(" The material is %s." % strct_name)
    for c in cutoff:
        for i in range(c):
            pos = np.mean(position_his[ll - ll//c*(i+1):ll - ll//c*(i),:,:],axis=0)
            latt = np.mean(lattice_his[ll - ll//c*(i+1):ll - ll//c*(i),:,:],axis=0)
            #print(pos)
            #print(latt)
            pos = np.matmul(pos,np.linalg.inv(latt))
            #print(pos)
            #print(pos)
            latt = [tuple(lat) for lat in latt]
            pos = [tuple(pos[i]) for i in range(natoms)]        
            #pos = [tuple(pos[i]) for i in range(int(natoms*0.6),natoms)]        
            strct = (
                latt,
                pos,
                [
                    1,
                ]
                * int(natoms),
            )
            
            #print("[get_spacegroup]")
            #print(lattice)
            #print(pos)
            #time.sleep(5)
            sym_iter = spglib.get_spacegroup(strct,symprec=5e-2)
            if sym_iter in sym.keys():
                sym[sym_iter] += 1
            else:
                sym[sym_iter] = 1
            sym_iter = spglib.get_spacegroup(strct,symprec=1e-1)
            if sym_iter in sym.keys():
                sym[sym_iter] += 1
            else:
                sym[sym_iter] = 1
            sym_iter = spglib.get_spacegroup(strct,symprec=3e-1)
            if sym_iter in sym.keys():
                sym[sym_iter] += 1
            else:
                sym[sym_iter] = 1
            sym_iter = spglib.get_spacegroup(strct,symprec=5e-1)
            if sym_iter in sym.keys():
                sym[sym_iter] += 1
            else:
                sym[sym_iter] = 1
            print(" Based on the %d(th) %d MD steps, \tthe spacegroup of the structure is \t %s (tol = 0.05) \t %s (tol = 0.1) \t %s (tol = 0.2) \t %s (tol = 0.5)." 
                  % ( c-i, position_his.shape[0]//c, spglib.get_spacegroup(strct,symprec=5e-2), spglib.get_spacegroup(strct,symprec=1e-1), spglib.get_spacegroup(strct,symprec=3e-1), spglib.get_spacegroup(strct,symprec=5e-1) ) )
    print("")
    # Sort the dictionary by value in descending order
    sorted_dict = dict(sorted(sym.items(), key=lambda item: item[1], reverse=True))
    
    # Print the sorted dictionary
    print('symmetry | occurrence')
    for key, value in sorted_dict.items():
        print(f'{key}: \t{value}')
    #print(sym)
    
#%%
spg_anal('')
