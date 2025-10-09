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
                #if counter%4000 == 0: print(counter)
                if counter % 4000 == 0:
                    if counter == 4000:
                        print("reading...", end="", flush=True)
                    else:
                        print(f" {counter}", end="", flush=True)

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
        if Line.startswith( ' energy of atom'):                              
            ntyp = Line.split( )                                             
            ntype = int(ntyp[3])   
        if Line.startswith( '   ions per type ='):
            natmtyp = Line.split( ) 
            natms = natmtyp[-ntype:]
        

        #else:
        #    continue
        #if not flag_simple: file_s.write(Line)
        
            
    file.close()

    #if not flag_simple: file_s.close()
    position_his = np.asarray(position_his)
    lattice_his = np.asarray(lattice_his)
    
    ##%%
    print()
    print(natms)
    #print(ntype)
    #ll = position_his.shape[0]
    #position_tmp = position_his[ll//2:ll,:,:]
    #position_tmp = position_his
    #pos = np.mean(position_his[position_his.shape[0]//2:,:,:],axis=0)
    #print(pos)
    #print(position_his[position_his.shape[0]//2:ll,:,:].shape)
    #lattice_his = lattice_his[lattice_his.shape[0]//2:,:,:]
    #latt = np.mean(lattice_his[lattice_his.shape[0]//2:,:,:],axis=0)
    #cart = frac*latt; frac= cart * inv(latt)
    #pos = np.matmul(pos,np.linalg.inv(latt))
    #print(pos)
    #lattice = [tuple(lat) for lat in latt]
    l = position_his.shape[1]
    #pos = [tuple(pos[i]) for i in range(int(l*0.6),l)]
    
    return [position_his,lattice_his,l, ntype,natms]
#%%
import time

def spg_anal(strct_name, frac):
    #return 0

#strct_name = 'Nd_H_2300'
    #cutoff = [1,2,4,8]
    cutoff = [2,4,8]
    [position_his,lattice_his,natoms,ntype,natms] = get_pos(strct_name)
    ll = position_his.shape[0]

    print(" The material is %s." % strct_name)
    for c in cutoff:
        pos = np.mean(position_his[ll - ll//c:,:,:],axis=0)
        latt = np.mean(lattice_his[ll - ll//c:,:,:],axis=0)
        pos = np.matmul(pos,np.linalg.inv(latt))
        #print(latt)
        #print(pos)
        latt = [tuple(lat) for lat in latt]
        pos = [tuple(pos[i]) for i in list(range(natoms))]        
        #pos = [tuple(pos)]  
#        atomind = [1] * int(natms[0]) + [2] * int(natms[1])      
        atomind = [item for sublist in [[str(1+i)] * int(natms[i]) for i in range(ntype)] for item in sublist]
        strct = (
            latt,
            pos,
            #[
            atomind,
            #],
        #    * natoms,
        )
        
        #print(strct)
        #print(pos)
        #print("[get_spacegroup]")
        #print(lattice)
        #print(pos)
        #time.sleep(5)
        print(" Based on the last %d MD steps, \tthe spacegroup of the structure is %s (tol = 0.02) \t %s (tol = 0.05) \t %s (tol = 0.1) \t %s (tol = 0.2) \t %s (tol = 0.5) \t %s (tol = 1.0) \t %s (tol = 2.0) \t %s (tol = 3.0)." 
              % ( position_his.shape[0]//c, spglib.get_spacegroup(strct,symprec=2e-2), spglib.get_spacegroup(strct,symprec=5e-2), spglib.get_spacegroup(strct,symprec=1e-1), spglib.get_spacegroup(strct,symprec=3e-1), spglib.get_spacegroup(strct,symprec=5e-1), spglib.get_spacegroup(strct,symprec=10e-1), spglib.get_spacegroup(strct,symprec=20e-1), spglib.get_spacegroup(strct,symprec=30e-1) ) )
        # Symmetrize the structure
        standardized_cell = spglib.standardize_cell(
            strct, symprec=10e-1,
            to_primitive=False,  # Set to True for a primitive cell, False for a conventional cell
            no_idealize=False    # Set to True to avoid idealizing lattice parameters
        )
        
        # Extract symmetrized lattice, positions, and atomic numbers
        sym_lattice, sym_positions, sym_atomic_numbers = standardized_cell
        
        print("Symmetrized Lattice:")
        print(sym_lattice)
        print("Symmetrized Positions:")
        print(sym_positions)
        print("Symmetrized Atomic Numbers:")
        print(sym_atomic_numbers)

    print("")
    with open('POSCAR_SLUSCHI', 'w') as out:
        out.write("POSCAR by SLUSCHI\n")
        out.write("1.0\n")
        # lattice vectors: each line is vector components (VASP expects 3 lines)
        for vec in sym_lattice:
            out.write(str(vec[0])+' '+str(vec[1])+' ' +str(vec[2]) + "\n")
        out.write("X\n")
        out.write(str(len(sym_atomic_numbers)) + "\n")
        out.write("Direct\n")
        for p in sym_positions:
            out.write("{: .16f} {: .16f} {: .16f}\n".format(p[0], p[1], p[2]))

    print(f"Wrote POSCAR")
    

    natms = [int(x) for x in natms]
    for nsep in range(ntype):                                          
       print(nsep + 1, natms[nsep])                                                               
       print(" The sublattice's symmetry is %s." % strct_name)
       if nsep == 0:                                                        
          start_index = 0                                                    
       else:                                                                
          start_index = sum(natms[:nsep])                               
                                                                                
       end_index = start_index + natms[nsep] 
       #start_index = 0
       #end_index = natms[0]+natms[1]

       for c in cutoff:                                                            
          pos = np.mean(position_his[ll - ll//c:,:,:],axis=0)                     
          latt = np.mean(lattice_his[ll - ll//c:,:,:],axis=0)                     
          pos = np.matmul(pos,np.linalg.inv(latt))                                
          #print(latt)                                                            
          #print(pos)                                                             
          latt = [tuple(lat) for lat in latt]                                     
          pos = [tuple(pos[i]) for i in range(start_index, end_index)]
          atomind = [1] * (end_index - start_index)
          strct = (                                                               
                 latt,                                                               
                 pos,                                                                
            #[                                                                  
                 atomind,                                                            
            #],                                                                 
        #    * natoms,                                                          
          )                               
        #print(strct)                                                           
        #print(pos)                                                             
        #print("[get_spacegroup]")                                              
        #print(lattice)                                                         
        #print(pos)                                                             
        #time.sleep(5)                                                          
          print(" Based on the last %d MD steps, \tthe spacegroup of the structure is %s (tol = 0.02) \t %s (tol = 0.05) \t %s (tol = 0.1) \t %s (tol = 0.2) \t %s (tol = 0.5) \t %s (tol = 1.0) \t %s (tol = 2.0) \t %s (tol = 3.0)." 
              % ( position_his.shape[0]//c, spglib.get_spacegroup(strct,symprec=2e-2), spglib.get_spacegroup(strct,symprec=5e-2), spglib.get_spacegroup(strct,symprec=1e-1), spglib.get_spacegroup(strct,symprec=3e-1), spglib.get_spacegroup(strct,symprec=5e-1), spglib.get_spacegroup(strct,symprec=10e-1), spglib.get_spacegroup(strct,symprec=20e-1), spglib.get_spacegroup(strct,symprec=30e-1) ) )
          # Symmetrize the structure
          standardized_cell = spglib.standardize_cell(
              strct, symprec=10e-1,
              to_primitive=False,  # Set to True for a primitive cell, False for a conventional cell
              no_idealize=False    # Set to True to avoid idealizing lattice parameters
          )
          
          # Extract symmetrized lattice, positions, and atomic numbers
          sym_lattice, sym_positions, sym_atomic_numbers = standardized_cell
          
          print("Symmetrized Lattice:")
          print(sym_lattice)
          print("Symmetrized Positions:")
          print(sym_positions)
          print("Symmetrized Atomic Numbers:")
          print(sym_atomic_numbers)
          print("")  

#          spg_anal('',1.0)
#    print("")

#%%

spg_anal('',1.0)
#spg_anlsep('',1.0)
