#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 14:10:47 2023

@author: qhong7
"""
import re
import numpy as np
from matplotlib import pyplot as plt

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

def read_natoms():
    file = open('param','r')
    Line = file.readline()
    Line = file.readline()
    #print(Line)
    tmp = Line.split('  ')
    #print(tmp)
    natoms = [int(item.replace(" ", "").replace("\n", "")) for item in tmp if item != ""]
    print(natoms)
    file.close()
    return natoms

def diff(filename):
    #strct = 'Er_L_3100'
    natoms = read_natoms()
    file = open(filename, 'r')
    
    i = 0
    flag_lattice, flag_position = False, False
    counter = 0
    dist_his = []
    time = 0
    time_his = []
    stepsize = 10.0
    while i==0 or Line:# and counter < 17000:
    #while i < 100000:
        i += 1
        try:
            Line = file.readline()
            #print(Line)
            
            if flag_lattice:
                lattice.append(line2vec(Line))
                remain_lattice -= 1
                if remain_lattice == 0: 
                    flag_lattice = False
                    lattice = np.asarray(lattice)
            
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
                        time += stepsize
                        #if time > 1000:
                        if time > 0:
                            dist_his.append(dist)
                            time_his.append(time)
                    counter += 1
                    if counter%4000 == 0: print(counter)
                    #if counter == 20000:
                    #    for i in range(dist_mat.shape[0]):
                    #        print(dist_mat[i])
                    position_last = position
    
            if Line.startswith( '      direct lattice vectors'): 
                flag_lattice = True
                remain_lattice = 3
                lattice = []
            if Line.startswith( ' POSITION'):  
                flag_position = True
                remain_position = natom+1  # 120 80
                position = []
            if Line.startswith( '   POTIM'): 
                tmp = Line.split(' ')
                if len(tmp) > 6:
                    stepsize = float(Line.split(' ')[6])
                else:
                    if stepsize == 10.0 and len(time_his) > 0:
                        print('WARNING: Error reading POTIM. It may cause failure to calculate diffusion coefficient.')
                        print(Line)
                #print(Line.split(' '))
                #stepsize = float(Line.split(' ')[5])
            if Line.startswith( '   number of dos'):            
                natom = int(Line.split('NIONS =    ')[1])
                #print("natom: ","{:.0f}".format(natom))
        except:
            file.close()
    print("natom: ","{:.0f}".format(natom))
    
    import os
    #os.environ["PATH"] = "/Library/TeX/texbin" + os.pathsep + os.getenv("PATH")
    #plt.rcParams["font.family"] = "Times New Roman"
    #plt.rcParams['text.usetex'] = True
    
    dist_his = np.asarray(dist_his)
    time_his = np.asarray(time_his)
    plt.close()
    clr = ['r','g','c','b','k']
    #print(dist_his.shape)
    for i in range(dist_his.shape[1]):
        natom_sum = 0
        for j in range(len(natoms)):
            natom_sum = natom_sum + natoms[j]
            if i < natom_sum:
                plt.plot(time_his,dist_his[:,i],clr[j])
                break
    max_y = np.max(dist_his)
    plt.text(0,max_y*0.95,'Oxygen',color='r')
    plt.text(0,max_y*0.9,'Iron',color='k')
    plt.xlabel('Time [fs]')
    plt.ylabel('Distance [\AA]')
    plt.grid(which='both')
    plt.savefig('diffusion.pdf')
      
    plt.close()
    dist_his_mean = []
    for i in range(dist_his.shape[0]):
        #mean = np.mean(dist_his[i,:145])
        #std = np.std(dist_his[i,:145])/np.sqrt(144.)
        natom_sum = 0
        for j in range(len(natoms)):
            mean = np.sqrt( np.mean( np.square(dist_his[i,natom_sum:natom_sum+natoms[j]]) ) )
            #std = np.sqrt( np.std( np.square(dist_his[i,:145]) )/np.sqrt(144.) )
            plt.plot([time_his[i],time_his[i]],[mean,mean],'.'+clr[j])
            dist_his_mean.append(mean)
            #mean = np.mean(dist_his[i,145:])
            #std = np.std(dist_his[i,145:])/np.sqrt(96.)
            #mean = np.sqrt( np.mean( np.square(dist_his[i,int(natom*0.5):]) ) )
            ##std = np.sqrt( np.std( np.square(dist_his[i,145:]) )/np.sqrt(96.) )
            #plt.plot([time_his[i],time_his[i]],[mean,mean],'.k')
            #dist_his_mean = np.asarray(dist_his_mean)
    max_y = np.max(dist_his_mean)
    plt.text(1000,max_y*0.95,'Oxygen',color='r')
    plt.text(1000,max_y*0.9,'Iron',color='k')
    plt.xlabel('Time [fs]')
    plt.ylabel('Distance [\AA]')
    plt.grid(which='both')
    plt.savefig('diffusion_avg.pdf')
    
    plt.close()
    #os.environ["PATH"] = "/Library/TeX/texbin" + os.pathsep + os.getenv("PATH")
    #plt.rcParams["font.family"] = "Times New Roman"                                 
    #plt.rcParams['text.usetex'] = True
    font = {'size'   : 18}
    plt.rc('font', **font)
    plt.figure(figsize=(9,6))
    #import matplotlib as mpl
    #mpl.rcParams.update(mpl.rcParamsDefault)
    natom_sum = 0
    for j in range(len(natoms)):
        print('Element # '+str(j+1))
        mean_his=[]
        for i in range(dist_his.shape[0]):
            mean = np.mean( np.square(dist_his[i,natom_sum:natom_sum+natoms[j]]) )
            std = np.std( np.square(dist_his[i,natom_sum:natom_sum+natoms[j]]) )/np.sqrt((natoms[j]))
            mean_his.append(mean)
            plt.plot([time_his[i]/1000,time_his[i]/1000],[(mean-std),(mean+std)],clr[j])
            #std = np.std(dist_his[i,:145])/np.sqrt(144.)
            #plt.plot([time_his[i],time_his[i]],[(mean-std)**2.,(mean+std)**2.],'r')
            #mean = np.mean(dist_his[i,145:])
            #std = np.std(dist_his[i,145:])/np.sqrt(96.)
            #plt.plot([time_his[i],time_his[i]],[(mean-std)**2.,(mean+std)**2.],'k')
        #print('here')
        #print(dist_his.shape)
        #for i in range(int(natom*0.6)):
            #print(np.square(dist_his[-1,i]))
        #print(np.square(time_his[-1]))
        mean_his = np.asarray(mean_his)
        #print(np.polyfit(time_his,mean_his,1, full=True))
        #print(np.floor(len(time_his)*0.1))
        #print(int(np.floor(len(time_his)*0.1)))
        #print(time_his[int(np.floor(len(time_his)*0.1)):])
        #[c,res,rank,sig,rcod] = np.polyfit(time_his[int(np.floor(len(time_his)*0.1)):],mean_his[int(np.floor(len(time_his)*0.1)):],1, full=True)
        [c,res,rank,sig,rcod] = np.polyfit(time_his[int(np.floor(len(time_his)*0.1)):],mean_his[int(np.floor(len(time_his)*0.1)):],1, full=True,)
        var = np.var(mean_his)*len(mean_his)
        R2 = 1.-res/var
        print(c,res,var, R2)
        print('Total MD length: \t\t'+str("{:.2f}".format(time_his[-1])+' fs'))
        print('Diffusion coefficient is: \t'+str("{:.2e}".format(c[0]/6.)) + ' Ang^2/fs or '+str("{:.2e}".format(c[0]/60.)) + ' cm^2/s')
        print('R2 of the linear fitting is: \t'+str("{:.2f}".format(R2[0])))
        plt.plot(time_his/1000,time_his*c[0]+c[1],color='m')
        max_mean_his = max(mean_his)
        #plt.text(max(time_his/1000)/100+1,max(mean_his)*0.95,"$\sigma^2$ = "+"{:.1e}".format(c[0]/10.)+"$\cdot t$+"+"{:.1e}".format(c[1]/60.)+", $D$: "+"{:.1e}".format(c[0]/60.)+" cm$^2$/s, R$^2$: "+"{:.2f}".format(R2[0]))
        plt.text(max(time_his/1000)/100,max(mean_his)*0.95,"$\sigma^2$ = "+"{:.1e}".format(c[0]/10.)+"$\cdot t$+"+"{:.1e}".format(c[1]/60.)+", $D$: "+"{:.1e}".format(c[0]/60.)+" cm$^2$/s, R$^2$: "+"{:.2f}".format(R2[0]))
        coeffs, cov = np.polyfit(time_his[int(np.floor(len(time_his)*0.1)):],mean_his[int(np.floor(len(time_his)*0.1)):], 1, cov=True)
        slope = coeffs[0]
        intercept = coeffs[1]
        slope_uncertainty = np.sqrt(cov[0, 0])  # Standard error of the slope
        intercept_uncertainty = np.sqrt(cov[1, 1])  # Standard error of the intercept
        print('The uncertainty of the slope is: \t'+str("{:.2e}".format(np.sqrt(cov[0, 0])/60.)))
        #print(mean_his)

        #plt.text(max(time_his/1000)/100,max(mean_his)*0.88,"R2: "+"{:.2f}".format(R2[0]))
        natom_sum = natom_sum + natoms[j]
    #mean_his=[]
    #for i in range(dist_his.shape[0]):
    #    mean = np.mean( np.square(dist_his[i,int(natom*0.5):]) )
    #    std = np.std( np.square(dist_his[i,int(natom*0.5):]) )/np.sqrt((natom*0.5))
    #    mean_his.append(mean)
    #    plt.plot([time_his[i]/1000,time_his[i]/1000],[(mean-std),(mean+std)],'k')
    #mean_his = np.asarray(mean_his)
    ##print(np.polyfit(time_his,mean_his,1, full=True))
    #[c,res,rank,sig,rcod] = np.polyfit(time_his[int(np.floor(len(time_his)*0.1)):],mean_his[int(np.floor(len(time_his)*0.1)):],1, full=True)
    #var = np.var(mean_his)*len(mean_his)
    #R2 = 1.-res/var
    #print(c,res,var, R2)
    #plt.plot(time_his/1000,time_his*c[0]+c[1],color='b')
    #plt.text(max(time_his/1000)/100,max_mean_his*0.75,"D2 = "+"{:.2e}".format(c[0])+"t+"+"{:.2f}".format(c[1]))
    #plt.text(max(time_his/1000)/100,max_mean_his*0.68,"R2: "+"{:.2f}".format(R2[0]))
    max_y = np.max(mean_his)
    #if c[1] > 0:
    #    plt.text(3000,max_y*0.6,r'$r^2 = '+"{:10.3g}".format(c[0])+' \cdot t + '+"{:10.3g}".format(c[1])+'$')
    #else:
    #    plt.text(3000,max_y*0.6,r'$r^2 = '+"{:10.3g}".format(c[0])+' \cdot t - '+"{:10.3g}".format(c[1]*-1.)+'$')
    #plt.text(3000,max_y*0.53,r'$\sigma_r^2 = 6Dt$')
    #plt.text(3000,max_y*0.46,r'$D = $'+"{:10.3g}".format(c[0]/6/10)+' cm$^2$/s')
    #plt.text(1000,max_y*0.95,'Oxygen',color='r')
    #plt.text(1000,max_y*0.9,'Iron',color='k')
    plt.xlabel('Time [ps]')
    plt.ylabel('Distance$^2$ [$\AA^2$]')
    plt.grid(which='both')
    #plt.xlim([0,16])
    #plt.ylim([0,50])
    plt.savefig('diffusion_coef.png',dpi=600)
    
diff('OUTCAR_collect')

#diff('Er_L_3300')
#diff('Er_L_3100')
#diff('Er_L_2900')
#diff('Er_L_2700')
#diff('Er_H_2900')
#diff('Er_H_2700')
#diff('Er_H_2500')
#diff('Er_C_2700')
#diff('Er_C_2500')
