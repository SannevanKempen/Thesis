#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 13:13:44 2022

@author: Ingrid
"""
import numpy as np
from scipy import io

network = 56; # 47 or 56

if network == 47:
    matlabfile = 'ChargingProbSCE47.mat'
    tikztemplate = 'tikz_template47.txt'
    outputfile = 'SCE47.txt'
elif network == 56:
    matlabfile = 'ChargingProbSCE56.mat'
    tikztemplate = 'tikz_template56.txt'
    outputfile = 'SCE56.txt'

mat1 = io.loadmat(matlabfile) 
v = mat1['p'].flatten()
pmax = 0.4
pmin = 0
phalf = (pmax+pmin)/2
# x = 355/np.amax(np.array([v1,v2]))
top = np.array([0.3*255, 0, 0])
middle = np.array([255, 0, 0])
bottom = np.array([255, 255, 0])


# with open('tikz_template.txt','r') as f:
with open(tikztemplate,'r') as f:
    with open(outputfile,'w') as n:
        contents = f.readlines()
        n.write("% legend  \n")
        for line in contents[1:2]:
            line = line.replace("TOP1",np.array2string(np.round(top[0]/255,2)))
            line = line.replace("TOP2",np.array2string(np.round(top[1]/255,2)))
            line = line.replace("TOP3",np.array2string(np.round(top[2]/255,2)))
            line = line.replace("MID1",np.array2string(np.round(middle[0]/255,2)))
            line = line.replace("MID2",np.array2string(np.round(middle[1]/255,2)))
            line = line.replace("MID3",np.array2string(np.round(middle[2]/255,2)))
            line = line.replace("BOT1",np.array2string(np.round(bottom[0]/255,2)))
            line = line.replace("BOT2",np.array2string(np.round(bottom[1]/255,2)))
            line = line.replace("BOT3",np.array2string(np.round(bottom[2]/255,2)))
            n.write(line)
        n.write("\draw (0cm,0cm) node {\pgfuseshading{legend}};  \n")
        n.write("% tikz on legend \n")
        index=0
        tik = pmax
        for line in contents[4:9]:
            n.write(line.replace("X",np.array2string(np.around(tik,1))))
            tik = tik - (pmax-pmin)/4
        index = 0
        n.write("% vertices \n")
        n.write(contents[10])
        for line in contents[11::]:
            p = v[index]
            if p == 0:
                col = np.array([255,255,255])
            elif p >= phalf:
                col = top - (top-middle)*((pmax-p)/(pmax-phalf))
            else:
                col = middle - (middle-bottom)*((phalf - p)/(phalf-pmin))
            # red1 = 255-np.maximum(0,np.floor(x*v1[index])-255)
            # green1 = np.maximum(0,255-np.floor(x*v1[index]))
            # line = line.replace("X",np.array2string(red1))
            # n1.write(line.replace("Z",np.array2string(green1)))
            
            line = line.replace("X",np.array2string(np.around(col[0]).astype(int)))
            line = line.replace("Y",np.array2string(np.around(col[1]).astype(int)))
            n.write(line.replace("Z",np.array2string(np.around(col[2]).astype(int))))
            index += 1
            
