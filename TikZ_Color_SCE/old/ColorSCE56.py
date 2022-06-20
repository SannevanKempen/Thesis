#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 13:13:44 2022

@author: Ingrid
"""
import numpy as np
from scipy import io

mat = io.loadmat('SCE56.mat') 
v1 = mat['probs'][0]
v2 = mat['probs'][1]
v2 = .25*v1+.75*v2
nu = mat['nu'].flatten()[0]
mu = mat['mu'].flatten()[1]
pmax = (1/nu)/(1/nu + 1/mu)
ar = np.array([v1,v2])
pmax = np.ceil(np.amax(ar)*10)/10
pmin = np.floor(np.amin(np.where(ar >0, ar, 100))*10)/10
phalf = (pmax-pmin)/2
# x = 355/np.amax(np.array([v1,v2]))
top = np.array([0.3*255, 0, 0])
middle = np.array([255, 0, 0])
bottom = np.array([255, 255, 0])

# performance measures
print(.25*mat['sumPower'][0,0]+.75*mat['sumPower'][0,1])
print(.25*mat['sumUcars'][0,0]+.75*mat['sumUcars'][0,1])
print(.25*mat['sumProb'][0,0]+.75*mat['sumProb'][0,1])



# with open('tikz_template.txt','r') as f:
with open('tikz_template56.txt','r') as f:
    with open('SCE56_NO_RES.txt','w') as n1:
        contents = f.readlines()
        n1.write("% legend  \n")
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
            n1.write(line)
        n1.write("\draw (0cm,0cm) node {\pgfuseshading{legend}};  \n")
        n1.write("% tikz on legend \n")
        index=0
        tik = np.around(pmax,2)
        for line in contents[4:9]:
            if tik < 0.01: tik = 0
            n1.write(line.replace("X",np.array2string(np.around(tik,2))))
            tik = tik - (pmax-pmin)/4
        index = 0
        n1.write("% vertices \n")
        for line in contents[10::]:
            p1 = v1[index]
            if p1 == 0:
                col = np.array([255,255,255])
            elif p1 >= phalf:
                col = top - (top-middle)*((pmax-p1)/(pmax-phalf))
            else:
                col = middle - (middle-bottom)*((phalf - p1)/(phalf-pmin))
            # red1 = 255-np.maximum(0,np.floor(x*v1[index])-255)
            # green1 = np.maximum(0,255-np.floor(x*v1[index]))
            # line = line.replace("X",np.array2string(red1))
            # n1.write(line.replace("Z",np.array2string(green1)))
            
            line = line.replace("X",np.array2string(np.around(col[0]).astype(int)))
            line = line.replace("Y",np.array2string(np.around(col[1]).astype(int)))
            n1.write(line.replace("Z",np.array2string(np.around(col[2]).astype(int))))
            index += 1
            
with open('tikz_template56.txt','r') as f:            
    with open('SCE56_RES.txt','w') as n2: 
        contents = f.readlines()
        n2.write("% legend  \n")
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
            n2.write(line)
        n2.write("\draw (0cm,0cm) node {\pgfuseshading{legend}};  \n")
        n2.write("% tikz on legend \n")
        index=0
        tik = pmax
        for line in contents[4:9]:
            if tik < 0.01: tik = 0
            n2.write(line.replace("X",np.array2string(np.around(tik,2))))
            tik = tik - (pmax-pmin)/4
        index = 0
        n2.write("% vertices \n")
        for line in contents[10::]:
            p2 = v2[index]
            if p2 == 0:
                col = np.array([255,255,255])
            elif p2 >= phalf:
                col = top - (top-middle)*((pmax-p2)/(pmax-phalf))
            else:
                col = middle - (middle-bottom)*((phalf - p2)/(phalf-pmin))
            
            line = line.replace("X",np.array2string(np.around(col[0]).astype(int)))
            line = line.replace("Y",np.array2string(np.around(col[1]).astype(int)))
            n2.write(line.replace("Z",np.array2string(np.around(col[2]).astype(int))))
            index += 1