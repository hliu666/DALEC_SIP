# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 12:18:48 2021

@author: Haoran
"""
import numpy as np
import sys
import os

def SIP_Model(cab):
    cab    = cab    #叶绿素浓度
    car    = 8     #类胡萝卜素浓度
    cbrown = 0     #棕色/衰老色素 范围[0,1]
    cw     = 0.01  #叶片等效水
    cm     = 0.005 #干物质
    ant    = 0.0   #花青素含量  
    alpha  = 600   # constant for the the optimal size of the leaf scattering element   

    prospectpro = np.loadtxt(os.path.join(os.path.split(os.path.abspath(sys.argv[0]))[0].replace("src", "pars"),"dataSpec_PDB.txt")) 
    
    lambdas = prospectpro[:,0]
    nr      = prospectpro[:,1]
    kab     = prospectpro[:,2]    
    kcar    = prospectpro[:,3]
    kant    = prospectpro[:,4]    
    kbrown  = prospectpro[:,5]
    kw      = prospectpro[:,6]    
    km      = prospectpro[:,7]

    kall    = (cab*kab + car*kcar + ant*kant + cbrown*kbrown + cw*kw + cm*km)/(cm*alpha)
    w0      = np.exp(-kall)
    
    # spectral invariant parameters
    fLMA = 2765.0*cm
    gLMA = 102.8 *cm
    
    p = 1-(1 - np.exp(-fLMA))/fLMA
    q = 1- 2 * np.exp(-gLMA)
    qabs = np.sqrt(q**2)
    
    # leaf single scattering albedo
    w = w0*(1-p)/(1-p*w0)
    
    # leaf reflectance and leaf transmittance
    refl  = w*(1/2+q/2*(1-p*w0)/(1-qabs*p*w0))
    trans = w*(1/2-q/2*(1-p*w0)/(1-qabs*p*w0))
    # rho, array_like 叶片朗伯反射率。
    # tau,  array_like 叶片透过率。    
    return [refl, trans]

#refl, trans = SIP_Model()
#import matplotlib.pyplot as plt
#fig = plt.figure(figsize=(16,4))
#plt.plot(refl, c = "red") 