# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 20:00:42 2021

@author: Deng
"""
import numpy as np
from math import exp, radians, cos, sin, tan, pi
from BRDF_func import weighted_sum_over_lidf, CIxy, sunshade_initial

def A_BRFv2_initial(tts, tto, CIs, CIo, CIy1, CIy2, lidf):
    hemi_pars     = A_BRFv2_single_hemi_initial(tts, CIs, CIy1, CIy2, lidf)
    dif_pars      = A_BRFv2_single_dif_initial(tto, CIo, CIy1, CIy2, lidf)  
    hemi_dif_pars = A_BRFv2_single_hemi_dif_initial(CIy1, CIy2, lidf)
    return hemi_pars, dif_pars, hemi_dif_pars
    
def A_BRFv2_single_hemi_initial(tts, CIs, CIy1, CIy2, lidf):
    xx=np.array([0.9602898565, -0.9602898565, 0.7966664774, -0.7966664774, 0.5255324099, -0.5255324099, 0.1834346425, -0.1834346425])
    
    ww=np.array([0.1012285363,  0.1012285363, 0.2223810345,  0.2223810345, 0.3137066459,  0.3137066459, 0.3626837834,  0.3626837834])   
    
    # * define limits of integration and the convertion factors for integration
    # * over thetaL (note the tL suffix!)
    upperlimit_tL = np.pi/2.0
    lowerlimit_tL = 0.0
    conv1_tL = (upperlimit_tL-lowerlimit_tL)/2.0
    conv2_tL = (upperlimit_tL+lowerlimit_tL)/2.0    
    neword_tL = conv1_tL*xx + conv2_tL   
    
    # * define limits of integration and the convertion factors for integration
    # * over phiL (note the pL suffix!)
    upperlimit_pL = 2.0*pi
    lowerlimit_pL = 0.0
    conv1_pL = (upperlimit_pL-lowerlimit_pL)/2.0
    conv2_pL = (upperlimit_pL+lowerlimit_pL)/2.0  
    neword_pL  = conv1_pL*xx + conv2_pL

    tta_arr  = neword_tL*180/pi;              # observer zenith angle
    psi_arr  = neword_pL*180/pi;              # relative azimuth angle    
    
    pars = []
    for i in range(len(psi_arr)):
        [Gs, Go, ks, ko, bf, sob, sof] = weighted_sum_over_lidf(lidf, tts, tta_arr[i], psi_arr[i])
        CIo = CIxy(CIy1, CIy2, tta_arr[i])   
        [Ps_arr, Po_arr, int_res_arr, nl] = sunshade_initial(tts, tta_arr[i], psi_arr[i], ks, ko, CIs, CIo)
        pars.append([tts, tta_arr[i], psi_arr[i], ks, ko, sob, sof, CIs, CIo, Ps_arr, Po_arr, int_res_arr, nl])
    
    return pars

def A_BRFv2_single_dif_initial(tto, CIo, CIy1, CIy2, lidf):
    xx=np.array([0.9602898565, -0.9602898565, 0.7966664774, -0.7966664774, 0.5255324099, -0.5255324099, 0.1834346425, -0.1834346425])
    
    ww=np.array([0.1012285363,  0.1012285363, 0.2223810345,  0.2223810345, 0.3137066459,  0.3137066459, 0.3626837834,  0.3626837834])   
    
    # * define limits of integration and the convertion factors for integration
    # * over thetaL (note the tL suffix!)
    upperlimit_tL = np.pi/2.0
    lowerlimit_tL = 0.0
    conv1_tL = (upperlimit_tL-lowerlimit_tL)/2.0
    conv2_tL = (upperlimit_tL+lowerlimit_tL)/2.0    
    neword_tL = conv1_tL*xx + conv2_tL   
    
    # * define limits of integration and the convertion factors for integration
    # * over phiL (note the pL suffix!)
    upperlimit_pL = 2.0*pi
    lowerlimit_pL = 0.0
    conv1_pL = (upperlimit_pL-lowerlimit_pL)/2.0
    conv2_pL = (upperlimit_pL+lowerlimit_pL)/2.0  
    neword_pL  = conv1_pL*xx + conv2_pL

    tta_arr = neword_tL*180/pi;              # observer zenith angle
    psi_arr = neword_pL*180/pi;              # relative azimuth angle    


    pars = []     
    for i in range(len(psi_arr)):                
        [Ga, Go, ks, ko, bf, sob, sof] = weighted_sum_over_lidf(lidf, tta_arr[i], tto, psi_arr[i])
        CIa = CIxy(CIy1, CIy2, tta_arr[i])          
        [Ps_arr, Po_arr, int_res_arr, nl] = sunshade_initial(tta_arr[i], tto, psi_arr[i], ks, ko, CIa, CIo)
        pars.append([tta_arr[i], tto, psi_arr[i], ks, ko, sob, sof, CIa, CIo, Ps_arr, Po_arr, int_res_arr, nl])                    
    
    return pars
        
def A_BRFv2_single_hemi_dif_initial(CIy1, CIy2, lidf):    
    xx=np.array([0.9602898565, -0.9602898565, 0.7966664774, -0.7966664774, 0.5255324099, -0.5255324099, 0.1834346425, -0.1834346425])
    
    ww=np.array([0.1012285363,  0.1012285363, 0.2223810345,  0.2223810345, 0.3137066459,  0.3137066459, 0.3626837834,  0.3626837834])   

    # * define limits of integration and the convertion factors for integration
    # * over thetaL (note the tL suffix!)
    upperlimit_mL = pi/2.0
    lowerlimit_mL = 0.0
    conv1_mL = (upperlimit_mL-lowerlimit_mL)/2.0
    conv2_mL = (upperlimit_mL+lowerlimit_mL)/2.0
    neword_mL = conv1_mL*xx + conv2_mL   

        
    #   * define limits of integration and the convertion factors for integration
    # * over phiL (note the pL suffix!)
    upperlimit_nL = 2.0*pi
    lowerlimit_nL = 0.0
    conv1_nL = (upperlimit_nL-lowerlimit_nL)/2.0
    conv2_nL = (upperlimit_nL+lowerlimit_nL)/2.0
    neword_nL = conv1_nL*xx + conv2_nL   

    
    # * define limits of integration and the convertion factors for integration
    # * over thetaL (note the tL suffix!)
    upperlimit_tL = np.pi/2.0
    lowerlimit_tL = 0.0
    conv1_tL = (upperlimit_tL-lowerlimit_tL)/2.0
    conv2_tL = (upperlimit_tL+lowerlimit_tL)/2.0    
    neword_tL = conv1_tL*xx + conv2_tL   

    
    # * define limits of integration and the convertion factors for integration
    # * over phiL (note the pL suffix!)
    upperlimit_pL = 2.0*pi
    lowerlimit_pL = 0.0
    conv1_pL = (upperlimit_pL-lowerlimit_pL)/2.0
    conv2_pL = (upperlimit_pL+lowerlimit_pL)/2.0  
    neword_pL  = conv1_pL*xx + conv2_pL

    tts_arr = neword_mL*180/pi
    tto_arr = neword_tL*180/pi
    psi_arr = abs(neword_nL*180/pi-neword_pL*180/pi)
    #psi_arr = abs(psi_arr-360*round(psi_arr/360))

    pars = []            
    for i in range(len(psi_arr)):        
        [Ga, Go, ks, ko, bf, sob, sof] = weighted_sum_over_lidf(lidf, tts_arr[i], tto_arr[i], psi_arr[i])
        CIs = CIxy(CIy1, CIy2, tts_arr[i])
        CIo = CIxy(CIy1, CIy2, tto_arr[i])     
        [Ps_arr, Po_arr, int_res_arr, nl] = sunshade_initial(tts_arr[i], tto_arr[i], psi_arr[i], ks, ko, CIs, CIo)
        pars.append([tts_arr[i], tto_arr[i], psi_arr[i], ks, ko, sob, sof, CIs, CIo, Ps_arr, Po_arr, int_res_arr, nl])  
        
    return pars