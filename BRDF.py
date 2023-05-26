# -*- coding: utf-8 -*-
"""
Created on Sat Apr 17 14:52:04 2021

@author: Haoran
"""
import sys
import os
import numpy as np
from BRDF_func import verhoef_bimodal, campbell, CIxy, weighted_sum_over_lidf, sunshade_initial, sunshade, A_BRFv2_single_hemi, A_BRFv2_single_dif,A_BRFv2_single_hemi_dif,i_hemi
from BRF import A_BRFv2_initial

def Soil_spectrum(soil_spectrum, min_gcc, band_Pars):
    [[rs, re], [gs, ge], [bs, be]] = band_Pars
    min_soil = min_gcc*(np.nanmean(soil_spectrum[rs:re]) + np.nanmean(soil_spectrum[bs:be]))/(1 - min_gcc)
    avg_delta = min_soil - np.nanmean(soil_spectrum[gs:ge])
    if avg_delta > 0:
        soil_spectrum[gs:ge] += avg_delta
    return soil_spectrum 

def BRDF_initial(min_gcc, band_Pars):
    
    lidfa = 30 # float Leaf Inclination Distribution at regular angle steps. 
    lidfb = -0.15 # float Leaf Inclination Distribution at regular angle steps. 
    lidftype = 2 # float Leaf Inclination Distribution at regular angle steps.
    
    tts = 30.0 # float Sun Zenith Angle 
    tto = 45.0 # float View(sensor) Zenith Angl
    psi = 90 # float Relative Sensor-Sun Azimuth Angle

    if psi > 180:
        psi = abs(psi - 360)
   
    rsoil = 1.0 #土壤标量1（亮度）
    psoil = 1.0  #土壤标量2（水分）
    
    #Read soil and prospect parameters
    soil = np.loadtxt(os.path.join(os.path.split(os.path.abspath(sys.argv[0]))[0].replace("src", "pars"),"soil_reflectance.txt"))     
    soil_spectrum1 = soil[:,0]
    soil_spectrum2 = soil[:,1]
    soil_spectrum1 = Soil_spectrum(soil_spectrum1, min_gcc, band_Pars)
    
    rg = rsoil * (
        psoil * soil_spectrum1 + (1.0 - psoil) * soil_spectrum2
    ) #soil reflectance
    
    # Calcualte leaf angle distribution
    if lidftype == 1:
        lidf = verhoef_bimodal(lidfa, lidfb, n_elements=13)
    elif lidftype == 2:
        lidf = campbell(lidfa, n_elements=13)
    else:
        raise ValueError(
            "lidftype can only be 1 (Campbell) or 2 (ellipsoidal)"
        )
        
    CIy1 = 1
    CIy2 = 1

    CIs = CIxy(CIy1, CIy2, tts)
    CIo = CIxy(CIy1, CIy2, tto)      

    [Gs, Go, ks, ko, bf, sob, sof] = weighted_sum_over_lidf(lidf, tts, tto, psi)

    [Ps_arr, Po_arr, int_res_arr, nl] = sunshade_initial(tts, tto, psi, ks, ko, CIs, CIo)
    
    [hemi_pars, dif_pars, hemi_dif_pars] = A_BRFv2_initial(tts, tto, CIs, CIo, CIy1, CIy2, lidf)
    
    return [CIs, CIo, CIy1, CIy2, ks, ko, sob, sof, tts, tto, psi, lidf, rg], [Ps_arr, Po_arr, int_res_arr, nl], [hemi_pars, dif_pars, hemi_dif_pars]


def BRDF(lai, year, julian_day, SIP_Pars, BRDF_Pars,SUN_Pars,BRF_Pars, band_Pars):
    [rho, tau] = SIP_Pars
    
    [CIs, CIo, CIy1, CIy2, ks, ko, sob, sof, tts, tto, psi, lidf, rg] = BRDF_Pars
    
    [Ps_arr, Po_arr, int_res_arr, nl] = SUN_Pars
    
    [hemi_pars, dif_pars, hemi_dif_pars] = BRF_Pars
    #soil and canopy properties
    w = rho + tau   #leaf single scattering albedo
    
    #计算lai    
    i0 = 1 - np.exp(-ks * lai * CIs)
    iv = 1 - np.exp(-ko * lai * CIo)
    
    t0 = 1 - i0
    tv = 1 - iv
    
  
    [kc, kg]    =   sunshade(tts, tto, psi, ks, ko, CIs, CIo, Ps_arr, Po_arr, int_res_arr, nl, lai)       
    
    [sob_vsla,          sof_vsla,          kgd]     = A_BRFv2_single_hemi(hemi_pars, lai)       
    
    [sob_vsla_dif,      sof_vsla_dif,      kg_dif]  = A_BRFv2_single_dif(dif_pars,   lai)
    
    [sob_vsla_hemi_dif, sof_vsla_hemi_dif, kgd_dif] = A_BRFv2_single_hemi_dif(hemi_dif_pars, lai)    

    
    rho2 = iv/2/lai
    
    iD = i_hemi(CIy1,CIy2,lai,lidf)  
    td = 1 - iD
    
    p  = 1 - iD/lai  

    rho_hemi     = iD/2/lai        
    rho_dif      = iv/2/lai        
    rho_dif_hemi = iD/2/lai  
 
    wso  = sob*rho + sof*tau

    Tdn   = t0+i0*w*rho_hemi/(1-p*w)
    Tup_o = tv+iD*w*rho2/(1-p*w)
    Rdn   = iD*w*rho_hemi/(1-p*w)
    
    BRFv = wso*kc/ko + i0*w*w*p*rho2/(1-p*w)      
    BRFs = kg*rg
    BRFm = rg*Tdn*Tup_o/(1-rg*Rdn)-t0*rg*tv       
    BRF  = BRFv + BRFs + BRFm

    Tup_hemi = td + iD*w*rho_hemi/(1-p*w)
    
    Rv = sob_vsla*rho + sof_vsla*tau+i0*w**2*p*rho_hemi/(1-p*w) 
    Rs = kgd*rg;    
    Rm = rg*(Tdn)*(Tup_hemi)/(1-rg*(Rdn))-t0*rg*td    
    R  = Rv + Rs + Rm   #albedo

    #absorption
    Av  = i0*(1-w)/(1-p*w)
    Aup = iD*(1-w)/(1-p*w)
    Am  = rg*(Tdn)*(Aup)/(1-rg*(Rdn))
    A   = Av + Am    #absorption

    Tdn_dif  = td+iD*w*rho_dif_hemi/(1-p*w)
    Tup_difo = tv+iD*w*rho_dif/(1-p*w)
    Rdn_dif  = iD*w*rho_dif_hemi/(1-p*w)
    
    BRF_difv = sob_vsla_dif*rho + sof_vsla_dif*tau+iD*w**2*p*rho_dif/(1-p*w)          
    BRF_difs = kg_dif*rg    
    BRF_difm = rg*(Tdn_dif)*(Tup_difo)/(1-rg*(Rdn_dif))-td*rg*tv    
    BRF_dif  = BRF_difv + BRF_difs + BRF_difm


    Tup_dif_hemi = td+iD*w*rho_dif_hemi/(1-p*w)
    
    R_difv = sob_vsla_hemi_dif*rho + sof_vsla_hemi_dif*tau+iD*w**2*p*rho_dif_hemi/(1-p*w)    
    R_difs = kgd_dif*rg    
    R_difm = rg.dot(Tdn_dif)*(Tup_dif_hemi)/(1-rg*(Rdn_dif))-td*rg*td
    R_dif  = R_difv + R_difs + R_difm

    #absorption
    Aup_dif = iD*(1-w)/(1-p*w)
    A_difv  = iD*(1-w)/(1-p*w)
    A_difm  = rg*(Tdn_dif)*(Aup_dif)/(1-rg*(Rdn_dif))
    A_dif   = A_difv + A_difm 

    fPAR = sum(A[0:301])/301
    
    [[rs, re], [gs, ge], [bs, be]] = band_Pars
    Red = float(np.nanmean(BRF[rs:re]))
    Green = float(np.nanmean(BRF[gs:ge]))
    Blue = float(np.nanmean(BRF[bs:be])) 
    """
    if year == 2012 and julian_day == 200:
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(16,4))
        plt.plot(BRF[0:310])  
    """ 
    return Green/(Red+Green+Blue),fPAR

    