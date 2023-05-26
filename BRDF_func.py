# -*- coding: utf-8 -*-
"""
Created on Sat Apr 17 17:57:31 2021

@author: Haoran

Improve the model efficiency 
"""
import numpy as np
from math import exp, radians, cos, sin, tan, pi

def verhoef_bimodal(a, b, n_elements=18):
    """Calculate the Leaf Inclination Distribution Function based on the
    Verhoef's bimodal LIDF distribution.
    Parameters
    ----------
    a : float
        controls the average leaf slope.
    b : float
        controls the distribution's bimodality.
            * LIDF type     [a,b].
            * Planophile    [1,0].
            * Erectophile   [-1,0].
            * Plagiophile   [0,-1].
            * Extremophile  [0,1].
            * Spherical     [-0.35,-0.15].
            * Uniform       [0,0].
            * requirement: |LIDFa| + |LIDFb| < 1.
    n_elements : int
        Total number of equally spaced inclination angles.
    Returns
    -------
    lidf : list
        Leaf Inclination Distribution Function at equally spaced angles.
    References
    ----------
    .. [Verhoef1998] Verhoef, Wout. Theory of radiative transfer models applied
        in optical remote sensing of vegetation canopies.
        Nationaal Lucht en Ruimtevaartlaboratorium, 1998.
        http://library.wur.nl/WebQuery/clc/945481.
    """

    freq = 1.0
    step = 90.0 / n_elements
    lidf = np.zeros(n_elements) * 1.0
    angles = (np.arange(n_elements) * step)[::-1]
    i = 0
    for angle in angles:
        tl1 = np.radians(angle)
        if a > 1.0:
            f = 1.0 - np.cos(tl1)
        else:
            eps = 1e-8
            delx = 1.0
            x = 2.0 * tl1
            p = float(x)
            while delx >= eps:
                y = a * np.sin(x) + 0.5 * b * np.sin(2.0 * x)
                dx = 0.5 * (y - x + p)
                x = x + dx
                delx = abs(dx)
            f = (2.0 * y + p) / np.pi
        freq = freq - f
        lidf[i] = freq
        freq = float(f)
        i += 1
    lidf = lidf[::-1]
    return lidf

def campbell(alpha, n_elements=13):
    """Calculate the Leaf Inclination Distribution Function based on the
    mean angle of [Campbell1990] ellipsoidal LIDF distribution.
    Parameters
    ----------
    alpha : float
        Mean leaf angle (degrees) use 57 for a spherical LIDF.
    n_elements : int
        Total number of equally spaced inclination angles .
    Returns
    -------
    lidf : list
        Leaf Inclination Distribution Function for 18 equally spaced angles.
    References
    ----------
    .. [Campbell1986] G.S. Campbell, Extinction coefficients for radiation in
        plant canopies calculated using an ellipsoidal inclination angle distribution,
        Agricultural and Forest Meteorology, Volume 36, Issue 4, 1986, Pages 317-321,
        ISSN 0168-1923, http://dx.doi.org/10.1016/0168-1923(86)90010-9.
    .. [Campbell1990] G.S Campbell, Derivation of an angle density function for
        canopies with ellipsoidal leaf angle distributions,
        Agricultural and Forest Meteorology, Volume 49, Issue 3, 1990, Pages 173-176,
        ISSN 0168-1923, http://dx.doi.org/10.1016/0168-1923(90)90030-A.
    """

    alpha = float(alpha)
    excent = exp(
        -1.6184e-5 * alpha ** 3.0
        + 2.1145e-3 * alpha ** 2.0
        - 1.2390e-1 * alpha
        + 3.2491
    )
    sum0 = 0.0
    freq = np.zeros(n_elements)
    step = 90.0 / n_elements
    for i in range(n_elements):
        tl1 = radians(i * step)
        tl2 = radians((i + 1.0) * step)
        x1 = excent / (np.sqrt(1.0 + excent ** 2.0 * np.tan(tl1) ** 2.0))
        x2 = excent / (np.sqrt(1.0 + excent ** 2.0 * np.tan(tl2) ** 2.0))
        if excent == 1.0:
            freq[i] = abs(np.cos(tl1) - np.cos(tl2))
        else:
            alph = excent / np.sqrt(abs(1.0 - excent ** 2.0))
            alph2 = alph ** 2.0
            x12 = x1 ** 2.0
            x22 = x2 ** 2.0
            if excent > 1.0:
                alpx1 = np.sqrt(alph2 + x12)
                alpx2 = np.sqrt(alph2 + x22)
                dum = x1 * alpx1 + alph2 * np.log(x1 + alpx1)
                freq[i] = abs(dum - (x2 * alpx2 + alph2 * np.log(x2 + alpx2)))
            else:
                almx1 = np.sqrt(alph2 - x12)
                almx2 = np.sqrt(alph2 - x22)
                dum = x1 * almx1 + alph2 * np.arcsin(x1 / alph)
                freq[i] = abs(
                    dum - (x2 * almx2 + alph2 * np.arcsin(x2 / alph))
                )
    sum0 = np.sum(freq)
    lidf = np.zeros(n_elements)
    for i in range(n_elements):
        lidf[i] = freq[i] / sum0

    return lidf

def CIxy(CIy1,CIy2,tts):
    CIs=(CIy2 - CIy1)/(75 - 0)*(tts - 0) + CIy1;
    return CIs

def weighted_sum_over_lidf(lidf, tts, tto, psi):
    ks = 0.0
    ko = 0.0
    bf = 0.0
    sob = 0.0
    sof = 0.0
    cts = np.cos(np.radians(tts))
    cto = np.cos(np.radians(tto))
    ctscto = cts * cto

    n_angles = len(lidf)
    angle_step = float(90.0 / n_angles)
    litab = np.arange(n_angles) * angle_step + (angle_step * 0.5)

    for i, ili in enumerate(litab):
        ttl = 1.0 * ili
        cttl = np.cos(np.radians(ttl))
        # SAIL volume scattering phase function gives interception and portions to be multiplied by rho and tau
        [chi_s, chi_o, frho, ftau] = volscatt(tts, tto, psi, ttl)
        # Extinction coefficients
        ksli = chi_s / cts
        koli = chi_o / cto
        # Area scattering coefficient fractions
        sobli = frho * np.pi / ctscto
        sofli = ftau * np.pi / ctscto
        bfli = cttl ** 2.0
        ks += ksli * float(lidf[i])
        ko += koli * float(lidf[i])
        bf += bfli * float(lidf[i])
        sob += sobli * float(lidf[i])
        sof += sofli * float(lidf[i])

    Gs = ks * cts
    Go = ko * cto   
     
    return Gs, Go, ks, ko, bf, sob, sof   

def volscatt(tts, tto, psi, ttl):
    """Compute volume scattering functions and interception coefficients
    for given solar zenith, viewing zenith, azimuth and leaf inclination angle.
    Parameters
    ----------
    tts : float
        Solar Zenith Angle (degrees).
    tto : float
        View Zenight Angle (degrees).
    psi : float
        View-Sun reliative azimuth angle (degrees).
    ttl : float
        leaf inclination angle (degrees).
    Returns
    -------
    chi_s : float
        Interception function  in the solar path.
    chi_o : float
        Interception function  in the view path.
    frho : float
        Function to be multiplied by leaf reflectance to obtain the volume scattering.
    ftau : float
        Function to be multiplied by leaf transmittance to obtain the volume scattering.
    References
    ----------
    Wout Verhoef, april 2001, for CROMA.
    """

    cts = np.cos(np.radians(tts))
    cto = np.cos(np.radians(tto))
    sts = np.sin(np.radians(tts))
    sto = np.sin(np.radians(tto))
    cospsi = np.cos(np.radians(psi))
    psir = np.radians(psi)
    cttl = np.cos(np.radians(ttl))
    sttl = np.sin(np.radians(ttl))
    cs = cttl * cts
    co = cttl * cto
    ss = sttl * sts
    so = sttl * sto
    cosbts = 5.0
    if np.abs(ss) > 1e-6:
        cosbts = -cs / ss
    cosbto = 5.0
    if np.abs(so) > 1e-6:
        cosbto = -co / so
    if np.abs(cosbts) < 1.0:
        bts = np.arccos(cosbts)
        ds = ss
    else:
        bts = np.pi
        ds = cs
    chi_s = 2.0 / np.pi * ((bts - np.pi * 0.5) * cs + np.sin(bts) * ss)
    if abs(cosbto) < 1.0:
        bto = np.arccos(cosbto)
        do_ = so
    else:
        if tto < 90.0:
            bto = np.pi
            do_ = co
        else:
            bto = 0.0
            do_ = -co
    chi_o = 2.0 / np.pi * ((bto - np.pi * 0.5) * co + np.sin(bto) * so)
    btran1 = np.abs(bts - bto)
    btran2 = np.pi - np.abs(bts + bto - np.pi)
    if psir <= btran1:
        bt1 = psir
        bt2 = btran1
        bt3 = btran2
    else:
        bt1 = btran1
        if psir <= btran2:
            bt2 = psir
            bt3 = btran2
        else:
            bt2 = btran2
            bt3 = psir
    t1 = 2.0 * cs * co + ss * so * cospsi
    t2 = 0.0
    if bt2 > 0.0:
        t2 = np.sin(bt2) * (
            2.0 * ds * do_ + ss * so * np.cos(bt1) * np.cos(bt3)
        )
    denom = 2.0 * np.pi ** 2
    frho = ((np.pi - bt2) * t1 + t2) / denom
    ftau = (-bt2 * t1 + t2) / denom
    if frho < 0.0:
        frho = 0.0
    if ftau < 0.0:
        ftau = 0.0

    return (chi_s, chi_o, frho, ftau)
  
def sunshade_initial(tts, tto, psi, ks, ko, CIs, CIo):

    tts0 = np.tan(np.radians(tts)) #tan_tts     
    tto0 = np.tan(np.radians(tto)) #tan_tto   

    nl   = 20
    xl   = np.linspace(0, -1, num = nl + 1).T
    #iLai = lai/nl
    dx   = 1/nl
    d    = 0.05
    H    = 1
    q    = d/H
    
    dso  = np.sqrt(np.power(tts0, 2) + np.power(tto0, 2) - 2*tts0*tto0*np.cos(np.radians(psi))) 
    """
    math.log() #只能对单个数求对数        
    math.exp() #只能对单个数求指数        
    numpy.log() #既能对单个数求对数，也可以对整个数组里的每一个数求对数（批量）        
    numpy.exp() #既能对单个数求指数，也可以对整个数组里的每一个数求指数（批量）
    """
    Ps_arr   = ks*xl*CIs
    Po_arr   = ko*xl*CIo

    #Ps0[1:nl] *= (1-np.exp(-ks*CIs*lai*dx))/(ks*CIs*lai*dx)  # Correct Ps/Po for finite dx
    #Po0[1:nl] *= (1-np.exp(-ko*CIo*lai*dx))/(ko*CIo*lai*dx)  # Correct Ps/Po for finite dx
    for i in range(len(xl)):
        int_res_sub  = np.array([])        
        #xl_sub = np.linspace(xl[i]-dx, xl[i], num = ni)
        a = xl[i]-dx
        b = xl[i]
        #http://liao.cpython.org/scipy18/
        int_res_sub = np.append(int_res_sub, Psofunction(a, ks, ko, CIs, CIo, q, dso))
        int_res_sub = np.append(int_res_sub, Psofunction(b, ks, ko, CIs, CIo, q, dso))
        int_res_sub = np.append(int_res_sub, Psofunction((a + b)/2, ks, ko, CIs, CIo, q, dso))        
        if i == 0:
           int_res_arr = int_res_sub.reshape(-1, 1)
        else:
           int_res_arr = np.hstack((int_res_arr, int_res_sub.reshape(-1, 1)))
         
    return Ps_arr, Po_arr, int_res_arr, nl     

def sunshade(tts, tto, psi, ks, ko, CIs, CIo, Ps_arr, Po_arr, int_res_arr, nl, lai):

    if tts == tto and psi == 0:
        kc = 1 - np.exp(-ks*CIs * lai)
        kg = np.exp(-ks*CIs * lai)        
     
    else:
        iLai = lai/nl
        dx   = 1/nl

        Ps   = np.exp(Ps_arr*lai)
        Po   = np.exp(Po_arr*lai)

        Ps[0:nl] *= (1-np.exp(-ks*CIs*lai*dx))/(ks*CIs*lai*dx)  # Correct Ps/Po for finite dx
        Po[0:nl] *= (1-np.exp(-ko*CIo*lai*dx))/(ko*CIo*lai*dx)  # Correct Ps/Po for finite dx
        
        int_res_arr = np.exp(int_res_arr*lai)/dx
        int_res = (int_res_arr[0,:] + int_res_arr[1,:] + 4*int_res_arr[2,:])*(dx/6)
        _arr = np.hstack((int_res.reshape(-1, 1), Ps.reshape(-1, 1), Po.reshape(-1, 1)))   
        Pso  = np.min(_arr, 1)
        
        kc = iLai * CIo * ko * sum(Pso[:nl])  #visable sunlit leaf
        kg = Pso[nl]                   #visable sunlit soil        
    return kc, kg 
    
def Psofunction(xl, ks, ko, CIs, CIo, q, dso):
    if dso != 0:
        alf = (dso / q) *2/(ks + ko)
        pso0 = (ks*CIs + ko*CIo)*xl + np.sqrt(ks*CIs*ko*CIo)/(alf)*(1-exp(xl*(alf))) 
    else:
        pso0 = (ks*CIs + ko*CIo)*xl - np.sqrt(ks*CIs*ko*CIo)*xl    # [nl+1]  factor for correlation of Ps and Po 
    return pso0 

def A_BRFv2_single_hemi(pars, lai):
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
    
    sum_pL   = 0
    sum_pL_f = 0
    sum_pL_g = 0
    
    for i in range(len(pars)):
        [tts, tta, psi, ks, ko, sob, sof, CIs, CIo, Ps_arr, Po_arr, int_res_arr, nl] = pars[i]
        [kc, kg] = sunshade(tts, tta, psi, ks, ko, CIs, CIo, Ps_arr, Po_arr, int_res_arr, nl, lai)
        
        sum_pL   += ww[i]*sob*kc/ko/pi
        sum_pL_f += ww[i]*sof*kc/ko/pi
        sum_pL_g += ww[i]*kg/pi              
        
        
    mu_tL  = np.cos(neword_tL)
    sin_tL = np.sin(neword_tL)
    
    sum_tL   = sum(ww* mu_tL*sin_tL*sum_pL)  *conv1_pL
    sum_tL_f = sum(ww* mu_tL*sin_tL*sum_pL_f)*conv1_pL
    sum_tL_g = sum(ww* mu_tL*sin_tL*sum_pL_g)*conv1_pL
    
    sob_vsla = sum_tL
    sof_vsla = sum_tL_f
    kgd      = sum_tL_g
    
    return sob_vsla, sof_vsla, kgd

def A_BRFv2_single_dif(pars, lai):
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

    sum_pL   = 0
    sum_pL_f = 0
    sum_pL_g = 0
   
    for i in range(len(pars)):                
        [tta, tto, psi, ks, ko, sob, sof, CIa, CIo, Ps_arr, Po_arr, int_res_arr, nl] = pars[i]
        [kc, kg] = sunshade(tta, tto, psi, ks, ko, CIa, CIo, Ps_arr, Po_arr, int_res_arr, nl, lai)
       
        sum_pL   += ww[i]*sob*kc/ko/pi
        sum_pL_f += ww[i]*sof*kc/ko/pi
        sum_pL_g += ww[i]*kg/pi                        
        
    mu_tL  = np.cos(neword_tL)
    sin_tL = np.sin(neword_tL)
 
    sum_tL   = sum(ww* mu_tL*sin_tL*sum_pL)  *conv1_pL
    sum_tL_f = sum(ww* mu_tL*sin_tL*sum_pL_f)*conv1_pL
    sum_tL_g = sum(ww* mu_tL*sin_tL*sum_pL_g)*conv1_pL
    
    sob_vsla = sum_tL
    sof_vsla = sum_tL_f
    kgd      = sum_tL_g
    
    return sob_vsla, sof_vsla, kgd
        
def A_BRFv2_single_hemi_dif(pars, lai):    
    xx=np.array([0.9602898565, -0.9602898565, 0.7966664774, -0.7966664774, 0.5255324099, -0.5255324099, 0.1834346425, -0.1834346425])
    
    ww=np.array([0.1012285363,  0.1012285363, 0.2223810345,  0.2223810345, 0.3137066459,  0.3137066459, 0.3626837834,  0.3626837834])   

    # * define limits of integration and the convertion factors for integration
    # * over thetaL (note the tL suffix!)
    upperlimit_mL = pi/2.0
    lowerlimit_mL = 0.0
    conv1_mL = (upperlimit_mL-lowerlimit_mL)/2.0
    conv2_mL = (upperlimit_mL+lowerlimit_mL)/2.0
    neword_mL = conv1_mL*xx + conv2_mL   
    mu_mL  = np.cos(neword_mL)
    sin_mL = np.sin(neword_mL)
    
    sum_mL = 0.0
    sum_mL_f = 0.0
    sum_mL_g = 0.0
        
    #   * define limits of integration and the convertion factors for integration
    # * over phiL (note the pL suffix!)
    upperlimit_nL = 2.0*pi
    lowerlimit_nL = 0.0
    conv1_nL = (upperlimit_nL-lowerlimit_nL)/2.0
    conv2_nL = (upperlimit_nL+lowerlimit_nL)/2.0
    neword_nL = conv1_nL*xx + conv2_nL   

    sum_nL = 0.0
    sum_nL_f = 0.0
    sum_nL_g = 0.0
    
    # * define limits of integration and the convertion factors for integration
    # * over thetaL (note the tL suffix!)
    upperlimit_tL = np.pi/2.0
    lowerlimit_tL = 0.0
    conv1_tL = (upperlimit_tL-lowerlimit_tL)/2.0
    conv2_tL = (upperlimit_tL+lowerlimit_tL)/2.0    
    neword_tL = conv1_tL*xx + conv2_tL   
    mu_tL  = np.cos(neword_tL)
    sin_tL = np.sin(neword_tL)

    sum_tL = 0.0
    sum_tL_f = 0.0
    sum_tL_g = 0.0
    
    # * define limits of integration and the convertion factors for integration
    # * over phiL (note the pL suffix!)
    upperlimit_pL = 2.0*pi
    lowerlimit_pL = 0.0
    conv1_pL = (upperlimit_pL-lowerlimit_pL)/2.0
    conv2_pL = (upperlimit_pL+lowerlimit_pL)/2.0  
    neword_pL  = conv1_pL*xx + conv2_pL

    sum_pL = 0.0
    sum_pL_f = 0.0
    sum_pL_g = 0.0
      
    for i in range(len(pars)):        
        [tts, tto, psi, ks, ko, sob, sof, CIs, CIo, Ps_arr, Po_arr, int_res_arr, nl] = pars[i]       
        [kca, kga] = sunshade(tts, tto, psi, ks, ko, CIs, CIo, Ps_arr, Po_arr, int_res_arr, nl, lai)
        
        sum_pL   += ww[i]*sob*kca/ko/pi
        sum_pL_f += ww[i]*sof*kca/ko/pi
        sum_pL_g += ww[i]*kga/pi           
            
    sum_pL = sum_pL*conv1_pL
    sum_tL = sum_tL + sum(ww* mu_tL*sin_tL*sum_pL)

    sum_pL_f = sum_pL_f*conv1_pL;
    sum_tL_f = sum_tL_f + sum(ww* mu_tL*sin_tL*sum_pL_f)
            
    sum_pL_g = sum_pL_g*conv1_pL
    sum_tL_g = sum_tL_g + sum(ww* mu_tL*sin_tL*sum_pL_g)
    
    sum_tL = sum_tL*conv1_tL;
    sum_nL = sum_nL + sum(ww*sum_tL/pi)

    sum_tL_f = sum_tL_f*conv1_tL    
    sum_nL_f = sum_nL_f + sum(ww*sum_tL_f/pi)
    
    sum_tL_g = sum_tL_g*conv1_tL    
    sum_nL_g = sum_nL_g + sum(ww*sum_tL_g/pi)
    
    sum_nL = sum_nL*conv1_nL;    
    sum_mL = sum_mL + sum(ww*mu_mL*sin_mL*sum_nL)

    sum_nL_f = sum_nL_f*conv1_nL;    
    sum_mL_f = sum_mL_f + sum(ww* mu_mL*sin_mL*sum_nL_f)
    
    sum_nL_g = sum_nL_g*conv1_nL;    
    sum_mL_g = sum_mL_g + sum(ww* mu_mL*sin_mL*sum_nL_g)    
    
    sum_mL = sum_mL*conv1_mL
    sum_mL_f = sum_mL_f*conv1_mL    
    sum_mL_g = sum_mL_g*conv1_mL

    sob_vsla = sum_mL
    sof_vsla = sum_mL_f
    kgd_dif  = sum_mL_g
    
    return sob_vsla, sof_vsla, kgd_dif

def i_hemi(CIy1,CIy2,lai,lidf):
    xx=np.array([0.9602898565, -0.9602898565, 0.7966664774, -0.7966664774, 0.5255324099, -0.5255324099, 0.1834346425, -0.1834346425])
    
    ww=np.array([0.1012285363,  0.1012285363, 0.2223810345,  0.2223810345, 0.3137066459,  0.3137066459, 0.3626837834,  0.3626837834])   
    
    # * define limits of integration and the convertion factors for integration
    # * over thetaL (note the tL suffix!)
    upperlimit_tL = np.pi/2.0
    lowerlimit_tL = 0.0
    conv1_tL = (upperlimit_tL-lowerlimit_tL)/2.0
    conv2_tL = (upperlimit_tL+lowerlimit_tL)/2.0        

    sum_tL = 0
    
    for i in range(len(ww)):
        
        neword_tL = conv1_tL*xx[i] + conv2_tL
        mu_tL     = np.cos(neword_tL)
        sin_tL    = np.sin(neword_tL)

        tta  =  neword_tL*180/pi    # observer zenith angle
            
        Ga,ka  = weighted_sum_over_lidf2(tta,lidf)
        
        CIa = CIxy(CIy1,CIy2,tta)
        
        ia=1-np.exp(-ka*lai*CIa)

        sum_tL = sum_tL + ww[i]* mu_tL*sin_tL*ia*2

    sum_tL = sum_tL*conv1_tL    
    return sum_tL
    
def weighted_sum_over_lidf2(tts, lidf):  
    
    litab   = np.array([5.,15.,25.,35.,45.,55.,65.,75.,81.,83.,85.,87.,89.]).T   
    cos_tts = np.cos(np.radians(tts))    
    chi_s   = volscat2(tts,litab)    
    ksli    = chi_s/cos_tts
    k       = np.dot(lidf,ksli)
    Gs      = k*cos_tts 
    return Gs,k

def volscat2(tts, ttli):
    #tts    [1]         Sun            zenith angle in degrees
    #tto    [1]         Observation    zenith angle in degrees
    #psi    [1]         Difference of  azimuth angle between solar and viewing position
    #ttli   [ttli]      leaf inclination array
    
    cos_ttli = np.cos(np.radians(ttli))                #   cosine of normal of upperside of leaf
    sin_ttli = np.sin(np.radians(ttli))                #   sine   of normal of upperside of leaf
    
    cos_tts = np.cos(np.radians(tts))                 #   cosine of sun zenith angle
    sin_tts = np.sin(np.radians(tts))                 #   sine   of sun zenith angle

    Cs = cos_ttli*cos_tts                 
    Ss = sin_ttli*sin_tts

    #As = max([Ss,Cs],[],2)                                   
    As = np.max(np.hstack((Cs.reshape(-1,1), Ss.reshape(-1,1))), axis=1)
    bts = np.arccos(-Cs/As)                   
     
    chi_s = 2/pi*((bts-pi/2)*Cs + np.sin(bts)*Ss) 
    return chi_s
    
    
    
    
    
    
    
    