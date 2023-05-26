import numpy as np

    
def acm(cf, clma, ceff, t_max, t_range, i, ca, julian_day, lat, acmc, fPAR):
    """ Aggregated canopy model (ACM) function
    ------------------------------------------
    Takes a foliar carbon (cf) value, leaf mass per area (clma) and canopy
    efficiency (ceff) and returns the estimated value for Gross Primary
    Productivity (gpp) of the forest at that time.
    :param cf: foliar carbon (g C m-2)
    :param clma: leaf mass area (g C m-2)
    :param ceff: canopy efficiency parameter
    :return: GPP value
    """
    
    hydro_resist = 1.0
    wp_diff = 2.5
    
    t_range = 0.5 * t_range
    L = cf / clma #LAI
    q = acmc[1] - acmc[2]
    gc = (abs(wp_diff))**acmc[8] / \
         (t_range + acmc[4]*hydro_resist)
    p = ((ceff*L) / gc)*np.exp(acmc[6]*t_max)
    ci = 0.5*(ca + q - p + np.sqrt((ca + q - p)**2 - 4*(ca*q - p*acmc[1])))
    #surface bidirectional reflectance factor    
    #E0 = (acmc[5]*L**2) / (L**2 + acmc[7]) #Canopy level quantum yield
    E0 = fPAR
    i = i * 0.46
    delta = -23.4*np.cos((360.*(julian_day + 10) / 365.) *
                         (np.pi/180.))*(np.pi/180.)
    s = 24*np.arccos((- np.tan(lat*np.pi/180.)*np.tan(delta))) / np.pi
    if s >= 24.:
        s = 24.
    elif s <= 0.:
        s = 0.
    else:
        s = s
        
    gpp = (E0*i*gc*(ca - ci))*(acmc[0]*s + acmc[3]) / (E0*i + gc*(ca - ci))
    return gpp
    
def fit_polynomial(ep, mult_fac):
    """ Polynomial used to find phi_f and phi (offset terms used in
    phi_onset and phi_fall), given an evaluation point for the polynomial
    and a multiplication term.
    :param ep: evaluation point
    :param mult_fac: multiplication term
    :return: fitted polynomial value ---- https://github.com/Ewan82/dalec2/blob/master/src/model/mod_class.py
    """
    cf = [2.359978471e-05, 0.000332730053021, 0.000901865258885,
          -0.005437736864888, -0.020836027517787, 0.126972018064287,
          -0.188459767342504]
    poly_val = cf[0]*ep**6 + cf[1]*ep**5 + cf[2]*ep**4 + cf[3]*ep**3 + cf[4]*ep**2 + \
        cf[5]*ep**1 + cf[6]*ep**0
    phi = poly_val*mult_fac
    return phi

"""
def onset(t, d_onset, c_ronset):
    s = 365.25 / np.pi
    release_coeff = np.sqrt(2.)*c_ronset / 2.
    mag_coeff = (np.log(1.+1e-3) - np.log(1e-3)) / 2.
    offset = fit_polynomial(1+1e-3, release_coeff)
    phi_onset = (2. / np.sqrt(np.pi))*(mag_coeff / release_coeff) * \
        np.exp(-(np.sin((t - d_onset + offset) /s)*(s / release_coeff))**2)
    return phi_onset
"""

def onset(year, julian_day, c_ronset, t_mean, pheOn_Pars):
    """Leaf onset function (controls labile to foliar carbon transfer)
    takes d_onset value, cronset value and returns a value for phi_onset.
    """
    [cu, d_onset0, gdd, onset_k, onset_b] = pheOn_Pars 
    Tc = -3.32
    if t_mean < Tc and julian_day >= 74:
        cu = cu  + 1
    a = 207.87
    b = 244.72
    c = -0.013
    GDDcrit = a + b * np.exp(c * cu)
    if len(gdd[gdd > GDDcrit]) == 0 and year == 2021:
        d_onset = 124
    else:
        d_onset = np.where(gdd == min(gdd[gdd > GDDcrit]))[0][-1] + 1 
    """
    d_onset = np.where(gdd == min(gdd[gdd > GDDcrit]))[0][-1] + 1 
    print(year, julian_day, cu, d_onset)
    """
    
    #d_onset = d_onset + int(c_ronset*0.3924 + 20.64) #0.39 + 20
    d_onset = d_onset + int(c_ronset*onset_k + onset_b)    
    s = 365.25 / np.pi
    release_coeff = np.sqrt(2.)*c_ronset / 2.
    mag_coeff = (np.log(1.+1e-3) - np.log(1e-3)) / 2.
    offset = fit_polynomial(1+1e-3, release_coeff)
    phi_onset = (2. / np.sqrt(np.pi))*(mag_coeff / release_coeff) * \
        np.exp(-(np.sin((julian_day - d_onset + offset) /s)*(s / release_coeff))**2)
    print(phi_onset)
    return phi_onset, cu, d_onset

    
def fall(t, d_fall, crfall, clspan):
    
    """Leaf fall function (controls foliar to litter carbon transfer) takes
    d_fall value, crfall value, clspan value and returns a value for
    phi_fall. ----https://github.com/Ewan82/dalec2/blob/master/src/model/mod_class.py
    """
    s = 365.25 / np.pi
    release_coeff = np.sqrt(2.)*crfall / 2.
    mag_coeff = (np.log(clspan) - np.log(clspan - 1.)) / 2.
    offset = fit_polynomial(clspan, release_coeff)
    phi_fall = (2. / np.sqrt(np.pi))*(mag_coeff / release_coeff) * \
        np.exp(-(np.sin((t - d_fall + offset) / s) * s / release_coeff)**2)
    
    return phi_fall
