from datetime import date

from dalec import dalec
from SIP import SIP_Model
from BRDF import BRDF_initial

def dalecmod(lat, leaf_Pars, time_Pars, meteo_Pars, RT_Pars, acm_constants):

    [lma, clab, onset_k, onset_b] = leaf_Pars
    [cab, min_gcc] = RT_Pars  
    [year, julianday] = time_Pars        
    [tmean, tmax, tmin, sw, gdd] = meteo_Pars 
    [cfol, clab, croo, cwoo, clit, csom] =  [0, clab, 14.897635, 3116.582, 61.11097, 28998.29]
    
    paraDict = {
        'auto_frac': 0.45,
        'fol_frac': 0.01,
        'roo_frac': 0.457,
        'theta_min':1.1e-5,
        'leaf_loss':1/3,
        'theta_woo':4.8e-5,
        'theta_roo':6.72e-3,
        'theta_lit':0.024,
        'theta_som':2.4e-5,
        'temp_expf':0.0193,
        'canopy_eff':90,
        'onset_day':140,
        'lab_frac':0.9288,
        'onset_per':27,
        'fall_day':308,
        'fall_per':35,
        'leaf_mass':lma,
        'acm_constants':acm_constants  
    }
    
    cu = 0
    d_onset = 0
    steps = len(tmean)     
    gcc_List, laiList = [], []

    model = dalec(paraDict)
    initial_cpool = [cfol, clab, croo, cwoo, clit, csom]
    
    band_Pars = [[170,350],[140,170],[60,100]] #Red, Green, Blue
    
    SIP_Pars = SIP_Model(cab)
    BRDF_Pars,SUN_Pars,BRF_Pars = BRDF_initial(min_gcc, band_Pars)
 
    for i in range(steps):
        
        d0 = date(year[0], 1, 1)
        d1 = date(year[i], 1, 1)
        d2 = date(year[i], 12, 31)
        
        sd = (d1 - d0).days
        ed = (d2 - d0).days
          
        time_i_Pars = [julianday[i], year[i]]
        pheOn_Pars = [cu, d_onset, gdd[sd:ed+1], onset_k, onset_b]
        meteo_i_Pars = [tmean[i], tmax[i], tmin[i], sw[i]]
        initial_cpool, cu, d_onset, gcc, lai = model.dalecv2(lat, initial_cpool, time_i_Pars, meteo_i_Pars, pheOn_Pars, SIP_Pars, BRDF_Pars, SUN_Pars, BRF_Pars, band_Pars)
        laiList.append(lai)
        gcc_List.append(gcc)
        
    return initial_cpool, gcc_List, laiList
