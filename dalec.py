import dalec_func as func
import numpy as np
from BRDF import BRDF

class dalec(object):

    def __init__(self, parameters:dict=None)->None:
      
        self.autof  = None
        self.labf   = None
        self.folf   = None
        self.roof   = None
        self.woot   = None
        self.mint   = None
        self.root   = None
        self.somt   = None
        self.litt   = None
        self.tempf  = None
        self.don    = None
        self.dfall  = None
        self.ce     = None
        self.pon    = None
        self.pfall  = None
        self.lma    = None
        self.loss   = None
        self.acmc   = None
        
        self.alias = {
            'auto_frac'     : 'autof',
            'lab_frac'      : 'labf',
            'fol_frac'      : 'folf',
            'roo_frac'      : 'roof',
            'theta_woo'     : 'woot',
            'theta_min'     : 'mint',
            'theta_roo'     : 'root',
            'theta_som'     : 'somt',
            'theta_lit'     : 'litt',
            'temp_expf'     : 'tempf',
            'onset_day'     : 'don',
            'fall_day'      : 'dfall',
            'canopy_eff'    : 'ce',
            'leaf_mass'     : 'lma',
            'leaf_loss'     : 'loss',
            'onset_per'     : 'pon',
            'fall_per'      : 'pfall',
            'acm_constants' : 'acmc'
        }
              
        assert len(parameters) == 18,'The input dictionary has more or less items than the module needed'
        for key, value in parameters.items():
            if self.alias.__contains__(key):
                setattr(self, self.alias[key], value)

    def dalecv2(self, lat, inital_state, time_Pars, meteo_Pars, pheOn_Pars, SIP_Pars, BRDF_Pars, SUN_Pars, BRF_Pars, band_Pars, water_stress=1):
        """
        inital_state, 六个碳库的初始状态
        julian_day,   一年中的第几天
        soil_moisture,土壤含水量
        water_threshold,水分胁迫时的土壤含水量
        t_mean，平均气温
        t_max，最高气温
        acm_constants，acm模型的9个常数
        wp_diff，冠层和土壤的最大水势差
        i，太阳辐射
        ca，大气二氧化碳浓度
        hydro_resist，水分从土壤到冠层遇到的阻力
        lat，纬度
        paraList，DALEC模型的参数
        """  
        [julian_day, year] = time_Pars        
        [t_mean, t_max, t_min, i] = meteo_Pars 
         
        ca = 390.  # Carbon dioxide concentration
        
        #prepare inputs to the calculation function       
        cfol, clab, croo, cwoo, clit, csom = inital_state
        i = i * 86400 / 1E6
        
        lai  = cfol / self.lma
        if lai <= 0.001:
            lai = 0.001
            
        gcc, fPAR = BRDF(lai, year, julian_day, SIP_Pars, BRDF_Pars, SUN_Pars, BRF_Pars, band_Pars)
        gpp = func.acm(cfol, self.lma, self.ce, t_max, t_max-t_min, i, ca, julian_day, lat, self.acmc, fPAR)*water_stress
        #phi_onset = func.onset(julian_day, self.don, self.pon)
        phi_onset, cu, d_onset = func.onset(year, julian_day, self.pon, t_mean, pheOn_Pars)
        
        leap = 0        
        if (year % 4 == 0) and (year % 100 != 0):
            leap = 1
        elif year % 400 == 0:
            leap = 1

        if leap == 0 and julian_day == 365:
            cu = 0
        if leap == 1 and julian_day == 366:
            cu = 0
            
        phi_fall  = func.fall(julian_day, self.dfall, self.pfall, 1/self.loss)
        temp      = np.exp(t_mean * self.tempf)

        # calculate the six carbom pools after a time step (one day)
        clab_next = (1 - phi_onset) * clab + (1 - self.autof) * (1 - self.folf) * self.labf * gpp
        cfol_next = (1 - phi_fall)  * cfol + (1 - self.autof) * self.folf * gpp + phi_onset * clab
        croo_next = (1 - self.root) * croo + (1 - self.autof) * (1 - self.folf) * (1 - self.labf) * self.roof * gpp
        cwoo_next = (1 - self.woot) * cwoo + (1 - self.autof) * (1 - self.folf) * (1 - self.labf) * (1 - self.roof) * gpp
        clit_next = (1 - (self.litt+self.mint)*temp)*clit + self.root*croo + phi_fall*cfol
        csom_next = (1 - self.somt*temp) *csom + self.woot*cwoo + self.litt*temp*clit
        
        if julian_day > 300 and julian_day < 366 and cfol_next != 0:
            cfol_next = 0 
           
        # variables that can help to calibrate the dalec mode
        rtot = self.autof *gpp + (self.litt*clit + self.somt*csom)*temp
        nee  = -(1. - self.autof) * gpp + rtot
        
        return [cfol_next,clab_next,croo_next,cwoo_next,clit_next,csom_next], cu, d_onset, gcc, lai