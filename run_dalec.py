# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 11:55:27 2021

@author: Haoran
"""

from dalec_model import dalecmod
from datetime import datetime

import xarray
import numpy as np
import pandas as pd
import os

def Date2Juday(date:tuple)->int:   
    # nonlear_year = [31,28,31,30,31,30,31,31,30,31,30,31]
    nonleap_year = [0,31,59,90,120,151,181,212,243,273,304,334]
    leap_year    = [0,31,60,91,121,152,182,213,244,274,305,335] 
    return leap_year[date[1]-1]+date[2] if (date[0]%4)==0 else nonleap_year[date[1]-1]+date[2]

def cal_gdd(year, julianday, tmean):
    Tc = -3.32 #Single base temperature        
    gdd = np.array([0])
    #初始化一个数组计算GDD
    for i in range(0, len(year)):
        if julianday[i] >= 74:
            if i == len(year) - 1:
                gdd = np.append(gdd, gdd[i] + max(tmean[i] - Tc, 0))             
            elif year[i] == year[i + 1]:
                gdd = np.append(gdd, gdd[i] + max(tmean[i] - Tc, 0))    
            else:
                gdd = np.append(gdd, max(tmean[i] - Tc, 0))
        else:
            gdd = np.append(gdd, 0)
    gdd = np.delete(gdd, 0)
    return gdd
    
def run_dalec(root, date):    
    #acmc constant 
    acm_constants = [0.0156935, 4.22273, 208.868, 0.0453194, 0.37836, 7.19298, 0.011136, 2.1001, 0.789798] 
    
    site_pd = pd.read_csv(os.path.join(root, "DALEC_F\pars\SiteInfo_F.csv"))  
    siteLists = []
    #对每个站点进行逐一的遍历
    for s in range(0, len(site_pd)):
        siteEns = os.path.join(root, "DALEC_F/nc/{0}/{1}/00".format(date, site_pd.iloc[s,1]))
        #if not site_pd.iloc[s,1] == "DELA":
        #    continue
        ensLists = [] #   
        ens = os.listdir(siteEns)
        for e in range(0, len(ens)):
            file = os.path.join(siteEns, ens[e])

            ds = xarray.open_dataset(file)
            lat = ds['latitude'].data
            
            """
            time-Parameters:
            """
            dates = zip(ds['year'].data.astype(int),ds['month'].data.astype(int),ds['day'].data.astype(int))
            julianday = []
            dateList = []
            for d in dates:
                julianday.append(Date2Juday(d))
                dateList.append(datetime(d[0], d[1], d[2]))
            year  = ds['year'].data  
            time_Pars = [year, julianday]
                            
            """
            leaf-Parameters:
            """        
            lma = ds['leaf_mass'].data
            clab = ds['clab'].data
            onset_k = ds['onset_k'].data
            onset_b = ds['onset_b'].data           
            leaf_Pars = [lma, clab, onset_k, onset_b]
            
            """
            radiative-transfer-Parameters:
            """          
            cab = ds['cab'].data        
            gcc = ds['Gcc'].data
            if len(gcc) > 365:
                min_gcc = np.nanmin(gcc[380:390])
            else:
                min_gcc = np.nanmin(gcc[20:30])          
            RT_Pars  = [cab, min_gcc]
            
            """
            meteorological-Parameters:
            """         
            tmean = ds['TA'].data
            tmax  = ds['TA_DAY'].data
            tmin  = ds['TA_NIGHT'].data
            sw    = ds['SW_IN'].data
            gdd   = cal_gdd(year, julianday, tmean) # gdd: array  
            meteo_Pars = [tmean, tmax, tmin, sw, gdd]

            """
            calculate the average values of ensembles 
            """
            merge = np.full((len(tmean)), 0.0)             
            count = np.full((len(tmean)), 0.0) 
               
            initial_cpool, gccList, laiList = dalecmod(lat, leaf_Pars, time_Pars, meteo_Pars, RT_Pars, acm_constants)
                        
            if ens[e].split(".")[0][-2:] == "00":
                ensLists.append(gccList[-16:]+ [np.nan]*19) 
                merge = merge + np.array(gccList+ [0]*19)
                count = count + np.array([1]*len(gccList)+ [0]*19)
                
            else:
                ensLists.append(gccList[-35:])
                merge = merge + np.array(gccList)
                count = count + np.array([1]*len(gccList))               

        from plot import plot_gcc, plot_forecast
        gccMerge = merge/count
        plot_gcc(site_pd.iloc[s,1], dateList, gcc, gccMerge)          
        plot_forecast(site_pd.iloc[s,1], dateList[-35:], gccMerge[-35:])
                    
        siteLists.append([site_pd.iloc[s,1], ensLists])
    
    #输出标准文件
    from forc_file import forcFile
    forcFile(date, siteLists)

root = r"D:\DELAC_Forests"
today = datetime(2021, 7, 1)
run_dalec(root, today)  
            