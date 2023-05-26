# -*- coding: utf-8 -*-
"""
Created on Tue May 18 16:14:29 2021

@author: Haoran
"""
import numpy as np
import pandas as pd
import datetime
import os, sys

def forcFile(date, siteLists):
    time = []
    
    siteID = []
    obs_flag = []
    forecast = []    
    data_assimilation = []

    statistic = []
    gcc_90 = []
    
    year, month, day = date.year, date.month, date.day
    for i in range(len(siteLists[0][1][0])):      #time
        for j in range(len(siteLists)):           #site
            ensl = []
            for k in range(len(siteLists[j][1])): #ensembles
                ensl.append(siteLists[j][1][k][i])
            
            time.extend([datetime.date(year,month,day)+ datetime.timedelta(days = i)]*2)
            
            siteID.extend([siteLists[j][0]]*2)
            obs_flag.extend([2]*2)
            forecast.extend([1]*2) 
            data_assimilation.extend([1]*2)
            
            gcc_90.append(np.mean(ensl))
            statistic.append('mean')
            
            gcc_90.append(np.std(ensl))
            statistic.append('sd')
            
    data = {'time':time,
            'siteID':siteID,
            'obs_flag':obs_flag,
            'forecast':forecast,
            'data_assimilation':data_assimilation,
            'statistic':statistic,
            'gcc_90':gcc_90            
            }
    df = pd.DataFrame(data)
    root = os.path.split(os.path.abspath(sys.argv[0]))[0][:-4]
    outPath = os.path.join(root, 'res\phenology-{0}-DALEC_SIP.csv'.format(date))
    df.to_csv(outPath, sep=',', index = False, header = True)   
