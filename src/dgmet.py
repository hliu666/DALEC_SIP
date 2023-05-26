# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 15:49:26 2021

@author: Administrator
"""

import os
import netCDF4 as nc
import numpy as np
import xarray as xr
import pandas as pd
import datetime

def daily_nc(root, date):    
    year, month, day = date.year, date.month, date.day
    site_pd = pd.read_csv(os.path.join(root, "DALEC_F\pars\SiteInfo_F.csv"))  
    #对每个站点进行逐一的遍历
    for s in range(0, len(site_pd)):
        siteName = site_pd.iloc[s,0] #获取站点名称
        siteID = site_pd.iloc[s,1]  
        
        start = site_pd.iloc[s,2]   
        end = site_pd.iloc[s,3]
        
        lat = site_pd.iloc[s,4]
        lon = site_pd.iloc[s,5]
        
        lma = site_pd.iloc[s,6]
        cab = site_pd.iloc[s,7]
        clab = site_pd.iloc[s,8]
        onset_k = site_pd.iloc[s,9]
        onset_b = site_pd.iloc[s,10]

        flag = 0 
        for y in range(start, end+1):
            if y < 2021:                
                histPath = os.path.join(root, "AWS\gridMET\hist\\00_20")                
                tmaxPath = os.path.join(histPath, "TMAX")
                tminPath = os.path.join(histPath, "TMIN")
                parPath = os.path.join(histPath, "PAR")
            else:
                histPath = os.path.join(root, "AWS\gridMET\hist\\21\{0}".format(date))
                tmaxPath = histPath
                tminPath = histPath
                parPath = histPath              
         
            gccPath = os.path.join(root, "AWS\MERRA-2\EC\DailyRGB") 
            
            tmaxFile = os.path.join(tmaxPath, "tmmx_{0}.nc".format(y))
            tminFile = os.path.join(tminPath, "tmmn_{0}.nc".format(y))
            parFile = os.path.join(parPath, "srad_{0}.nc".format(y))
            
            tmaxObj = xr.open_dataset(tmaxFile)  
            tminObj = xr.open_dataset(tminFile)  
            parObj = xr.open_dataset(parFile)  
    
            tmaxData = tmaxObj['air_temperature'].sel(lon=lon, lat=lat, method='nearest')
            tminData = tminObj['air_temperature'].sel(lon=lon, lat=lat, method='nearest')
            parData = parObj['surface_downwelling_shortwave_flux_in_air'].sel(lon=lon, lat=lat, method='nearest')
            
            tmax_df = tmaxData.to_dataframe()
            tmin_df = tminData.to_dataframe()
            par_df = parData.to_dataframe()
            
            tmax_df.index = tmax_df.index.to_series().dt.date
            tmin_df.index = tmin_df.index.to_series().dt.date
            par_df.index = par_df.index.to_series().dt.date
            
            tmax_df = tmax_df[tmax_df.index < date]  
            tmin_df = tmin_df[tmin_df.index < date]            
            par_df = par_df[par_df.index < date]            
            
            tmax_df.index = pd.to_datetime(tmax_df.index)
            tmin_df.index = pd.to_datetime(tmin_df.index)
            par_df.index = pd.to_datetime(par_df.index)
            
            tmax_df.loc[:, 'year'] = tmax_df.index.year
            tmax_df.loc[:, 'month'] = tmax_df.index.month
            tmax_df.loc[:, 'day'] = tmax_df.index.day
            
            tmax_arr = tmax_df.drop(['lon', 'lat'], axis=1).to_numpy()
            tmin_arr = tmin_df.drop(['lon', 'lat'], axis=1).to_numpy()    
            par_arr = par_df.drop(['lon', 'lat'], axis=1).to_numpy()
            
            length = min(len(tmax_arr), len(tmin_arr), len(par_arr))
            
            tmax_arr = tmax_arr[:length,:]
            tmin_arr = tmin_arr[:length,:]
            par_arr  = par_arr[:length,:]
            
            if y < 2021:
                gccFile =  os.path.join(gccPath, "{0}/{1}_{2}_RGB.csv".format(siteName, y, siteName))   
                gcc_df = pd.read_csv(gccFile, sep = ',', encoding = 'utf-8',header=None)
                gcc_arr = np.array(gcc_df.iloc[:,1])
            else:
                gcc_arr = np.full((len(par_arr)), np.nan)
                
            hist_arr_sub = tmax_arr[:, [1,2,3,0]]                                                     
            hist_arr_sub = np.hstack((hist_arr_sub, tmin_arr.reshape(-1, 1)))                             
            hist_arr_sub = np.hstack((hist_arr_sub, ((hist_arr_sub[:,3] + hist_arr_sub[:,4])/2).reshape(-1, 1)))  
            hist_arr_sub[:,3:6] = hist_arr_sub[:,3:6] - 273.0
            hist_arr_sub = np.hstack((hist_arr_sub, par_arr.reshape(-1, 1)))                               
            hist_arr_sub = np.hstack((hist_arr_sub, gcc_arr.reshape(-1, 1)))                               
    
            if flag == 0:
                hist_arr = hist_arr_sub
                flag = 1
            else:
                hist_arr = np.vstack((hist_arr, hist_arr_sub))
       
        forcPath = os.path.join(root, "AWS/gridMET/forc/{0}/{1}/00".format(date, siteID))
        forcFiles = os.listdir(forcPath)
        
        for file in forcFiles:
            forcDs = nc.Dataset(os.path.join(forcPath, file))  
    
            taData = np.array(forcDs['air_temperature'][:])[0][0][:-1]
            parData = np.array(forcDs['surface_downwelling_shortwave_flux_in_air'][:])[0][0][:-1]
                  
            if len(taData)%24 == 0:
                forc_arr = np.zeros((int(len(taData)/24), 8), dtype=float)
                for d in range(len(forc_arr)):
                    #计算预测时间
                    _date = datetime.datetime(year, month, day)+ datetime.timedelta(days = d)
                    
                    forc_arr[d, 0] = _date.year
                    forc_arr[d, 1] = _date.month
                    forc_arr[d, 2] = _date.day
                    
                    #计算温度
                    forc_arr[d, 3] = np.max(taData[(d*24):((d+1)*24)]) - 273.0
                    forc_arr[d, 4] = np.min(taData[(d*24):((d+1)*24)]) - 273.0
                    forc_arr[d, 5] = np.mean(taData[(d*24):((d+1)*24)]) - 273.0             
                    
                    #计算辐射
                    forc_arr[d, 6] = np.mean(parData[(d*24):((d+1)*24)])
    
                    #gcc空值
                    forc_arr[d, 7] = np.nan
                    
                #_arr = np.vstack((hist_arr, forc_arr[1:,:]))
                _arr = np.vstack((hist_arr, forc_arr))
                #生成nc文件
                Lat = np.array([lat])	
                Lon = np.array([lon])	
                LMA = np.array([lma])	
                Cab = np.array([cab])	
                Clab = np.array([clab])	
                Onset_k = np.array([onset_k])	
                Onset_b = np.array([onset_b])	
                
                out = os.path.join(root, "DALEC_F/nc/{0}/{1}/00".format(date, siteID))
                if not os.path.exists(out):
                    os.makedirs(out)
    
                outfile = os.path.join(out, "{0}.{1}_{2}".format(siteName.split("_")[0], siteName.split("_")[1], file.split("_")[-1]))
                if os.path.exists(outfile):       
                    os.remove(outfile)
                nc_w = nc.Dataset(outfile,'w',format = 'NETCDF4')   #创建一个格式为.nc文件                   
    
                interval = len(_arr)      # 两日期差距
                #确定基础变量的维度信息 相对与坐标系的各个轴(x,y,z)
                nc_w.createDimension('latitude', 1)   
                nc_w.createDimension('longitude', 1)  
                nc_w.createDimension('leaf_mass', 1)   
                nc_w.createDimension('cab', 1)  
                nc_w.createDimension('clab', 1)   
                nc_w.createDimension('onset_k', 1)  
                nc_w.createDimension('onset_b', 1) 
                
                nc_w.createDimension('year', interval)
                nc_w.createDimension('month', interval)
                nc_w.createDimension('day', interval)
                            
                nc_w.createDimension('TA', interval)
                nc_w.createDimension('TA_DAY', interval)
                nc_w.createDimension('TA_NIGHT', interval)
                nc_w.createDimension('SW_IN', interval)
                
                nc_w.createDimension('Gcc', interval)
    
                
                ##创建变量 参数依次为:'变量名称','数据类型','基础维度信息'
                nc_w.createVariable('latitude',np.float32,('latitude'))  
                nc_w.createVariable('longitude',np.float32,('longitude'))
                nc_w.createVariable('leaf_mass',np.float32,('leaf_mass'))  
                nc_w.createVariable('cab',np.float32,('cab'))
                nc_w.createVariable('clab',np.float32,('clab'))
                nc_w.createVariable('onset_k',np.float32,('onset_k'))  
                nc_w.createVariable('onset_b',np.float32,('onset_b'))
    
                nc_w.createVariable('year',np.int,('year'))  
                nc_w.createVariable('month',np.int,('month'))  
                nc_w.createVariable('day',np.int,('day'))  
                
                nc_w.createVariable('TA',np.float32,('TA'))  
                nc_w.createVariable('TA_DAY',np.float32,('TA_DAY'))
                nc_w.createVariable('TA_NIGHT',np.float32,('TA_NIGHT'))  
                nc_w.createVariable('SW_IN',np.float32,('SW_IN'))
    
                nc_w.createVariable('Gcc',np.float32,('Gcc'))
    
                #写入变量的数据 维度必须与定义的一致。
                nc_w.variables['latitude'][:] = Lat  
                nc_w.variables['longitude'][:] = Lon  
                nc_w.variables['leaf_mass'][:] = LMA  
                nc_w.variables['cab'][:] = Cab  
                nc_w.variables['clab'][:] = Clab  
                nc_w.variables['onset_k'][:] = Onset_k  
                nc_w.variables['onset_b'][:] = Onset_b  
    
                nc_w.variables['year'][:] = _arr[:, 0]  
                nc_w.variables['month'][:] = _arr[:, 1]    
                nc_w.variables['day'][:] = _arr[:, 2]  
                
                nc_w.variables['TA_DAY'][:] = _arr[:, 3]    
                nc_w.variables['TA_NIGHT'][:] = _arr[:, 4]    
                nc_w.variables['TA'][:] = _arr[:, 5]    
                nc_w.variables['SW_IN'][:] = _arr[:, 6]                            
    
                nc_w.variables['Gcc'][:] = _arr[:, 7]                            
                
                nc_w.close()   
                    
            else:
                print("time length Not divisible by 24, the length is {0}".format(len(taData)))
                