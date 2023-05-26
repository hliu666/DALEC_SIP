# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 01:41:24 2021

@author: hliu
"""
import urllib
import requests
import datetime
import os

import warnings
warnings.filterwarnings('ignore')

from dgmet import daily_nc
from run_dalec import run_dalec

"""  
sdt = datetime.date(2021, 6, 1) 
edt = datetime.date(2021, 6, 20) 
for i in range((edt - sdt).days+1):
    today = sdt + datetime.timedelta(days=i)
    print (today)
    root =  os.path.dirname(os.path.dirname(os.getcwd()))
    scodes = ["BART", "CLBJ", "DELA", "GRSM", "HARV", "STEI", "UKFS"]
    for s in range(len(scodes)):
        outforc = os.path.join(root, "AWS/gridMET//forc/{0}/{1}/00".format(today, scodes[s]))   
        if not os.path.exists(outforc):
            os.makedirs(outforc)         
        for en in range(0, 31):
            if en == 0:
                future = today + datetime.timedelta(16)
            else:
                future = today + datetime.timedelta(35)
             
            url = "https://data.ecoforecast.org/minio/download/drivers/noaa/NOAAGEFS_1hr/{0}/{1}/00/NOAAGEFS_1hr_{2}_{3}T00_{4}T00_ens{5}.nc?token=".format(scodes[s], today, scodes[s], today, future, str(en).zfill(2))
            forcfile = os.path.join(outforc, "NOAAGEFS_1hr_{0}_{1}T00_{2}T00_ens{3}.nc".format(scodes[s], today, future, en))
            opener=urllib.request.build_opener()
            opener.addheaders=[('User-Agent','Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/36.0.1941.0 Safari/537.36')]
            urllib.request.install_opener(opener)
            urllib.request.urlretrieve(url, forcfile)
    
    outhist = os.path.join(root, "AWS\gridMET\hist\\21\{0}".format(today))
    if not os.path.exists(outhist):
        os.makedirs(outhist) 
    meteos = ['pr', 'srad', 'tmmn', 'tmmx']  
    for m in range(len(meteos)):   
        url = "https://www.northwestknowledge.net/metdata/data/{0}_2021.nc".format(meteos[m])
        histfile = os.path.join(outhist, "{0}_2021.nc".format(meteos[m]))    
        r = requests.get(url)
        with open (histfile, "wb") as code:
            code.write(r.content)
    daily_nc(root, today)                 
    run_dalec(root, today) 
    
"""    
root =  os.path.dirname(os.path.dirname(os.getcwd()))
today = datetime.date.today() + datetime.timedelta(-1)
print(today)
scodes = ["BART", "CLBJ", "DELA", "GRSM", "HARV", "STEI", "UKFS"]
for s in range(len(scodes)):
    outforc = os.path.join(root, "AWS/gridMET//forc/{0}/{1}/00".format(today, scodes[s]))   
    if not os.path.exists(outforc):
        os.makedirs(outforc)         
    for en in range(0, 31):
        if en == 0:
            future = today + datetime.timedelta(16)
        else:
            future = today + datetime.timedelta(35)
         
        url = "https://data.ecoforecast.org/minio/download/drivers/noaa/NOAAGEFS_1hr/{0}/{1}/00/NOAAGEFS_1hr_{2}_{3}T00_{4}T00_ens{5}.nc?token=".format(scodes[s], today, scodes[s], today, future, str(en).zfill(2))
        forcfile = os.path.join(outforc, "NOAAGEFS_1hr_{0}_{1}T00_{2}T00_ens{3}.nc".format(scodes[s], today, future, en))
        opener=urllib.request.build_opener()
        opener.addheaders=[('User-Agent','Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/36.0.1941.0 Safari/537.36')]
        urllib.request.install_opener(opener)
        urllib.request.urlretrieve(url, forcfile)

outhist = os.path.join(root, "AWS\gridMET\hist\\21\{0}".format(today))
if not os.path.exists(outhist):
    os.makedirs(outhist) 
meteos = ['pr', 'srad', 'tmmn', 'tmmx']  
for m in range(len(meteos)):   
    url = "https://www.northwestknowledge.net/metdata/data/{0}_2021.nc".format(meteos[m])
    histfile = os.path.join(outhist, "{0}_2021.nc".format(meteos[m]))    
    r = requests.get(url)
    with open (histfile, "wb") as code:
        code.write(r.content)
       
daily_nc(root, today)
run_dalec(root, today)        
