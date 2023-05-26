# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 18:57:57 2021

@author: Haoran
"""
import numpy as np
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt

plt.rc('font',family='Times New Roman')

def plot_forecast(siteName, dateList, gccSim):
    Gcc_s = np.array(gccSim)
    Gcc_df = pd.DataFrame(Gcc_s, index = dateList)  
    Gcc_df = Gcc_df.dropna()

    fig = plt.figure(figsize=(8,6))
    ax = plt.subplot(1,1,1)
    ax.plot(Gcc_df.index, Gcc_df.iloc[:,0],color='black',linestyle='solid')
    
    ax.spines['left'].set_linewidth(0.7)
    ax.spines['right'].set_linewidth(0.7)
    ax.spines['top'].set_linewidth(0.7)
    ax.spines['bottom'].set_linewidth(0.7)
    ax.tick_params(axis='both', length = 4, width = 1, labelsize = 4.5)    

    ax.set_title("{0}".format(siteName), fontsize = 20, family="Times New Roman") 
    ax.tick_params(axis='x', labelrotation = 30)  
    #ax.set_xlabel('Number of Observation Points', fontsize = 20, family="Times New Roman")    
    ax.set_ylabel('Reflectance of GCC', fontsize = 20, family="Times New Roman")
    
    
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    #plt.locator_params(axis='x', nbins=5)
    
    sinPath = os.path.join(os.path.split(os.path.abspath(sys.argv[0]))[0][:-4], "res\plot")
    fig.tight_layout() 
    fig.savefig(os.path.join(sinPath, '{0} of GCC_forecast.jpg'.format(siteName)),dpi=300) 
    
def plot_gcc(siteName, dateList, gccObs, gccSim):

    Gcc_o = np.array(gccObs)
    Gcc_s = np.array(gccSim)

    Gcc = np.hstack((Gcc_s.reshape(-1,1), Gcc_o.reshape(-1,1)))
    #Gcc = np.delete(Gcc, np.where(np.isnan(Gcc)), axis = 0)

    #print(np.asarray(dateList).shape, Gcc.shape)
    
    Gcc_df = pd.DataFrame(Gcc, index = dateList)  
    #Gcc_df = Gcc_df.dropna()
    
    fig = plt.figure(figsize=(8,6))
    ax = plt.subplot(1,1,1)
    ax.plot(Gcc_df.index, Gcc_df.iloc[:,0],color='black',linestyle='solid')
    ax.plot(Gcc_df.index, Gcc_df.iloc[:,1],'r.',markersize=3, label='Observation data')
    
    ax.spines['left'].set_linewidth(0.7)
    ax.spines['right'].set_linewidth(0.7)
    ax.spines['top'].set_linewidth(0.7)
    ax.spines['bottom'].set_linewidth(0.7)
    ax.tick_params(axis='both', length = 1.5, width = 0.8, labelsize = 4.5)    

    ax.set_title("{0}".format(siteName), fontsize = 20, family="Times New Roman") 
    ax.tick_params(axis='x', labelrotation = 30)  
    #ax.set_xlabel('Number of Observation Points', fontsize = 20, family="Times New Roman")    
    ax.set_ylabel('Reflectance of GCC', fontsize = 20, family="Times New Roman")
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(loc='upper right')
    sinPath = os.path.join(os.path.split(os.path.abspath(sys.argv[0]))[0][:-4], "res\plot")
    fig.tight_layout() 
    fig.savefig(os.path.join(sinPath, '{0} of GCC.jpg'.format(siteName)),dpi=300)   
