#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 11:31:56 2023

@author: eli
"""

# =============================================================================
# %% Modify in this cell
# =============================================================================
# how many model runs to plot? (1-4)
nruns = 1

# run names?
'''
Model output files have standard nomenclature, which the code will use to find
all variables automatically. Format is as follows:
RES_SPC.dat
where RES is the temporal resolution of the data (1Hz or 2min) and SPC is the
species included in the model run (all or SNA or noNVC)
'''
run1 = '2min_SNAorgMeta'
# run1 = '2min_SNAorgLevPVOLPMeta'

model1 = 'EAIM'

# which flights to plot?
'''
Each flight is plotted individually in a for loop, so this will tell the code
which flights to loop through and plot.
Must be a list of strings that match the Flight column in the obs data or 'all'
'''
# flights = 'all'
flights = ['RF02','RF03']

# plot raw 1Hz output?
'''
For model runs at 1Hz, data is averaged back to the 2-minute PILS timesteps by
default. To override this and plot the raw 1Hz data (it will be a noisy plot),
set raw_1Hz to True
'''
raw_1Hz = False

make_gif = False

# =============================================================================


## %% Import Modules
#_--------------------------------------------
# import statements
#---------------------------------------------
import pandas as pd
import numpy as np
import glob
from datetime import datetime
from scipy import stats
import subprocess # to copy text onto my clipboard (v cool!)
import os
import time

#---------------------------------------------
# plotting commands
#---------------------------------------------
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import dates
from pytz import timezone
import plotly.graph_objs as go
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap

# set figure defaults
mpl.rcParams['figure.dpi'] = 150
plt.rcParams['figure.figsize'] = (14.0/2, 8.0/2)
mpl.rcParams['lines.linewidth'] = 1.0

#%%
M_NH3  = 17.031e6 # ammonia (https://pubchem.ncbi.nlm.nih.gov/compound/222)
M_HNO3 = 63.013e6 # nitric acid (https://pubchem.ncbi.nlm.nih.gov/compound/944)
M_H2SO4= 98.08e6  # sulfuric acid (https://pubchem.ncbi.nlm.nih.gov/compound/1118)
M_Cl   = 35.45e6  # chloride ion (https://pubchem.ncbi.nlm.nih.gov/compound/312)
M_Na   = 22.989e6 # sodium ion (https://pubchem.ncbi.nlm.nih.gov/compound/5360545)
M_K    = 39.098e6 # Potassium (https://pubchem.ncbi.nlm.nih.gov/compound/813)
M_Mg   = 24.305e6 # Magnesium (https://pubchem.ncbi.nlm.nih.gov/compound/888)
M_Ca   = 40.08e6  # Calcium (https://pubchem.ncbi.nlm.nih.gov/compound/271)

#%%
## Get obs
df = pd.read_csv( './Data/TRANS2AM_Merge_1Hz.csv' )
for RF in np.unique(df['Flight']):# For each flight (because code relies on UTC since midnight)
    df.loc[(df['Flight']==RF),'datetime'] = pd.to_datetime(df.loc[(df['Flight']==RF),'Datetime_UTC'], unit='s',
                                                           origin=datetime.strptime('20%i'%np.unique(df.loc[(df['Flight']==RF),'yymmdd'])[0], '%Y%m%d'),
                                                           utc=True)    # Add datetime in UTC
df['datetime'] = pd.to_datetime(df['datetime'])# Convert object to datetime
df['datetime'] = df.datetime.dt.tz_convert('America/Denver')# Convert UTC to MT

df2 = pd.read_csv( './Data/TRANS2AM_Merge_2min.csv' )
for RF in np.unique(df2['Flight']):# For each flight (because code relies on UTC since midnight)
    df2.loc[(df2['Flight']==RF),'datetime'] = pd.to_datetime(df2.loc[(df2['Flight']==RF),'Datetime_UTC'], unit='s',
                                                           origin=datetime.strptime('20%i'%np.unique(df2.loc[(df2['Flight']==RF),'yymmdd'])[0], '%Y%m%d'),
                                                           utc=True)    # Add datetime in UTC
df2['datetime'] = pd.to_datetime(df2['datetime'])# Convert object to datetime
df2['datetime'] = df2.datetime.dt.tz_convert('America/Denver')# Convert UTC to MT
# Rewrite RH artifact with nan's
I_rh = []
for i in np.where(df.rh>=100)[0]:
    I = np.where( (df2.Flight==df.Flight[i]) &
                  (df2.Datetime_UTC_start < df.Datetime_UTC[i]) &
                  (df2.Datetime_UTC_stop >= df.Datetime_UTC[i]) )[0]
    if I not in I_rh and len(I)>0:
        I_rh.append(I[0])
        # print(I_rh)

## Get model output
x, y = [], []
for i in range(1,nruns+1):
    exec("x.append(run%i)"%i)
    exec("y.append(model%i)"%i)
counter = 0
for run_name, model in zip(x, y):
    counter+=1

    # Get model runs parameters
    res,spc = run_name.split('_')
    if res=='1Hz': dt=1
    if res=='2min': dt=120
    
    # Get obs
    df_obs = pd.read_csv( glob.glob('./Data/TRANS2AM_Merge_%s.csv'%(res))[0] )
    # Calculate missing organics based on PILS-PCASP volume concentration
    if 'PVOLP' in run_name:
        # assume sna, inorg (NVC+Cl), org densities.8*4
        p_org = 1.2
        p_nvc =  2.1
        p_inorg = 1.7

        # connversion coefficient from STP to ambient
        C = (273.15*df_obs['ps_hads_a']) / (1013.25*(df_obs['temp']+273.15))

        # for each PILS species
        for spcs in [x for x in df_obs.columns if 'ug_m3' in x]:
            # get mass conc
            exec("mass_conc = (df_obs.%s)"%spcs)

            # convert mass conc from STP to ambient
            mass_conc = mass_conc*C
            
            # get density of spcs
            if any(x in spcs for x in ['Ammonium', 'Sulfate','Nitrate']):
                density = p_inorg
            elif any(x in spcs for x in ['Potassium', 'Sodium','Magnesium','Calcium','Chloride']):
                density = p_nvc
            else:
                density = p_org
            
            # calc volume conc
            vol_conc = mass_conc / density
            
            # save volume conc
            exec( "%s = vol_conc"%spcs.replace('ug_m3','um3_cm3') )

        # sum volume conc of PILS species
        vol_conc = np.nansum([Ammonium_um3_cm3, Sulfate_um3_cm3, Nitrate_um3_cm3,
                              Potassium_um3_cm3, Propionate_um3_cm3,
                              Sodium_um3_cm3, Succinate_um3_cm3,
                              Levoglucosan_um3_cm3, Magnesium_um3_cm3, MSA_um3_cm3,
                              Nitrite_um3_cm3, Oxalate_um3_cm3, Formate_um3_cm3,
                              Glutarate_um3_cm3, Acetate_um3_cm3, Calcium_um3_cm3,
                              Chloride_um3_cm3],axis=0)
        # if volume conc is zero, set to nan (bc nansum tuns 0 into nans)
        vol_conc[np.where(vol_conc==0)]=np.nan
        
        # get PCASP PM3 volume concs
        obs = df_obs.PVOLP_OBR

        # calculate missing volume
        misvol = obs - vol_conc
        # where misvol is negative (i.e., PILS higher than PCASP), set to nan
        misvol[misvol<0] = np.nan
        
        # convert missing volume to mass conc
        mismas = misvol * p_org
        # convert missing mass conc from ambient to STP
        mismas = mismas / C
        
        # add to dataframe
        df_obs['PVOLP_ug_m3'] = mismas
    
    # Get model output
    if model=='ISORROPIA':
        headers=['NATOT','SO4TOT','NH4TOT','NO3TOT','CLTOT','CATOT','KTOT',
               'MGTOT','RH','TEMP','GNH3','GHCL','GHNO3','CNACL','CNANO3',
               'CNA2SO4','CNAHSO4','CNH4CL','CNH4NO3','CNH42S4','CNH4HS4',
               'CLC','CCASO4','CCANO32','CCACL2','CK2SO4','CKHSO4','CKNO3',
               'CKCL','CMGSO4','CMGNO32','CMGCL2','HLIQ','NALIQ','NH4LIQ',
               'CLLIQ','SO4LIQ','HSO4LIQ','NO3LIQ','CaLIQ','KLIQ','MgLIQ',
               'NH4AER','CLAER','NO3AER','WATER','LMASS','SMASS']#,'CASE']
        df_mod = pd.read_table('./ISORROPIA/outputs/%s.dat'%(run_name),
                               header=0,sep='\s+',names=headers,comment='"')

    if model=='EAIM':
        # Get model run parameters
        runs = len(df_obs) # how many total runs
        batches = int(np.ceil(runs/60)) # in how many batches (since online E-AIM can only run 100 max at a time)
        runs_in_batch = [60]*(batches-1) # how many runs per batch
        runs_in_batch.append(runs%60)

        # Get data
        # input
        for batch in range(batches):
            # print('%i of %i batches'%(batch+1, batches))
            fid = './EAIM/input/%s/input_%iof%i_%s.dat'%(run_name,batch+1,batches,run_name)
            tmp = pd.read_table(fid,header=None,skiprows=1,sep='\s+')
            if batch==0:
                df_in = tmp
            else:
                df_in = pd.concat([df_in,tmp]).reset_index(drop=True)
        if 'org' not in run_name:
            if 'SNA' in run_name:
                df_in = df_in[[5, 7, 9, 10]]
                df_in.columns = ['RH', 'NH4+', 'SO4=', 'NO3-']
            else:
                df_in = df_in[[5, 7, 8, 9, 10, 11]]
                df_in.columns = ['RH', 'NH4+', 'Na+', 'SO4=', 'NO3-','Cl-']
                df_in['Na+']*=M_Na
                df_in['Cl-']*=M_Cl
        else:
            if 'SNA' in run_name:
                df_in = df_in[[5, 7, 9, 10]]
                df_in.columns = ['RH', 'NH4+', 'SO4=', 'NO3-']
            else:
                df_in = df_in[[5, 7, 8, 9, 10, 11]]
                df_in.columns = ['RH', 'NH4+', 'Na+', 'SO4=', 'NO3-','Cl-']
                df_in['Na+']*=M_Na
                df_in['Cl-']*=M_Cl
        df_in['RH']*=100
        df_in['NH4+']*=(M_NH3)
        df_in['SO4=']*=(M_H2SO4)
        df_in['NO3-']*=(M_HNO3)

        # output
        for batch in range(batches):
            # print('%i of %i batches'%(batch+1, batches))
            fid = './EAIM/output/%s/output_%iof%i_%s.dat'%(run_name,batch+1,batches,run_name)
            if 'org' in run_name: nskip = 17+2*runs_in_batch[batch]
            if 'org' not in run_name: nskip = 10+runs_in_batch[batch]
            tmp = pd.read_fwf(fid, skiprows=nskip, delimiter=' ', infer_nrows=runs_in_batch[batch])[:runs_in_batch[batch]]
            if batch>0:
                df_mod = pd.concat([df_mod,tmp]).reset_index(drop=True)
            else:
                df_mod = tmp
        pH = -np.log10(df_mod['m_H(aq)'].astype(float))

        
        for batch in range(batches):
            # print('%i of %i batches'%(batch+1, batches))
            fid = './EAIM/output/%s/output_%iof%i_%s.dat'%(run_name,batch+1,batches,run_name)
            tmp = pd.read_fwf(fid, skiprows=3, delimiter=' ', infer_nrows=runs_in_batch[batch])[:runs_in_batch[batch]]
            if batch>0:
                df_mod = pd.concat([df_mod,tmp]).reset_index(drop=True)
            else:
                df_mod = tmp

        # Get moles
        vol_aq = df_mod['Volume(aq)'].astype(float) # vol of aqeous phase in cm3
        n_H_aq = df_mod['n_H(aq)'].astype(float)
        n_H2O_aq = df_mod['n_H2O(aq)'].astype(float)
        n_NH4_aq = df_mod['n_NH4(aq)'].astype(float)
        n_HSO4_aq = df_mod['n_HSO4(aq)'].astype(float)
        n_SO4_aq = df_mod['n_SO4(aq)'].astype(float)
        n_NO3_aq = df_mod['n_NO3(aq)'].astype(float)
        n_NH3_aq = df_mod['n_NH3(aq)'].astype(float)
        n_HNO3_g = df_mod['n_HNO3(g)'].astype(float)
        n_NH3_g = df_mod['n_NH3(g)'].astype(float)
        n_NH3_aq = df_mod['n_NH3(aq)'].astype(float)
        num_s = int(len([x for x in df_mod.columns if x.startswith('id_')])) # max number of solids formed at once
        
        if 'SNA' in run_name.split('_')[1]:
            ids = [10,11,12,13] # id's of solids
            n_s010 = np.array([np.nan]*len(df_mod))
            n_s011 = np.array([np.nan]*len(df_mod))
            n_s012 = np.array([np.nan]*len(df_mod))
            n_s013 = np.array([np.nan]*len(df_mod))
        else:
            n_Na_aq = df_mod['n_Na(aq)'].astype(float)
            n_Cl_aq = df_mod['n_Cl(aq)'].astype(float)
            
            ids = [10,11,12,13,17,18,20,22,25,27]
            n_s010 = np.array([np.nan]*len(df_mod))
            n_s011 = np.array([np.nan]*len(df_mod))
            n_s012 = np.array([np.nan]*len(df_mod))
            n_s013 = np.array([np.nan]*len(df_mod))
            n_s017 = np.array([np.nan]*len(df_mod))
            n_s018 = np.array([np.nan]*len(df_mod))
            n_s020 = np.array([np.nan]*len(df_mod))
            n_s022 = np.array([np.nan]*len(df_mod))
            n_s025 = np.array([np.nan]*len(df_mod))
            n_s027 = np.array([np.nan]*len(df_mod))
        if num_s>0:
            for id_s in ids:
                for num in range(1,num_s+1):
                    exec("I = np.where(df_mod.id_s%02i=='%03i')[0]"%(num,id_s))
                    exec("n_s%03i[I] = df_mod.moles_s%02i[I].astype(float)"%(id_s,num))

        # Add moles
        if 'SNA' in run_name:
            n_SO4TOT = np.nansum([ n_HSO4_aq, n_SO4_aq, n_s010, 2*n_s011, n_s012 ], axis=0)
            n_NH4AER = np.nansum([ n_NH4_aq, n_NH3_aq, 2*n_s010, 3*n_s011, n_s012, n_s013 ], axis=0)
            n_NO3AER = np.nansum([ n_NO3_aq, n_s013 ], axis=0)
        else:
            n_SO4TOT = np.nansum([ n_HSO4_aq, n_SO4_aq, n_s010, 2*n_s011, n_s012, n_s018, 2*n_s020, n_s022 ], axis=0)
            n_NH4AER = np.nansum([ n_NH4_aq, n_NH3_aq, 2*n_s010, 3*n_s011, n_s012, n_s013, n_s017 ], axis=0)
            n_NO3AER = np.nansum([ n_NO3_aq, n_s013, n_s025 ], axis=0)
            n_NATOT = np.nansum([ n_Na_aq, 2*n_s018, 3*n_s020, n_s022, n_s025, n_s027 ], axis=0)
            n_CLTOT = np.nansum([ n_Cl_aq, n_s017, n_s027 ], axis=0)

        n_NH4TOT = np.nansum([ n_NH3_g, n_NH4AER], axis=0)
        n_NO3TOT = np.nansum([ n_HNO3_g, n_NO3AER], axis=0)

        # Get moles organics
        if 'org' in run_name:
            
            n_Oxalic = df_mod['n_Oxalic(aq)'].astype(float)
            n_HOxal  = df_mod['n_HOxal-(aq)'].astype(float)
            n_Oxal2  = df_mod['n_Oxal2-(aq)'].astype(float)

            n_Acetic = df_mod['n_Acetic(aq)'].astype(float)
            n_Acet   = df_mod['n_Acet-(aq)'].astype(float)

            n_Formic = df_mod['n_Formic(aq)'].astype(float)            
            n_Form   = df_mod['n_Form-(aq)'].astype(float)

            M_OXAL = 88.02e6 # oxalate (https://pubchem.ncbi.nlm.nih.gov/compound/71081)
            M_ACET = 59.04e6 # acetate (https://pubchem.ncbi.nlm.nih.gov/compound/175)
            M_FORM = 45.017e6 # formate (https://pubchem.ncbi.nlm.nih.gov/compound/283)
            
            # OXAL_DIS = np.nansum([n_HOxal,n_Oxal2],axis=0)#*M_OXAL
            HOXAL_DIS = n_HOxal#*M_OXAL
            OXAL2_DIS = n_Oxal2#*M_OXAL
            ACET_DIS = n_Acet#*M_ACET
            FORM_DIS = n_Form#*M_FORM

            OXAL_NONDIS = n_Oxalic#*M_OXAL
            ACET_NONDIS = n_Acetic#*M_ACET
            FORM_NONDIS = n_Formic#*M_FORM
        
        if 'Lev' in run_name:
            
            n_Lev = df_mod['n_Levogl(aq)'].astype(float)
            
            M_LEV = 162.14e6 # levoglucosan (https://pubchem.ncbi.nlm.nih.gov/compound/2724705)
            
            LEV = n_Lev#*M_LEV
        
        # save where inputs were effectively NaNs
        I = np.where(np.array([float('%.1g'%x) for x in n_SO4TOT]) == 1e-20)[0]
        
        # Convert mol to ug (already on a per m3 basis)
        M_HNO3 = 63.013e6 # nitric acid (https://pubchem.ncbi.nlm.nih.gov/compound/944)
        M_NH3  = 17.031e6 # ammonia (https://pubchem.ncbi.nlm.nih.gov/compound/222)
        M_H2SO4= 98.08e6  # sulfuric acid (https://pubchem.ncbi.nlm.nih.gov/compound/1118)
        M_NA   = 22.989e6 # sodium ion (https://pubchem.ncbi.nlm.nih.gov/compound/5360545)
        M_CL   = 35.45e6  # chlroide (https://pubchem.ncbi.nlm.nih.gov/compound/312)
        M_H2O  = 18.015e6 # water (https://pubchem.ncbi.nlm.nih.gov/compound/962)
        SO4TOT = n_SO4TOT#*M_H2SO4
        NH4TOT = n_NH4TOT#*M_NH3
        NO3TOT = n_NO3TOT#*M_HNO3
        NH4AER = n_NH4AER#*M_NH3
        NO3AER = n_NO3AER#*M_HNO3
        WATER  = n_H2O_aq#*M_H2O
        if 'SNA' not in run_name:
            NATOT = n_NATOT#*M_NA
            CLTOT = n_CLTOT#*M_CL
            
        # Calculate pH
        # mH = 
        # pH = -np.log10( (n_H_aq) / (vol_aq*1e-3) )

        # Get it into a dataframe
        df_mod = pd.DataFrame({'WATER':WATER,'pH':pH,'SO4TOT':SO4TOT,'NH4TOT':NH4TOT,'NO3TOT':NO3TOT,'NH4AER':NH4AER,'NO3AER':NO3AER})

        if 'SNA' not in run_name:
            df_mod = pd.concat([df_mod, pd.DataFrame({'NATOT':NATOT,'CLTOT':CLTOT})],axis=1)
        if 'org' in run_name:
            df_mod = pd.concat([df_mod, pd.DataFrame({'HOXAL_DIS':HOXAL_DIS, 'OXAL2_DIS':OXAL2_DIS, 'OXAL_NONDIS':OXAL_NONDIS,
                                                      'ACET_DIS':ACET_DIS, 'ACET_NONDIS':ACET_NONDIS,
                                                      'FORM_DIS':FORM_DIS, 'FORM_NONDIS':FORM_NONDIS })],axis=1)
        if 'Lev' in run_name:
            df_mod = pd.concat([df_mod, pd.DataFrame({'LEV':LEV})],axis=1)

        # Write NaNs into dataframe
        df_mod.loc[I] = np.nan
        df_mod.loc[I_rh] = np.nan

    # Average 1Hz to 2-min PILS timesteps
    if res=='1Hz':
        # Get 2-min data
        df_2min = pd.read_csv('./Data/TRANS2AM_Merge_2min.csv')
        # Loop through 2-min PILS timesteps
        for i in range(len(df_2min)):
            # Average model output in each 2-min PILS timestep
            df_mod[i:i+1] = np.nanmean(df_mod[(df_obs['Flight']==df_2min['Flight'][i]) &
                                              (df_obs['Datetime_UTC']>df_2min['Datetime_UTC_start'][i]) &
                                              (df_obs['Datetime_UTC']<=df_2min['Datetime_UTC_stop'][i])],axis=0)
        # Overwrite dataframes
        df_mod = df_mod[:len(df_2min)]
        df_obs = df_2min

    ## Stitch DataFrames together
    df = pd.concat([df_obs,df_mod,df_in],axis=1)
    # Redo Datetime_UTC there are nans
    df['Datetime_UTC'] = df[['Datetime_UTC_start','Datetime_UTC_stop']].mean(axis=1).astype(int)
    
    ## Add datetime
    # df['datetime']=np.nan
    # For each flight (because code relies on UTC since midnight)
    for RF in np.unique(df['Flight']):
        # Add datetime in UTC
        df.loc[(df['Flight']==RF),'datetime'] = pd.to_datetime(df.loc[(df['Flight']==RF),'Datetime_UTC'], unit='s',
                                                               origin=datetime.strptime('20%i'%np.unique(df.loc[(df['Flight']==RF),'yymmdd'])[0], '%Y%m%d'),
                                                               utc=True)
    # Convert object to datetime
    df['datetime'] = pd.to_datetime(df['datetime'])
    # Convert UTC to MT
    df['datetime'] = df.datetime.dt.tz_convert('America/Denver')

    ## add STP_factor conversion factor
    df['STP_factor'] = (273.15*df['ps_hads_a']) / (1013.25*(df['temp']+273.15))
    
    # Overwrite RH with nans
    df['rh'][I_rh]=np.nan
    
    # Add smoke flag
    df['type'] = 'smoke'
    df.loc[(df['Flight'] == 'RF10'),'type'] = 'smoke+dust'
    df.loc[(df['Flight'] == 'RF13'),'type'] = 'clear'
    df.loc[(df['Flight'] == 'RF14'),'type'] = 'clear'
    
    exec('df%i = df'%(counter))
   
df = pd.read_csv( './Data/TRANS2AM_Merge_1Hz.csv' )
for RF in np.unique(df['Flight']):# For each flight (because code relies on UTC since midnight)
    df.loc[(df['Flight']==RF),'datetime'] = pd.to_datetime(df.loc[(df['Flight']==RF),'Datetime_UTC'], unit='s',
                                                           origin=datetime.strptime('20%i'%np.unique(df.loc[(df['Flight']==RF),'yymmdd'])[0], '%Y%m%d'),
                                                           utc=True)    # Add datetime in UTC
df['datetime'] = pd.to_datetime(df['datetime'])# Convert object to datetime
df['datetime'] = df.datetime.dt.tz_convert('America/Denver')# Convert UTC to MT

#%%
def insert_nan(data,res):
    if res=='2min': rundt=120
    if res=='1Hz': rundt=1
    
    # # insert nans where obs is nan
    # data.loc[data.Ammonium_ug_m3.isna(), ['NH4AER','NO3AER']]  = np.nan
    
    # insert nans where PILS cartrige change
    dt = np.diff(data.Datetime_UTC)[:-1]
    if np.any(dt!=rundt):
        nan_df = pd.DataFrame(np.nan, index=[0], columns=data.columns)
        # print(np.where(dt!=rundt))
        # print(np.where(dt!=rundt)[0])
        for i in np.where(dt!=rundt)[0]:
            # print(dt[i])
            data = pd.concat([data[:i+1],nan_df,nan_df,data[i+1:]]).reset_index(drop=True)
            data.datetime.values[i+1] = data.datetime.values[i-1]+np.abs(data.datetime.values[i]-data.datetime.values[i-1])
            data.datetime.values[i+2] = data.datetime.values[i+3]-np.abs(data.datetime.values[i]-data.datetime.values[i-1])
    return(data)

#%%
# for RF in ['RF02','RF03','RF08']:
for RF in ['RF13']:
# for RF in np.roll(np.unique(df1.Flight),1):
    if RF not in ['TF04','RF07','RF10']:
        data = df[df.Flight==RF]
        data = data[data.AVlat.values<=41]
        data = data[data.AVzmsl-data.topo<=1500]
        data = insert_nan(data,'1Hz')

    
        data1 = df1[df1.Flight==RF]
        data1 = data1[data1.AVlat.values<=41]
        data1 = data1[data1.AVzmsl-data1.topo<=1500]
        data1 = insert_nan(data1,'2min')
        if RF=='RF12':
            data1 = data1[:-1]

        fig, ax= plt.subplots(3,1,sharex=True,figsize=(6,5),dpi=300,height_ratios=[1/3, 1, 1])
    
        fmtr = dates.DateFormatter('%H:%M')
        fmtr.set_tzinfo(timezone('America/Denver'))
        
        y=0
        ax[1].fill_between(data1.datetime,y,y+data1.NH4AER*1e9,step='mid',fc='#ebaf3d',alpha=0.5,label='NH$_4^+$')
        ax[1].step(data1.datetime,y+data1.Ammonium_ug_m3*data1.STP_factor/M_NH3*1e9,where='mid',c='#ebaf3d',lw=1,label='NH$_4^+$')
    
        y=0
        ax[2].fill_between(data1.datetime,y,y+data1.NO3AER*1e9,step='mid',fc='#2a58a0',alpha=0.5,label='NO$_3^-$')
        ax[2].step(data1.datetime,y+data1.Nitrate_ug_m3*data1.STP_factor/M_HNO3*1e9,where='mid',c='#2a58a0',lw=1,label='NO$_3^-$')
        y+=data1.NO3AER*1e9
        ax[2].fill_between(data1.datetime,y,y+data1.SO4TOT*1e9*2,step='mid',fc='#e6382f',alpha=0.5,label='SO$_4^{2-}$')
        y+=data1.SO4TOT*2*1e9
        ax[2].fill_between(data1.datetime,y,y+data1.ACET_DIS*1e9,step='mid',fc='#39814a',alpha=0.5,label='Acetate')
        y+=data1.ACET_DIS*1e9
        ax[2].fill_between(data1.datetime,y,y+data1.HOXAL_DIS*1e9,step='mid',fc='#299357',alpha=0.5,label='H-oxalate')
        y+=data1.HOXAL_DIS*1e9
        ax[2].fill_between(data1.datetime,y,y+data1.FORM_DIS*1e9,step='mid',fc='#7abc51',alpha=0.5,label='Formate')
        y+=data1.FORM_DIS*1e9
        ax[2].fill_between(data1.datetime,y,y+data1.OXAL2_DIS*1e9*2,step='mid',fc='#c8de85',alpha=0.5,label='Oxalate')
        y+=data1.OXAL2_DIS*1e9*2
        ymin, ymax = ax[2].get_ylim()
        dytick = np.diff(ax[2].get_yticks())[0]
        ymax = dytick*np.ceil(ymax/dytick)
        ax[2].set_ylim(0,ymax)
        ax[2].xaxis.set_major_formatter(fmtr)
        ax[2].invert_yaxis()
    
        ymax, ymin = ax[2].get_ylim()
        ax[1].set_ylim(ymin,ymax)
        ax[0].step(data1.datetime,data1.Ammonium_ug_m3*data1.STP_factor/M_NH3*1e9,where='mid',c='#ebaf3d',lw=1)
        ax[0].set_ylim(ymax)
        ax[0].set_title(RF+'\n',loc='left')
        ax[0].set_title(data1.datetime.dt.strftime('%b %d %Y').values[0].upper()+'\n',loc='right')
        dytick = np.diff(ax[1].get_yticks())[0]
        ax[0].set_yticks(np.arange(ymax,ax[0].get_ylim()[1],dytick),minor=True)
        ax[0].tick_params(axis='y', which='minor', length=3.5, color='k')
    
        ax[0].autoscale(enable=True, axis='x', tight=True)
        ax[1].autoscale(enable=True, axis='x', tight=True)
        ax[2].autoscale(enable=True, axis='x', tight=True)
    
        xmin,xmax = plt.xlim()
        ax0 = ax[0].twinx()
        ax0.plot(data.datetime,data.NH3_ppbv,c='grey',lw=1,alpha=0.5)
        ax0.set_ylim(100,1000)
        ax0.set_yscale('log')
        ax0.set_xlim(xmin,xmax)
        ax0.spines['right'].set_color('grey')
        ax0.tick_params(axis='y', colors='grey')
        ax0.set_yticklabels([],minor=True)
        ax0.set_yticks([1000])
        ax0.set_yticklabels([1000])
        # ax0.set_yticklabels(['', '', '', '', '', '', '', '', '','1000'])
        # ax0.set_yticks(np.arange(100,ax0.get_ylim()[1],20), minor=True)
        # ax0.tick_params(axis='y', which='minor', length=3.5, color='grey')
    
        ax1 = ax[1].twinx()
        ax1.plot(data.datetime,data.NH3_ppbv,c='grey',lw=1,zorder=0,alpha=0.5)
        ax1.set_ylim(0,100)
        ax1.set_xlim(xmin,xmax)
        ax1.set_ylabel('NH$_3$ (ppbv)', color='grey')#, rotation = 270, labelpad = 15)
        ax1.spines['right'].set_color('grey')
        ax1.tick_params(axis='y', colors='grey')
        
        ax2 = ax[2].twinx()
        ax2.plot(data.datetime,data.HNO3_ppbv,c='grey',lw=1,zorder=0,alpha=0.5)
        ax2.set_ylim(0,10)
        ax2.set_xlim(xmin,xmax)
        ax2.set_ylabel('HNO$_3$ (ppbv)', color='grey')#, rotation = 270, labelpad = 15)
        ax2.spines['right'].set_color('grey')
        ax2.tick_params(axis='y', colors='grey')
        ax2.invert_yaxis()
    
        lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
        lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
        lines_mod = [lines[x] for x in [0,2,4,5,6,7,8]]
        labels_mod = [labels[x] for x in [0,2,4,5,6,7,8]]
        leg = fig.legend(lines_mod, labels_mod, title='\n', loc='lower left', frameon=False,
                          ncol=7,bbox_to_anchor=(0.13,0.81+0.06),labelcolor='linecolor',handlelength=0,columnspacing=0.3)
    
        for t in leg.texts: t.set_alpha(1)
            
        legend_elements = [Line2D([0], [0], color='k', lw=1, label='Obs'),
                            Patch(facecolor='k', alpha=0.3, label='Model')]
        fig.legend(handles=legend_elements, loc='lower center', ncol=2,bbox_to_anchor=(0.5,0.87+0.05),frameon=False)
    
        
        fig.supylabel('Aerosol Composition (neq m$^{-3}$)')
        plt.tight_layout()
        plt.subplots_adjust(wspace=0, hspace=0)
    
        ax[0].spines[['top', 'right','bottom']].set_visible(False)
        ax0.spines[['top','bottom']].set_visible(False)
        ax[1].spines[['top', 'right','bottom']].set_visible(False)
        ax1.spines[['top','bottom']].set_visible(False)
        ax[2].spines[['right','top']].set_visible(False)
        ax2.spines[['top','bottom']].set_visible(False)
        
        ax[0].grid(lw=1,alpha=0.5,axis='x')
        ax[1].grid(lw=1,alpha=0.5,axis='x')
        ax[2].grid(lw=1,alpha=0.5,axis='x')
    
        if ax[0].get_ylim()[1]>=100: ins = plt.axes([0.143,-0.35,0.737,0.349])
        else: ins = plt.axes([0.129,-0.35,0.751,0.349])
        ins.step(data1.datetime,data1.rh,c='#469b9b',lw='2')
        ins.set_xlim(xmin,xmax)
        ins.set_ylim(0,100)
        ins.set_ylabel('RH (%)', color='#469b9b')
        ins.xaxis.set_major_formatter(fmtr)
        ins2 = ins.twinx()
        interval = np.hstack([np.linspace(0, 0.35), np.linspace(1-0.35, 1)])
        colors = plt.cm.coolwarm(interval)
        cmap = LinearSegmentedColormap.from_list('name', colors)
        # im = ins2.scatter(data.datetime,(data.AVzmsl-data.topo)/1000,c=data.temp,cmap=cmap,s=3,vmin=10,vmax=30)
        im = ins2.scatter(data.datetime,(data.AVzmsl-data.topo)/1000,c=data.temp,cmap='coolwarm',s=3,vmin=15,vmax=30)
        ins2.set_xlim(xmin,xmax)
        ins2.set_ylim(0,1.5)
        ins2.set_yticks([0,0.5,1.0,1.5])
        ins2.set_ylabel('Altitude (km)')#, rotation = 270, labelpad = 15)
        ins.spines['left'].set_color('#469b9b')
        ins.tick_params(axis='y', colors='#469b9b')
        ins.grid(lw=1,alpha=0.5,axis='x')
        ins.spines[['top']].set_visible(False)
        ins2.spines[['top','left']].set_visible(False)
        '''
        cax = plt.axes([1.01,-0.35,0.025,0.349])
        cbar = plt.colorbar(im,cax=cax,extend='both')
        cbar.set_label('Temp (˚C)')
        cbar.set_ticks([15,20,25,30])
        '''
        cax = plt.axes([0.635,-0.001-0.03,0.2,0.03])
        cbar = plt.colorbar(im,cax=cax,extend='both',orientation='horizontal')
        cbar.set_ticks([15,22.5,30])
        cbar.ax.set_xticklabels(['15˚C','22.5˚C','30˚C'])
        ins.set_zorder(1)  # default zorder is 0 for ax1 and ax2
        ins.patch.set_visible(False)  # prevents ax1 from hiding ax2
    
        ins = plt.axes([1.225-0.15,0.6,0.3,0.3])
        x,y = df1.Ammonium_ug_m3*df1.STP_factor,df1.NH4AER*M_NH3
        x,y = x[df1.Flight!=RF],y[df1.Flight!=RF]
        x,y = x[df1.AVlat<=41], y[df1.AVlat<=41]
        x,y = x[df1.AVzmsl-df1.topo<=1500], y[df1.AVzmsl-df1.topo<=1500]
        I = np.where( (np.isfinite(x)) & (np.isfinite(y)))
        r = stats.linregress(x.iloc[I],y.iloc[I]).rvalue
        ins.plot(x,y,'o',c='#ebaf3d',markersize=0.5)
        ins.text(0.5, 0.95, 'r$_{\mathrm{other}}$=%.2f'%r,c='#ebaf3d', ha='left', va='top', transform=ins.transAxes)
        x,y = data1.Ammonium_ug_m3*data1.STP_factor,data1.NH4AER*M_NH3
        x,y = x[data1.Flight==RF],y[data1.Flight==RF],
        I = np.where( (np.isfinite(x)) & (np.isfinite(y)))
        r = stats.linregress(x.iloc[I],y.iloc[I]).rvalue
        ins.plot(x,y,'o',c='r',markersize=0.5)
        ins.text(0.5, 0.85, 'r$_{\mathrm{%s}}$=%.2f'%(RF,r),c='r', ha='left', va='top', transform=ins.transAxes)
        ins.axline((0,0),slope=1,lw=0.5,c='k')
        ins.axis('square')
        axmax = np.max([plt.gca().get_xlim()[1],plt.gca().get_ylim()[1]])
        ins.set_xlim(0,axmax)
        ins.set_ylim(0,axmax)
        ins.set_xticks(np.arange(0,3+1,1))
        ins.set_xticks(np.arange(0,3+0.5,0.5),minor=True)
        ins.set_yticks(np.arange(0,3+1,1))
        ins.set_yticks(np.arange(0,3+0.5,0.5),minor=True)
        ins.set_xlabel('Obs NH$_4^+$(µg m$^{-3}$)')
        ins.set_ylabel('Model NH$_4^+$ (µg m$^{-3}$)')
        
        ins = plt.axes([1.225-0.15,0.15,0.3,0.3])
        x,y = df1.Nitrate_ug_m3*df1.STP_factor,df1.NO3AER*M_HNO3
        x,y = x[df1.Flight!=RF],y[df1.Flight!=RF]
        x,y = x[df1.AVlat<=41], y[df1.AVlat<=41]
        x,y = x[df1.AVzmsl-df1.topo<=1500], y[df1.AVzmsl-df1.topo<=1500]
        I = np.where( (np.isfinite(x)) & (np.isfinite(y)))
        r = stats.linregress(x.iloc[I],y.iloc[I]).rvalue
        ins.plot(x,y,'o',c='#2a58a0',markersize=0.5)
        ins.text(0.5, 0.95, 'r$_{\mathrm{other}}$=%.2f'%r,c='#2a58a0', ha='left', va='top', transform=ins.transAxes)
        x,y = data1.Nitrate_ug_m3*data1.STP_factor,data1.NO3AER*M_HNO3
        x,y = x[data1.Flight==RF],y[data1.Flight==RF],
        I = np.where( (np.isfinite(x)) & (np.isfinite(y)))
        r = stats.linregress(x.iloc[I],y.iloc[I]).rvalue
        ins.plot(x,y,'o',c='r',markersize=0.5)
        ins.text(0.5, 0.85, 'r$_{\mathrm{%s}}$=%.2f'%(RF,r),c='r', ha='left', va='top', transform=ins.transAxes)
        ins.axline((0,0),slope=1,lw=0.5,c='k')
        ins.axis('square')
        axmax = np.max([plt.gca().get_xlim()[1],plt.gca().get_ylim()[1]])
        ins.set_xscale('log')
        ins.set_yscale('log')
        # ins.set_xlim(0,axmax)
        # ins.set_ylim(0,axmax)
        # ins.set_xticks(np.arange(0,5+1,2))
        # ins.set_xticks(np.arange(0,5+1,1),minor=True)
        # ins.set_yticks(np.arange(0,5+1,2))
        # ins.set_yticks(np.arange(0,5+1,1),minor=True)
        ins.set_xlabel('Obs NO$_3^-$ (µg m$^{-3}$)')
        ins.set_ylabel('Model NO$_3^-$ (µg m$^{-3}$)')
        
        ins = plt.axes([1.25-0.15,-0.3,0.247,0.3])
        x,y = df1.RH,df1.WATER*M_H2O
        x,y = x[df1.AVlat<=41], y[df1.AVlat<=41]
        x,y = x[df1.AVzmsl-df1.topo<=1500], y[df1.AVzmsl-df1.topo<=1500]
        ins.scatter(x[df1.Flight!=RF],y[df1.Flight!=RF],s=0.5,c='#469b9b')
        ins.scatter(x[df1.Flight==RF],y[df1.Flight==RF],s=0.5,c='red')
        ins.set_xlim(0,100)
        ins.set_ylim(0)
        ins.set_xticks(np.arange(0,101,20))
        ins.set_xticks(np.arange(0,101,10),minor=True)
        ins.set_yticks(np.arange(0,15+5,5))
        ins.set_yticks(np.arange(0,15+2.5,2.5),minor=True)
        ins.set_xlabel('Obs RH (%)')
        ins.set_ylabel('Model AWC (µg m$^{-3}$)')
        
        ins = plt.axes([1.305-0.15,-0.1,0.135,0.095])
        x,y = df1.temp, (df1.AVzmsl-df1.topo)/1000
        x,y = x[df1.AVlat<=41], y[df1.AVlat<=41]
        x,y = x[df1.AVzmsl-df1.topo<=1500], y[df1.AVzmsl-df1.topo<=1500]
        ins.plot(x[df1.Flight!=RF],y[df1.Flight!=RF],'o',c='grey',markersize=0.5)
        ins.plot(x[df1.Flight==RF],y[df1.Flight==RF],'o',c='r',markersize=0.5)
        ins.set_xlabel('T (˚C)',size=7)
        ins.set_ylabel('z (km)',size=7)
        ins.set_yticks([0.5,1.5],minor=True)
        ins.set_ylim(0)
        
        ins.tick_params(axis='y', which='minor', length=2.5, color='k')
        ins.tick_params(axis='both', which='major', labelsize=7)
        
        plt.show()
        
        #%%
        plt.plot(data.datetime,data.NH3_ppbv,c='grey',lw=1,alpha=0.5)
        # plt.ylim(100,1000)
        plt.yscale('log')
        plt.ylim(top=1000)
