#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 14:58:05 2023

@author: eli
"""

import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import matplotlib.patheffects as pe
from matplotlib.font_manager import FontProperties
import matplotlib.image as mpimg
import cmasher as cmr

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy
#%% Import data
# Get TRANS2Am data
df = pd.read_csv( './Data/TRANS2AM_Merge_1Hz.csv' )

# Get map data
mp = pd.read_csv( './Data/CO_map.csv' )

# Get CAFO data
CAFO = pd.read_csv( './Data/CO_CAFOs.csv' )
CAFO.Longitude = -np.abs(CAFO.Longitude)
# Get CAFO data
OG = pd.read_csv( './Data/CO_OG.csv' )

# Get CAFO data
coal = pd.read_csv( './Data/CO_coal.csv' )

## Get obs
df1 = pd.read_csv( './Data/TRANS2AM_Merge_1Hz.csv' )
for RF in np.unique(df1['Flight']):# For each flight (because code relies on UTC since midnight)
    df1.loc[(df1['Flight']==RF),'datetime'] = pd.to_datetime(df1.loc[(df1['Flight']==RF),'Datetime_UTC'], unit='s',
                                                            origin=datetime.strptime('20%i'%np.unique(df1.loc[(df1['Flight']==RF),'yymmdd'])[0], '%Y%m%d'),
                                                            utc=True)    # Add datetime in UTC
df1['datetime'] = pd.to_datetime(df1['datetime'])# Convert object to datetime
df1['datetime'] = df1.datetime.dt.tz_convert('America/Denver')# Convert UTC to MT
## add STP_factor conversion factor
df1['STP_factor'] = (273.15*df1['ps_hads_a']) / (1013.25*(df1['temp']+273.15))

df2 = pd.read_csv( './Data/TRANS2AM_Merge_2min.csv' )
for RF in np.unique(df2['Flight']):# For each flight (because code relies on UTC since midnight)
    df2.loc[(df2['Flight']==RF),'datetime'] = pd.to_datetime(df2.loc[(df2['Flight']==RF),'Datetime_UTC'], unit='s',
                                                            origin=datetime.strptime('20%i'%np.unique(df2.loc[(df2['Flight']==RF),'yymmdd'])[0], '%Y%m%d'),
                                                            utc=True)    # Add datetime in UTC
df2['datetime'] = pd.to_datetime(df2['datetime'])# Convert object to datetime
df2['datetime'] = df2.datetime.dt.tz_convert('America/Denver')# Convert UTC to MT
## add STP_factor conversion factor
df2['STP_factor'] = (273.15*df2['ps_hads_a']) / (1013.25*(df2['temp']+273.15))


#%%


label_loc=True
xmin, xmax = -105.5,-102.4
ymin, ymax = 39.7,41.0

data1 = df1[~df1['Flight'].isin(['TF04','RF07','RF10'])]
data1 = data1[data1.AVzmsl - data1.topo <= 1500]
data1 = data1[data1.AVlat <= 41]

fig = plt.figure(figsize=(9,3.5),dpi=600)

ax = plt.subplot(111)
# urban region
plt.fill(mp.CO_urban_lon,mp.CO_urban_lat, linewidth=0,c='#9b6a0a',label='Urban areas', alpha=0.5)

# label towns
if label_loc:
    plt.text(-104.9903, 39.7392,'Denver',
              weight='heavy', ha='center',
              path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])
    plt.text(-105.2705, 40.0150,'Boulder',
              weight='heavy', ha='center', size=7,
              path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])
    plt.text(-105.1019, 40.1672,'Longmont',
              weight='heavy', ha='center', size=7,
              path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])
    plt.text(-104.7091, 40.4233,'Greeley',
              ha='center', va='top', weight='heavy', size=7,
              path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])
    plt.text(-105.0750, 40.3978,'Loveland',
              ha='center', weight='heavy', size=7,
              path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])
    plt.text(-105.0844, 40.5853,'Fort Collins',
              ha='center', weight='heavy', size=7,
              path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])
    plt.text(-103.8000, 40.2503,'Fort Morgan',
              weight='heavy', va='bottom', ha='center', size=7,
              path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])
    plt.text(-103.2077, 40.6255,'Sterling',
              weight='heavy', va='bottom', ha='center', size=7,
              path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])

# flight track colored by ammonia
im = plt.scatter(data1.AVlon,data1.AVlat,c=data1.NH3_ppbv,s=0.1,vmin=0,vmax=20,cmap='Spectral_r')
cbar = plt.colorbar(im,pad=0,extend='max',ticks=[0,5,10,15,20])
cbar.set_label('NH$_3$ (ppbv)')
cbar.ax.yaxis.set_ticks([2.5,7.5,12.5,17.5], minor=True)

# highways
# plt.plot(mp.CO_urban_lon,mp.CO_urban_lat,
#           linewidth=1,c='#9b6a0a',alpha=0.5)
plt.plot(mp.CO_highways_lon,mp.CO_highways_lat,
          linewidth=1,c='#a0a0a0',label='Major roads')

# Power plant
plt.scatter(coal.pawnee_PP_lon,coal.pawnee_PP_lat,
            c='#f6be47', edgecolors='k', marker='^', label='Coal plant')


# AFOs
plt.scatter(CAFO.Longitude,CAFO.Cattle_LAT,edgecolors='k', linewidth=0.4,
            c='k',s=10,label='CAFO')
plt.scatter(CAFO.Longitude,CAFO.Dairy_all_LAT,edgecolors='k', linewidth=0.4,
            c='k',s=10)
plt.scatter(CAFO.Longitude,CAFO.Sheep_LAT,edgecolors='k', linewidth=0.4,
            c='k',s=10)
plt.scatter(CAFO.Longitude,CAFO.Swine_LAT,edgecolors='k', linewidth=0.4,
            c='k',s=10)
plt.scatter(CAFO.Longitude,CAFO.Poultry_LAT,edgecolors='k', linewidth=0.4,
            c='k',s=10)


# set lat lon limits
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)

# x- and y-axis ticks and labels
xticks = plt.xticks()[0][:-1]
plt.xticks(xticks,np.abs(xticks))
plt.xlabel('Longitude (˚W)')
plt.ylabel('Latitude (˚N)')

# set aspect ratio (to preserve lat:lon)
f = 1.0/np.cos(40*np.pi/180)
plt.gca().set_aspect(f)

plt.legend(loc='lower right',ncol=2,bbox_to_anchor=(1.0,0.0),columnspacing=0.5,
           handletextpad=0.5, markerscale=0.8, framealpha=1)

# add NH3 spread
ax1 = plt.axes([0.91,0.12,0.13,0.5])
flierprops = dict(marker='o', markerfacecolor='k', markersize=2,
                  linestyle='none', markeredgecolor='None')
ax1.boxplot([df1.NH3_ppbv[np.isfinite(df1.NH3_ppbv)],
              df2.NH3_ppbv[np.isfinite(df2.NH3_ppbv)],
              df2.Ammonium_ug_m3[np.isfinite(df2.Ammonium_ug_m3)]],
            flierprops=flierprops,widths=0.7)
ax1.set_xticks([1,2,3],['NH$_3$\n(1 Hz)','NH$_3$\n          '+r'($\frac{1}{120}$ Hz)','NH$_4^+$'])
ax1.set_yscale('log')
ax1.set_ylim(1e-1,1e3)

# Add topographic US map for context of CO in US
m = plt.axes([0.25,0.75,0.17,0.5],projection=ccrs.PlateCarree(),frameon=False)
img = mpimg.imread('./Data/SR_LR/SR_LR.tif') # https://www.naturalearthdata.com/downloads/10m-shaded-relief/10m-shaded-relief-basic/
img_extent = (-180, 180, -90, 90)
m.imshow(img, origin='upper', extent=img_extent, transform=ccrs.PlateCarree(),cmap='Greys_r',vmin=120,vmax=200)
m.add_feature(cfeature.COASTLINE,lw=0.5)
m.add_feature(cfeature.STATES,lw=0.1)
m.add_feature(cfeature.OCEAN,lw=0.5,zorder=999,fc='white')
m.plot([xmin,xmin,xmax,xmax,xmin],[ymin,ymax,ymax,ymin,ymin],c='r',lw=0.1)
m.fill([xmin,xmin,xmax,xmax,xmin],[ymin,ymax,ymax,ymin,ymin],c='r',lw=0.1,alpha=0.5)
m.set_xlim(-125, -69)
m.set_ylim(24, 49.1)
mark_inset(m, ax, loc1=2, loc2=1, fc="none", ec="0.5",zorder=999)
m.set_zorder(0)


plt.show()
#%%


label_loc=True
xmin, xmax = -105.5,-102.4
ymin, ymax = 39.7,41.0

data1 = df1[~df1['Flight'].isin(['TF04','RF07','RF10'])]
data1 = data1[data1.AVzmsl - data1.topo <= 1500]
data1 = data1[data1.AVlat <= 41]

fig = plt.figure(figsize=(9,3.5),dpi=600)

ax = plt.subplot(111)
# urban region
plt.fill(mp.CO_urban_lon,mp.CO_urban_lat, linewidth=0,c='#9b6a0a',label='Urban areas', alpha=0.5)

# label towns
if label_loc:
    plt.text(-104.9903, 39.7392,'Denver',
              weight='heavy', ha='center',
              path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])
    plt.text(-105.2705, 40.0150,'Boulder',
              weight='heavy', ha='center', size=7,
              path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])
    plt.text(-105.1019, 40.1672,'Longmont',
              weight='heavy', ha='center', size=7,
              path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])
    plt.text(-104.7091, 40.4233,'Greeley',
              ha='center', va='top', weight='heavy', size=7,
              path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])
    plt.text(-105.0750, 40.3978,'Loveland',
              ha='center', weight='heavy', size=7,
              path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])
    plt.text(-105.0844, 40.5853,'Fort Collins',
              ha='center', weight='heavy', size=7,
              path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])
    plt.text(-103.8000, 40.2503,'Fort Morgan',
              weight='heavy', va='bottom', ha='center', size=7,
              path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])
    plt.text(-103.2077, 40.6255,'Sterling',
              weight='heavy', va='bottom', ha='center', size=7,
              path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])

# flight track colored by ammonia
im = plt.scatter(data1.AVlon,data1.AVlat,c=data1.NH3_ppbv,s=0.1,vmin=0,vmax=20,cmap='Spectral_r')
cbar = plt.colorbar(im,pad=0,extend='max',ticks=[0,5,10,15,20])
cbar.set_label('NH$_3$ (ppbv)')
cbar.ax.yaxis.set_ticks([2.5,7.5,12.5,17.5], minor=True)

# highways
# plt.plot(mp.CO_urban_lon,mp.CO_urban_lat,
#           linewidth=1,c='#9b6a0a',alpha=0.5)
plt.plot(mp.CO_highways_lon,mp.CO_highways_lat,
          linewidth=1,c='#a0a0a0',label='Major roads')

# Power plant
plt.scatter(coal.pawnee_PP_lon,coal.pawnee_PP_lat,
            c='#f6be47', edgecolors='k', marker='^', label='Coal plant')


# AFOs
plt.scatter(CAFO.Longitude,CAFO.Cattle_LAT,edgecolors='k', linewidth=0.4,
            c='k',s=10,label='CAFO')
plt.scatter(CAFO.Longitude,CAFO.Dairy_all_LAT,edgecolors='k', linewidth=0.4,
            c='k',s=10)
plt.scatter(CAFO.Longitude,CAFO.Sheep_LAT,edgecolors='k', linewidth=0.4,
            c='k',s=10)
plt.scatter(CAFO.Longitude,CAFO.Swine_LAT,edgecolors='k', linewidth=0.4,
            c='k',s=10)
plt.scatter(CAFO.Longitude,CAFO.Poultry_LAT,edgecolors='k', linewidth=0.4,
            c='k',s=10)


# set lat lon limits
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)

# x- and y-axis ticks and labels
xticks = plt.xticks()[0][:-1]
plt.xticks(xticks,np.abs(xticks))
plt.xlabel('Longitude (˚W)')
plt.ylabel('Latitude (˚N)')

# set aspect ratio (to preserve lat:lon)
f = 1.0/np.cos(40*np.pi/180)
plt.gca().set_aspect(f)

plt.legend(loc='lower right',ncol=2,bbox_to_anchor=(1.0,0.0),columnspacing=0.5,
           handletextpad=0.5, markerscale=0.8, framealpha=1)

# Add topographic US map for context of CO in US
m = plt.axes([0.88,0.5,0.17,0.5],projection=ccrs.PlateCarree(),frameon=False)
img = mpimg.imread('./Data/SR_LR/SR_LR.tif') # https://www.naturalearthdata.com/downloads/10m-shaded-relief/10m-shaded-relief-basic/
img_extent = (-180, 180, -90, 90)
m.imshow(img, origin='upper', extent=img_extent, transform=ccrs.PlateCarree(),cmap='Greys_r',vmin=120,vmax=200)
m.add_feature(cfeature.COASTLINE,lw=0.5)
m.add_feature(cfeature.STATES,lw=0.1)
m.add_feature(cfeature.OCEAN,lw=0.5,zorder=999,fc='white')
m.plot([xmin,xmin,xmax,xmax,xmin],[ymin,ymax,ymax,ymin,ymin],c='r',lw=0.1)
m.fill([xmin,xmin,xmax,xmax,xmin],[ymin,ymax,ymax,ymin,ymin],c='r',lw=0.1,alpha=0.5)
m.set_xlim(-125, -69)
m.set_ylim(24, 49.1)

# add NH3 spread
ax = plt.axes([0.91,0.12,0.13,0.5])
flierprops = dict(marker='o', markerfacecolor='k', markersize=2,
                  linestyle='none', markeredgecolor='None')
ax.boxplot([df1.NH3_ppbv[np.isfinite(df1.NH3_ppbv)],
              df2.NH3_ppbv[np.isfinite(df2.NH3_ppbv)],
              df2.Ammonium_ug_m3[np.isfinite(df2.Ammonium_ug_m3)]],
            flierprops=flierprops,widths=0.7)
ax.set_xticks([1,2,3],['NH$_3$\n(1 Hz)','NH$_3$\n          '+r'($\frac{1}{120}$ Hz)','NH$_4^+$'])
plt.yscale('log')
plt.ylim(1e-1,1e3)

plt.show()
#%%

label_loc=True
xmin, xmax = -105.5,-102.4
ymin, ymax = 39.7,41.0

data1 = df1[~df1['Flight'].isin(['TF04','RF07','RF10'])]
data1 = data1[data1.AVzmsl - data1.topo <= 1500]
data1 = data1[data1.AVlat <= 41]

fig = plt.figure(figsize=(9,3.5),dpi=600)

ax = plt.subplot(111)
# urban region
plt.fill(mp.CO_urban_lon,mp.CO_urban_lat, linewidth=0,c='#9b6a0a',label='Urban areas', alpha=0.5)

# label towns
if label_loc:
    plt.text(-104.9903, 39.7392,'Denver',
              weight='heavy', ha='center',
              path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])
    plt.text(-105.2705, 40.0150,'Boulder',
              weight='heavy', ha='center', size=7,
              path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])
    plt.text(-105.1019, 40.1672,'Longmont',
              weight='heavy', ha='center', size=7,
              path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])
    plt.text(-104.7091, 40.4233,'Greeley',
              ha='center', va='top', weight='heavy', size=7,
              path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])
    plt.text(-105.0750, 40.3978,'Loveland',
              ha='center', weight='heavy', size=7,
              path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])
    plt.text(-105.0844, 40.5853,'Fort Collins',
              ha='center', weight='heavy', size=7,
              path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])
    plt.text(-103.8000, 40.2503,'Fort Morgan',
              weight='heavy', va='bottom', ha='center', size=7,
              path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])
    plt.text(-103.2077, 40.6255,'Sterling',
              weight='heavy', va='bottom', ha='center', size=7,
              path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])

# flight track colored by ammonia
im = plt.scatter(data1.AVlon,data1.AVlat,c=data1.NH3_ppbv,s=0.1,vmin=0,vmax=20,cmap='Spectral_r')
cbar = plt.colorbar(im,pad=0,extend='max',ticks=[0,5,10,15,20])
cbar.set_label('NH$_3$ (ppbv)')
cbar.ax.yaxis.set_ticks([2.5,7.5,12.5,17.5], minor=True)

# highways
plt.plot(mp.CO_highways_lon,mp.CO_highways_lat,
          linewidth=1,c='#a0a0a0',label='Major roads')

# facilities
plt.scatter(coal.pawnee_PP_lon,coal.pawnee_PP_lat,
            c='#f6be47', edgecolors='k', marker='^', label='Coal plant')

plt.scatter(CAFO.Longitude,CAFO.Cattle_LAT,edgecolors='k', linewidth=0.4,
            c='k',s=10,label='CAFO')
plt.scatter(CAFO.Longitude,CAFO.Dairy_all_LAT,edgecolors='k', linewidth=0.4,
            c='k',s=10)
plt.scatter(CAFO.Longitude,CAFO.Sheep_LAT,edgecolors='k', linewidth=0.4,
            c='k',s=10)
plt.scatter(CAFO.Longitude,CAFO.Swine_LAT,edgecolors='k', linewidth=0.4,
            c='k',s=10)
plt.scatter(CAFO.Longitude,CAFO.Poultry_LAT,edgecolors='k', linewidth=0.4,
            c='k',s=10)

# set lat lon limits
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)

# x- and y-axis ticks and labels
xticks = plt.xticks()[0][:-1]
plt.xticks(xticks,np.abs(xticks))
plt.xlabel('Longitude (˚W)')
plt.ylabel('Latitude (˚N)')

# set aspect ratio (to preserve lat:lon)
f = 1.0/np.cos(40*np.pi/180)
plt.gca().set_aspect(f)

plt.legend(loc='lower left',bbox_to_anchor=(1.15,0.3))

# Add topographic US map for context of CO in US
m = plt.axes([0.85,0.5,0.17,0.5],projection=ccrs.PlateCarree(),frameon=False)
img = mpimg.imread('./Data/SR_LR/SR_LR.tif') # https://www.naturalearthdata.com/downloads/10m-shaded-relief/10m-shaded-relief-basic/
img_extent = (-180, 180, -90, 90)
m.imshow(img, origin='upper', extent=img_extent, transform=ccrs.PlateCarree(),cmap='Greys_r',vmin=120,vmax=200)
m.add_feature(cfeature.COASTLINE,lw=0.5)
m.add_feature(cfeature.STATES,lw=0.1)
m.add_feature(cfeature.OCEAN,lw=0.5,zorder=999,fc='white')
m.plot([xmin,xmin,xmax,xmax,xmin],[ymin,ymax,ymax,ymin,ymin],c='r',lw=0.1)
m.fill([xmin,xmin,xmax,xmax,xmin],[ymin,ymax,ymax,ymin,ymin],c='r',lw=0.1,alpha=0.3)
m.set_xlim(-125, -69)
m.set_ylim(24, 49.1)

# add NH3 spread
ax = plt.axes([0.87,0.12,0.15,0.15])#,frameon=False)
flierprops = dict(marker='o', markerfacecolor='k', markersize=2,
                  linestyle='none', markeredgecolor='None')
ax.boxplot([df1.NH3_ppbv[np.isfinite(df1.NH3_ppbv)],
              df2.NH3_ppbv[np.isfinite(df2.NH3_ppbv)],
              df2.Ammonium_ug_m3[np.isfinite(df2.Ammonium_ug_m3)]],
            flierprops=flierprops)
# ax.set_xticks([1,2,3],['(g)\n(1 Hz)','(g)\n          (120$^{-1}$ Hz)','(p)'])
ax.set_xticks([1,2,3],['NH$_3$\n(1 Hz)','NH$_3$\n          (120$^{-1}$ Hz)','NH$_4^+$'])
plt.title('NH$_x$ (ppbv)',size=10)
plt.yscale('log')
ax.set_yticks(ticks = [1e-1,1e0,1e1,1e2,1e3],labels=['0.1','1','10','100','1,000'])
plt.ylim(1e-1,1e3)

plt.show()