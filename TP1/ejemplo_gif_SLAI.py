#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 18:31:27 2021

@author: laurare
"""
# EJEMPLO DE COMO HACER UN GIF CON DATOS DE ALTIMETRIA
import numpy as np
import xarray as xr
from geopy.distance import geodesic
import matplotlib.pyplot as plt
import imageio
import os
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

pru = xr.open_dataset('/Users/laurare/OneDrive - cima.fcen.uba.ar/frente_talud/dataset-duacs-rep-global-merged-allsat-phy-l4-monthly_1604530037460.nc')                 

# Vemos la metadata
pru

# cargo las variables mensuales
lat = pru.latitude.values # data
lon = (pru.longitude.values) - 360 #data
time = pru.time.values # ESTO ES IMPORTANTE QUE SE CARGUE ASI
sla = pru.sla.values #data

# En el caso de tener datos diarios y quiero datos mensuales
#time2 = pru.time.resample(time='M').mean()
#time1 = time2.dt.strftime('%m/%y').values
#time2 = time2.data
#msla = pru.sla.resample(time = 'M').mean()
#msla = msla.data


#### HAGO GIF MENSUAL DE SEA LEVEL ANOMALY
lo,la = np.meshgrid(lon, lat)

# GRAFICO UN MES PARA ULTIMAR DETALLES DE ESTETICA QUE LUEGO REPITO EN EL GIFF
fig = plt.figure(figsize=(10,8))
######## FIG DEF
# ax1 = plt.axes(projection=ccrs.PlateCarree())  #ejes

# crs_latlon = ccrs.PlateCarree()  # proyeccion
# #defino extension ejes
# lonmin = np.floor(lon[0]) 
# lonmax = -45
# latmin = np.floor(lat[0])
# latmax = np.ceil(lat[-1])

# ax1.set_extent([lonmin,lonmax, latmin, latmax], crs = crs_latlon)

# # defino la cantidad de isolineas y el minimo-maximo de los datos
# levels = MaxNLocator(nbins=40).tick_values(-1, 1)
# # pick the desired colormap, sensible levels, and define a normalization
# # instance which takes data values and translates those into levels.
# cmap = plt.get_cmap('jet')
# norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

# plt.pcolormesh(lo, la, sla[0,:,:],cmap = cmap,norm=norm,transform=ccrs.PlateCarree())
# ax1.add_feature(cfeature.LAND,facecolor='silver')
# ax1.add_feature(cfeature.BORDERS, linestyle=':')
# ax1.add_feature(cfeature.COASTLINE)

# gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.8, linestyle='dotted')
# gl.xlabels_top = False
# gl.ylabels_right = False
# gl.ylocator = mticker.FixedLocator(np.arange(latmin,latmax+0.1,1))
# gl.xlocator = mticker.FixedLocator(np.arange(lonmin,lonmax+0.1,3))
# gl.xformatter = LONGITUDE_FORMATTER
# gl.yformatter = LATITUDE_FORMATTER
# gl.xlabel_style = {'size': 14}
# gl.ylabel_style = {'size': 14}

# lon_formatter = LongitudeFormatter(zero_direction_label=True)
# lat_formatter = LatitudeFormatter()

# ax1.xaxis.set_major_formatter(lon_formatter)
# ax1.yaxis.set_major_formatter(lat_formatter)

# ax1.set_aspect('auto', adjustable=None)
# cb = plt.colorbar( orientation='vertical',ticks=[np.arange(-1,1,0.2)])
# plt.title('SLA mensual (m)', fontsize=14)
# fecha = plt.text(-60,-37,str(time[0])[0:10],fontsize=10,transform=ccrs.PlateCarree())



########### HAGO EL GIFF DE DATOS MENSUALES DE SLA ######################
# solo hago un gif de 3 anios
lonmin = np.floor(lon[0]) 
lonmax = -45
latmin = np.floor(lat[0])
latmax = np.ceil(lat[-1])

fig = plt.figure(figsize=(10,8))
filenames = []

with imageio.get_writer('/Users/laurare/Desktop/sla_mensual_Altimetria.gif', mode='I') as writer:
    for j in range(36):
        ax1 = plt.axes(projection=ccrs.PlateCarree())  #ejes

        crs_latlon = ccrs.PlateCarree()  # proyeccion
        #defino extension ejes
        ax1.set_extent([lonmin,lonmax, latmin, latmax], crs = crs_latlon)
        levels = MaxNLocator(nbins=40).tick_values(-1, 1)
        # pick the desired colormap, sensible levels, and define a normalization
        # instance which takes data values and translates those into levels.
        cmap = plt.get_cmap('jet')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        plt.pcolormesh(lo, la, sla[j,:,:],cmap = cmap,norm=norm,transform=ccrs.PlateCarree())
        ax1.add_feature(cfeature.LAND,facecolor='silver')
        ax1.add_feature(cfeature.BORDERS, linestyle=':')
        ax1.add_feature(cfeature.COASTLINE)

        gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.8, linestyle='dotted')
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.ylocator = mticker.FixedLocator(np.arange(latmin,latmax+0.1,1))
        gl.xlocator = mticker.FixedLocator(np.arange(lonmin,lonmax+0.1,3))
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 14}
        gl.ylabel_style = {'size': 14}
    
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        
        ax1.xaxis.set_major_formatter(lon_formatter)
        ax1.yaxis.set_major_formatter(lat_formatter)
        
        ax1.set_aspect('auto', adjustable=None)
        cb = plt.colorbar( orientation='vertical',ticks=[np.arange(-1,1,0.2)])
        plt.title('SLA mensual (m)', fontsize=14)
        fecha = plt.text(-60,-37,str(time[j])[0:10],fontsize=10,transform=ccrs.PlateCarree())
        
        # create file name and append it to a list
        filename = f'/Users/laurare/Desktop/{j}.png'
        #filenames.append(filename)
        plt.savefig(filename)
        image = imageio.imread(filename)
        writer.append_data(image)
        os.remove(filename)
        fecha.remove()
        fig.clf()    
        
         
