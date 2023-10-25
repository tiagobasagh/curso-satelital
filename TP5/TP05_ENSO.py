#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 14:19:32 2023

@author: laura
"""
#import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import path
#import julian as jd
import xarray as xr


#%pip uninstall shapely
#%pip install shapely --no-binary shapely # HAY QUE RESPONDER Y (YES)
#%pip install --upgrade cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

#### EJERCICIO 2 ####
# Cargo SST para fase negativa
#sst1 = xr.open_dataset('/content/g4.timeAvgMap.MODISA_L3m_NSST_Monthly_9km_R2019_0_sst.20150801-20151031.134W_45S_33W_17N.nc',decode_times=False)
sst1 = xr.open_dataset('/Users/laura/OneDrive - cima.fcen.uba.ar/Oceanografia_satelital/2023/TP5/g4.timeAvgMap.MODISA_L3m_NSST_Monthly_9km_R2019_0_sst.20150801-20151031.134W_45S_33W_17N.nc')

# miro metadata
sst1

# chequeo unidades
sst1.MODISA_L3m_NSST_Monthly_9km_R2019_0_sst.units

# cargo variables
lon = sst1.lon.data
lat = sst1.lat.data
sst_pos = sst1.MODISA_L3m_NSST_Monthly_9km_R2019_0_sst.data

# grafico
fig = plt.figure(figsize=(10, 8), facecolor='white')
ax = plt.axes(projection=ccrs.Mercator())
lon_min = -180; lon_max = -45
lat_min = -45; lat_max = 20
ax.set_extent([lon_min, lon_max, lat_min, lat_max],crs = ccrs.PlateCarree()) # min lon, max lon, min lat, max lat
ax.add_feature(cfeature.LAND, color='silver')
ax.add_feature(cfeature.LAKES, color='lightcyan')
ax.add_feature(cfeature.RIVERS, edgecolor='black')
ax.coastlines(resolution='50m', color='black', linewidth=1)
plt.title('SST, fase positiva',color='black', size=14)
cm = plt.pcolormesh(lon,lat,sst_pos,vmin = np.nanmin(sst_pos),vmax = np.nanmax(sst_pos),cmap = 'jet',transform = ccrs.PlateCarree())
divider = make_axes_locatable(ax)
ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)

fig.add_axes(ax_cb)
plt.colorbar(cm,cax=ax_cb)

# elegir el rango de valores de sst que destaque mejor lo que queremos ver
cm.set_clim(10,28)

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.8, linestyle='dotted')
#gl.top_labels = False; gl.right_labels = False
gl.ylabels_right = False; gl.xlabels_top = False
gl.xlocator = mticker.FixedLocator([-180,-160,-140,-120,-100,-80,-60,-40])
gl.ylocator = mticker.FixedLocator([-50,-40,-30,-20,-10,0,10,20])
gl.xformatter = LONGITUDE_FORMATTER; gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 13}
gl.ylabel_style = {'size': 13}
plt.show()


# Cargo SST para fase negativa
#sst2 = xr.open_dataset('/content/g4.timeAvgMap.MODISA_L3m_NSST_Monthly_9km_R2019_0_sst.20220801-20221031.134W_45S_33W_17N.nc',decode_times=False)
sst2 = xr.open_dataset('/Users/laura/OneDrive - cima.fcen.uba.ar/Oceanografia_satelital/2023/TP5/g4.timeAvgMap.MODISA_L3m_NSST_Monthly_9km_R2019_0_sst.20220801-20221031.134W_45S_33W_17N.nc')

# miro metadata
sst2

# chequeo unidades
sst2.MODISA_L3m_NSST_Monthly_9km_R2019_0_sst.units

# cargo variables
#lon = sst2.lon.data
#lat = sst2.lat.data
sst_neg = sst2.MODISA_L3m_NSST_Monthly_9km_R2019_0_sst.data

# grafico
fig = plt.figure(figsize=(10, 8), facecolor='white')
ax = plt.axes(projection=ccrs.Mercator())
lon_min = -180; lon_max = -45
lat_min = -45; lat_max = 20
ax.set_extent([lon_min, lon_max, lat_min, lat_max],crs = ccrs.PlateCarree()) # min lon, max lon, min lat, max lat
ax.add_feature(cfeature.LAND, color='silver')
ax.add_feature(cfeature.LAKES, color='lightcyan')
ax.add_feature(cfeature.RIVERS, edgecolor='black')
ax.coastlines(resolution='50m', color='black', linewidth=1)
plt.title('SST, fase negativa',color='black', size=14)

# Plot the wind vectors
cm = plt.pcolormesh(lon,lat,sst_neg,vmin = np.nanmin(sst_neg),vmax = np.nanmax(sst_neg),cmap = 'jet',transform = ccrs.PlateCarree())
divider = make_axes_locatable(ax)
ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)

fig.add_axes(ax_cb)
plt.colorbar(cm,cax=ax_cb)

# elegir el rango de valores de sst que destaque mejor lo que queremos ver
cm.set_clim(10,28)

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.8, linestyle='dotted')
#gl.top_labels = False; gl.right_labels = False
gl.ylabels_right = False; gl.xlabels_top = False
gl.xlocator = mticker.FixedLocator([-180,-160,-140,-120,-100,-80,-60,-40])
gl.ylocator = mticker.FixedLocator([-50,-40,-30,-20,-10,0,10,20])
gl.xformatter = LONGITUDE_FORMATTER; gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 13}
gl.ylabel_style = {'size': 13}
plt.show()


# Hago la diferencia entre fase positiva y negativa
# grafico la diferencia total entre el Ninio y la Ninia. 
# Tambien se puede hacer la diferencia entre el Ninio con respecto a la fase neutra 
# y la diferencia entre la Ninia con respcto a la fase neutra

fig = plt.figure(figsize=(10, 8), facecolor='white')
ax = plt.axes(projection=ccrs.Mercator())
lon_min = -180; lon_max = -45
lat_min = -45; lat_max = 20
ax.set_extent([lon_min, lon_max, lat_min, lat_max],crs = ccrs.PlateCarree()) # min lon, max lon, min lat, max lat
ax.add_feature(cfeature.LAND, color='silver')
ax.add_feature(cfeature.LAKES, color='lightcyan')
ax.add_feature(cfeature.RIVERS, edgecolor='black')
ax.coastlines(resolution='50m', color='black', linewidth=1)
plt.title('Anomalia de SST, fase positiva - fase negativa',color='black', size=14)
cm = plt.pcolormesh(lon,lat,sst_pos-sst_neg,vmin = np.nanmin(sst_pos-sst_neg),vmax = np.nanmax(sst_pos-sst_neg),cmap = 'seismic',transform = ccrs.PlateCarree())
divider = make_axes_locatable(ax)
ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)

fig.add_axes(ax_cb)
plt.colorbar(cm,cax=ax_cb)

# elegir el rango de valores de sst que destaque mejor lo que queremos ver
cm.set_clim(-5,5)

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.8, linestyle='dotted')
#gl.top_labels = False; gl.right_labels = False
gl.ylabels_right = False; gl.xlabels_top = False
gl.xlocator = mticker.FixedLocator([-180,-160,-140,-120,-100,-80,-60,-40])
gl.ylocator = mticker.FixedLocator([-50,-40,-30,-20,-10,0,10,20])
gl.xformatter = LONGITUDE_FORMATTER; gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 13}
gl.ylabel_style = {'size': 13}
plt.show()



#### EJERCICIO 4 ####
# cargo SLA fase positiva
alt = xr.open_dataset('/Users/Laura/OneDrive - cima.fcen.uba.ar/Oceanografia_satelital/2023/TP5/SLA_fase_positiva_ENSO.nc')

# Vemos la metadata del xarray
alt

lon = alt.longitude.data-360 
lat = alt.latitude.data
sla_pos = alt.sla.data
time = alt.time.data

# Grafico 
fig = plt.figure(figsize=(10, 8), facecolor='white')
ax = plt.axes(projection=ccrs.Mercator())
#ax = plt.axes(projection=ccrs.Orthographic(central_longitude=-120.0, central_latitude=-30.0))
lon_min = -180; lon_max = -45
lat_min = -45; lat_max = 20
ax.set_extent([lon_min, lon_max, lat_min, lat_max],crs = ccrs.PlateCarree()) # min lon, max lon, min lat, max lat
ax.add_feature(cfeature.LAND, color='silver')
ax.add_feature(cfeature.LAKES, color='lightcyan')
ax.add_feature(cfeature.RIVERS, edgecolor='black')
ax.coastlines(resolution='50m', color='black', linewidth=1)
plt.title('SLA, fase positiva',color='black', size=14)
cm = plt.pcolormesh(lon,lat,np.nanmean(sla_pos,axis=0),cmap = 'seismic',transform = ccrs.PlateCarree())
divider = make_axes_locatable(ax)
ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)

fig.add_axes(ax_cb)
plt.colorbar(cm,cax=ax_cb)
cm.set_clim(-.5,.5)

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.8, linestyle='dotted')
#gl.top_labels = False; gl.right_labels = False
gl.ylabels_right = False; gl.xlabels_top = False
gl.xlocator = mticker.FixedLocator([-180,-160,-140,-120,-100,-80,-60,-40])
gl.ylocator = mticker.FixedLocator([-40,-30,-20,-10,0,10,20])
gl.xformatter = LONGITUDE_FORMATTER; gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 13}
gl.ylabel_style = {'size': 13}
plt.show()


# cargo SLA fase negativa
alt = xr.open_dataset('/Users/Laura/OneDrive - cima.fcen.uba.ar/Oceanografia_satelital/2023/TP5/SLA_fase_negativa_ENSO.nc')

# Vemos la metadata del xarray
alt

lon = alt.longitude.data-360 
lat = alt.latitude.data
sla_neg = alt.sla.data
time = alt.time.data

# Grafico 
fig = plt.figure(figsize=(10, 8), facecolor='white')
ax = plt.axes(projection=ccrs.Mercator())
#ax = plt.axes(projection=ccrs.Orthographic(central_longitude=-120.0, central_latitude=-30.0))
lon_min = -180; lon_max = -45
lat_min = -45; lat_max = 20
ax.set_extent([lon_min, lon_max, lat_min, lat_max],crs = ccrs.PlateCarree()) # min lon, max lon, min lat, max lat
ax.add_feature(cfeature.LAND, color='silver')
ax.add_feature(cfeature.LAKES, color='lightcyan')
ax.add_feature(cfeature.RIVERS, edgecolor='black')
ax.coastlines(resolution='50m', color='black', linewidth=1)
plt.title('SLA, fase negativa',color='black', size=14)
cm = plt.pcolormesh(lon,lat,np.nanmean(sla_neg,axis=0),cmap = 'seismic',transform = ccrs.PlateCarree())
divider = make_axes_locatable(ax)
ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)

fig.add_axes(ax_cb)
plt.colorbar(cm,cax=ax_cb)
cm.set_clim(-.5,.5)

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.8, linestyle='dotted')
#gl.top_labels = False; gl.right_labels = False
gl.ylabels_right = False; gl.xlabels_top = False
gl.xlocator = mticker.FixedLocator([-180,-160,-140,-120,-100,-80,-60,-40])
gl.ylocator = mticker.FixedLocator([-40,-30,-20,-10,0,10,20])
gl.xformatter = LONGITUDE_FORMATTER; gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 13}
gl.ylabel_style = {'size': 13}
plt.show()


# cargo SLA fase neutra
alt = xr.open_dataset('/Users/Laura/OneDrive - cima.fcen.uba.ar/Oceanografia_satelital/2023/TP5/SLA_fase_neutra.nc')

# Vemos la metadata del xarray
alt

lon = alt.longitude.data-360 
lat = alt.latitude.data
sla_n = alt.sla.data
time = alt.time.data

# Grafico 
fig = plt.figure(figsize=(10, 8), facecolor='white')
ax = plt.axes(projection=ccrs.Mercator())
lon_min = -180; lon_max = -45
lat_min = -45; lat_max = 20
ax.set_extent([lon_min, lon_max, lat_min, lat_max],crs = ccrs.PlateCarree()) # min lon, max lon, min lat, max lat
ax.add_feature(cfeature.LAND, color='silver')
ax.add_feature(cfeature.LAKES, color='lightcyan')
ax.add_feature(cfeature.RIVERS, edgecolor='black')
ax.coastlines(resolution='50m', color='black', linewidth=1)
plt.title('SLA, fase neutra',color='black', size=14)

cm = plt.pcolormesh(lon,lat,np.nanmean(sla_n,axis=0),cmap = 'seismic',transform = ccrs.PlateCarree())
divider = make_axes_locatable(ax)
ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)

fig.add_axes(ax_cb)
plt.colorbar(cm,cax=ax_cb)
cm.set_clim(-.5,.5)

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.8, linestyle='dotted')
#gl.top_labels = False; gl.right_labels = False
gl.ylabels_right = False; gl.xlabels_top = False
gl.xlocator = mticker.FixedLocator([-180,-160,-140,-120,-100,-80,-60,-40])
gl.ylocator = mticker.FixedLocator([-40,-30,-20,-10,0,10,20])
gl.xformatter = LONGITUDE_FORMATTER; gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 13}
gl.ylabel_style = {'size': 13}
plt.show()


## Hago la diferencia entre fase negativa y fase positiva
# Grafico SLA fase negativa - SLA fase neutra 
fig = plt.figure(figsize=(10, 8), facecolor='white')
ax = plt.axes(projection=ccrs.Mercator())
lon_min = -180; lon_max = -45
lat_min = -45; lat_max = 20
ax.set_extent([lon_min, lon_max, lat_min, lat_max],crs = ccrs.PlateCarree()) # min lon, max lon, min lat, max lat
ax.add_feature(cfeature.LAND, color='silver')
ax.add_feature(cfeature.LAKES, color='lightcyan')
ax.add_feature(cfeature.RIVERS, edgecolor='black')
ax.coastlines(resolution='50m', color='black', linewidth=1)
plt.title('SLA, negativa - neutra',color='black', size=14)

cm = plt.pcolormesh(lon,lat,np.nanmean(sla_neg,axis=0) - np.nanmean(sla_n,axis=0),cmap = 'seismic',transform = ccrs.PlateCarree())
divider = make_axes_locatable(ax)
ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)

fig.add_axes(ax_cb)
plt.colorbar(cm,cax=ax_cb)
cm.set_clim(-.5,.5)

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.8, linestyle='dotted')
#gl.top_labels = False; gl.right_labels = False
gl.ylabels_right = False; gl.xlabels_top = False
gl.xlocator = mticker.FixedLocator([-180,-160,-140,-120,-100,-80,-60,-40])
gl.ylocator = mticker.FixedLocator([-40,-30,-20,-10,0,10,20])
gl.xformatter = LONGITUDE_FORMATTER; gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 13}
gl.ylabel_style = {'size': 13}

# Add a title to the plot
plt.show()


# Grafico SLA fase positiva - SLA fase neutra
fig = plt.figure(figsize=(10, 8), facecolor='white')
ax = plt.axes(projection=ccrs.Mercator())
#ax = plt.axes(projection=ccrs.Orthographic(central_longitude=-120.0, central_latitude=-30.0))
lon_min = -180; lon_max = -45
lat_min = -45; lat_max = 20
ax.set_extent([lon_min, lon_max, lat_min, lat_max],crs = ccrs.PlateCarree()) # min lon, max lon, min lat, max lat
ax.add_feature(cfeature.LAND, color='silver')
ax.add_feature(cfeature.LAKES, color='lightcyan')
ax.add_feature(cfeature.RIVERS, edgecolor='black')
ax.coastlines(resolution='50m', color='black', linewidth=1)
plt.title('SLA, positiva - neutra',color='black', size=14)

cm = plt.pcolormesh(lon,lat,np.nanmean(sla_pos,axis=0) - np.nanmean(sla_n,axis=0),cmap = 'seismic',transform = ccrs.PlateCarree())
divider = make_axes_locatable(ax)
ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)
plt.colorbar(cm,cax=ax_cb)
cm.set_clim(-.5,.5)

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.8, linestyle='dotted')
#gl.top_labels = False; gl.right_labels = False
gl.ylabels_right = False; gl.xlabels_top = False
gl.xlocator = mticker.FixedLocator([-180,-160,-140,-120,-100,-80,-60,-40])
gl.ylocator = mticker.FixedLocator([-40,-30,-20,-10,0,10,20])
gl.xformatter = LONGITUDE_FORMATTER; gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 13}
gl.ylabel_style = {'size': 13}
plt.show()


# Cargo los datos de percentil de humedad de suelo
# Fase negativa
grace = xr.open_dataset('/Users/Laura/OneDrive - cima.fcen.uba.ar/Oceanografia_satelital/2023/TP5/g4.timeAvgMap.GRACEDADM_CLSM025GL_7D_3_0_sfsm_inst.20221001-20221231.134W_45S_33W_17N.nc')

# Vemos la metadata del xarray
grace

lon = grace.lon.data-360 
lat = grace.lat.data
hum_neg = grace.GRACEDADM_CLSM025GL_7D_3_0_sfsm_inst.data


# Grafico 
fig = plt.figure(figsize=(8, 9), facecolor='white')
ax = plt.axes(projection=ccrs.Mercator())
#ax = plt.axes(projection=ccrs.Orthographic(central_longitude=-120.0, central_latitude=-30.0))
lon_min = -90; lon_max = -30
lat_min = -55; lat_max = 20
ax.set_extent([lon_min, lon_max, lat_min, lat_max],crs = ccrs.PlateCarree()) # min lon, max lon, min lat, max lat
ax.add_feature(cfeature.LAND, color='silver')
ax.add_feature(cfeature.LAKES, color='lightcyan')
ax.add_feature(cfeature.RIVERS, edgecolor='black')
ax.coastlines(resolution='50m', color='black', linewidth=1)
plt.title('% Humedad de Suelo, fase negativa',color='black', size=14)
cm = plt.pcolormesh(lon,lat,hum_neg,cmap = 'jet',transform = ccrs.PlateCarree())
divider = make_axes_locatable(ax)
ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)

fig.add_axes(ax_cb)
plt.colorbar(cm,cax=ax_cb)
cm.set_clim(0,100)

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.8, linestyle='dotted')
#gl.top_labels = False; gl.right_labels = False
gl.ylabels_right = False; gl.xlabels_top = False
gl.xlocator = mticker.FixedLocator([-90,-80,-70,-60,-50,-40,-30])
gl.ylocator = mticker.FixedLocator([-60,-50,-40,-30,-20,-10,0,10,20])
gl.xformatter = LONGITUDE_FORMATTER; gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 13}
gl.ylabel_style = {'size': 13}
plt.show()

# repetir el grafico de la fase positiva, de fase (positiva - neutra) y fase (negativa - neutra)

# A tener en cuenta: Realizamos un analisis sencillo de como impacta el fenomeno 
# climatico ENSO en la temperatura, altura del mar y humedad de suelo. Sin embargo, 
# para analizarlo correctamente deberiamos eliminar la tendencia a largo plazo de todas las variables, 
# eliminar el ciclo estacional (ya sea restando el anio climatologico o filtrando las variabilidades menores a un anio)
# Lo de la tendencia a largo plazo es evidente en la SLA. Mirar el siguiente video de 
# la variacion de la anomalia del nivel del mar (SLA) y de temperatura
# https://svs.gsfc.nasa.gov/30756/
# https://svs.gsfc.nasa.gov/30489
# https://www.youtube.com/watch?v=48flIFPAL78