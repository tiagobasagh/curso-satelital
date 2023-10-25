#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 10:17:13 2019

@author: Ramiro
"""
#importamos las librerias que vamos a utilizar
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np


from mpl_toolkits.basemap import Basemap # si no funciona hacer el paso siguiente
#import os
#import conda
#conda_file_dir = conda.__file__
#conda_dir = conda_file_dir.split('lib')[0]
#proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
#os.environ["PROJ_LIB"] = proj_lib
#from mpl_toolkits.basemap import Basemap


import warnings
warnings.filterwarnings('ignore')

'''cargamos los datos'''

#nombre del archivo netCDF
filename = 'SST_noaa_OI_2014_2017.nc'

#utilizamos la libreria xarray para cargar los datos del archivo seleccionado
#modificar la ruta según donde hayan guardado el netCDF
data = xr.open_dataset('/Users/RamiroFerrari/OneDrive - cima.fcen.uba.ar/Documents/Cursos_y_clases/Satelital/2021/TP1 lectura datos netcdf-20210816/'+filename)

#visualisamos el contenido del netCDF
print(data)

#vemos que el archivo cuenta con una variable = sst
#esta variable cuenta con 3 coordenadas = tiempo, latitud y longitud.
# Por lo tanto, es una matriz de 3 dimensiones

#creamos una variable y le asignamos la matriz de sst
SST = data.sst
#vemos el contenido de sst, sus coordenadas y atributos
print(SST)
#esto nos devuelve en la primera linea : <xarray.DataArray 'sst' (time: 1461, lat: 92, lon: 120)>
#vemos que el orden de las coordenadas en la matriz es time,lat,lon..por lo tanto es una matriz de 1461x92x120
#otra forma de ver las dimensiones de la matriz es utilizando el atributo 'shape' :
print(SST.shape)

#es importante ver los atributos de la variable..por ejemplo, las unidades utilizadas :
print(SST.units)

'''graficamos el campo espacial de SST para un determinado tiempo'''

#creamos las variables lats y lons con los datos de las latitudes y longitudes de los datos descargados
lons = -1*(360-SST.lon) #SST.lon viene en formato [0°,360°] por lo tanto debo pasarlo al formato [-180°,180°] para graficar
lats = SST.lat

#seleccionamos un tiempo determinado para graficar la SST
#vemos el formato del tiempo
print(SST.time)
#en los atributos vemos que el delta_t es de 1 día..por lo tanto tenemos 1461 datos diarios de sst
#Por ejemplo para ver cual es el enesimo dia ejecutamos: SST.time[n].values
print(SST.time[0].values) # que vemos que corresponde al 1 de enero de 2014
print(SST.time[50].values) # que vemos que corresponde al 20 de febrero de 2014

#selecciono el primer dia de la serie y grafico el mapa de contornos
#necesito seleccionar el primer dato de tiempo, todos las latitudes y todas las longitudes:
t_sst = SST[0,:,:] #los dos puntos ':' significa que seleccionamos todos el array

#cargo los datos del mapa de la libreria 'basemap'
fig=plt.figure(figsize=(12,12), dpi= 60)
#selecciono los límites del mapa con los extremos de los datos de lat-lon
map = Basemap(llcrnrlon=lons[0],llcrnrlat=lats[0],urcrnrlon=lons[-1],urcrnrlat=lats[-1],projection='cyl',resolution='l')
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.drawmeridians(np.arange(0,360,3),labels=[True,False,False,True],linewidth=0.5,fontsize=20)
map.drawparallels(np.arange(-80,80,3),labels=[True,False,True,False],linewidth=0.5,fontsize=20)
map.fillcontinents(color='#ddaa66', lake_color='#0000ff')

#creamos la grilla de lat-lon con los valores obtenidos
xx, yy = np.meshgrid(lons, lats)
#selecciono el rango de valores de temperatura a graficar
clvls = np.arange(4,28,0.5)
#ploteo el mapa de contornos            
ax=map.contourf(xx,yy,t_sst,clvls,cmap='jet',extend='both')
#creo la barra de colores
cbar = map.colorbar(ax,location='bottom',pad="9%",format="%.1f")
cbar.set_label(u'(°C)',fontsize=20)
cbar.ax.tick_params(labelsize=20)
#le pongo título al gráfico
plt.title('SST',fontsize=20)
plt.show()

'''graficamos la serie temporal de SST para una determinada latitud y longitud'''

#seleccinamos una latitud y longitud dentro de los valores de lat-lon del netCDF
#SST.lat[n].values .... por ejemplo:
n_lat = 10 # seleccionamos entre [0,91]. Recordemos para ver la dimension de latitudes =  lats.shape
n_lon = 14 # seleccionamos entre [0,119]. Recordemos para ver la dimension de latitudes = lons.shape
_lat = lats[n_lat].values #corresponde a la latitud -52.375
_lon = lons[n_lon].values #corresponde a la longitud -66.375

#creo la variable con la serie temporal de sst para esa ubicación
_sst = SST[:,n_lat,n_lon]

#creo una variable con las fechas para el eje x
dates = _sst.time.values
fechas = dates[0:len(dates):int(len(dates)/10)]
xticks = np.arange(0,len(dates),int(len(dates)/10))

#graficamos la serie temporal
plt.figure(figsize=(12,8), dpi= 60)
plt.plot(_sst)
plt.xticks(xticks,[str(f)[5:10] for f in fechas],rotation=45,fontsize=20)
plt.xticks(fontsize=20)
plt.xlabel('fecha',fontsize=20)
plt.ylabel('°C',fontsize=20)
plt.title('SST --- Lat = '+str(_lat)+' --- Lon = '+str(_lon),fontsize=20)
plt.show()

'''graficamos dos series en un mismo panel'''
n_lat2 = 80
n_lon2 = 94
_lat2 = lats[n_lat2].values
_lon2 = lons[n_lon2].values
_sst2 = SST[:,n_lat2,n_lon2]

#graficamos la serie temporal
plt.figure(figsize=(12,8), dpi= 60)
plt.plot(_sst, label='SST1')
plt.plot(_sst2, label='SST2')
plt.xticks(xticks,[str(f)[5:10] for f in fechas],rotation=45,fontsize=20)
plt.xticks(fontsize=20)
plt.xlabel('fecha',fontsize=20)
plt.ylabel('°C',fontsize=20)
plt.title('SST1 --- Lat = '+str(_lat)+' --- Lon = '+str(_lon)+'\n SST2 --- Lat = '+str(_lat2)+' --- Lon = '+str(_lon2),fontsize=20)
plt.grid()
plt.legend()
plt.show()