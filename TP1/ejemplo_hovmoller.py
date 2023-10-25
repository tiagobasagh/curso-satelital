#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 12:03:28 2020

@author: Laura 
@co-author: Ramiro 
"""
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
#pip install xarray
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

#creo una variable con las fechas para el eje x
dates = SST.time.values
fechas = dates[0:len(dates):int(len(dates)/10)]

#Hovmoller: Lon versus tiempo
#Opcion 1:
'''graficamos hovmoler a lo largo de la longitud'''
# Get times and make array of datetime objects
vtimes = data.time.values.astype('datetime64[ms]').astype('O') # array de tiempo en formato fecha
n_lat2 = 80 # corresponde a la latitud -34.875
_lat2 = lats[n_lat2].values
_sst3 = SST[:,n_lat2,:]

lo,tt = np.meshgrid(lons.values,vtimes)

#graficamos la serie temporal
plt.figure(figsize=(12,8), dpi= 60)
plt.pcolormesh(lo,tt,_sst3,cmap='jet',vmin= 10,vmax = 25)
plt.yticks(rotation=45,fontsize=20)
plt.yticks(fontsize=20)
plt.xlim(-55,-40)
plt.xlabel('Longitude',fontsize=20)
plt.ylabel('Time',fontsize=20)
plt.title('Hovmoller SST at  Lat = '+str(_lat2),fontsize=20)
plt.grid()
plt.legend()
plt.show()

#plt.savefig("TP1_Hov_a.jpg") # guarda el grafico en jpg

#Opción 2)
'''graficamos hovmoler a lo largo de la longitud'''
# Get times and make array of datetime objects
n_lat2 = 80
_sst3 = SST[:,n_lat2,:]

time = np.array(np.arange(1,1462,1), dtype=np.float32)

lo,tt=np.meshgrid(lons.values,time)
yticks = np.arange(0,len(dates),int(len(dates)/10))

#graficamos la serie temporal
plt.figure(figsize=(12,8), dpi= 60)
plt.pcolormesh(lo,tt,_sst3,cmap='jet',vmin= 10,vmax = 25)
plt.yticks(yticks,[str(f)[5:10] for f in fechas],rotation=45,fontsize=20)
plt.yticks(fontsize=20)
plt.xlim(-55,-40)
plt.xlabel('Longitude',fontsize=20)
plt.ylabel('Time',fontsize=20)
plt.title('Hovmoller SST at  Lat = '+str(_lat2),fontsize=20)
plt.grid()
plt.legend()
plt.show()



#Opcion 1)
'''graficamos hovmoler a lo largo de la latitud'''
# Get times and make array of datetime objects
n_lon2 = 94
_lon2 = lons[n_lon2].values
_sst3 = SST[:,:,n_lon2]

time = np.array(np.arange(1,1462,1), dtype=np.float32)

la,tt=np.meshgrid(lats.values,time)
xticks = np.arange(0,len(dates),int(len(dates)/10))

#graficamos la serie temporal
plt.figure(figsize=(12,8), dpi= 60)
plt.pcolormesh(tt,la,_sst3,cmap='jet',vmin= 5,vmax = 25)
plt.xticks(xticks,[str(f)[5:10] for f in fechas],rotation=45,fontsize=20)
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
plt.ylim(-55,-30)
plt.xlabel('Latitude',fontsize=20)
plt.ylabel('Time',fontsize=20)
plt.title('Hovmoller SST at  Lon = '+str(_lat2),fontsize=20)
plt.grid()
plt.legend()
plt.show()

#Opcion 2:
'''graficamos hovmoler a lo largo de la longitud'''
# Get times and make array of datetime objects
vtimes = data.time.values.astype('datetime64[ms]').astype('O')
n_lon2 = 94 # corresponde a 
_lon2 = lons[n_lon2].values
_sst3 = SST[:,:,n_lon2]

#time = np.array(np.arange(1,1462,1), dtype=np.float32)

la,tt=np.meshgrid(lats.values,vtimes)

#graficamos la serie temporal
plt.figure(figsize=(12,8), dpi= 60)
plt.pcolormesh(tt,la,_sst3,cmap='jet',vmin= 5,vmax = 25)
plt.xticks(rotation=45,fontsize=20)
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
plt.ylim(-55,-30)
plt.xlabel('Latitude',fontsize=20)
plt.ylabel('Time',fontsize=20)
plt.title('Hovmoller SST at  Lon = '+str(_lat2),fontsize=20)
plt.grid()
plt.legend()
plt.show()
