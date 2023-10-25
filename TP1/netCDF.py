
# coding: utf-8

# In[1]:

import numpy as np
import netCDF4
import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')


# In[6]:

ncfile = netCDF4.Dataset('/Users/RamiroFerrari/OneDrive - cima.fcen.uba.ar/Documents/Cursos_y_clases/Satelital/2021/TP1 lectura datos netcdf-20210816/air.sig995.2012.nc')  # Open file

ncfile #esto me lista las cosas que tiene el archivo

#leo las distintas variables:
#air = ncfile.variables["air"]   # Get Variable object for "temperature"
#print(type(air))
#air_data = air[:]          # Extract NumPy data array from Variable

#leo la variable directamente:
air = ncfile.variables["air"][:]   # Get Variable object for "temperature"

lon = ncfile.variables["lon"][:]   # Get Variable object for "temperature"
print(type(lon))
#lon_data = lon[:]          # Extract NumPy data array from Variable
#otra forma mas corta:
#lon = ncfile.variables["lon"][:]

lat = ncfile.variables["lat"][:]   # Get Variable object for "temperature"
print(type(lat))
#lat_data = lat[:]          # Extract NumPy data array from Variable

time = ncfile.variables["time"][:]   # Get Variable object for "temperature"
print(type(time))
#time_data = time[:]          # Extract NumPy data array from Variable

#si quiero ver uno de los atributos:
#T.units                # Units attribute

ncfile.close()

air.shape

#probemos una figura:
plt.contourf(lon, lat, air[1,:,:], 10, cmap=plt.cm.rainbow)
plt.axis([0, 360, -90, 90])
plt.colorbar(orientation="horizontal")

cont = plt.contour(lon, lat, air[1,:,:], 10)   # Save object
plt.clabel(cont)
plt.show()

import mpl_toolkits.basemap as bm

mapproj = bm.Basemap(projection='cyl',
                     llcrnrlat=-60.0, llcrnrlon=0.0,
                     urcrnrlat=60.0, urcrnrlon=360.0)
mapproj.drawcoastlines()
mapproj.drawparallels(np.array([-60, -40, -20, 0, 20, 40, 60]), labels=[1,0,0,0])
mapproj.drawmeridians(np.array([-180, -90, 0, 90, 180]), labels=[0,0,0,1])

lonall, latall = np.meshgrid(lon, lat)
lonproj, latproj = mapproj(lonall, latall)
plt.contour(lonproj, latproj, air[1,:,:], 10)
plt.show()













# In[7]:

alt = ncfile.variables['altitude']
temp = ncfile.variables['temperature']
plt.plot(temp, alt)
plt.grid()
plt.title("Temperature as a function of Altitude")
plt.ylabel("Altitude [" + alt.units + "]")
plt.xlabel("Temperature [" + temp.units + "]")
plt.show()

ncfile.close()


# In[8]:

ncfile = netCDF4.Dataset('../Hands-On/data/air.mon.ltm.nc')  # Open file
air_temp = ncfile.variables["air"]   # Get Variable object for "air"

print(type(air_temp))

temp_data = air_temp[:]       # Extract NumPy data array from Variable

print(type(temp_data))
print(air_temp.ncattrs())     # List all variable attributes

air_temp.units               # Units attribute


# In[9]:

ncf = netCDF4.Dataset("../Hands-On/data/air.mon.ltm.nc", "r")
lon = ncf.variables["lon"][:]
lat = ncf.variables["lat"][:]
temp_data = ncf.variables["air"][0,0,:,:]
temp_data.shape

plt.contourf(lon, lat, temp_data, 10, cmap=plt.cm.rainbow)
plt.axis([0, 360, -90, 90])
plt.colorbar(orientation="horizontal")

cont = plt.contour(lon, lat, temp_data, 10)   # Save object
plt.clabel(cont)
plt.show()


# In[10]:

ncf = netCDF4.Dataset("../Hands-On/data/air.mon.ltm.nc", "r")
level_data = ncf.variables["level"][:]
lat_data = ncf.variables["lat"][:]

# Extract data for month_index 0, all levels, all lats, lon_index 0
temp_data = ncf.variables["air"][0,:,:,0]
temp_data.shape

plt.contourf(lat_data, level_data, temp_data, 10, cmap=plt.cm.rainbow)
plt.ylim(1000, 0)
plt.colorbar()

plt.contour(lat_data, level_data, temp_data, 10)


# In[11]:

import mpl_toolkits.basemap as bm

ncf = netCDF4.Dataset("../Hands-On/data/air.mon.ltm.nc", "r")
lon = ncf.variables["lon"][:]
lat = ncf.variables["lat"][:]
temp_data = ncf.variables["air"][0,0,:,:]


mapproj = bm.Basemap(projection='cyl',
                     llcrnrlat=-60.0, llcrnrlon=0.0,
                     urcrnrlat=60.0, urcrnrlon=360.0)
mapproj.drawcoastlines()
mapproj.drawparallels(np.array([-60, -40, -20, 0, 20, 40, 60]), labels=[1,0,0,0])
mapproj.drawmeridians(np.array([-180, -90, 0, 90, 180]), labels=[0,0,0,1])

lonall, latall = np.meshgrid(lon, lat)
lonproj, latproj = mapproj(lonall, latall)
plt.contour(lonproj, latproj, temp_data)


# In[9]:

def us_plot(month):

    ncf = netCDF4.Dataset("../Hands-On/data/air.mon.ltm.nc", "r")
    lon = ncf.variables["lon"][:]
    lat = ncf.variables["lat"][:]
    temp_data = ncf.variables["air"][month-1,0,:,:]

    mapproj = bm.Basemap(projection='lcc', lat_0=30, lon_0=250,
                     llcrnrlat=20.0, llcrnrlon=230,
                     urcrnrlat=55.0, urcrnrlon=310)
    mapproj.drawcoastlines()
    lonall, latall = np.meshgrid(lon, lat)
    lonproj, latproj = mapproj(lonall, latall)
    plt.contourf(lonproj, latproj, temp_data, 20)
    plt.colorbar()
    plt.title("Month "+str(month))

us_plot(12)


# In[10]:

us_plot(5)


# Documentación oficial de netcdf4-python:
# 
# https://unidata.github.io/netcdf4-python/
# 
# Galeria de ejemplos de basemaps:
# 
# http://matplotlib.org/basemap/users/examples.html

# In[11]:

#Este css es trabajo de @LorenaABarba y su grupo
from IPython.core.display import HTML
css_file = 'css/personal.css'
HTML(open(css_file, "r").read())

