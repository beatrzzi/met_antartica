#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 10:21:15 2022
@author: beatrzzi
"""

#importando bibliotecas
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.colors 
import metpy.calc as mpcalc
from metpy.units import units
import numpy as np
import xarray as xr
import cartopy.io.shapereader as shpreader # Import shapefiles
from datetime import datetime, timedelta  # basicas datas e tipos de tempo
import cmocean

#dataset
file_1 = xr.open_dataset(
    r'C:\Users\User\OneDrive\Documentos\Meteorologia\bia\antartica\dados\antartica-pressure-levels.nc'
    ).metpy.parse_cf()

file_1 = file_1.assign_coords(dict(
    longitude = (((file_1.longitude.values + 180) % 360) - 180))
    ).sortby('longitude')

file_2 = xr.open_dataset(
    r'C:\Users\User\OneDrive\Documentos\Meteorologia\bia\antartica\dados\antartica-surface-levels2.nc'
    ).metpy.parse_cf()

file_2 = file_2.assign_coords(dict(
    longitude = (((file_2.longitude.values + 180) % 360) - 180))
    ).sortby('longitude')

#extent
lon_slice = slice(-180., -50.)
lat_slice = slice(-15., -80.)

#pega as lat/lon
lats = file_1.latitude.sel(latitude=lat_slice).values
lons = file_1.longitude.sel(longitude=lon_slice).values


#seta as variaveis
level_1 = 850
level_2 = 500
level_3 = 1000 

# cria uma escala de cores:
colors = ["#0000CD","#FF0000"]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
cmap.set_over("#FF0000")
cmap.set_under("#0000CD")


for i in range(len(file_1.variables['time'])):
    
    geopotencial_500 = file_1.z.metpy.sel(
        time = file_1.time[i], 
        level=level_2, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()*0.1
    
    u = file_1['u'].metpy.sel(
        time = file_1.time[i], 
        level=level_1, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()
    
    v = file_1['v'].metpy.sel(
        time = file_1.time[i], 
        level=level_1, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()
    
    t = file_1.t.metpy.sel(
        time = file_1.time[i],  
        level=level_1, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()

    
    #data
    vtime = file_1.time.data[i].astype('datetime64[ms]').astype('O')
    
    dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)
    
    adv_temp = mpcalc.advection(t, u=u, v=v, dx=dx, dy=dy, x_dim=- 1, y_dim=- 2)*1000
    
    # escolha o tamanho do plot em polegadas (largura x altura)
    plt.figure(figsize=(25,25))
    
    # usando a proje????o da coordenada cilindrica equidistante 
    ax = plt.axes(projection=ccrs.PlateCarree())
    gl = ax.gridlines(crs=ccrs.PlateCarree(),
                      color='gray',
                      alpha=1.0, 
                      linestyle='--', 
                      linewidth=0.5, 
                      xlocs=np.arange(-180, 180, 10), 
                      ylocs=np.arange(-90, 90, 10), 
                      draw_labels=True
                      )
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 29, 'color': 'black'}
    gl.ylabel_style = {'size': 29, 'color': 'black'}
    
 
    # intevalos da adv de temp
    intervalo_min3 = -1.5
    intervalo_max3 = 1.6
    interval_3 = 0.1              # de quanto em quanto voce quer que varie
    levels_3 = np.arange(intervalo_min3, intervalo_max3, interval_3)
    
    # intevalos da geopotencial
    intervalo_min1 = np.amin(np.array(geopotencial_500))
    intervalo_max1 = np.amax(np.array(geopotencial_500))
    interval_1 = 50              # de quanto em quanto voce quer que varie
    levels_1 = np.arange(intervalo_min1, intervalo_max1, interval_1)

    
    # plota a imagem adv de temp
    sombreado = ax.contourf(lons, 
                            lats, 
                            adv_temp, 
                            cmap='seismic', 
                            levels = levels_3, 
                            extend = 'neither'
                            )
    
    # plota a imagem geopotencial
    contorno_2 = ax.contour(lons,
                          lats, 
                          geopotencial_500,
                          cmap=cmap, 
                          linewidths=3, 
                          linestyles='dashed',
                          levels=levels_1
                          )
    
    ax.clabel(contorno_2, 
              inline = 1, 
              inline_spacing = 1, 
              fontsize=24,
              fmt = '%3.0f', 
              colors= 'black'
             )
    
    #adicionando shapefile
    shapefile = list(
        shpreader.Reader(
        r'C:\Users\User\OneDrive\Documentos\Meteorologia\es2\shapefiles\BR_UF_2021\BR_UF_2021.shp'
        ).geometries()
        )
    
    ax.add_geometries(
        shapefile, ccrs.PlateCarree(), 
        edgecolor = 'black', 
        facecolor='none', 
        linewidth=0.5
        )
    
    # adiciona continente e bordas
    ax.add_feature(cfeature.LAND)
    ax.coastlines(resolution='10m', color='black', linewidth=3)
    ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=3)
    
    # adiciona legenda 
    barra_de_cores = plt.colorbar(sombreado, 
                                  orientation = 'horizontal', 
                                  pad=0.04, 
                                  fraction=0.04
                                  )
    font_size = 20 # Adjust as appropriate.
    barra_de_cores.ax.tick_params(labelsize=font_size)
    
   
    # Add a title
    plt.title('Adv de temperatura (K/s) em 850 hPa e\nGeopotencial em 500hPa',
              fontweight='bold', 
              fontsize=30, 
              loc='left'
              )
    
    #previsao
    plt.title('Time: {}'.format(vtime), fontsize=25, loc='right')
    #analise
    #plt.title('An??lise: {}'.format(vtime), fontsize=30, loc='right')
    
    #--------------------------------------------------------------------------
    # Salva imagem
    #plt.savefig(f'/adveccao_de_temperatura_{vtime}.png', bbox_inches='tight')
