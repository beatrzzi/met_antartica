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

#pega as lat\lon
lats = file_1.latitude.sel(latitude=lat_slice).values
lons = file_1.longitude.sel(longitude=lon_slice).values

#seta as variaveis
level_1 = 1000
level_2 = 500

 
for i in range(len(file_1.variables['time'])):
    
    geopotencial_1000 = file_1.z.metpy.sel(
        time = file_1.time[i], 
        level=level_1, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()*0.1
    
    geopotencial_500 = file_1.z.sel(
        time = file_1.time[i], 
        level=level_2, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()*0.1
    
    pnmm = file_2.msl.metpy.sel(
        time = file_1.time[i], 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()*0.01
    
    #data
    vtime = file_1.time.data[i].astype('datetime64[ms]').astype('O')
    
    dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)
    
    espessura = geopotencial_500 - geopotencial_1000
    
    # escolha o tamanho do plot em polegadas (largura x altura)
    plt.figure(figsize=(25,25))
    
    # usando a projeção da coordenada cilindrica equidistante 
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
    
    # cria uma escala de cores:
    colors = ["#2d001c", "#5b0351", "#780777", "#480a5e", "#1e1552", 
              "#1f337d", "#214c9f", "#2776c6", "#2fa5f1", "#1bad1d", 
              "#8ad900", "#ffec00", "#ffab00", "#f46300", "#de3b00", 
              "#ab1900", "#6b0200", '#3c0000']
    
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
    cmap.set_over('#3c0000')
    cmap.set_under('#28000a')
    
    # intevalos da espessura
    intervalo_min2 = 4900
    intervalo_max2 = 5900
    interval_2 = 25             # de quanto em quanto voce quer que varie
    levels_2 = np.arange(intervalo_min2, intervalo_max2, interval_2)
    
    # intevalos da pnmm
    intervalo_min1 = np.amin(np.array(pnmm))
    intervalo_max1 = np.amax(np.array(pnmm))
    interval_1 = 3              # de quanto em quanto voce quer que varie
    levels_1 = np.arange(intervalo_min1, intervalo_max1, interval_1)
    
    # plota a imagem espessura
    sombreado = ax.contourf(lons, 
                            lats, 
                            espessura, 
                            cmap=cmap, 
                            levels = levels_2, 
                            extend = 'both'
                            )
    
    # plota a imagem pressao
    contorno_1 = ax.contour(lons,
                          lats, 
                          pnmm, 
                          colors='black', 
                          linewidths=2, 
                          levels=levels_1
                          )
    
    ax.clabel(contorno_1, 
              inline = 1, 
              inline_spacing = 1, 
              fontsize=20, 
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
    
    # adiciona mascara de terra
    ax.add_feature(cfeature.LAND)
    # adiciona continente e bordas
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
    plt.title('Espessura (500-1000) hPa',
              fontweight='bold', 
              fontsize=35, 
              loc='left'
              )
    
    #previsao
    #plt.title('Valid time: {}'.format(vtime), fontsize=35, loc='right')
    #analise
    plt.title('Time: {}'.format(vtime), fontsize=35, loc='right')
    
    #--------------------------------------------------------------------------
    # Salva imagem
    #plt.savefig(f'\espessura_.png', bbox_inches='tight')
    
