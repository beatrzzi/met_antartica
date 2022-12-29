# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 12:27:07 2022

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
    r'C:\Users\User\OneDrive\Documentos\Meteorologia\bia\antartica\dados\antartica-surface-levels2.nc'
    ).metpy.parse_cf()

file_1 = file_1.assign_coords(dict(
    longitude = (((file_1.longitude.values + 180) % 360) - 180))
    ).sortby('longitude')

#extent
lon_slice = slice(-180., -50.)
lat_slice = slice(-15., -80.)

#pega as lat\lon
lats = file_1.latitude.sel(latitude=lat_slice).values
lons = file_1.longitude.sel(longitude=lon_slice).values

for i in range(len(file_1.variables['time'])):
    
    #fluxo de calor latente
    latente = file_1['slhf'].metpy.sel(
           time = file_1.time[i],
           latitude=lat_slice, 
           longitude=lon_slice
           ).metpy.unit_array.squeeze()*0.00001
    
    #fluxo de calor sensivel
    sensivel = file_1['sshf'].metpy.sel(
           time = file_1.time[i],
           latitude=lat_slice, 
           longitude=lon_slice
           ).metpy.unit_array.squeeze()*0.00001
    
    #pressao reduzida ao nivel medio do mar
    pnmm = file_1.msl.metpy.sel(
        time = file_1.time[i], 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()*0.01
    
    #data
    vtime = file_1.time.data[i].astype('datetime64[ms]').astype('O')
    dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)
    
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
    
    # intevalos da pnmm
    intervalo_min2 = np.amin(np.array(pnmm))
    intervalo_max2 = np.amax(np.array(pnmm))
    interval_2 = 3              # de quanto em quanto voce quer que varie
    levels_2 = np.arange(intervalo_min2, intervalo_max2, interval_2)
    
    # intevalos do fluxo de calor lat
    intervalo_min3 = -8
    intervalo_max3 = 8
    interval_3 = 1 # de quanto em quanto voce quer que varie
    levels_3 = np.arange(intervalo_min3, intervalo_max3, interval_3)
    
    # intevalos do fluxo de calor sensivel
    intervalo_min4 = -19200
    intervalo_max4 = 4000
    interval_4 = 1000 # de quanto em quanto voce quer que varie
    levels_4 = np.arange(intervalo_min4, intervalo_max4, interval_4)
    
    # plota o fluxo de calor lat
    sombreado = ax.contourf(lons, 
                            lats, 
                            latente, 
                            cmap='PuOr', 
                            levels = levels_3, 
                            extend = 'both')
    
    # plota a imagem pressao
    contorno_1 = ax.contour(lons,
                          lats, 
                          pnmm, 
                          colors='black', 
                          linewidths=2, 
                          levels=levels_2)
    
    ax.clabel(contorno_1, 
              inline = 1, 
              inline_spacing = 1, 
              fontsize=20, 
              fmt = '%3.0f', 
              colors= 'black')
    
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
    plt.title('Fluxo de calor latente (J/m²).10⁵',
              fontweight='bold', 
              fontsize=35, 
              loc='left'
              )
    #previsao
    #plt.title('Valid time: {}'.format(vtime), fontsize=35, loc='right')
    #analise
    plt.title('Time: {}'.format(vtime), fontsize=35, loc='right')