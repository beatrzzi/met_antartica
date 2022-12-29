#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 12:03:04 2022

@author: bmiranda
"""
#importando bibliotecas
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.units import units
import numpy as np
import xarray as xr
import cartopy.io.shapereader as shpreader # Import shapefiles
from datetime import datetime, timedelta  # basicas datas e tipos de tempo
import cmocean

file_1 = xr.open_dataset(
    '/home/bmiranda/Desktop/ES2/bia-isa/dados/GFS_Global_0p25deg_ana_20221005_1200.grib2.nc4'
    ).metpy.parse_cf()

file_1 = file_1.assign_coords(dict(
    longitude = (((file_1.longitude.values + 180) % 360) - 180))
    ).sortby('longitude')

#
#extent
lon_slice = slice(-90., -10.,7)
lat_slice = slice(10., -70.,7)

#pega as lat/lon
lats = file_1.latitude.sel(latitude=lat_slice).values
lons = file_1.longitude.sel(longitude=lon_slice).values

#seta as variaveis
level_1 = 850 * units('hPa')

for i in range(len(file_1.variables['time'])):
    
    u = file_1['u-component_of_wind_isobaric'].metpy.sel(
        time = file_1.time[i], 
        vertical=level_1, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze().to('kt')
    
    v = file_1['v-component_of_wind_isobaric'].metpy.sel(
        time = file_1.time[i], 
        vertical=level_1, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze().to('kt')
    
    pnmm = file_1.Pressure_reduced_to_MSL_msl.metpy.sel(
        time = file_1.time[i], 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()* 0.01 * units.hPa/units.Pa
    
    t = file_1.Temperature_isobaric.metpy.sel(
        time = file_1.time[i],  
        vertical=level_1, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze().to('degC')
    
    vtime = file_1.time.data[i].astype('datetime64[ms]').astype('O')
    
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
    interval_2 = 4              # de quanto em quanto voce quer que varie
    levels_2 = np.arange(intervalo_min2, intervalo_max2, interval_2)
    
    # intevalos de temp
    intervalo_min3 = -25
    intervalo_max3 = 40
    interval_3 = 5 # de quanto em quanto voce quer que varie
    levels_3 = np.arange(intervalo_min3, intervalo_max3, interval_3)

    sombreado = ax.contourf(lons, 
                            lats, 
                            t, 
                            cmap=cmocean.cm.balance, 
                            levels = levels_3, 
                            extend = 'neither'
                            )
    
    plt.barbs(
        lons,
        lats,
        u,
        v,
        fill_empty=True,
        length=7,
        sizes=dict(emptybarb=0, height=0.6),
        barbcolor="black",
        barb_increments=dict(flag=50),)
    
    
    # plota a imagem pressao
    contorno = ax.contour(lons,
                          lats, 
                          pnmm, 
                          colors='black', 
                          linewidths=0.8, 
                          levels=levels_2
                          )
    
    # ax.clabel(contorno, 
    #           inline = 1, 
    #           inline_spacing = 1, 
    #           fontsize=20, 
    #           fmt = '%3.0f', 
    #           colors= 'black'
    #           )
    
    #adicionando shapefile
    shapefile = list(
        shpreader.Reader(
        '/home/bmiranda/Desktop/ES2/bia-isa/shapefiles/BR_UF_2021/BR_UF_2021.shp'
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
    plt.title('Temperatura (°C) em 850 hPa',
              fontweight='bold', 
              fontsize=35, 
              loc='left'
              )
    
    #previsao
    #plt.title('Valid Time: {}'.format(vtime), fontsize=35, loc='right')
    #analise
    plt.title('Análise: {}'.format(vtime), fontsize=30, loc='right')
    
    #--------------------------------------------------------------------------
    # Salva imagem
    plt.savefig(f'/home/bmiranda/Desktop/ES2/bia-isa/temperatura_850-{vtime}.png', bbox_inches='tight')