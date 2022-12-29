# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 13:07:10 2022
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
lon_slice = slice(-180., -50.,7)
lat_slice = slice(-15., -80.,7)

#pega as lat\lon
lats = file_1.latitude.sel(latitude=lat_slice).values
lons = file_1.longitude.sel(longitude=lon_slice).values

for i in range(len(file_1.variables['time'])):

    altura_sig_onda = file_1['swh'].metpy.sel(
           time = file_1.time[i],
           latitude=lat_slice, 
           longitude=lon_slice
           ).metpy.unit_array.squeeze()
    
    u10 = file_1['u10'].metpy.sel(
        time = file_1.time[i], 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze().to('kt')
    
    v10 = file_1['v10'].metpy.sel(
        time = file_1.time[i],
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze().to('kt')
    
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
    
    # intevalos da altura significativa de onda
    intervalo_min2 = 0.5
    intervalo_max2 = 15
    interval_2 = 1.0
    levels_2 = np.arange(intervalo_min2, intervalo_max2, interval_2)
    
    sombreado = ax.contourf(lons, 
                            lats, 
                            altura_sig_onda, 
                            cmap= 'magma_r', 
                            levels = levels_2, 
                            extend = 'neither'
                            )
    
    plt.barbs(
        lons,
        lats,
        u10,
        v10,
        fill_empty=True,
        length=7,
        sizes=dict(emptybarb=0, height=0.6),
        barbcolor="black",
        barb_increments=dict(flag=50),)
    
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
    plt.title('Altura significativa de onda (m) e vento a 10m',
              fontweight='bold', 
              fontsize=32, 
              loc='left'
              )
    
    #previsao
    #plt.title('Valid Time: {}'.format(vtime), fontsize=35, loc='right')
    #analise
    plt.title('Time: {}'.format(vtime), fontsize=32, loc='right')
    
    #--------------------------------------------------------------------------
    # Salva imagem
    #plt.savefig(f'/altura_sig_onda_{vtime}.png', bbox_inches='tight')
    
    
