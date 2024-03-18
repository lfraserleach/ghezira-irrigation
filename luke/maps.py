
#%%
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import xarray as xr
import rioxarray as rxr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import sys
os.chdir('/Users/lfl/google_drive/phd/utcdw_hackathon/ghezira-irrigation')
from gheziralib import *

path_to_data = '/Users/lfl/utcdw_data/'
path_to_obs = path_to_data
path_to_mod = path_to_data
path_to_savefig = f'{path_to_data}/figures/'

#%%
# load obs
start = 1980
end = 2011
f_precip = xr.open_dataset(path_to_obs+f"precip.e5.daily.highres.{start}-{end-1}.nc")
f_rhumid = xr.open_dataset(path_to_obs+f"RH.e5.daily.highres.{start}-{end-1}.nc")
f_tmax = xr.open_dataset(path_to_obs+f"tmax.e5.daily.highres.{start}-{end-1}.nc")
f_tmin = xr.open_dataset(path_to_obs+f"tmin.e5.daily.highres.{start}-{end-1}.nc")

#%%
# load historical cesm
print('loading cesm...')
f_precip = xr.open_dataset(
    path_to_mod+f"precip.cesm.daily.historical.{start}-{end-1}.nc")
f_rhumid = xr.open_dataset(
    path_to_mod+f"RH.cesm.daily.historical.{start}-{end-1}.nc")
f_temp = xr.open_dataset(
    path_to_mod+f"tas.cesm.daily.historical.{start}-{end-1}.nc")

# precip_mod.attrs['units'] = 'mm/day' # update the unit attributes
precip_mod=mask_region_cesm(
    f_precip.pr) * 3600 * 24

rhumid_mod=mask_region_cesm(
    f_rhumid.hur) * 100

temp_mod=mask_region_cesm(
    f_temp.tas) - 273.15

#%%
# load historcial downscaled data
print('loading downscaled...')
f_precip_ds = xr.open_dataset(path_to_mod+
    f"pr_his_dbcca.nc")
f_temp_ds = xr.open_dataset(path_to_mod+
    f"tas_his_dbcca.nc")

precip_ds=mask_region(f_precip_ds.pr,
    shapefile_name=f"{path_to_data}Gezira.shp")

temp_ds=mask_region(f_temp_ds.tas,
    shapefile_name=f"{path_to_data}Gezira.shp")

#%%
# load future cesm data
start_fut = 2070
end_fut = 2101
# load cesm
f_precip_370 = xr.open_dataset(
    path_to_mod+f"pr.cesm.daily.ssp370.{start_fut}-{end_fut-1}.nc")
f_rhumid_370 = f_rhumid
f_temp_370 = xr.open_dataset(
    path_to_mod+f"tas.cesm.daily.ssp370.{start_fut}-{end_fut-1}.nc")

precip_370=mask_region_cesm(
    f_precip_370.pr) * 3600 * 24
temp_370=mask_region_cesm(
    f_temp_370.tas) * 3600 * 24

#%%
# load future downscaled
print('loading downscaled...')
f_precip_370_ds = xr.open_dataset(
    path_to_mod+f"pr_370_dbcca.nc")
f_rhumid_370_ds = f_rhumid
f_temp_370_ds = xr.open_dataset(
    path_to_mod+f"tas_370_dbcca.nc")

precip_370_ds=mask_region(f_precip_370_ds.pr,
    shapefile_name=f"{path_to_data}Gezira.shp")

temp_370_ds=mask_region(f_temp_370_ds.tas,
    shapefile_name=f"{path_to_data}Gezira.shp")


#%%
# plot future downscaled precip
# precip_mod.attrs['units'] = 'mm/day' # update the unit attributes
cmap_pr = plt.get_cmap('BrBG')
min_pr = -1
max_pr = 1
step_pr = 0.2
lon_bnds = [32,34]
lat_bnds = [13.5,15.5]
lvl_pr = np.arange(min_pr,max_pr+step_pr,step_pr)
extend_pr = 'both'
precip_370_ds=mask_region(
    f_precip_370_ds.pr, shapefile_name=f"{path_to_data}Gezira.shp")*3600*24
precip_ds['time'] = precip_370_ds.time
lat_ds = precip_370_ds.lat
lon_ds = precip_370_ds.lon

proj = ccrs.PlateCarree()
f, ax = plt.subplots(subplot_kw={'projection': proj},
    constrained_layout=True,figsize=(5,5))
cb = ax.contourf(
    lon_ds, lat_ds, (precip_370_ds-precip_ds).mean(dim='time')/86400,
    cmap=cmap_pr,levels=lvl_pr,extend=extend_pr)
f.colorbar(cb,label='mm/day')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_xticks(np.arange(lon_bnds[0],lon_bnds[1]+1,0.5))
ax.set_yticks(np.arange(lat_bnds[0],lat_bnds[1]+1,0.5))
ax.set_extent([lon_bnds[0], lon_bnds[1], lat_bnds[0], lat_bnds[1]], 
    crs=ccrs.PlateCarree())
f.savefig(f"{path_to_savefig}pr_ds.pdf")

#%%
# plot future model precip
# precip_mod.attrs['units'] = 'mm/day' # update the unit attributes
cmap_pr = plt.get_cmap('BrBG')
min_pr = -1
max_pr = 1
step_pr = 0.2
lon_bnds = [32,34]
lat_bnds = [13.5,15.5]
lvl_pr = np.arange(min_pr,max_pr+step_pr,step_pr)
extend_pr = 'both'
# precip_370=mask_region(
#     f_precip_370.pr, shapefile_name=f"{path_to_data}Gezira.shp")*3600*24
precip_mod['time'] = precip_370.time
lat = precip_mod.lat
lon = precip_mod.lon

proj = ccrs.PlateCarree()
f, ax = plt.subplots(subplot_kw={'projection': proj},
    constrained_layout=True,figsize=(5,5))
cb = ax.contourf(
    lon, lat, (precip_370-precip_mod).mean(dim='time')/86400,
    cmap=cmap_pr,levels=lvl_pr,extend=extend_pr)
f.colorbar(cb,label='mm/day')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_xticks(np.arange(lon_bnds[0],lon_bnds[1]+1,0.5))
ax.set_yticks(np.arange(lat_bnds[0],lat_bnds[1]+1,0.5))
ax.set_extent([lon_bnds[0], lon_bnds[1], lat_bnds[0], lat_bnds[1]], 
    crs=ccrs.PlateCarree())
f.savefig(f"{path_to_savefig}pr_coarse.pdf")

#%%
# plot future downscaled precip
# precip_mod.attrs['units'] = 'mm/day' # update the unit attributes
cmap_tas = plt.get_cmap('BrBG')
min_tas = -1
max_tas = 1
step_tas = 0.2
lon_bnds = [32,34]
lat_bnds = [13.5,15.5]
lvl_tas = np.arange(min_tas,max_tas+step_tas,step_tas)
extend_tas = 'both'
temp_370_ds=mask_region(
    f_temp_370_ds.tas, shapefile_name=f"{path_to_data}Gezira.shp")*3600*24
temp_ds['time'] = temp_370_ds.time
lat_ds = temp_370_ds.lat
lon_ds = temp_370_ds.lon

proj = ccrs.PlateCarree()
f, ax = plt.subplots(subplot_kw={'projection': proj},
    constrained_layout=True,figsize=(5,5))
cb = ax.contourf(
    lon_ds, lat_ds, (temp_370_ds-temp_ds).mean(dim='time')/86400,
    cmap=cmap_pr,levels=lvl_pr,extend=extend_pr)
f.colorbar(cb,label='mm/day')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_xticks(np.arange(lon_bnds[0],lon_bnds[1]+1,0.5))
ax.set_yticks(np.arange(lat_bnds[0],lat_bnds[1]+1,0.5))
ax.set_extent([lon_bnds[0], lon_bnds[1], lat_bnds[0], lat_bnds[1]], 
    crs=ccrs.PlateCarree())
f.savefig(f"{path_to_savefig}pr_ds.pdf")

#%%
# plot future model precip
# precip_mod.attrs['units'] = 'mm/day' # update the unit attributes
cmap_pr = plt.get_cmap('BrBG')
min_pr = -1
max_pr = 1
step_pr = 0.2
lon_bnds = [32,34]
lat_bnds = [13.5,15.5]
lvl_pr = np.arange(min_pr,max_pr+step_pr,step_pr)
extend_pr = 'both'
# precip_370=mask_region(
#     f_precip_370.pr, shapefile_name=f"{path_to_data}Gezira.shp")*3600*24
precip_mod['time'] = precip_370.time
lat = precip_mod.lat
lon = precip_mod.lon

proj = ccrs.PlateCarree()
f, ax = plt.subplots(subplot_kw={'projection': proj},
    constrained_layout=True,figsize=(5,5))
cb = ax.contourf(
    lon, lat, (precip_370-precip_mod).mean(dim='time')/86400,
    cmap=cmap_pr,levels=lvl_pr,extend=extend_pr)
f.colorbar(cb,label='mm/day')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_xticks(np.arange(lon_bnds[0],lon_bnds[1]+1,0.5))
ax.set_yticks(np.arange(lat_bnds[0],lat_bnds[1]+1,0.5))
ax.set_extent([lon_bnds[0], lon_bnds[1], lat_bnds[0], lat_bnds[1]], 
    crs=ccrs.PlateCarree())
f.savefig(f"{path_to_savefig}pr_coarse.pdf")
