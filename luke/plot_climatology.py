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

#### PATHS ###
path_to_data = '/Users/lfl/utcdw_data/'
path_to_obs = path_to_data
path_to_mod = path_to_data
path_to_savefig = f'{path_to_data}/figures/'
# path_to_obs = 'C:/Users/prate/Desktop/Climate Impacts Hackathon/data/era5_Geriza_Highres/'
# path_to_mod = 'C:/Users/prate/Desktop/Climate Impacts Hackathon/data/cesm/'
# path_to_savefig = 'C:/Users/prate/Desktop/Climate Impacts Hackathon/figures/'

#### FUNCTIONS ####
def plot(obs, mod, ds, title, unit):
    print(f'plotting {title} obs...')
    plt.figure(figsize=(6,3.75))
    clim = obs.groupby('time.dayofyear').mean(dim='time')
    max=obs.groupby('time.dayofyear').max(dim='time')
    min=obs.groupby('time.dayofyear').min(dim='time')
    time=np.arange(366)
    plt.plot(time, clim, label='TerraClim')
    plt.fill_between(time,min,max,alpha=0.5)
    
    print(f'plotting {title} model...')
    clim = mod.groupby('time.dayofyear').mean(dim='time')
    max=mod.groupby('time.dayofyear').max(dim='time')
    min=mod.groupby('time.dayofyear').min(dim='time')
    time=np.arange(365)
    plt.plot(time,clim,label='CESM Historical Data')
    plt.fill_between(time,min,max,alpha=0.5)

    print(f'plotting {title} downscaled...')
    clim = ds.groupby('time.dayofyear').mean(dim='time')
    max=ds.groupby('time.dayofyear').max(dim='time')
    min=ds.groupby('time.dayofyear').min(dim='time')
    time=np.arange(365)
    plt.plot(time,clim,label='CESM DBCCA')
    plt.fill_between(time,min,max,alpha=0.5)

    plt.grid()
    plt.xlim(0,365)
    plt.legend()
    plt.xlabel('Days since January 1st')
    plt.ylabel(f'{title} ({unit})')
    plt.title(f'Mean Gezira Scheme Climatology {start}-{end-1} ({unit})')
    plt.savefig(path_to_savefig+f'{title}_climatology.pdf')

#%%
# load obs
# initialize arrays for terraclim data
start, end = 1980, 2011
yrs = end - start
days = 365

print('loading terraclim...')
f_precip = xr.open_dataset(path_to_obs+f"precip.e5.daily.highres.{start}-{end-1}.nc")
f_rhumid = xr.open_dataset(path_to_obs+f"RH.e5.daily.highres.{start}-{end-1}.nc")
f_tmax = xr.open_dataset(path_to_obs+f"tmax.e5.daily.highres.{start}-{end-1}.nc")
f_tmin = xr.open_dataset(path_to_obs+f"tmin.e5.daily.highres.{start}-{end-1}.nc")

precip_obs=f_precip.tp #* 1000 # do unit conversion
precip_obs.attrs['units'] = 'mm/day' # update the unit attributes
precip_obs = mask_region(precip_obs,x_dim='lon',y_dim='lat',shapefile_name=f"{path_to_data}/Gezira.shp")
precip_obs=precip_obs.mean(dim=('lat','lon'))

rhumid_obs= mask_region(f_rhumid.relative_humidity,x_dim='lon',y_dim='lat',shapefile_name=f"{path_to_data}/Gezira.shp")
rhumid_obs= rhumid_obs.mean(dim=('lat','lon'))

tmax= mask_region(f_tmax.t2m,x_dim='lon',y_dim='lat',shapefile_name=f"{path_to_data}/Gezira.shp")
tmax=tmax.mean(dim=('lat','lon'))

tmin= mask_region(f_tmin.t2m,x_dim='lon',y_dim='lat',shapefile_name=f"{path_to_data}/Gezira.shp")
tmin=tmin.mean(dim=('lat','lon'))

temp_obs = (tmax + tmin) / 2

#%%
# load cesm
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
precip_mod=precip_mod.mean(dim=('lat','lon'))
print(np.shape(precip_mod))
print('meaned precip...')

rhumid_mod=mask_region_cesm(
    f_rhumid.hur) * 100
rhumid_mod=rhumid_mod.mean(dim=('lat','lon'))
print('meaned humidity...')

temp_mod=mask_region_cesm(
    f_temp.tas) - 273.15
temp_mod=temp_mod.mean(dim=('lat','lon'))
print('meaned temp...')

revap_obs=reference_crop_evapotranspiration(temp_obs,rhumid_obs)
revap_mod=reference_crop_evapotranspiration(temp_mod,rhumid_mod)

#%%
# load downscaled data
print('loading downscaled...')
f_precip_ds = xr.open_dataset(path_to_mod+
    f"pr_his_dbcca.nc")
f_temp_ds = xr.open_dataset(path_to_mod+
    f"tas_his_dbcca.nc")

precip_ds=mask_region(f_precip_ds.pr,
    shapefile_name=f"{path_to_data}Gezira.shp")
precip_ds=precip_ds.mean(dim=('lat','lon'))
print(np.shape(precip_ds))
print('meaned precip...')

temp_ds=mask_region(f_temp_ds.tas,
    shapefile_name=f"{path_to_data}Gezira.shp")
temp_ds=temp_ds.mean(dim=('lat','lon'))
print('meaned temp...')

revap_ds=reference_crop_evapotranspiration(temp_ds,rhumid_obs)


#%%
print('plotting...')
plot(precip_obs, precip_mod, precip_ds, 'Daily Precipitation', 'mm/day')
plot(temp_obs, temp_mod, temp_ds, 'Surface Air Temperature', 'C')
plot(rhumid_obs, rhumid_mod, rhumid_obs.convert_calendar('noleap'), 'Relative Humidity', '%')
plot(revap_obs, revap_mod, revap_obs.convert_calendar('noleap'), 'Reference Evapotransporation', 'mm/day')
plt.show()
