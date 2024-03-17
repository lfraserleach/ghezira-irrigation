import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import xarray as xr
import rioxarray as rxr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import sys
sys.path.insert(1, 'C:/Users/prate/Desktop/Climate Impacts Hackathon/ghezira-irrigation/')
from gheziralib import *

#### PATHS ###
path_to_downscaled = 'C:/Users/prate/Desktop/Climate Impacts Hackathon/data/downscaled/'
path_to_obs = 'C:/Users/prate/Desktop/Climate Impacts Hackathon/data/era5_Geriza_HighRes/'

# SETTINGS
forcing = '370' # or '370'
bias_correction='dbcca' # or 'bcca'
efficiency=0.65 #0.55

f_precip = xr.open_dataset(path_to_downscaled+f"pr_{forcing}_{bias_correction}.nc")
# f_rhumid = xr.open_dataset(path_to_downscaled+f"hur.{forcing}.{bias_correction}.nc")
f_rhumid= xr.open_dataset(path_to_obs+f"RH.e5.daily.highres.1980-2010.nc")
f_temp = xr.open_dataset(path_to_downscaled+f"tas_{forcing}_{bias_correction}.nc")
                         
precip_mod=mask_region(f_precip.pr) 

rhumid_mod=mask_region(f_rhumid.relative_humidity) 
rhumid_mod=rhumid_mod.convert_calendar("noleap")
rhumid_mod['time']=precip_mod.time

temp_mod=mask_region(f_temp.tas) #+273.15
# print(temp_mod.dims, temp_mod.shape)
print('Mean daily precipitation:', precip_mod.mean().data, 'mm/day')
revap_mod=reference_crop_evapotranspiration(temp_mod,rhumid_mod)
print('Mean daily reference evapotranspiration:', revap_mod.mean().data, 'mm/day')
CROPS = [WHEAT, SORGHUM, GROUNDNUTS, COTTON]
LANDUSE = [0.39, 0.37, 0.15, 0.09]
area = 750000 * 10000 # convert ha to m^2
daily_water_use = np.zeros_like(temp_mod)
for i in range(4):
    cwd = crop_water_demand(CROPS[i], temp_mod, rhumid_mod)
    print(f'Mean daily {CROPS[i]} crop water deman:', cwd.mean().data, 'mm/day')
    iwd = irrigation_water_demand(CROPS[i], temp_mod, rhumid_mod, precip_mod, efficiency)
    print(f'Mean daily {CROPS[i]} irrigation water demand:', iwd.mean().data, 'mm/day')
    daily_crop_water_use = iwd * area * LANDUSE[i] /1000 # convert from mm/day m^3/day
    print(f'Mean daily {CROPS[i]} water volume demand:', daily_crop_water_use.mean().data * 1000, 'L/day')
    daily_water_use += daily_crop_water_use
    
daily_water_use=daily_water_use.mean(dim=('lat','lon')) 

mean_annual_water_use = daily_water_use.groupby('time.year').sum(dim='time').mean(dim='year') * 1000 # convert cubic m to cubic L
print(f'Mean annual total water usage:', mean_annual_water_use.data, 'L/year')
    


# print(daily_water_use.groupby('time.year').dims) #
# annual_water_use = daily_water_use.groupby('time.year').sum(dim='time') # cubic m / year
# print(annual_water_use.dims)
# mean_daily_water_use = daily_water_use.mean(dim='year')
# print('Mean annual water requirement:', mean_annual_water_use.data / (1000)**3, 'cubic km / year') # cubic km / year

# bn_flow_rate=16500 #1548 ##cubic m / s
# print('Days of blue nile flow required:', mean_annual_water_use.data / (bn_flow_rate * 3600 * 24), 'days') # cubic km / year