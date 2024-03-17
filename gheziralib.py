#%%
# import modules
import xarray as xr
from xclim import sdba
from xclim.core.calendar import convert_calendar
import xclim.indices as xci
import xclim.ensembles as xce
import numpy as np
import pandas as pd
import geopandas as gpd
import scipy.stats as stats
import seaborn as sns
import matplotlib.pyplot as plt
import ec3
import gcsfs
import zarr
import cartopy.crs as ccrs

# names of crops
SORGHUM='sorghum'
WHEAT='wheat'
GROUNDNUTS='groundnuts'
COTTON='cotton'

monx=[31,59,90,120,151,181,212,243,273,304,334,365]
KC = {}
KC[SORGHUM]=np.zeros(365)
KC[WHEAT]=np.zeros(365)
KC[GROUNDNUTS]=np.zeros(365)
KC[COTTON]=np.zeros(365)
# setup sorghum values 
KC[SORGHUM][:monx[5]], KC[SORGHUM][monx[10]:monx[11]]=0, 0
KC[SORGHUM][monx[5]:monx[6]]=0.39
KC[SORGHUM][monx[6]:monx[7]]=0.95
KC[SORGHUM][monx[7]:monx[8]]=1.10
KC[SORGHUM][monx[8]:monx[9]]=0.81
KC[SORGHUM][monx[9]:monx[10]]=1.96

# set up wheat values
KC[WHEAT][:monx[0]]=1.15
KC[WHEAT][monx[0]:monx[1]]=1.04
KC[WHEAT][monx[1]:monx[2]]=0.43
KC[WHEAT][monx[2]:monx[9]]=0
KC[WHEAT][monx[9]:monx[10]]=0.23
KC[WHEAT][monx[10]:monx[11]]=0.75

# set up ground nuts values
KC[GROUNDNUTS][:monx[3]], KC[GROUNDNUTS][monx[8]:monx[11]] = 0, 0
KC[GROUNDNUTS][monx[3]:monx[4]]=0.42
KC[GROUNDNUTS][monx[4]:monx[5]]=0.96
KC[GROUNDNUTS][monx[5]:monx[6]]=1.44
KC[GROUNDNUTS][monx[6]:monx[7]]=1.08
KC[GROUNDNUTS][monx[7]:monx[8]]=1.00

# set up cotton values
KC[COTTON][:monx[1]], KC[COTTON][monx[9]:monx[11]] = 0, 0
KC[COTTON][monx[1]:monx[2]]=0.35
KC[COTTON][monx[2]:monx[3]]=0.63
KC[COTTON][monx[3]:monx[4]]=0.91
KC[COTTON][monx[4]:monx[5]]=1.20
KC[COTTON][monx[5]:monx[6]]=1.20
KC[COTTON][monx[6]:monx[7]]=1.00
KC[COTTON][monx[7]:monx[8]]=0.80
KC[COTTON][monx[8]:monx[9]]=0.60

# KC[SORGHUM]=[np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,0.39,0.95,1.10,0.81,1.96,np.nan]
# KC[WHEAT]=[1.15,1.04,0.43,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,0.23,0.75]
# KC[GROUNDNUTS]=[np.nan,np.nan,np.nan,np.nan,0.42,0.96,1.44,1.08,1.00,np.nan,np.nan,np.nan]
# KC[COTTON]=[np.nan,np.nan,0.35,0.63,0.91,1.20,1.20,1.00,0.80,0.60,np.nan,np.nan]

def reference_crop_evapotranspiration(temp, rh):
    return 16 * temp / rh 
#rh.isel(time=slice(0,timen))
    
def crop_water_demand(crop, temp, rh):
    """
    Return the total water demand for a given crop,
    temperature (temp), and relative humidity (rh) during 
    a given time of the year (time as datetime)
    """
    kc = np.resize(KC[crop], np.shape(temp))
    return kc * reference_crop_evapotranspiration(temp, rh)

def irrigation_water_demand(crop, temp, rh, er, eff):
    """
    Returns the amount of water required from supplementary
    sources for irrigation of a given crop, at temperature temp,
    relative humidity rh, effective rainfall er, and irrigation 
    efficiency eff. Returns 0 if crop is not irrigated that day
    """
    cwd = crop_water_demand(crop, temp, rh)
    iwd_raw=(cwd-er)/eff
    iwd = iwd_raw.where(cwd!=0, other=0)
    return iwd

def mask_region(arr,x_dim='lon',y_dim='lat',shapefile_name="C:/Users/prate/Desktop/Climate Impacts Hackathon/data/Gezira_shapefile/Gezira.shp"):
    df_shapefile = gpd.read_file(shapefile_name, crs="epsg:4326")
    arr = arr.rio.set_spatial_dims(x_dim="lon", y_dim="lat")
    arr = arr.rio.write_crs("epsg:4326")
    return arr.rio.clip(df_shapefile.geometry.values, df_shapefile.crs, drop = False, invert = False)

def mask_region_cesm(arr):
    return arr.isel(lat=[3,4]).isel(lon=[2,3])

def load_var(
    var,member_id='r1i1p1f1',experiment_id='ssp370',lon_bnds=[30,40],
    lat_bnds=[10,20],years=np.arange(2070, 2100, 1),
    source_id='CESM2',table_id='day'
    ):
    # download the CESM2 model data from Google Cloud - this code is mostly the same as 6.1

    # URL for the Google Cloud model data catalog
    url_gcsfs_catalog = 'https://storage.googleapis.com/cmip6/cmip6-zarr-consolidated-stores.csv'
    # lon_bnds360 = [l + 360 for l in lon_bnds]

    # open the Google Cloud model data catalog with pandas
    df_catalog = pd.read_csv(url_gcsfs_catalog)

    # search for our selected model, both historical and SSP3-7.0 scenarios
    search_string = f"table_id == '{table_id}' & source_id == '{source_id}' & variable_id == '{var}'" # continue on the next line
    search_string += f" & experiment_id == '{experiment_id}'"
    df_search = df_catalog.query(search_string)

    # filter the search results further for the ensemble member we want to use
    df_search_member = df_search.query(f"member_id == '{member_id}'")

    # authenticate access to Google Cloud
    gcs = gcsfs.GCSFileSystem(token='anon')

    # get the URLs for each dataset and turn into zarr store mappers
    # url_hist = df_search_member[df_search_member.experiment_id == 'historical'].zstore.values[0]
    # mapper_hist = gcs.get_mapper(url_hist)
    # url_ssp3 = df_search_member[df_search_member.experiment_id == f'{experiment_id}'].zstore.values[0]
    url = df_search_member[df_search_member.experiment_id == f'{experiment_id}'].zstore.values[0]
    mapper = gcs.get_mapper(url)

    # future
    ds_raw = xr.open_zarr(mapper, consolidated = True)
    var_raw = ds_raw[f'{var}'].sel(lat = slice(*lat_bnds), lon = slice(*lon_bnds))

    # select time periods
    # var_hist_raw = var_hist_raw.sel(time = var_hist_raw.time.dt.year.isin(years_hist))
    var_raw = var_raw.sel(time = var_raw.time.dt.year.isin(years))
    return var_raw