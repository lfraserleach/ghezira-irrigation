import xarray as xr
import numpy as np

class Simulation:
    def __init__(self,runid,runmode):
        self.runid = runid
        self.runmode = runmode

class VariableCtlPrt:
    """
    Class to store variable data for control and perturbation experiments
    """
    def __init__(self, var_ctl=None, var_prt=None):
        self.ctl = var_ctl
        self.prt = var_prt

class Variable2Simulations: 
    """
    Class to store variable data for bulk and pam control and perturbation 
    experiments
    """
    def __init__(self, var_sim1=None, var_sim2=None):
        self.sim1 = var_sim1
        self.sim2 = var_sim2
    
    def __add__(self, other):
        var_sim1 = self.sim1 + other.sim1
        var_sim2 = self.sim2 + other.sim2
        return Variable2Simulations(var_sim1, var_sim2)
    
    def __sub__(self, other):
        var_sim1 = self.sim1 - other.sim1
        var_sim2 = self.sim2 - other.sim2
        return Variable2Simulations(var_sim1, var_sim2)

    def __neg__(self):
        var_sim1 = -self.sim1
        var_sim2 = -self.sim2
        return Variable2Simulations(var_sim1, var_sim2)
    
    def load_sim1(
            self, var_name, var_dir, var_file, yearly=False, yri=None,
            yrf=None):
        """If yearly is False, var_file is the full filename. If yearly is
        True, var_file has format 
        {var_name}_{table}_{model}_{runid}_{ripf}_gn, with the years to be
        filled in with the provided yri and yrf."""
        if yearly:
            self.sim1 = load_var_yearly(var_name, var_dir, var_file, yri, yrf)
        else:
            self.sim1 = xr.open_dataset(f'{var_dir}/{var_file}')[var_name]

    def load_sim2(
            self, var_name, var_dir, var_file, yearly=False, yri=None,
            yrf=None):
        """If yearly is False, var_file is the full filename. If yearly is
        True, var_file has format 
        {var_name}_{table}_{model}_{runid}_{ripf}_gn, with the years to be
        filled in with the provided yri and yrf."""
        if yearly:
            self.sim2 = load_var_yearly(var_name, var_dir, var_file, yri, yrf)
        else:
            self.sim2 = xr.open_dataset(f'{var_dir}/{var_file}')[var_name]

    def load_sim1_zonalmean(
            self, var_name, var_dir, var_file, yearly=False, yri=None,
            yrf=None):
        if yearly:
            da = load_var_yearly(var_name, var_dir, var_file, yri, yrf)
            self.sim1 = da.mean(dim='lon')
            da.close()
        else:
            da = xr.open_dataset(f'{var_dir}/{var_file}')[var_name]
            self.sim1 = da.mean(dim='lon')
            da.close()


    def load_sim2_zonalmean(
            self, var_name, var_dir, var_file, yearly=False, yri=None,
            yrf=None):
        if yearly:
            da = load_var_yearly(var_name, var_dir, var_file, yri, yrf)
            self.sim2 = da.mean(dim='lon')
            da.close()
        else:
            da = xr.open_dataset(f'{var_dir}/{var_file}')[var_name]
            self.sim2 = da.mean(dim='lon')
            da.close()

class PlotParameters:
    """"
    Class to store parameters for plotting
    """
    def __init__(
        self, cmap, extend, min, max, levels, label=None, mult_factor=1):
        self.cmap = cmap
        self.extend = extend
        self.min = min
        self.max = max
        self.levels = levels
        self.mult_factor = mult_factor
        self.label = label

def output_dir_user(user, runid, model, config):
    output_dir = (f'/space/hall5/sitestore/eccc/crd/ccrn/users/{user}' 
        + f'/canesm_runs/{runid}/data/nc_output/CMIP6/CCCma/CCCma/'
        + f'{model}-{runid}/{config}')
    return output_dir

def output_dir_cmip(mip, runid, model):
    output_dir = ('/space/hall5/sitestore/eccc/crd/ccrn/model_output/'
        + f'CMIP6/final/CMIP6/{mip}/CCCma/{model}/{runid}')
    return output_dir

def output_dir_niagara(runid):
    output_dir = ('/project/p/pjk/lfl/pam_data/hpfx.collab.science.gc.ca/'
    + f'~rud001/4Luke/{runid}/')
    return output_dir

def output_dir_niagara_cmip(runid, model):
    output_dir = ('/project/p/pjk/lfl/pam_data/hpfx.collab.science.gc.ca/'
    + f'~rud001/4Luke/{model}/{runid}/')
    return output_dir

def subdir_path_cmor(ripf, table, var_name, version='gn/v20190429'):
    """Returns path from data directory to output files assuming
    CMOR directory structure"""
    path = f'{ripf}/{table}/{var_name}/{version}'
    return path

def subdir_path_niagara(table, var_name, version='gn/v20190429'):
    """Returns path from data directory to output files assuming
    CMOR directory structure"""
    path = f'{table}/{var_name}/{version}/'
    return path

def subdir_path_niagara_cmip(ripf, table, var_name, version='gn/v20190429'):
    """Returns path from data directory to output files assuming
    CMOR directory structure"""
    path = f'{ripf}/{table}/{var_name}/{version}/'
    return path

def fname_user(var_name, table, model, runid, config, ripf, yri, yrf):
    fname = (
        f'{var_name}_{table}_{model}-{runid}_{config}_{ripf}_gn_' 
        + f'{yri}01-{yrf}12.nc')
    return fname

def fname_noyears(var_name, table, model, runid, config, ripf, yri, yrf):
    """For loading yearly files"""
    fname = (
        f'{var_name}_{table}_{model}-{runid}_{ripf}_{config}_gn')
    return fname

def fname_cmip(var_name, table, model, runid, ripf, yri, yrf):
    fname = (
        f'{var_name}_{table}_{model}_{runid}_{ripf}_gn_{yri}01-{yrf}12.nc')
    return fname

def load_var_yearly(varname, var_dir, var_file, yri, yrf):
    """Load variable data from yearly files."""
    var_list = []
    for yr in range(yri, yrf+1):
        var_1yr = (xr.open_dataset(f'{var_dir}/{var_file}_{yr}01-{yr}12.nc')
            [varname])
        var_list.append(var_1yr)
    var_da = xr.concat(var_list, dim='time')
    return var_da
