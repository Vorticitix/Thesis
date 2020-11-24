import pandas as pd
import numpy as np

plot_dic = {
    'envelope':{'data_variable':'v','title':'Envelope','label':'E (m/s)','filename':'envelope'},
    'phasespeed':{'data_variable':'v','title':'Phase Speed','label':'Cp (m/s)','filename':'phasespeed'},
    'u_wind':{'data_variable':'u','title':'Zonal Wind','label':'u (m/s)','filename':'u_wind'},
    'v_wind':{'data_variable':'v','title':'Meridional Wind','label':'v (m/s)','filename':'v_wind'},
    'groupspeed':{'data_variable':'v','title':'Group Speed','label':'Cg (m/s)','filename':'groupspeed'}
}

file_dic = {
    'envelope':{'ERA5':'era51_mars_env_wledit2000-10000_latavg_v300_79-19_6hourly_smoothed.nc',
               'GFS':'gefsrf2_env_wledit2000-10000_latavg_v300_control0-252h_6hourly_2x2_dec84-nov19.nc',
               'CFS':'cfsr-cfsv2_env_wledit2000-10000_latavg_v300_79-19_6hourly_anom_from_smoothed04_clim_smoothed.nc'},
    'phasespeed':{'ERA5':'era51_mars_phasevel_wledit2000-10000_latavg_v300_envgt15_79-19_6hourly_setvrange_-100to100.nc',
               'GFS':'gefsrf2_phasevel_wledit2000-10000_latavg_v300_envgt15_control0-252h_6hourly_2x2_dec84-nov19_setvrange_-100to100.nc',
               'CFS':'cfsr-cfsv2_phasevel_wledit2000-10000_latavg_v300_envgt15_79-19_6hourly_anom_from_smoothed04_clim_setvrange_-100to100.nc'},
    'u_wind':{'ERA5':'era5_u300_79-19_6hourly.nc',
               'GFS':'gefsrf2_u300_control0-252h_6hourly_2x2_dec84-nov19.nc'}
}
#The function below is used to convert real datetimes to hours sicne 1 Jan 1800 because GEFS data are in these units
def convert_date_gefs(actual_time):
    #set timestamp at 1 Jan 1800, because GEFS data is in hours since this data
    init_time = pd.Timestamp(1800,1,1)
    #Calculate hours since 1 Jan 1800
    diff = (actual_time-init_time)/ np.timedelta64(1, 'h') 
    return diff

#The function below is used to convert real datetimes to hours sicne 1 Jan 1979 because ERA5 data are in these units
def convert_date_era(actual_time):
    #set timestamp at 1 Jan 1979, because GEFS data is in hours since this data
    init_time = pd.Timestamp(1979,1,1)
    #Calculate hours since 1 Jan 1979
    diff = (actual_time-init_time)/ np.timedelta64(1, 'h') 
    return diff

#The function below converts hours since 1 Jan 1800 to real datetime values
def convert_date_gefs_r(hours):
    init_time = pd.Timestamp(1800,1,1)
    #Differentiate between 1 value and 1D array
    if (np.shape(hours)!=()):
        actual_time = [init_time + pd.Timedelta(hour,unit='h') for hour in hours]
    else:
        actual_time = init_time + pd.Timedelta(hours,unit='h')
    return actual_time

#The function below converts hours since 1 Jan 11979 to real datetime values
def convert_date_era_r(hours):
    init_time = pd.Timestamp(1979,1,1)
    #Differentiate between 1 value and 1D array
    if np.shape(hours)!=():
        actual_time = [init_time + pd.Timedelta(hour,unit='h') for hour in hours]
    else:
        actual_time = init_time + pd.Timedelta(hours,unit='h')
    return actual_time

#The function below receives a 2D xarray dataset and returns one value for latitude weighted mean
def weighted_average_area(Dataset):
    #create 2D grid with latitude weights. 
	lons,lats = np.meshgrid(Dataset.lon,Dataset.lat)
    #list datavariable name
	var_name = list(Dataset.keys())[0]
    #Calculate latitude weighted mean
	numerator = np.sum(Dataset[var_name]*np.cos(np.deg2rad(lats)))
	denominator = np.sum(np.cos(np.deg2rad(lats)))
	weighted_area = numerator/denominator
	return weighted_area

