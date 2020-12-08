import pandas as pd
import numpy as np
import xarray as xr

plot_dic = {
    'envelope':{'data_variable':'v','title':'Envelope','label':'E (m/s)','filename':'envelope'},
    'phasespeed':{'data_variable':'v','title':'Phase Speed','label':'Cp (m/s)','filename':'phasespeed'},
    'u_wind':{'data_variable':'u','title':'Zonal Wind','label':'u (m/s)','filename':'u_wind'},
    'v_wind':{'data_variable':'v','title':'Meridional Wind','label':'v (m/s)','filename':'v_wind'},
    'groupspeed':{'data_variable':'v','title':'Group Speed','label':'Cg (m/s)','filename':'groupspeed'},
    'wavenumber':{'data_variable':'v','title':'Wave Number','label':'k (rad/m)','filename':'wavenumber'}
}

file_dic = {
    'envelope':{'ERA5':'era51_mars_env_wledit2000-10000_latavg_v300_79-19_6hourly_smoothed.nc',
               'GFS':'gefsrf2_env_wledit2000-10000_latavg_v300_control0-252h_6hourly_2x2_dec84-nov19.nc',
               'CFS':'cfsr-cfsv2_env_wledit2000-10000_latavg_v300_79-19_6hourly_anom_from_smoothed04_clim_smoothed.nc'},
    'phasespeed':{'ERA5':'era51_mars_phasevel_wledit2000-10000_latavg_v300_envgt15_79-19_6hourly_setvrange_-100to100.nc',
               'GFS':'gefsrf2_phasevel_wledit2000-10000_latavg_v300_envgt15_control0-252h_6hourly_2x2_dec84-nov19_setvrange_-100to100.nc',
               'CFS':'cfsr-cfsv2_phasevel_wledit2000-10000_latavg_v300_envgt15_79-19_6hourly_anom_from_smoothed04_clim_setvrange_-100to100.nc'},
    'u_wind':{'ERA5':'era5_u300_79-19_6hourly.nc',
               'GFS':'gefsrf2_u300_control0-252h_6hourly_2x2_dec84-nov19.nc'},
    'wavenumber':{'ERA5':'era51_mars_k_wledit2000-10000_latavg_v300_envgt15_79-19_6hourly_setvrange_0to1.nc',
                  'GFS':'gefsrf2_k_wledit2000-10000_latavg_v300_envgt15_control0-252h_6hourly_2x2_dec84-nov19_setvrange_0to1.nc'}
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
def weighted_average_area_2D(Dataset):
    #create 2D grid with latitude weights. 
	lons,lats = np.meshgrid(Dataset.lon,Dataset.lat)
    #list datavariable name
	var_name = list(Dataset.keys())[0]
    #Calculate latitude weighted mean
	numerator = np.sum(Dataset[var_name]*np.cos(np.deg2rad(lats)))
	denominator = np.sum(np.cos(np.deg2rad(lats)))
	weighted_area = numerator/denominator
	return weighted_area

def weighted_average_area_3D(Dataset,variable):
    #create 2D grid with latitude weights. 
	lons,lats = np.meshgrid(Dataset.lon,Dataset.lat)
	
	var_name = list(Dataset.keys())[0]	
	#make 3D array of lats
	lats_3D = np.repeat(lats[np.newaxis,:,:],len(Dataset.time),axis=0)
	#For wave speed values we only want to get an average if we have 25% of values that are not nan
	if variable=='phasespeed':
		#Count number of not_nans in dataset for each timestep
		counts_non_nan = np.count_nonzero(np.invert(np.isnan(Dataset[var_name].values)),axis=(1,2))
		boolean = counts_non_nan > (len(lons)*len(lats))*0.2
		#Calculate latitude weighted mean only for timesteps where more than 25% of values are defined
		numerator = np.nansum(Dataset[var_name][boolean,:,:]*np.cos(np.deg2rad(lats_3D[boolean,:,:])),axis=(1,2))
		denominator = np.nansum(np.cos(np.deg2rad(lats_3D[boolean,:,:])),axis=(1,2))
		weighted_area_np = numerator/denominator
		weighted_area = xr.Dataset(data_vars=dict(
			v=(['time'], weighted_area_np)),coords=dict(
			time=Dataset.time.values[boolean]))

	else:
		numerator = np.sum(Dataset[var_name]*np.cos(np.deg2rad(lats_3D)),axis=(1,2))
		denominator = np.sum(np.cos(np.deg2rad(lats_3D)),axis=(1,2))		
		weighted_area = numerator/denominator
	return weighted_area
    
    