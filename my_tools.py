import pandas as pd
import numpy as np
import xarray as xr

plot_dic = {
    'envelope':{'data_variable':'v','title':'Envelope','label':'E (m/s)','filename':'envelope'},
    'phasespeed':{'data_variable':'v','title':'Phase Speed','label':'Cp (m/s)','filename':'phasespeed'},
    'u_wind':{'data_variable':'u','title':'Zonal Wind','label':'u (m/s)','filename':'u_wind'},
    'v_wind':{'data_variable':'v','title':'Meridional Wind','label':'v (m/s)','filename':'v_wind'},
    'groupspeed':{'data_variable':'v','title':'Group Speed','label':'Cg (m/s)','filename':'groupspeed'},
    'wavenumber':{'data_variable':'v','title':'Wave Number','label':'k (rad/m)','filename':'wavenumber'},
    'T850':{'data_variable':'t','title':'850 hPa Temperature','label':'T (°C)','filename':'T850'}
}

file_dic = {
    'envelope':{'ERA5':'era51_mars_env_wledit2000-10000_latavg_v300_79-19_6hourly_smoothed.nc',
               'GFS':'gefsrf2_env_wledit2000-10000_latavg_v300_control0-252h_6hourly_2x2_dec84-nov19.nc',
               'CFS':'cfsr-cfsv2_env_wledit2000-10000_latavg_v300_79-19_6hourly_anom_from_smoothed04_clim_smoothed.nc',
               'ERA5RF':'era5rf_env_wledit2000-10000_latavg_v300_0-240h_12hourly_2x2nh_jan79-dec19.nc'},
    'phasespeed':{'ERA5':'era51_mars_phasevel_wledit2000-10000_latavg_v300_envgt15_79-19_6hourly_setvrange_-100to100.nc',
               'GFS':'gefsrf2_phasevel_wledit2000-10000_latavg_v300_envgt15_control0-252h_6hourly_2x2_dec84-nov19_setvrange_-100to100.nc',
               'CFS':'cfsr-cfsv2_phasevel_wledit2000-10000_latavg_v300_envgt15_79-19_6hourly_anom_from_smoothed04_clim_setvrange_-100to100.nc',
               'ERA5RF':'era5rf_phasevel_wledit2000-10000_latavg_v300_envgt15_0-120h_6hourly_2x2nh_jan79-dec19_setvrange_-100to100.nc'},
    'u_wind':{'ERA5':'era51_mars_u300_79-19_6hourly.nc',
               'GFS':'gefsrf2_u300_control0-252h_6hourly_2x2_dec84-nov19.nc',
             'ERA5RF':'era5rf_u300_0-240h_12hourly_2x2nh_jan79-dec19.nc'},
    'v_wind':{'ERA5':'era51_mars_v300_79-19_6hourly.nc',
               'GFS':'gefsrf2_v300_control0-252h_6hourly_2x2_dec84-nov19.nc',
             'ERA5RF':'era5rf_v300_0-240h_12hourly_2x2nh_jan79-dec19.nc'},
    'wavenumber':{'ERA5':'era51_mars_k_wledit2000-10000_latavg_v300_envgt15_79-19_6hourly_setvrange_0to1.nc',
                  'GFS':'gefsrf2_k_wledit2000-10000_latavg_v300_envgt15_control0-252h_6hourly_2x2_dec84-nov19_setvrange_0to1.nc'},
    'T850':{'ERA5':'era51_mars_t850_79-19_6hourly.nc',
           'GFS':'gefsrf2_t850_control0-252h_6hourly_2x2_dec84-nov19.nc',
           'ERA5RF':'era5rf_t850_0-240h_12hourly_2x2nh_jan79-dec19.nc'},
    'T850_grid':{'ERA5':'era51_fldmean_mars_t850_79-19_24hourly_lon_{}_{}_lat_{}_{}.nc',
                'GFS':'gefsrf2_fldmean_t850_control0-252h_24hourly_lon_{}_{}_lat_{}_{}.nc',
                'ERA5RF':'era5rf_fldmean_t850_0-240h_24hourly_2x2nh_dec84-nov19_lon_{}_{}_lat_{}_{}.nc'},
    'Z500':{'ERA5':'era51_mars_phi500_79-19_6hourly.nc'},
    'MSLP':{'ERA5':'era5_mslp_79-19_6hourly_remapbil2x2.nc'}
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

#The function below converts hours since 1 Jan 1800 to real datetime values
def convert_date_era5rf_r(hours):
    init_time = pd.Timestamp(1900,1,1)
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
	if Dataset[var_name].count()/len(Dataset[var_name])**2<0.2:
		weighted_area = np.nan
	return weighted_area

def weighted_average_area_3D(Dataset,variable,multiplier=0.2):
    #create 2D grid with latitude weights. 
	lons,lats = np.meshgrid(Dataset.lon,Dataset.lat)
	
	var_name = list(Dataset.keys())[0]	
	#make 3D array of lats
	lats_3D = np.repeat(lats[np.newaxis,:,:],len(Dataset.time),axis=0)
	#For wave speed values we only want to get an average if we have 25% of values that are not nan
	if variable=='phasespeed':
		#Count number of not_nans in dataset for each timestep
		counts_non_nan = np.count_nonzero(np.invert(np.isnan(Dataset[var_name].values)),axis=(1,2))
		boolean = counts_non_nan >= (len(lons)*len(lats))*multiplier
		#Calculate latitude weighted mean only for timesteps where more than 25% of values are defined
		numerator = np.nansum(Dataset[var_name][:,:,:]*np.cos(np.deg2rad(lats_3D[:,:,:])),axis=(1,2))
		denominator = np.nansum(np.cos(np.deg2rad(lats_3D[:,:,:])),axis=(1,2))
		weighted_area_np = numerator/denominator
		weighted_area_np[np.invert(boolean)]=np.nan
		weighted_area = xr.DataArray(data=weighted_area_np,coords=dict(
			time=Dataset.time.values),dims='time')

	else:
		numerator = np.sum(Dataset[var_name]*np.cos(np.deg2rad(lats_3D)),axis=(1,2))
		denominator = np.sum(np.cos(np.deg2rad(lats_3D)),axis=(1,2))		
		weighted_area = numerator/denominator
	return weighted_area

def weighted_average_area_3D_lead(Dataset,variable,multiplier=0.2):
    #create 2D grid with latitude weights. 
	lons,lats = np.meshgrid(Dataset.lon,Dataset.lat)
	
	var_name = list(Dataset.keys())[0]	
	#make 3D array of lats
	lats_3D = np.repeat(lats[np.newaxis,:,:],len(Dataset.lead),axis=0)
	#For wave speed values we only want to get an average if we have 25% of values that are not nan
	if variable=='phasespeed':
		#Count number of not_nans in dataset for each timestep
		counts_non_nan = np.count_nonzero(np.invert(np.isnan(Dataset[var_name].values)),axis=(1,2))
		boolean = counts_non_nan >= (len(lons)*len(lats))*multiplier
		#Calculate latitude weighted mean only for timesteps where more than 25% of values are defined
		numerator = np.nansum(Dataset[var_name][:,:,:]*np.cos(np.deg2rad(lats_3D[:,:,:])),axis=(1,2))
		denominator = np.nansum(np.cos(np.deg2rad(lats_3D[:,:,:])),axis=(1,2))
		weighted_area_np = numerator/denominator
		weighted_area_np[np.invert(boolean)]=np.nan
		weighted_area = xr.DataArray(data=weighted_area_np,coords=dict(
			lead=Dataset.lead.values),dims='lead')

	else:
		numerator = np.sum(Dataset[var_name]*np.cos(np.deg2rad(lats_3D)),axis=(1,2))
		denominator = np.sum(np.cos(np.deg2rad(lats_3D)),axis=(1,2))		
		weighted_area = numerator/denominator
	return weighted_area
    

def detect_heatwaves(ds):
    anom = ds.t
    #create boolean array and return true if anomaly is positve (higher than 90th percentile)
    anom_bool = anom>0
    #convert boolean array to 0 and 1
    iszero = np.concatenate(([0], np.equal(anom_bool.values, True).view(np.int8), [0]))
    #create array where where start and end of anomaly is indicated by a 1 instead of 0
    abs_diff = np.abs(np.diff(iszero))
    idxz = np.where(abs_diff == 1)[0].reshape(-1, 2)
    #Persistent extremes
    idx_bool_persistent = (idxz[:,1]-idxz[:,0])>=4
    pot_hw_idxz_persistent = idxz[idx_bool_persistent,:]
    #subtract 1 day to set end date at last day of temp extreme
    persistent_datez = np.vstack((ds.time.values[pot_hw_idxz_persistent[:,0]],ds.time.values[pot_hw_idxz_persistent[:,1]]-pd.Timedelta(1,'days'))).T
    #short lived extremes
    idx_bool_short = (idxz[:,1]-idxz[:,0])<=2
    pot_hw_idxz_short = idxz[idx_bool_short,:]
    #subtract 1 day to set end date at last day of temp extreme
    short_datez = np.vstack((ds.time.values[pot_hw_idxz_short[:,0]],ds.time.values[pot_hw_idxz_short[:,1]]-pd.Timedelta(1,'days'))).T
    return persistent_datez, short_datez

def detect_coldwaves(ds):
    anom = ds.t
    #create boolean array and return true if anomaly is negative (lower than 10th percentile)
    anom_bool = anom<0
    #convert boolean array to 0 and 1
    iszero = np.concatenate(([0], np.equal(anom_bool.values, True).view(np.int8), [0]))
    #create array where where start and end of anomaly is indicated by a 1 instead of 0
    abs_diff = np.abs(np.diff(iszero))
    idxz = np.where(abs_diff == 1)[0].reshape(-1, 2)
    #Differentiate between persistent and short lived extremes
    idx_bool_persistent = (idxz[:,1]-idxz[:,0])>=4
    pot_hw_idxz_persistent = idxz[idx_bool_persistent,:]
    persistent_datez = np.vstack((ds.time.values[pot_hw_idxz_persistent[:,0]],ds.time.values[pot_hw_idxz_persistent[:,1]]-pd.Timedelta(1,'days'))).T
    idx_bool_short = (idxz[:,1]-idxz[:,0])<=2
    pot_hw_idxz_short = idxz[idx_bool_short,:]
    short_datez = np.vstack((ds.time.values[pot_hw_idxz_short[:,0]],ds.time.values[pot_hw_idxz_short[:,1]]-pd.Timedelta(1,'days'))).T
    return persistent_datez, short_datez

def gilbert_skill_score(hits,miss,false_alarm,non_event):
    hitsc = (hits + miss)*(hits+false_alarm)/(hits+miss+non_event+false_alarm)
    GSS = (hits-hitsc)/(hits+miss+false_alarm-hitsc)
    return GSS

def hansen_kuiper_discrimant(hits,miss,false_alarm,non_event):
	HK = hits/(hits+miss)-false_alarm/(false_alarm+non_event)
	return HK

def gilbert_skill_score_dataframe(df):
    hits = df['hits']
    miss = df['miss']
    false_alarm = df['false_alarm']
    non_event = df['non_event']
    hitsc = (hits + miss)*(hits+false_alarm)/(hits+miss+non_event+false_alarm)
    GSS = (hits-hitsc)/(hits+miss+false_alarm-hitsc)
    return(GSS)
    

#the function below yields the same results as ydaysub from the cdo framework
def calculate_anomalies(multi_year,climatology):
    #make a copy of forecast/reanalysis and climatology timeseries to alternate there time data
    multi_year_alt = multi_year.load().copy()
    clim_new_time = climatology.load().copy()
    #convert time data of climatology to dayofyear (1 = 1 Jan, 366 = 31 December)
    clim_new_time['time'] = clim_new_time.load().time.dt.dayofyear
    #create empty list to add every individual year
    multi_year_total = []
    #loop over all years in data and select individual years from data
    for year in range(multi_year.time.dt.year.values.min(),multi_year.time.dt.year.values.max()+1):
        #select individual years from data
        year_sub = multi_year.sel(time=slice(pd.Timestamp('{}-01-01'.format(year)),
                                pd.Timestamp('{}-12-31'.format(year))))
        #copy subsetted data to alternate time data
        year_sub_new_time = year_sub.copy()
        #Convert time data to from timestamp to dayofyear
        if year%4!=0:
            # if year is not a leap year the dayofyear numbers do not exactly match the climatology because Feb 29th is not present
            # To account for this we skip th 59th day (which would be Feb 29 in a leap year) and add 1 day to day number thereafter
            year_sub_new_time['time'] = [i+1 if i>59 else i for i in year_sub_new_time.time.dt.dayofyear.values]
        else:
            #During a leap year the method is straightforward to convert to daynumber
            year_sub_new_time['time'] = year_sub_new_time.time.dt.dayofyear
        #calculate the anomalies
        #every single day from the year subset is substracted with the value in the climatology with the same day number
        #So in non leap years Feb 29th in the climatology data is skipped
        anom_sub = year_sub_new_time - clim_new_time
        #replace absolute temperatures with anomalies in the year subset data. Now we also have normal dates as time data again
        year_sub.t.values = anom_sub.t.values
        #add year subset to list
        multi_year_total.append(year_sub)
    #concatenate all objects in list along datetime axis
    multi_year_anom = xr.concat(multi_year_total,dim='time')
    return multi_year_anom

#For a given boolean array with anomalies and non anomalous values detect length of longest consecutive True values
#I assume here that maximum length of consecutive true values resembles extreme duration
def max_len_boolean_array(array):
    #create boolean array where True inidcates non-zero values
    iszero = np.concatenate(([0], np.equal(array, True).view(np.int8), [0]))
    #indicate at which indices the value changes from True to False or vice versa. Convert to 2D array to se longest consecutive series of True values
    abs_diff = np.abs(np.diff(iszero))
    idx_2D = np.where(abs_diff == 1)[0].reshape(-1, 2)
    diff = idx_2D[:,1]-idx_2D[:,0]
    #caluclate longest consecutive series of true Values. Add one because the last day counts as well
    max_length = np.max(diff)+1
    return max_length
    