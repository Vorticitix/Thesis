import pandas as pd
import numpy as np

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

