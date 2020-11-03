import pandas as pd
import numpy as np
def convert_date_gefs(actual_time):
    init_time = pd.Timestamp(1800,1,1)
    diff = (actual_time-init_time)/ np.timedelta64(1, 'h') 
    return diff

def convert_date_era(actual_time):
    init_time = pd.Timestamp(1979,1,1)
    diff = (actual_time-init_time)/ np.timedelta64(1, 'h') 
    return diff
    
def convert_date_gefs_r(hours):
    init_time = pd.Timestamp(1800,1,1)
    if (np.shape(hours)!=()):
        actual_time = [init_time + pd.Timedelta(hour,unit='h') for hour in hours]
    else:
        actual_time = init_time + pd.Timedelta(hours,unit='h')
    return actual_time

def convert_date_era_r(hours):
    init_time = pd.Timestamp(1979,1,1)
    if np.shape(hours)!=():
        actual_time = [init_time + pd.Timedelta(hour,unit='h') for hour in hours]
    else:
        actual_time = init_time + pd.Timedelta(hours,unit='h')
    return actual_time

def weighted_average_area(Dataset):
	lons,lats = np.meshgrid(Dataset.lon,Dataset.lat)
	var_name = list(Dataset.keys())[0]
	numerator = np.sum(Dataset[var_name]*np.cos(np.deg2rad(lats)))
	denominator = np.sum(np.cos(np.deg2rad(lats)))
	weighted_area = numerator/denominator
	return weighted_area

