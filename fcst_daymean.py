import xarray as xr 
from glob import glob
import numpy as np 
import sys

path = '/media/onno/Algemeen/Thesis/'

def daymean(file):
	print(file)
	ds = xr.open_dataset(path+file,decode_times=False)
	ds['lead'] = ds.lead.values//24
	latz = [(90,45),(45,0),(0,-45),(-45,-90)]
	ds_list = []
	for lat in latz:
		print(lat)
		ds_lat = ds.sel(lat=slice(lat[0],lat[1]))
		ds_lat = ds_lat.groupby('lead').mean()
		ds_list.append(ds_lat)
	ds_new = xr.concat(ds_list,dim='lat')
	ds_new.to_netcdf(path+file[:-3]+'_daymean.nc')


filez_all = ['era5rf_env_wledit2000-10000_latavg_v300_0-240h_12hourly_2x2nh_jan79-dec19.nc','era5rf_phasevel_wledit2000-10000_latavg_v300_envgt15_0-120h_6hourly_2x2nh_jan79-dec19_setvrange_-100to100.nc',
			'gefsrf2_env_wledit2000-10000_latavg_v300_control0-252h_6hourly_2x2_dec84-nov19.nc','gefsrf2_phasevel_wledit2000-10000_latavg_v300_envgt15_control0-252h_6hourly_2x2_dec84-nov19_setvrange_-100to100.nc',
			'era5rf_t850_0-240h_12hourly_2x2nh_jan79-dec19.nc','gefsrf2_t850_control0-252h_6hourly_2x2_dec84-nov19.nc']


for file in filez_all:
	daymean(file)
	sys.exit()

