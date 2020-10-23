import pandas as pd
import numpy as np
def convert_date_gefs(actual_time):
    init_time = pd.Timestamp(1800,1,1)
    diff = (actual_time-init_time)/ np.timedelta64(1, 'h') 
    return diff

