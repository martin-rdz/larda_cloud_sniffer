
import ast
import datetime
import numpy as np


def ts_to_dt(ts):
    return datetime.datetime.utcfromtimestamp(ts)
def dt_to_ts(dt):
    return (dt - datetime.datetime(1970, 1, 1)).total_seconds()
    
toarray = lambda s: np.array(ast.literal_eval(s))
    

def standard(cloud):
    ffcloud = lambda s: float(cloud[s])

    dt_begin = ts_to_dt(ffcloud('Begin_Date_Unix'))
    dt_end = ts_to_dt(ffcloud('End_Date_Unix'))
    duration = dt_end-dt_begin

    conds = [ffcloud('Date_Year')>=2014, f"{ffcloud('Date_Year')}",
             duration > datetime.timedelta(seconds=20*60), f"{duration}",
             ffcloud('N_Profiles')*30/(duration.seconds) > 0.8, f"fill {ffcloud('N_Profiles')*30/duration.seconds:.3f}",
             cloud['Cloud_Run'] in ['layered'], f"{cloud['Cloud_Run']}",
             ffcloud('Cloud_Thickness_MED') < 350, f"thick {ffcloud('Cloud_Thickness_MED'):.1f}",
             ffcloud('CTH_STD') < 150, f"cth_std {ffcloud('CTH_STD'):.1f}",
             ffcloud('CTT') < 273.5, f"ctt {ffcloud('CTT'):.1f}",
             ffcloud('CTT') > 235.15, f"ctt {ffcloud('CTT'):.1f}",
            ]
            
    return conds
    
    
def with_gt_cbh(cloud, thres):
    ffcloud = lambda s: float(cloud[s])

    conds = standard(cloud)
    conds += [ffcloud('CBH') > thres, f"{ffcloud('CBH'):.1f}",]
    #conds += [ffcloud('CTH') > thres, f"{ffcloud('CTH'):.1f}",]
    
    return conds
    
#         conds = [duration > datetime.timedelta(seconds=20*60), f"{duration}",
#                  ffcloud('N_Profiles')/(duration.seconds/30.) > 0.8, f"{ffcloud('N_Profiles')*30/duration.seconds:.3f}",
#                  cloud['Cloud_Run'] in ['layered'], f"{cloud['Cloud_Run']}",
#                  ffcloud('Cloud_Thickness_MED') < 350, f"{ffcloud('Cloud_Thickness_MED'):.1f}",
#                  ffcloud('CTH_STD') < 150, f"{ffcloud('CTH_STD'):.1f}",
#                  ffcloud('CBH') > 2000, f"{ffcloud('CBH'):.1f}",
#                 ]

def with_lt_cbh(cloud, thres):
    ffcloud = lambda s: float(cloud[s])

    conds = standard(cloud)
    conds += [ffcloud('CBH') < thres, f"{ffcloud('CBH'):.1f}",]
    #conds += [ffcloud('CTH') < thres, f"{ffcloud('CTH'):.1f}",]
    
    return conds
    
    
def not_wave(cloud, autocorr_thres):
    ffcloud = lambda s: float(cloud[s])

    conds = standard(cloud)
    
    #if 'v_dl_autocor_time' in cloud and len(cloud['v_dl_autocor_time']) != '[nan]':
    if 'v_dl_autocor_time' in cloud and len(cloud['v_dl_autocor_time']) > 5:
        v_dl_autocor_time = cloud['v_dl_autocor_time']
        v_dl_autocor_coeff = cloud['v_dl_autocor_coeff']
        #print(v_dl_autocor_time[-10:])
        #print(v_dl_autocor_coeff[-10:])
        if not v_dl_autocor_time[-1] == ']':
            v_dl_autocor_time += ']'
        if not v_dl_autocor_coeff[-1] == ']':
            v_dl_autocor_coeff += ']'

        vel = max(float(cloud['VEL']),0.1)
        autocor_time = toarray(v_dl_autocor_time)
        autocor_coeff = toarray(v_dl_autocor_coeff)
        autocorr_lt_thres = np.where(autocor_coeff > 0.8)[0]
        i_above_thres = autocorr_lt_thres[-1] if len(autocorr_lt_thres) > 0 else 0
        autocorr_at_thres = autocor_time[i_above_thres]*vel if len(autocor_time) > 0 else 9999
        conds += [autocorr_at_thres < autocorr_thres, "autocorr_at_thres < autocorr_thres"]
    
    else:
        # doppler lidar obs not long enough
        conds += [False, "v_dl_autocor not in cloud/DL missing"]

    return conds
    
    
    
def conditions_ice_w_CTH(cloud):
    ffcloud = lambda s: float(cloud[s])

    conds = standard(cloud)
    conds += [ffcloud('CTH') > 2000, f"{ffcloud('CTH'):.1f}",
             ffcloud('IWC_TOP_N')/ffcloud('N_Profiles') > 0.05, f"{ffcloud('IWC_TOP_N')/ffcloud('N_Profiles'):.2f}",
            ]

    return conds


def conditions_ice_wo_CTH(cloud):
    ffcloud = lambda s: float(cloud[s])

    conds = standard(cloud)
    conds += [ffcloud('IWC_TOP_N')/ffcloud('N_Profiles') > 0.05, f"{ffcloud('IWC_TOP_N')/ffcloud('N_Profiles'):.2f}",
            ]

    return conds