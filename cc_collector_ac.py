#!/usr/bin/python

import argparse
import datetime
import numpy as np
import pickle
import os
import CLS_Clouds as CLS
import scipy.stats

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--date', type=str, help='date on which the algorithm should be run')
parser.add_argument('--campaign', type=str, help='Cloudnet station on which the algorithm should be run')
parser.add_argument('--type', metavar='type_str', type=str, help='Cloud type on which the algorithm should be run')
args = parser.parse_args()

campaign=args.campaign
date=datetime.datetime.strptime(args.date, "%Y%m%d").date()
cloud_type=args.type
print("starting collector {} {} {}".format(args.campaign, args.date, args.type))

def load_object(filename):
    ifile=open(filename,"rb")
    obj = pickle.load(ifile)
    ifile.close()
    return obj

def dict_to_str(source):
    outstr=""
    keys=sorted(source.keys())
    for key in keys:
        outstr += str(source[key])+";"
    return outstr


def keys_to_str(source):    
    outstr=""    
    keys=sorted(source.keys())    
    for key in keys:        
        outstr += str(key)+";"    
    return outstr

if cloud_type==None:
    cloud_type='all'

recognized_types=['pure_ice','pure_liquid','mixed-phase','liquid-based','tower','all']
assert cloud_type in recognized_types, "{} not in {}".format(cloud_type, recognized_types)

infile='../cloud_properties/cloud_properties_'+campaign+'/'+date.strftime("%Y_%m_%d")+'_clouds.dat'
outfile='../cloud_collections/'+'cloud_collection_'+campaign+'_'+cloud_type+'.csv'

if not(os.path.exists(outfile)):
    newfile=True
else:
    if os.path.getsize(outfile)==0:
        newfile=True
    else:
        newfile=False

clouds = load_object(infile)

ice_only=[4]
liquid=[1,3,5,7]
droplets_only=[1,5]
sep=6
output={}

for i in range(len(clouds)):

    #print clouds[i].cloud_type
    #print clouds[i].top_variation(),clouds[i].time_length(),clouds[i].fill_factor()

    if (clouds[i].cloud_type==cloud_type or cloud_type=='all') and clouds[i].n_profiles()>0:

       if clouds[i].top_variation()<600.0 and clouds[i].time_length()>900 and clouds[i].fill_factor()>0.7:
            
            begin=datetime.datetime.utcfromtimestamp(clouds[i].begin_time)
            end=datetime.datetime.utcfromtimestamp(clouds[i].end_time)
            middle=begin+datetime.timedelta(seconds=clouds[i].end_time-clouds[i].begin_time)

            clouds[i].validate_features()

            output["A_Unique_Identifier"]=begin.strftime("%Y%m%d%H%M%S")+str(int(clouds[i].tops[0]))
            output["Cloud_Run"]=clouds[i].cloud_type
            output["Cloud_Type"]=clouds[i].most_common_type()
            output["Date_Year"]=middle.year
            output["Date_Month"]=middle.month
            output["Date_Day"]=middle.day
            output["Date_Hour"]=middle.hour
            output["Begin_Date"]=begin.strftime("%Y_%m_%d_%H_%M_%S")
            output["Begin_Date_Unix"]=clouds[i].begin_time
            output["End_Date"]=end.strftime("%Y_%m_%d_%H_%M_%S")
            output["End_Date_Unix"]=clouds[i].end_time

            output["CTH"]=np.max(clouds[i].tops)
            output["CTH_AVG"]=np.average(clouds[i].tops)
            output["CTH_STD"]=clouds[i].top_variation()
            output["CTH_DIFFSUM"]=clouds[i].top_diffsum()
            output["LLH_STD"]=clouds[i].liquid_layer_variation()
            output["CBH"]=np.min(clouds[i].bases)
            output["CBT"],output["CTT_MED"],output["CTT"]=clouds[i].temperature_range()
            output["CBP"],output["CTP"]=clouds[i].pressure_range()

            output["VEL"],output["DIR"]=clouds[i].vel_and_dir()

            output["N_Profiles"]=clouds[i].n_profiles()
            output["N_Liquid"]=clouds[i].n_profiles(["pure liquid","liquid-based layer"])
            output["N_MixedPhase"]=clouds[i].n_profiles(["mixed-phase layer"])

            output["LWC_AVG"],output["LWC_MED"],output["LWC_STD"],output["LWC_N"]=clouds[i].average("LWC",CLS.droplets_only)
            output["PATH_LWP_AVG"],output["PATH_LWP_S_AVG"],output["PATH_IWP_AVG"],output["PATH_IWP_STD"]=clouds[i].average_paths()
            output["Cloud_Thickness_AVG"],output["Cloud_Thickness_MED"],output["Cloud_Thickness_STD"]=clouds[i].cloud_top_thickness()

            output["IWC_AVG"],output["IWC_MED"],output["IWC_STD"],output["IWC_N"]=clouds[i].average("IWC",ice_only)
            output["IWC_TOP_AVG"],output["IWC_TOP_MED"],output["IWC_TOP_STD"],output["IWC_TOP_N"]=clouds[i].separation_average("IWC",sep)
            output["ILCR_AVG"],output["ILCR_MED"],output["ILCR_N"],output["ILCR_values"]=clouds[i].ilr(sep)

            clouds[i].correct_LDR(-20,-33)
            output["LDR_AVG"],output["LDR_MED"],output["LDR_STD"],output["LDR_N"]=clouds[i].average("LDR",ice_only)
            output["LDR_values"]=list(np.histogram(clouds[i].return_values("LDR",ice_only),70,(-35.0,1))[0])
            output["LDR_TOP_AVG"],output["LDR_TOP_MED"],output["LDR_TOP_STD"],output["LDR_TOP_N"]=clouds[i].separation_average("LDR",sep)

            output["N_Melting"]=clouds[i].n_melting()
            output["N_Drizzle"]=clouds[i].n_drizzle()

            output["v_AVG"],output["v_MED"],output["v_STD"],output["v_N"]=clouds[i].average("v",ice_only)
            output["v_TOP_AVG"],output["v_TOP_MED"],output["v_TOP_STD"],output["v_TOP_N"]=clouds[i].separation_average("v",sep)
            output["v_values"]=list(np.histogram(clouds[i].return_values("v",ice_only),60,(-1.5,0.0))[0])

            #output["tfv_AVG"],output["tfv_MED"],output["tfv_STD"],output["tfv_N"]=clouds[i].average("tfv",ice_only)
            #output["tfv_TOP_AVG"],output["tfv_TOP_MED"],output["tfv_TOP_STD"],output["tfv_TOP_N"]=clouds[i].separation_average("tfv",sep)
            #output["vair_AVG"],output["vair_MED"],output["vair_STD"],output["vair_N"]=clouds[i].combined_average("vair","v",ice_only)

            #output["v_comb_tfv_AVG"],output["v_comb_tfv_MED"],output["v_comb_tfv_STD"],output["v_comb_tfv_N"]=clouds[i].combined_average("v","tfv",ice_only)

            output["Z_AVG"],output["Z_MED"],output["Z_STD"],output["Z_N"]=clouds[i].average("Z",ice_only)
            output["Z_AVG_drop"],output["Z_MED_drop"],output["Z_STD_drop"],output["Z_N_drop"]=clouds[i].average("Z",droplets_only)
            output["Z_TOP_AVG"],output["Z_TOP_MED"],output["Z_TOP_STD"],output["Z_TOP_N"]=clouds[i].separation_average("Z",sep)
            z_top_vals = clouds[i].return_values_separation("Z",sep)
            output["Z_TOP_values"]=list(np.histogram(10*np.log10(z_top_vals),90,(-70.0,20))[0])
            z_vals = clouds[i].return_values("Z",ice_only)
            #print(z_vals)
            #print(10*np.log10(z_vals))
            output["Z_values"]=list(np.histogram(10*np.log10(z_vals),90,(-70.0,20))[0])

            output["ZE_TOP_AVG"],output["ZE_TOP_MED"],output["ZE_TOP_STD"],output["ZE_TOP_N"]=clouds[i].separation_average("ratio_z_e",sep)
            #print('ZE_TOP_stat ', output["ZE_TOP_AVG"],output["ZE_TOP_MED"],output["ZE_TOP_STD"],output["ZE_TOP_N"])
            output["ZE_AVG"],output["ZE_MED"],output["ZE_STD"],output["ZE_N"]=clouds[i].average("ratio_z_e",ice_only)
            #print('ZE_stat ', output["ZE_AVG"],output["ZE_MED"],output["ZE_STD"],output["ZE_N"])

            output["SNR_TOP_AVG"],output["SNR_TOP_MED"],output["SNR_TOP_STD"],output["SNR_TOP_N"]=clouds[i].separation_average("SNR",sep)

            output["width_AVG"],output["width_MED"],output["width_STD"],output["width_N"]=clouds[i].average("width",ice_only)

            #output["SNR_AVG"],output["SNR_MED"],output["SNR_STD"],output["SNR_N"]=clouds[i].average("SNR",ice_only)
            #output["SNR_10pp"],output["SNR_90pp"]=clouds[i].pp90("SNR",ice_only)

            output["alpha_Hogan_AVG"],output["alpha_Hogan_MED"],output["alpha_Hogan_STD"],output["alpha_Hogan_N"]=clouds[i].average("alpha_hogan",ice_only)
            output["alpha_Hogan_TOP_AVG"],output["alpha_Hogan_TOP_MED"],output["alpha_Hogan_TOP_STD"],output["alpha_Hogan_TOP_N"]=clouds[i].separation_average("alpha_hogan",sep)
            output["beta_AVG"],output["beta_MED"],output["beta_STD"],output["beta_N"]=clouds[i].average("beta",ice_only)
            output["delta_AVG"],output["delta_MED"],output["delta_STD"],output["delta_N"]=clouds[i].average("delta",ice_only)
            output["voldepol_AVG"],output["voldepol_MED"],output["voldepol_STD"],output["voldepol_N"]=clouds[i].average("voldepol",ice_only)
            
            #output["v_lidar_AVG"],output["v_lidar_STD"],output["v_lidar_N"],output["v_lgt0_AVG"],output["v_lgt0_STD"],output["v_lgt0_N"],vv_values=clouds[i].velocities()
            #output["v_lidar_histo"]=list(np.histogram(vv_values,60,(-3.0,3.0))[0])

            # refactored vertical velocities in liquid layers
            output['v_dl_mean'], output['v_dl_std'], output['v_dl_n'], v_base, _ = clouds[i].velocities()
            if len(v_base) > 0:
                output['v_dl_perc'] = np.percentile(v_base, [10,25,50,75,90]).tolist()
                output['v_dl_skew'] = scipy.stats.skew(v_base)
            else:
                output['v_dl_perc'] = [-99,-99,-99,-99,-99] 
                output['v_dl_skew'] = -99

            output['v_cr_mean'], output['v_cr_std'], output['v_cr_n'], v_base, _ = clouds[i].velocities_liquid_radar('whole')
            if len(v_base) > 0:
                output['v_cr_perc'] = np.percentile(v_base, [10,25,50,75,90]).tolist()
                output['v_cr_skew'] = scipy.stats.skew(v_base)
            else:
                output['v_cr_perc'] = [-99,-99,-99,-99,-99] 
                output['v_cr_skew'] = -99

            #manually corrected LDR
            output["LDRcorr_TOP_AVG"], output["LDRcorr_TOP_MED"], output["LDRcorr_TOP_STD"], output["LDRcorr_TOP_N"] = clouds[i].separation_average("LDRcorr", sep)
            output["LDRcorr_AVG"], output["LDRcorr_MED"], output["LDRcorr_STD"], output["LDRcorr_N"] = clouds[i].average("LDRcorr", ice_only)

            #interpolation
            #il=16
            #vv_int=[]
            #if len(vv_values)>il+1:
            #    for n in range(len(vv_values)/il-il-1):
            #        vv_int.append(np.average(vv_values[il*n:il*n+il]))
            #vv_int=np.array(vv_int)
            #vv_values=vv_int

            #output["v_lidar_histo_INT"]=list(np.histogram(vv_values,60,(-3.0,3.0))[0])

            output["v_radar_AVG"],output["v_radar_STD"],output["v_radar_N"],output["v_radar_lgt0_AVG"],output["v_radar_lgt0_STD"],output["v_radar_lgt0_N"],output["v_radar_WIDTH"],vv_radar_values,output["Z_top"]=clouds[i].velocities_radar()
            output["v_radar_histo"]=list(np.histogram(vv_radar_values,60,(-3.0,3.0))[0])

            output["file_hist_class"] = clouds[i].get_class_file_history()

            print(clouds[i].cloud_type,output["Z_top"])

            if int(output["A_Unique_Identifier"])==201207172309455452:
                print(output["A_Unique_Identifier"])
                for n in range(len(vv_values)):
                    #print clouds[i].n_valid, clouds[i].n_invalid
                    print(n, vv_values[n])

            #output["adv_wind_profiler"],output["std_wind_profiler"],output["max_wind_profiler"],output["min_wind_profiler"],output["dvdh_wind_profiler"]=clouds[i].horizontal_wind(output["CTH"],4)
            #if len(clouds[i].features)>100:
            #    clouds[i].FA()
            #if output["Begin_Date"]=='2018_12_30_21_05_14':
            #    exit()

            with open(outfile, "a+") as f:
                if newfile:
                    keys=keys_to_str(output)
                    f.write(keys+"\n")
                    newfile=False

                output_str=dict_to_str(output)
                f.write(output_str+"\n")

