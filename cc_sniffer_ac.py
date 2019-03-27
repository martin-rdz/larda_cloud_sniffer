#!/usr/bin/python

import sys
sys.path.append('/home/larda3/larda/')

import argparse
import CLS_Clouds as clouds
import datetime
import numpy as np
import pickle
import copy
from scipy import interpolate

import pyLARDA
import pyLARDA.Transformations as lT
import pyLARDA.helpers as h

import logging
log = logging.getLogger('pyLARDA')
#log.setLevel(logging.DEBUG)
#log.setLevel(logging.INFO)
log.setLevel(logging.WARNING)
log.addHandler(logging.StreamHandler())

parser = argparse.ArgumentParser(description=' ')
parser.add_argument('--date', metavar='date_str', type=str, help='date on which the algorithm should be run')
parser.add_argument('--campaign', type=str, help='coudnet station for which the algorithm should be run')
args = parser.parse_args()
print("starting sniffer {} {}".format(args.campaign, args.date))

def unmask(cloudnet_profile):
    
    cloudnet_profile[cloudnet_profile==999]=0
    cloudnet_profile[cloudnet_profile==-999]=0
    cloudnet_profile[cloudnet_profile==-99]=0
    cloudnet_profile[cloudnet_profile==99]=0    
    cloudnet_profile[np.isnan(cloudnet_profile)]=0
    cloudnet_profile[np.isinf(cloudnet_profile)]=0
    
    return cloudnet_profile
            

def translate_bit(bit):
    
    cn_bit_names=dict(enumerate([
        "Clear Sky","Cloud droplets only","Drizzle or rain","Drizzle/rain & cloud droplets",
        "Ice","Ice & supercooled droplets", "Melting ice", "Melting ice & cloud droplets", 
        "Aerosol","Insects","Aerosol & insects"]))
    
    return cn_bit_names.get(bit,"out of range")


def bounds_of_runs(array, maxstepsize=1):
    """ [ 0,  1,  9, 10, 11, 14, 15, 16] => [[0, 1], [ 9, 11], [14, 16]]"""
    if len(array) > 0:
        bounds = [[e[0], e[-1]] for e in \
            np.split(array, np.where(np.diff(array) > maxstepsize)[0]+1)]
    else:
        bounds = []
    return bounds



def search(elements,target):

    hit=[]
    for e in np.ravel(elements):
        hit.append(np.any(e in target))
    print('elements, target, hit', elements, target, hit)
    return np.any(hit)

def save_object(dump_object,filename):
    
    dump=open(filename,'wb')
    pickle.dump(dump_object,dump,pickle.HIGHEST_PROTOCOL)
    dump.close()
    
def load_object(filename):
    
    dump=open(filename,'rb')
    
    return pickle.load(dump)

def connect_features(detected_features,h_threshold=4000.0,v_threshold=200.0,profile_time_length=30,advection_speed=10,temporal_search_range=10,cloud_type="none"):

    #connect mixed-phase features
    pl=profile_time_length #temporal profile length [s]
    c=advection_speed #advection speed [m/s]
    d=temporal_search_range #temporal search range

    cloud_counter=0
    found_clouds=[]

    found_clouds.append(clouds.cloud())
    found_clouds[-1].cloud_type=cloud_type
    # go through all features, but with kind of a -d,d sliding window
    for i in range(d,len(detected_features)-d): 
        print(i, "out of", d, len(detected_features)-d)
        #read current feature
        f0=detected_features[i]

        connections=[]
        
        #read features adjacent in the list
        for n in range(-d,d):
            if n==0:
                continue

            f1=detected_features[i+n]

            v_dist=np.abs(f0.top_range - f1.top_range)
            h_dist=np.abs(c*(f0.time-f1.time))
            print('v_dist', v_dist, "h_dist", h_dist, v_threshold, h_threshold)
            if v_dist<v_threshold and h_dist<h_threshold:
                print('found feature close to', i, i+n) 
                # none of the features is assign to a cloud
                if f0.cloud_system==-1 and f1.cloud_system==-1:
                    f0.cloud_system=cloud_counter
                    f1.cloud_system=cloud_counter

                    found_clouds[cloud_counter].append_feature(f0)
                    found_clouds[cloud_counter].append_feature(f1)

                    cloud_counter+=1
                    found_clouds.append(clouds.cloud())
                    found_clouds[-1].cloud_type=cloud_type
                    print("new cloud system, total:", cloud_counter, "run:", cloud_type)

                # first feature is assigned to a cloud, second not
                elif f0.cloud_system!=-1 and f1.cloud_system==-1:
                    f1.cloud_system=f0.cloud_system
                    found_clouds[f1.cloud_system].append_feature(f1)

                # second is assigned to a cloud, first not
                elif f0.cloud_system==-1 and f1.cloud_system!=-1:
                    f0.cloud_system=f1.cloud_system
                    found_clouds[f0.cloud_system].append_feature(f0)                

    return found_clouds



campaign=args.campaign
build_lists=True 
larda = pyLARDA.LARDA().connect(campaign, build_lists=build_lists)

print(larda.camp.info_dict)
station_altitude=larda.camp.info_dict["altitude"]

#hand over date
begin_dt=datetime.datetime.strptime(args.date, "%Y%m%d")
end_dt=begin_dt + datetime.timedelta(hours=23, minutes=59, seconds=59)
time_interval = [begin_dt, end_dt]

for i in range(0,11):
    print(i, translate_bit(i))

#advanced datasets
lidar_present=False
doppler_present=False
windprofiler_present=False
corrected_tfv_present=False


var_shortcuts = {"cc": "CLASS", "LWC_S": "LWC_S", "LWC": "LWC",
        "IWC": "IWC",
        "Z":"Z", "v": "VEL", "uwind":"UWIND", "vwind":"VWIND",
        "T": "T", "p": "P", "beta": "beta", "width": "WIDTH"}

data = {k:larda.read("CLOUDNET", v, time_interval, [0, 'max']) for k,v in var_shortcuts.items()}


def calc_snr(data):
    var = h.lin2z(data['var']) + 10*(-2.*np.log10(data['rg']) + 2.*np.log10(5000.) - np.log10(0.00254362123253))
    return var, data['mask']
data["SNR"] = lT.combine(calc_snr, data['Z'], {'name': "SNR"})
#h.pprint(data["SNR"], verbose=True)

try:
    data["LDR"] = larda.read("CLOUDNET","LDR",time_interval, [0,'max'])
    def calc_LDR_thres(data):
        snr_co = -20.
        ldr_limit = -33
        var = snr_co-data["var"]
        mask = data["mask"]
        return var, mask

    data["LDR_thr"] = lT.combine(calc_LDR_thres, data["SNR"], {})
    data["LDR"]["mask"] = np.logical_or(data["LDR"]["mask"], h.lin2z(data["LDR"]["var"])<-33)
    data["LDR"]["mask"] = np.logical_or(data["LDR"]["mask"], h.lin2z(data["LDR"]["var"])<data["LDR_thr"]["var"])

except Exception as e:
    print("No Radar depol found")
    #create an empty depol dataset
    data['LDR']=copy.deepcopy(data["Z"])
    data['LDR']['var'][:]=-999
    data['LDR']['mask'][:]=True
    data['LDR_thr']=copy.deepcopy(data["LDR"])

data["T"] = lT.interpolate2d(data["T"], new_range=data["Z"]["rg"])
data["p"] = lT.interpolate2d(data["p"], new_range=data["Z"]["rg"])
data["uwind"] = lT.interpolate2d(data["uwind"], new_range=data["Z"]["rg"])
data["vwind"] = lT.interpolate2d(data["vwind"], new_range=data["Z"]["rg"])


try:
    data["delta"] = larda.read("CLOUDNET","depol",time_interval, [0,'max'])
    lidar_present=True
except Exception as e:
    print("Error:", e)
    print("No lidar data found, continue with lidar_present=False")

#try:
#    v_lidar=larda.read("DL_Velocity",begin,end)
#    a_lidar=larda.read("DL_Amplitude",begin,end)
#
#    #v_lidar.abs_reference(a_lidar, 0.1e9)
#    
#    doppler_present=True
#
#except Exception as e:
#    print("Error:", e)
#    print("No WiLi data found, continue with doppler_present=False")

try:
    data["v_lidar"] = larda.read("SHAUN","VEL",time_interval, [0,'max'])
    data["a_lidar"] = larda.read("SHAUN","beta_raw",time_interval, [0,'max'])
    #v_lidar.abs_reference(a_lidar, 0.1e9)
    doppler_present=True
except Exception as e:
    print("Error:", e)
    print("No SHAUN data found, continue with doppler_present=False")

#try:
#    wp_vel=larda.read("WIPRO_ADV",begin-3600,end+3600)
#    wp_dir=larda.read("WIPRO_DIR",begin-3600,end+3600)
#
#    windprofiler_present=True
#
#except Exception as e:
#    print("Error:", e)
#    print("No Wind Profiler found, continue with windprofiler_present=False")
#
#try:
#
#    tfv=larda.read("CLOUDNET_TFV",begin,end)
#    vair=larda.read("CLOUDNET_VAIR",begin,end)
#    
#    #TODO: noise cut in original wind profiler converted files
#
#    tfv.data[tfv.data<-6.0]=-999
#    tfv.data[tfv.data>6.0]=-999
#
#    vair.data[vair.data<-6.0]=-999
#    vair.data[vair.data>6.0]=-999
#
#    corrected_tfv_present=True
#
#except Exception as e:
#    print("Error:", e)
#    print("No terminal fall velocity and air velocity information found")


features_in_timestep=[]

bridge=5 #Possible gap between a liquid cloud base and radar precipitation below [Cloudnet height steps]
min_tower_height=1000.0 #minimum height of a cloud feature to be classified as part of a cumulus congestus or nimbostratus cloud


verbose = False

for i in range(data["cc"]["ts"].shape[0]):
    #i = 2502
    #classfication profile
    profiles = {}
    profiles['cc'] = lT.slice_container(data['cc'], index={'time': [i]})
    #h.pprint(profiles['cc'])
    profiles['IWC'] = lT.slice_container(data['IWC'], index={'time': [i]})
    profiles['LWC_S'] = lT.slice_container(data['LWC_S'], index={'time': [i]})
    profiles['LWC'] = lT.slice_container(data['LWC'], index={'time': [i]})
    profiles['Z'] = lT.slice_container(data['Z'], index={'time': [i]})
    profiles['SNR'] = lT.slice_container(data['SNR'], index={'time': [i]})
    profiles['LDR'] = lT.slice_container(data['LDR'], index={'time': [i]})
    profiles['LDR_thr'] = lT.slice_container(data['LDR_thr'], index={'time': [i]})
    profiles['v'] = lT.slice_container(data['v'], index={'time': [i]})
    profiles['width'] = lT.slice_container(data['width'], index={'time': [i]})
    profiles['beta'] = lT.slice_container(data['beta'], index={'time': [i]})
   
    #if ~np.any(cc_profile==1) and ~np.any(cc_profile==3):
    #    continue

    if lidar_present:
        profiles['delta'] = lT.slice_container(data['delta'], index={'time': [i]})

    if doppler_present:
        #ih.pprint(data["v_lidar"])
        it_b_dl = h.argnearest(data["v_lidar"]['ts'], data["cc"]["ts"][i]-15)
        it_e_dl = h.argnearest(data["v_lidar"]['ts'], data["cc"]["ts"][i]+15)
        if not it_b_dl == it_e_dl:
            print("no doppler lidar for this profile", i)
            profiles['v_lidar'] = lT.slice_container(data['v_lidar'], 
                    index={'time': [it_b_dl, it_e_dl]})
            profiles['a_lidar'] = lT.slice_container(data['a_lidar'], 
                    index={'time': [it_b_dl, it_e_dl]})

#    if windprofiler_present:
#        wp_timebin=wp_vel.get_time_bin(Z.times[i])
#    
#    if corrected_tfv_present:
#        tfv_bin=tfv.get_time_bin(v.times[i])
#        if tfv.times[tfv_bin]-v.times[i]>75.0:
#            tfv_profile=copy.deepcopy(v.data[i])
#            tfv_profile[:]=0
#            vair_profile=copy.deepcopy(v.data[i])
#            vair_profile[:]=0
#        else:
#            tfv_profile=unmask(tfv.data[tfv_bin])
#            vair_profile=unmask(vair.data[tfv_bin])
#

    profiles['T'] = lT.slice_container(data['T'], index={'time': [i]})
    profiles['p'] = lT.slice_container(data['p'], index={'time': [i]})
    profiles['uwind'] = lT.slice_container(data['uwind'], index={'time': [i]})
    profiles['vwind'] = lT.slice_container(data['vwind'], index={'time': [i]})

    features_in_profile=[]
    bounds_imperative = []
    #print cc_profile
    print("Time:",i)
    ih = len(profiles["cc"]['rg']) - 2

    # this assignment is possible outside the loop?
    cc_profile = profiles["cc"]['var']
    print('cc_profile', cc_profile)
    
    indices_cloud = np.where(np.isin(profiles['cc']['var'], clouds.cloud_particles))[0]
    # 1. option: bridge everything up to 5 bins
    bounds_feature = bounds_of_runs(indices_cloud, maxstepsize=5)
    # 2. option: bridge up to 5 bins only if top is liquid
    bounds_feature = bounds_of_runs(indices_cloud, maxstepsize=1)
    merge_after = []
    # find bounds that fulfill the cond and save the index
    for i, [lower, upper] in enumerate(zip(bounds_feature, bounds_feature[1::])):
        cond = [profiles['cc']['var'][upper[0]] in clouds.liquid,
                np.any(np.isin(profiles['cc']['var'][upper[0]-6:upper[0]], clouds.ice))]
        merge_after.append(i) if all(cond) else None
    # and now merge at these bounds
    for i in merge_after[::-1]:
        bounds_feature[i][1] = bounds_feature[i+1][1]
        bounds_feature.pop(i+1)

    # topmost pixels are not considered
    bounds_feature = list(filter(lambda x: x[0] != len(profiles['cc']['rg'])-1, bounds_feature))
    print('bounds from function   ', len(bounds_feature), bounds_feature) if verbose else None

    #populate the feature with data
    features_in_profile = []
    for b in bounds_feature:
        #print('bound ', b)
        b_top = min(b[1]+1, cc_profile.shape[0]-1)
        #populate the feature with data
        feature = clouds.cloud_feature()
        feature.time = profiles['cc']['ts']
        feature.classifications = cc_profile[b[0]:b_top]
        feature.base_range = profiles["cc"]["rg"][b[0]]
        feature.top_range = profiles["cc"]["rg"][b_top]
        feature.ranges = profiles["cc"]["rg"][b[0]:b_top]
                        
        feature.has_melting=search(clouds.melting, feature.classifications)
        feature.has_drizzle=search(clouds.drizzle, feature.classifications)
        keys_to_feature = ["IWC", "LWC_S", "LWC", "Z", "v", "width", "T", "p", "SNR",
                "uwind", "vwind", "beta", "LDR"]
        if lidar_present:
            keys_to_feature += ["delta"]

        for k in keys_to_feature:
            feature.measurements[k] = lT.slice_container(profiles[k], index={"range": [b[0], b_top]})

        def calc_alpha_hogan(datalist):
            T_C = datalist[1]['var']
            Zdb = h.lin2z(datalist[0]['var'])
            var = 10**(0.000477*Zdb*T_C + 0.0683*Zdb - 0.0171*T_C - 3.11)
            mask = np.logical_or(datalist[0]["mask"], datalist[1]["mask"])
            return var, mask
        feature.measurements[k] = lT.combine(calc_alpha_hogan, 
                [feature.measurements['Z'], feature.measurements["T"]], 
                {'name': "alpha_hogan"})

        if "v_lidar" in profiles:
            ir_b_dl = h.argnearest(profiles["v_lidar"]["rg"], feature.base_range)
            ir_t_dl = h.argnearest(profiles["v_lidar"]["rg"], feature.top_range)
            lidar_range = [ir_b_dl, ir_t_dl+1]
            feature.measurements["v_lidar"] = lT.slice_container(
                    profiles['v_lidar'], index={"range": lidar_range})
            feature.measurements["a_lidar"] = lT.slice_container(
                    profiles['a_lidar'], index={"range": lidar_range})

        features_in_profile.append(feature)

    #if no features are detected, skip
    # empty list is skipped in for anyways
    #if len(features_in_profile)==0:
    #    continue
    #Classify detected features
    print("classify features") if verbose else None
    for f in features_in_profile:
        
        #search for the presence of ice crystals and liquid droplets
        ice_present = search(clouds.ice,f.classifications)
        liquid_present = search(clouds.liquid,f.classifications)
        melting_present = search(clouds.melting,f.classifications)
        print("ice, liquid, melting", ice_present, liquid_present, melting_present) 
        #discriminate feature types
        ignore_assert = False
        if ice_present and not(liquid_present) and not(melting_present):
            
            f.type="pure_ice"
            
        elif not(ice_present) and liquid_present:
            
            f.type="pure_liquid"
        
        elif melting_present and ice_present and not(liquid_present):
        
            if f.top_range-f.base_range<min_tower_height:
                f.type="pure_ice"
            else:
                f.type="tower"
            
        # tower not classified when
        # ice and liqud are present, only liquid and melting
        # TODO make tower the last/top question
        elif ice_present and liquid_present:
            print('ice and liquid present') if verbose else None
            #
            #search for end of liquid sub-layer
            print(f.base_range, f.top_range) if verbose else None
            ih=len(f.classifications)-1
            feature_cc = f.classifications.copy()
            print("ice and liquid", f.classifications) if verbose else None
            f.type = 'mixed-phase'
    
            # maybe the tower classification is connected to the 
            # overwriting of the classification
            indices_liquid = np.where(np.isin(feature_cc, clouds.liquid))[0]
            bounds_liquid = bounds_of_runs(indices_liquid, maxstepsize=1)
            print(bounds_liquid) if verbose else None

            if feature_cc[0] in clouds.liquid:
                f.type = 'liquid-based'
            #elif feature_cc[0] in clouds.melting:
            #    # should be named melting instead
            #    type = 'liquid-based'
            else:
                f.type = 'mixed-phase'

            # other option ice with liquid at top
            # case not covered liquid at base and liquid in between
            indices_ice = np.where(np.isin(feature_cc, clouds.ice_only))[0]
            if np.any(np.isin(indices_liquid, indices_ice+1)) \
                    and feature_cc[-1] in clouds.liquid:
                print("yeah found transition, overwrite type") if verbose else None
                f.type = 'mixed-phase'
            # and ice above and below liquid?
            # 1 1 1 1 5 4 4 4 4 4 5 5 5 5 4 4 4 4 
            # [1 1 1 1 1 1 1 5 4 4 4 4 5 5 5 5 5 5 5 4 4 4]
            f.liquid_layers = len(bounds_liquid)
            print("no liquid", f.liquid_layers) if verbose else None
            f.liquid_layer_base = [f.ranges[il[0]] for il in bounds_liquid]
            f.liquid_layer_top = [f.ranges[il[1]] for il in bounds_liquid]
            print('liquid layer base', f.liquid_layer_base) if verbose else None
            print('liquid layer top', f.liquid_layer_top) if verbose else None

            # index where ice precipiation starts
            f.precipitation_top = indices_liquid[0]
            print('precip_top ', f.precipitation_top) if verbose else None

            # tower is more important than the categories above
            if f.top_range-min(f.liquid_layer_base) > min_tower_height:
                f.type="tower"

            #if f.type != 'liquid-based': # for some reason there are as many liquid layers as bins in feature
            #    assert f.liquid_layers == len(bounds_liquid)
            #input()
        print("time", i)

    #if f.type == 'tower':
    #    input()
    features_in_timestep.append(features_in_profile)


#extract mixed-phase features and put them in linear array
detected_features_mixed=[]
detected_features_ice=[]
detected_features_tower=[]

for features_in_profile in features_in_timestep:
    for f in features_in_profile:
        if f.type=="mixed-phase" or f.type=="liquid-based" or f.type=="pure_liquid":
            detected_features_mixed.append(f)
        elif f.type=="pure_ice":
            detected_features_ice.append(f)
        elif f.type=="tower":
            detected_features_tower.append(f)

print("Searching for layered clouds")
clouds_mixed=connect_features(detected_features_mixed,h_threshold=4000.0,v_threshold=200.0,cloud_type="layered")
#clouds_mixed = []
print("Searching for cirrus (pure ice) clouds")
#
# 10000 as h_threshold seems a littlebit too much
clouds_ice=connect_features(detected_features_ice,h_threshold=5000.0,v_threshold=500,cloud_type="ice")
print("Searching for deep clouds")
clouds_tower=connect_features(detected_features_tower,h_threshold=5000.0,v_threshold=500.0,cloud_type="tower")
#clouds_tower = []

all_clouds = clouds_mixed + clouds_ice + clouds_tower

cloud_rectangles=[]
for i, cloud in enumerate(all_clouds):
    
    if all_clouds[i].n_profiles()==0:
        continue
    
    c_type=cloud.most_common_type()
    print("cloud type ", cloud.cloud_type)
    print("len of feature", len(cloud.features))
    #clouds[i].type=c_type
    if c_type=="pure_liquid":
        color='blue'
    elif c_type=="pure_ice":
        color='green'
    elif c_type=="liquid-based":
        color='blue'
    elif c_type=="mixed-phase":
        color='red'
    elif c_type=="tower":
        color="black"
    else:
        color='gray'
    
    cg = cloud.geometry()
    
    #if clouds[i].top_variation()<200.0 and clouds[i].time_length()>1800 and clouds[i].cloud_top_thickness()[0]<400.0 and  clouds[i].fill_factor()>0.75:
    #if c_type=="tower":
    if cloud.time_length()>900 and cloud.fill_factor()>0.60:
        cloud_rectangles.append((cg[0],cg[1],cg[2],cg[3],color))

print("Saving")
save_object(all_clouds,'../cloud_properties/cloud_properties_'+campaign+'/'+begin_dt.strftime("%Y_%m_%d")+'_clouds.dat')

fig, ax = lT.plot_timeheight(data['cc'], range_interval=[0, 12000])
import matplotlib.patches as patches
for cm in cloud_rectangles:
    print(cm)
    begin = h.ts_to_dt(cm[0])
    duration=datetime.timedelta(seconds=cm[2])
    rect = patches.Rectangle(
            (begin,cm[1]),duration,cm[3],linewidth=2,
            edgecolor=cm[4],facecolor=cm[4],alpha=0.2)

    # Add the patch to the Axes
    ax.add_patch(rect)

savename='../plots/cloud_overview_'+campaign+'/'+begin_dt.strftime("%Y_%m_%d")+'_cloud_detection.png'
fig.savefig(savename, dpi=250)


#cc is here the larda object with cc
#cc.plot(begin_dt,end_dt,0,12000,savename='../plots/cloud_overview_'+campaign+'/'+begin_dt.strftime("%Y_%m_%d")+'_cloud_detection.png',cloud_marks=cloud_rectangles)
#cc.plot(begin,end,0,12000,savename='cloud_overview_'+campaign+'/'+begin_dt.strftime("%Y_%m_%d")+'_cloud_detection.png')

