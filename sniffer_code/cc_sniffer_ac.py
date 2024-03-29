#!/usr/bin/python

import sys, os
sys.path.append('/home/larda3/larda/')
import traceback

import argparse
import CLS_Clouds as clouds
import datetime
import numpy as np
import pickle, csv
import copy
from scipy import interpolate
from collections import Counter

import pyLARDA
import pyLARDA.Transformations as lT
import pyLARDA.helpers as h

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('--date', metavar='date_str', type=str, help='date on which the algorithm should be run')
    parser.add_argument('--campaign', type=str, help='coudnet station for which the algorithm should be run')
    args = parser.parse_args()
    print("starting sniffer {} {}".format(args.campaign, args.date))


# default settings
# ====================================
MAX_VERT_GAP = 3

CONN_THRES_LAYERED = {
    # try with 3000 = 5min as threshold, otherwise
    # there will be connections over the scans
    'h_thres': 3000.0, 'v_thres': 200.0
}

CONN_THRES_ICE = {
    'h_thres': 5000.0, 'v_thres': 500.0
}

CONN_THRES_TOWER = {
    'h_thres': 5000.0, 'v_thres': 500.0
}

#minimum height of a cloud feature to be classified as part of a 
# cumulus congestus or nimbostratus cloud
MIN_TOWER_HEIGHT = 2000

THRES_SHALLOW_PREC = 4000
# ====================================


def unmask(cloudnet_profile):
    
    cloudnet_profile[cloudnet_profile==999]=0
    cloudnet_profile[cloudnet_profile==-999]=0
    cloudnet_profile[cloudnet_profile==-99]=0
    cloudnet_profile[cloudnet_profile==99]=0    
    cloudnet_profile[np.isnan(cloudnet_profile)]=0
    cloudnet_profile[np.isinf(cloudnet_profile)]=0
    
    return cloudnet_profile
            

def translate_bit(bit):
    """get the string from the bit"""
    
    cn_bit_names=dict(enumerate([
        "Clear Sky","Cloud droplets only","Drizzle or rain","Drizzle/rain & cloud droplets",
        "Ice","Ice & supercooled droplets", "Melting ice", "Melting ice & cloud droplets", 
        "Aerosol","Insects","Aerosol & insects"]))
    
    return cn_bit_names.get(bit,"out of range")


def bounds_of_runs(array, maxstepsize=1):
    """transfer a list of indices in to nested list of the boudaries of runs

    [ 0,  1,  9, 10, 11, 14, 15, 16] => [[0, 1], [ 9, 11], [14, 16]]
    """
    if len(array) > 0:
        bounds = [[e[0], e[-1]] for e in \
            np.split(array, np.where(np.diff(array) > maxstepsize)[0]+1)]
    else:
        bounds = []
    return bounds


def search(elements,target):
    """check if element contains any target"""
    hit=[]
    for e in np.ravel(elements):
        hit.append(np.any(e in target))
    #print('elements, target, hit', elements, target, hit)
    return np.any(hit)


def save_object(dump_object,filename):
    
    dump=open(filename,'wb')
    pickle.dump(dump_object,dump,pickle.HIGHEST_PROTOCOL)
    dump.close()

    
def load_object(filename):
    
    dump=open(filename,'rb')
    
    return pickle.load(dump)

def connect_features(detected_features, 
    h_thres=None, v_thres=None, 
    profile_time_length=30, advection_speed=10, temporal_search_range=10, 
    cloud_type="none", verbose=False):
    """connect features (chunks of single profiles based on the thresholds)

    Args:
        detected_features: list of detected features
        h_thres: horizontal threshold in [m] 
        v_thres: vertical threshold in [m]
        profile_time_length: default 30s 
        advection_speed: default 10 [m/s]
        temporal_search_range: default 10
        cloud_type: string identifying cloud type
    """

    #connect mixed-phase features
    d = temporal_search_range #temporal search range [no. profiles?]

    cloud_counter=0
    found_clouds=[]
    
    assert h_thres is not None
    assert v_thres is not None
    assert profile_time_length*advection_speed*2*temporal_search_range \
        > h_thres, "search window smaller than threshold"

    found_clouds.append(clouds.cloud())
    found_clouds[-1].cloud_type=cloud_type
    # go through all features, but with kind of a -d,d sliding window
    print('number features ', len(detected_features), d)
    for i in range(d,len(detected_features)-d): 
        print(i, "out of", d, len(detected_features)-d) if verbose else None
        #read current feature
        f0=detected_features[i]

        connections=[]
        #read features adjacent in the list
        for n in range(-d,d):
            if n==0:
                continue

            f1=detected_features[i+n]

            v_dist=np.abs(f0.top_range - f1.top_range)
            h_dist=np.abs(advection_speed*(f0.time-f1.time))
            #print('v_dist', v_dist, "h_dist", h_dist, v_threshold, h_thres)
            if v_dist<v_thres and h_dist<h_thres:
                #print('found feature close to', i, i+n) 
                # none of the features is assign to a cloud
                if f0.cloud_system==-1 and f1.cloud_system==-1:
                    f0.cloud_system=cloud_counter
                    f1.cloud_system=cloud_counter

                    found_clouds[cloud_counter].append_feature(f0)
                    found_clouds[cloud_counter].append_feature(f1)

                    cloud_counter+=1
                    found_clouds.append(clouds.cloud())
                    found_clouds[-1].cloud_type=cloud_type
                    #print("new cloud system, total:", cloud_counter, "run:", cloud_type)

                # first feature is assigned to a cloud, second not
                elif f0.cloud_system!=-1 and f1.cloud_system==-1:
                    f1.cloud_system=f0.cloud_system
                    found_clouds[f1.cloud_system].append_feature(f1)

                # second is assigned to a cloud, first not
                elif f0.cloud_system==-1 and f1.cloud_system!=-1:
                    f0.cloud_system=f1.cloud_system
                    found_clouds[f0.cloud_system].append_feature(f0)                

    # remove empty clouds
    found_clouds = list(filter(lambda x: x.begin_time != -1 and x.end_time != -1, found_clouds))
    return found_clouds



def find_features_in_profile(profiles, keys_to_feature, verbose=False):
    """ 
    
    Args:
        profiles (dict): range slices of data, cc with classification flag is required
        keys_to_features (list): keys that are added to features
        verbose (optional, bool): verbose print output
    
    """

    #classfication profile

    #if ~np.any(cc_profile==1) and ~np.any(cc_profile==3):
    #    continue

    #print cc_profile
    #ih = len(profiles["cc"]['rg']) - 2

    # this assignment is possible outside the loop?
    cc_profile = profiles["cc"]['var']
    #print('cc_profile', cc_profile)
    
    indices_cloud = np.where(np.isin(profiles['cc']['var'], clouds.cloud_particles))[0]
    # original comment: Possible gap between a liquid cloud base and radar precipitation below [Cloudnet height steps]
    # 1. option: bridge everything up to 5 bins
    bounds_feature = bounds_of_runs(indices_cloud, maxstepsize=MAX_VERT_GAP)
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
    #print('bounds from function   ', len(bounds_feature), bounds_feature) if verbose else None

    #populate the feature with data
    features_in_profile = []
    for b in bounds_feature:
        #print('bound ', b)
        b_top = min(b[1]+1, cc_profile.shape[0]-1)
        #populate the feature with data
        feature = clouds.cloud_feature()
        feature.time = profiles['cc']['ts']
        feature.classifications = cc_profile[b[0]:b_top]
        feature.classification_file_history = profiles['cc']['meta']['history']
        feature.base_range = profiles["cc"]["rg"][b[0]]
        feature.top_range = profiles["cc"]["rg"][b_top]
        feature.ranges = profiles["cc"]["rg"][b[0]:b_top]
        feature.dh = profiles['cc']['rg'][1] - profiles['cc']['rg'][0]
                        
        feature.has_melting=search(clouds.melting, feature.classifications)
        feature.has_drizzle=search(clouds.drizzle, feature.classifications)

        for k in keys_to_feature:
            feature.measurements[k] = lT.slice_container(profiles[k], index={"range": [b[0], b_top]})

        def calc_alpha_hogan(datalist):
            T_C = datalist[1]['var'] - 273.15
            Zdb = h.lin2z(datalist[0]['var'])
            var = 10**(0.000477*Zdb*T_C + 0.0683*Zdb - 0.0171*T_C - 3.11)
            mask = np.logical_or(datalist[0]["mask"], datalist[1]["mask"])
            return var, mask
        if 'Z' in profiles:
            feature.measurements['alpha_hogan'] = lT.combine(calc_alpha_hogan, 
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

        if 'pT_no' in profiles:
            ir_b_pT = h.argnearest(profiles["pT_no"]["rg"], feature.base_range)
            ir_t_pT = h.argnearest(profiles["pT_no"]["rg"], feature.top_range)
            feature.measurements["pT_no"] = lT.slice_container(
                    profiles['pT_no'], index={"range": [ir_b_pT, ir_t_pT+1]})
            print(feature.measurements["pT_no"]['var'])

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
        print("ice, liquid, melting", ice_present, liquid_present, melting_present) if verbose else None 
        #discriminate feature types
        ignore_assert = False
        if ice_present and not(liquid_present) and not(melting_present):
            
            f.type="pure_ice"
            
        elif not(ice_present) and liquid_present:
            
            f.type="pure_liquid"
            #
            # 2020-10-12: liquid only layers were missing the liquid layer base, which
            #             introduced errors in the statistics (i.e. when using base_range)
            #
            feature_cc = f.classifications.copy()
            # maybe the tower classification is connected to the 
            # overwriting of the classification
            indices_liquid = np.where(np.isin(feature_cc, clouds.liquid))[0]
            bounds_liquid = bounds_of_runs(indices_liquid, maxstepsize=1)
            f.liquid_layers = len(bounds_liquid)
            f.liquid_layer_base = [f.ranges[il[0]] for il in bounds_liquid]
            f.liquid_layer_top = [f.ranges[il[1]] for il in bounds_liquid]
        
        elif melting_present and ice_present and not(liquid_present):
        
            if f.top_range - f.base_range < MIN_TOWER_HEIGHT:
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
            #ih=len(f.classifications)-1
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

            # print('sniffer ', f.classifications, f.base_range, f.liquid_layer_base)

            # index where ice precipiation starts
            f.precipitation_top = indices_liquid[0]
            print('precip_top ', f.precipitation_top) if verbose else None

            # tower is more important than the categories above
            if f.top_range-min(f.liquid_layer_base) > MIN_TOWER_HEIGHT:
                f.type="tower"

            #if f.type != 'liquid-based': # for some reason there are as many liquid layers as bins in feature
            #    assert f.liquid_layers == len(bounds_liquid)
            #input()
        print(f"final class of feature {f.type}, ({h.ts_to_dt(f.time).strftime('%Y%m%d-%H%M%S')} {f.base_range:.2f})")

    #if f.type == 'tower':
    #    input()
    return features_in_profile




def classify_wo_tower(classifications):

    ice_present = search(clouds.ice, classifications)
    liquid_present = search(clouds.liquid, classifications)
    melting_present = search(clouds.melting, classifications)

    #print("ice, liquid, melting", ice_present, liquid_present, melting_present) 
    if ice_present and not(liquid_present):
        ctype="pure_ice"
    elif not(ice_present) and liquid_present:
        ctype="pure_liquid"
    elif ice_present and liquid_present:
        ctype = 'mixed-phase'
    return ctype


def features_to_simple_class(features_in_profile):

    if not features_in_profile:
        #print('clear profile')
        return 'clear'

    rain = list(filter(lambda f: np.any(np.isin(f.classifications, [2,3,6,7])), 
                         features_in_profile))
    #print('precip ', len(rain))
    if rain:
        #print([(f.base_range, f.top_range) for f in rain])
        if rain[-1].top_range - rain[-1].base_range >= THRES_SHALLOW_PREC:
            return 'rain_deep'
        else:
            return 'rain_shallow'

    if len(features_in_profile) > 1:
        ctype = [classify_wo_tower(f.classifications) for f in features_in_profile]
        #print([f.classifications for f in features_in_profile])
        #print('multi_layer present ', ctype)
        most_common = Counter(ctype).most_common(1)[0]
        if most_common[1]/len(features_in_profile) > 0.66:
            return 'multi_layer_mostly_' + most_common[0]
        else:
            return "multi_layer_various"
    else:
        #print(features_in_profile[0].type)
        return 'single_layer_' + classify_wo_tower(features_in_profile[0].classifications)


if __name__ == "__main__":
    import logging
    log = logging.getLogger('pyLARDA')
    #log.setLevel(logging.DEBUG)
    #log.setLevel(logging.INFO)
    log.setLevel(logging.WARNING)
    log.addHandler(logging.StreamHandler())

    campaign=args.campaign
    build_lists=False 
    larda = pyLARDA.LARDA().connect(campaign, build_lists=build_lists)


    #larda_rsd2 = pyLARDA.LARDA('remote', uri = 'http://larda3.tropos.de').connect(campaign, build_lists=build_lists)
    larda_rsd2 = pyLARDA.LARDA('local').connect(campaign, build_lists=build_lists)

    print(larda.camp.info_dict)
    larda_polly_ari = pyLARDA.LARDA('local').connect('polly_arielle_leipzig', build_lists=build_lists)
    larda_polly_ift = pyLARDA.LARDA('local').connect('polly_ift_leipzig', build_lists=build_lists)

    #hand over date
    begin_dt=datetime.datetime.strptime(args.date, "%Y%m%d")
    end_dt=begin_dt + datetime.timedelta(hours=23, minutes=59, seconds=59)
    time_interval = [begin_dt, end_dt]

    for i in range(0,11):
        print(i, translate_bit(i))

    #advanced datasets
    lidar_present=False
    doppler_present=False
    qbsc_present=False
    windprofiler_present=False
    corrected_tfv_present=False
    peakTree_present=False


    var_shortcuts = {"cc": "CLASS", "LWC_S": "LWC_S", "LWC": "LWC",
            "IWC": "IWC",
            "Z":"Z", "v": "VEL", "uwind":"UWIND", "vwind":"VWIND",
            "T": "T", "p": "P", "beta": "beta", "width": "WIDTH"}

    #data = {k:larda.read("CLOUDNET", v, time_interval, [0, 'max']) for k,v in var_shortcuts.items()}
    data = {}
    for k,v in var_shortcuts.items():
        print('loading ', k, v)
        data[k] = larda.read("CLOUDNET", v, time_interval, [0, 'max'])
        #data[k]['var'] = data[k]['var'].data.astype(np.float64)

    def calc_snr(data):
        var = h.lin2z(data['var']) + 10*(-2.*np.log10(data['rg']) + 2.*np.log10(5000.) - np.log10(0.00254362123253))
        return var, data['mask']
    data["SNR"] = lT.combine(calc_snr, data['Z'], {'name': "SNR"})
    data["SNR"]['var'] = data["SNR"]['var'].data.astype(np.float64)

    #try:
    data["LDR"] = larda.read("CLOUDNET","LDR",time_interval, [0,'max'])
    #data["LDR"]['var'] = data["LDR"]['var'].data.astype(np.float64)
    #data["LDR"]['mask'] = data["LDR"]['mask'].data
    def calc_Zcx(data):
        Z = data[0]
        LDR = data[1]
        var = Z['var']*LDR['var']
        mask = np.logical_or(Z['mask'], LDR['mask'])
        return var, mask
    data['Zcx'] = lT.combine(calc_Zcx, [data['Z'], data['LDR']], {'name': "Zcx", 'var_lims': [-47,-20]})

    sensitivity_cx = np.broadcast_to((data['Zcx']['rg']**2*h.z2lin(-34)/5000**2), data['Zcx']['var'].shape)
    new_ldr_mask = np.logical_or(data['Zcx']['var'] < sensitivity_cx, data['Zcx']['var'] < data['Z']['var']/h.z2lin(-30))
    sensitivity_mask = (data['Zcx']['var'] < sensitivity_cx)
    decoupling_mask = data['Zcx']['var'] < data['Z']['var']*h.z2lin(-30)
    new_ldr_mask = np.logical_or(sensitivity_mask, decoupling_mask)
    new_ldr_mask = np.logical_or(new_ldr_mask, data['LDR']['var'] > h.z2lin(-11))
    data['LDRcorr'] = {**data['LDR']}
    data['LDRcorr']['mask'] = np.logical_or(data['LDR']['mask'], new_ldr_mask)


    try:
        data["qbsc"] = larda_rsd2.read("POLLYNET","qbsc532",time_interval, [0, 'max'])
        #data["qbsc"]['var'] = data["qbsc"]['var'].data.astype(np.float64)
        qbsc_present = True
    except:
        traceback.print_exc()
        print('quasibackscatter not available from lacros polly')

    time_interval = [begin_dt-datetime.timedelta(minutes=5), end_dt]
    if campaign == 'lacros_leipzig' and not qbsc_present:
        try:
            data["qbsc"] = larda_polly_ari.read("POLLYNET","qbsc532",time_interval, [0, 'max'])
            #data["qbsc"]['var'] = data['qbsc']['var'].data.astype(np.float64)
            qbsc_present = True
        except:
            traceback.print_exc()
            print('quasibackscatter not available from arielle')
        if not qbsc_present:
            try:
                data["qbsc"] = larda_polly_ift.read("POLLYNET","qbsc532",time_interval, [0, 'max'])
                data["qbsc"]['var'] = data["qbsc"]['var'].data.astype(np.float64)
                qbsc_present = True
            except:
                traceback.print_exc()
                print('quasibackscatter not available from polly ift')

    #try:
    #    data["voldepol"] = larda.read("CLOUDNET","depol",time_interval, [0,'max'])
    #    lidar_present=True
    #except Exception as e:
    #    print("Error:", e)
    #    print("No lidar data found, continue with lidar_present=False")
    try:
        voldepol = larda_rsd2.read("POLLYNET","voldepol532",time_interval, [0, 'max'])
        #voldepol['var'] = voldepol['var'].data.astype(np.float64)
        lidar_present=True
    except:
        traceback.print_exc()
        print("No lidar data found pollyxt_lacros")

    if campaign == 'lacros_leipzig' and not lidar_present:
        try:
            voldepol = larda_polly_ari.read("POLLYNET","voldepol532",time_interval, [0, 'max'])
            #voldepol['var'] = voldepol['var'].data.astype(np.float64)
            lidar_present=True
        except:
            traceback.print_exc()
            print('voldepol not available from arielle')
        if not lidar_present:
            try:
                voldepol = larda_polly_ift.read("POLLYNET","voldepol532",time_interval, [0, 'max'])
                #voldepol['var'] = voldepol['var'].data.astype(np.float64)
                lidar_present=True
            except:
                traceback.print_exc()
                print('voldepol not available from polly ift')

    if lidar_present:
        voldepol['var'][~np.isfinite(voldepol['var'])] = 0.0
        voldepol['mask'] = np.logical_or(voldepol['mask'], ~np.isfinite(voldepol['var']))
        print(voldepol['ts'])
        data['voldepol'] = pyLARDA.Transformations.interpolate2d(
            voldepol, new_time=data['cc']['ts'], new_range=data['cc']['rg'])

    if qbsc_present:
        data['qbsc']['var'][~np.isfinite(data['qbsc']['var'])] = 0.0
        data['qbsc']['mask'] = np.logical_or(data['qbsc']['mask'], ~np.isfinite(data['qbsc']['var']))
        qbsc_interp = pyLARDA.Transformations.interpolate2d(
            data['qbsc'], new_time=data['cc']['ts'], new_range=data['cc']['rg'])
        def calc_ice_qext(datalist):
            qext = datalist[0]['var']*32
            mask = np.logical_or(qext < 1e-10, datalist[0]['mask'])
            return qext, mask
        data['qiceext'] = pyLARDA.Transformations.combine(
            calc_ice_qext, [qbsc_interp], 
            {'var_unit': 'm^-1', 'var_lims': [5e-6,1e-2], 'name': 'ice qext'})
        data['qbsc'] = qbsc_interp
    
        def calc_z_e(datalist):
            assert np.all(datalist[0]['ts'] == datalist[1]['ts'])
            assert np.all(datalist[0]['rg'] == datalist[1]['rg'])
            ratio = datalist[0]['var']/datalist[1]['var']
            mask = np.logical_or(datalist[0]['mask'], datalist[1]['mask'])
            mask = np.logical_or(mask, ~np.isfinite(ratio))
            ratio[mask] = 0.0
            return ratio, mask
    
        data["ratio_z_e"] = pyLARDA.Transformations.combine(
            calc_z_e, [data['Z'], data['qiceext']], 
            {'var_unit': 'm mm^-6', 'var_lims': [1e-0,3e3], 'name': 'Z/E'})


    #new_ldr_mask = np.logical_or(data['Zcx']['var'] < h.z2lin(-37), data['Zcx']['var'] > (data['Z']['var']-h.z2lin(-32)))
    #data['LDRcorr'] = {**data['LDR']}
    #data['LDRcorr']['mask'] = np.logical_or(data['LDR']['mask'], new_ldr_mask)

    #except Exception as e:
    #    print("No Radar depol found")
    #    #create an empty depol dataset
    #    data['LDR']=copy.deepcopy(data["Z"])
    #    data['LDR']['var'][:]=-999
    #    data['LDR']['mask'][:]=True
    #    data['LDRcorr']=copy.deepcopy(data["LDR"])

    # interpolate the model data
    data["T"] = lT.interpolate2d(data["T"], new_range=data["Z"]["rg"])
    data["p"] = lT.interpolate2d(data["p"], new_range=data["Z"]["rg"])
    data["uwind"] = lT.interpolate2d(data["uwind"], new_range=data["Z"]["rg"])
    data["vwind"] = lT.interpolate2d(data["vwind"], new_range=data["Z"]["rg"])

    try:
        data["v_lidar"] = larda.read("SHAUN","VEL",time_interval, [0,'max'])
        data["a_lidar"] = larda.read("SHAUN","beta_raw",time_interval, [0,'max'])
        #data["v_lidar"]['var'] = data['v_lidar']['var'].data.astype(np.float64)
        #data['a_lidar']['var'] = data['a_lidar']['var'].data.astype(np.float64)
        data["v_lidar"]['mask'] = np.logical_or(data["v_lidar"]['mask'],
                                                np.isclose(data['v_lidar']['var'], -999))
        data["a_lidar"]['mask'] = np.logical_or(data["a_lidar"]['mask'],
                                                np.isclose(data['a_lidar']['var'], -999))
        #v_lidar.abs_reference(a_lidar, 0.1e9)
        doppler_present = True
    except:
        traceback.print_exc()
        print("No SHAUN data found, continue with doppler_present=False")


    try:
        data["pT_no"] = larda.read("peakTree","no_nodes",time_interval, [0,'max'])
        peakTree_present = True
    except:
        traceback.print_exc()
        print("No peakTree data found, continue with peakTree_present=False")
        #input()

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

    for k in data.keys():
        print('checking ', k, type(data[k]['mask']), type(data[k]['var']))
        assert not isinstance(data[k]['mask'], np.ma.MaskedArray), data[k]['mask']
        assert not isinstance(data[k]['var'], np.ma.MaskedArray), data[k]['var']

    verbose = False

    features_in_timestep = []
    simple_class = []
    for i in range(data["cc"]["ts"].shape[0]):
        
        print("Time:",i)
        profiles = {}
        profiles['cc'] = lT.slice_container(data['cc'], index={'time': [i]})
        #h.pprint(profiles['cc'])
        profiles['IWC'] = lT.slice_container(data['IWC'], index={'time': [i]})
        profiles['LWC_S'] = lT.slice_container(data['LWC_S'], index={'time': [i]})
        profiles['LWC'] = lT.slice_container(data['LWC'], index={'time': [i]})
        profiles['Z'] = lT.slice_container(data['Z'], index={'time': [i]})
        profiles['SNR'] = lT.slice_container(data['SNR'], index={'time': [i]})
        profiles['LDR'] = lT.slice_container(data['LDR'], index={'time': [i]})
        profiles['LDRcorr'] = lT.slice_container(data['LDRcorr'], index={'time': [i]})
        profiles['v'] = lT.slice_container(data['v'], index={'time': [i]})
        profiles['width'] = lT.slice_container(data['width'], index={'time': [i]})
        profiles['beta'] = lT.slice_container(data['beta'], index={'time': [i]})
        
        if lidar_present:
            profiles['voldepol'] = lT.slice_container(data['voldepol'], index={'time': [i]})

        if doppler_present:
            it_b_dl = h.argnearest(data["v_lidar"]['ts'], data["cc"]["ts"][i]-15)
            it_e_dl = h.argnearest(data["v_lidar"]['ts'], data["cc"]["ts"][i]+15)
            if not it_b_dl == it_e_dl:
                profiles['v_lidar'] = lT.slice_container(data['v_lidar'], 
                        index={'time': [it_b_dl, it_e_dl]})
                profiles['a_lidar'] = lT.slice_container(data['a_lidar'], 
                        index={'time': [it_b_dl, it_e_dl]})
            else:
                print("no doppler lidar for this profile", i)

        if peakTree_present:
            it_b_pT = h.argnearest(data["pT_no"]['ts'], data["cc"]["ts"][i]-15)
            it_e_pT = h.argnearest(data["pT_no"]['ts'], data["cc"]["ts"][i]+15)
            if not it_b_pT == it_e_pT:
                profiles['pT_no'] = lT.slice_container(data['pT_no'], 
                        index={'time': [it_b_pT, it_e_pT]})
            else:
                print("peakTree not available for this profile", i)

        if qbsc_present:
            profiles['ratio_z_e'] = lT.slice_container(data['ratio_z_e'], index={'time': [i]})
            profiles['qbsc'] = lT.slice_container(data['qbsc'], index={'time': [i]})


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

        keys_to_feature = ["IWC", "LWC_S", "LWC", "Z", "v", "width", "T", "p", "SNR",
            "uwind", "vwind", "beta", "LDR", "LDRcorr"]
        if lidar_present:
            keys_to_feature += ["voldepol"]
        if qbsc_present:
            keys_to_feature += ['ratio_z_e']
            keys_to_feature += ['qbsc']
        #if peakTree_present:
        #    keys_to_feature += ['pT_no']

        features_in_profile = find_features_in_profile(profiles, keys_to_feature)
        
        # estimate this more naive classification
        if features_in_profile:
            simple_class.append(features_to_simple_class(features_in_profile))
            print('estimated class ', simple_class[-1])
        else:
            print('cloud free profile', data["cc"]["ts"][i])
            simple_class.append('clear')
        features_in_timestep.append(features_in_profile)


    #extract mixed-phase features and put them in a list
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
    clouds_mixed=connect_features(
        detected_features_mixed, cloud_type="layered", **CONN_THRES_LAYERED)
    #clouds_mixed = []
    print("Searching for cirrus (pure ice) clouds")
    #
    # 10000 as h_threshold seems a littlebit too much
    clouds_ice=connect_features(
        detected_features_ice, cloud_type="ice", **CONN_THRES_ICE)
    print("Searching for deep clouds")
    clouds_tower=connect_features(
        detected_features_tower, cloud_type="tower", **CONN_THRES_TOWER)
    #clouds_tower = []

    all_clouds = clouds_mixed + clouds_ice + clouds_tower
    print("Saving")
    save_object(all_clouds,
        '../cloud_properties/cloud_properties_'+campaign+'/'\
        +begin_dt.strftime("%Y_%m_%d")+'_clouds.dat')

    #save the profile wise classification
    class_counter = Counter(simple_class)
    profile_classes = ['clear', 'rain_deep', 'rain_shallow',
            'single_layer_pure_liquid', 'single_layer_pure_ice', 
            'single_layer_mixed-phase',
            'multi_layer_various', 'multi_layer_mostly_pure_liquid', 
            'multi_layer_mostly_pure_ice', 'multi_layer_mostly_mixed-phase' ]


    #print(simple_class)
    assert len(simple_class) == np.sum([class_counter[e] for e in profile_classes])
    outfile = '../cloud_collections/class_stat_'+campaign+'.csv'
    if not os.path.isfile(outfile):
        with open(outfile, 'w') as csvfile:
            w = csv.writer(csvfile, delimiter=',')
            w.writerow(['date', 'no_profiles']+profile_classes)

    with open(outfile, 'a') as csvfile:
        w = csv.writer(csvfile, delimiter=',')
        w.writerow([begin_dt.strftime("%Y%m%d"), len(simple_class)] \
            + [class_counter[e] for e in profile_classes])


    #plotting
    cloud_rectangles=[]
    for i, cloud in enumerate(all_clouds):
        
        if all_clouds[i].n_profiles()==0:
            continue
        
        c_type=cloud.most_common_type()
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

        print(f"cloud type {cloud.cloud_type:>10s} number of features {len(cloud.features):>5}  \
{h.ts_to_dt(cg[0]).strftime('%Y%m%d-%H%M')} + {cg[2]:7.2f} | {cg[1]:7.2f} + {cg[3]:7.2f}")
        #if clouds[i].top_variation()<200.0 and clouds[i].time_length()>1800 and clouds[i].cloud_top_thickness()[0]<400.0 and  clouds[i].fill_factor()>0.75:
        #if c_type=="tower":
        if cloud.time_length()>900 and cloud.fill_factor()>0.60:
            cloud_rectangles.append((cg[0],cg[1],cg[2],cg[3],color))

    fig, ax = lT.plot_timeheight(data['cc'], range_interval=[0, 12000])
    import matplotlib.patches as patches
    for cm in cloud_rectangles:
        print(f"plot {h.ts_to_dt(cm[0]).strftime('%Y%m%d-%H%M')} + {cm[2]:7.2f} | {cm[1]:7.2f} + {cm[3]:7.2f}")
        begin = h.ts_to_dt(cm[0])
        duration=datetime.timedelta(seconds=cm[2])
        rect = patches.Rectangle(
                (begin,cm[1]),duration,cm[3],linewidth=2,
                edgecolor=cm[4],facecolor=cm[4],alpha=0.2)

        # Add the patch to the Axes
        ax.add_patch(rect)

    savename='../plots/cloud_overview_'+campaign+'/'+begin_dt.strftime("%Y_%m_%d")+'_cloud_detection.png'
    fig.savefig(savename, dpi=250)

