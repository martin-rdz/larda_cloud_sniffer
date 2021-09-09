import datetime
import numpy as np
import pickle
from scipy import interpolate
import scipy.signal

from collections import Counter

import sys
sys.path.append('/home/larda3/larda/')
import pyLARDA.helpers as h

# flags for cloud particles
cloud_particles=[1,3,4,5,6,7]
liquid=[1,3,5,7]
droplets_only=[1,5]
ice=[4,5,6]
ice_only=4
melting=[6,7]
drizzle=[2,3]

class cloud_feature :
    
    def __init__(self):
        
        self.time=-1
        
        self.top=-1
        self.base=-1
        
        self.classifications=[]
        self.hasl=[]
        self.temperature=[]
        self.has_melting=False
        self.has_drizzle=False
        self.type="none"
        self.valid=False
        
        self.measurements={}
        self.measurements["IWC"]=[]
        self.measurements["LWC"]=[]
        self.measurements["Z"]=[]
        self.measurements["LDR"]=[]
        self.measurements["LDR_min"]=[]
        self.measurements["v"]=[]
        self.measurements["width"]=[]
        self.measurements["T"]=[]
        self.measurements["p"]=[]
        self.measurements["SNR"]=[]
        self.measurements["alpha_hogan"]=[]
        self.measurements["beta"]=[]
        self.measurements["voldepol"]=[]
        #self.measurements["tfv"]=[]
        #self.measurements["vair"]=[]

        self.liquid_layers=0
        self.liquid_layer_base=[]
        self.liquid_layer_top=[]
                
        self.precipitation_top=-1
        self.precipitation_base=-1
        
        self.cloud_system=-1
    
    def show(self):
        
        print("Top height:", self.top, "- Base height:", self.base, "- Type:", self.type, "- Time:", self.time)


def get_only_valid(data):
    return data['var'][~data['mask']]

def fill_with(data, fill):
    """fill an array where mask is true with fill value"""
    filled = data['var'].copy()
    filled[data['mask']] = fill
    return filled

def flatten(xs):
    """flatten inhomogeneous deep lists
       e.g. ``[[1,2,3],4,5,[6,[7,8],9],10]``
    """
    result = []
    if isinstance(xs, (list, tuple)):
        for x in xs:
            result.extend(flatten(x))
    else:
        result.append(xs)
    return result


def autocorr(x):
    result = np.correlate(x, x, mode='full')
    return result[result.size // 2:]



def time_analysis_from_vel(locations_of_vel, idx):
    """get the autocorrelation and the periodogram from the dl observations

    Args:
        locations_of_vel: list of tupels (ts, rg, vel)

    Returns:
        (f, Pxx_den), (time_shifts, v_autocorr)
    """

    locations_of_vel_s = sorted(locations_of_vel, key=lambda k: k[0])

    sep_time = np.array([e[0] for e in locations_of_vel_s])
    #sep_range = [e[1] for e in locations_of_vel_s]
    vel_line = np.array([e[idx] for e in locations_of_vel_s])

    #print(vel_line)

    delta_t = np.median(sep_time[1:] - sep_time[:-1])
    print('delta t', delta_t, vel_line[:10])

    if len(vel_line) > 0:
        f, Pxx_den = scipy.signal.welch(vel_line[vel_line != None],
                                        fs=1/delta_t)

        v_autocorr = autocorr(vel_line[vel_line != None])
        v_autocorr = v_autocorr/float(v_autocorr.max())
        time_shifts = np.arange(v_autocorr.shape[0])*delta_t
    else:
        f, Pxx_den = np.array([]), np.array([])
        time_shifts, v_autocorr = np.array([]), np.array([])

    return (f, Pxx_den), (time_shifts[:500], v_autocorr[:500])


class cloud():
    
    def __init__(self):
        
        self.cloud_type="none"
        
        self.begin_time=-1
        self.end_time=-1
                
        #self.top = -1
        #self.base = -1
        self.top_range = -1
        self.base_range = -1
        
        self.times=[]
        self.tops=[]
        self.bases=[]
        self.types=[]
        
        self.features=[]

        self.n_valid=0
        self.n_invalid=0
        
    def update_geometry(self):
        
        begin_t,base_h,duration_t,extent_h=self.geometry()
        self.begin_time=begin_t
        self.end_time=begin_t+duration_t
        self.base=base_h
        self.top=base_h+extent_h
        
    def append_feature(self,feature):
        
        self.features.append(feature)
        
        self.times.append(feature.time)
        self.tops.append(feature.top_range)
        self.bases.append(feature.base_range)
        self.types.append(feature.type)
        
        self.update_geometry()

    def validate_features(self):

        sep_height=[]

        for f in self.features:

            ll_base=-1
            if len(f.liquid_layer_base)>0:
                ll_base=f.liquid_layer_base[0]
            elif f.type=="pure liquid":
                ll_base=f.base
            sep_height.append(ll_base)
            
        med_sep_height=np.median(sep_height)

        print('validate list of liquid layer bases', med_sep_height, sep_height)
       
        avg_th, med_th, std_th, _ = self.cloud_top_thickness()

        valid,invalid=0,0

        for f in self.features:

            ll_base=-1
            if len(f.liquid_layer_base)>0:
                ll_base=f.liquid_layer_base[0]
            elif f.type=="pure liquid":
                ll_base=f.base

            print('time of feature ', h.ts_to_dt(f.time), ll_base)
            if np.abs(ll_base-med_sep_height)>150.0:
                f.valid=False
                invalid+=1
            else:
                f.valid=True
                valid+=1

        self.n_valid=valid
        self.n_invalid=invalid

    def n_profiles(self,types=[]):
        
        n=0
        if len(types)==0:
            n=len(self.features)
        else:
            for f in self.features:
                if f.type in types:
                    n+=1

        return n
        
    def geometry(self):
        """return the geometrical boundaries of the cloud""" 
        cloud_min_t = np.min(self.times)
        cloud_min_h = np.min(self.bases)
        cloud_width = self.time_length()
        cloud_height = np.max(self.tops)-np.min(self.bases)
        
        #print('current geometry ', datetime.datetime.utcfromtimestamp(cloud_min_t),
        #        " dur ", cloud_width, " range ", cloud_min_h, cloud_height)
        return cloud_min_t,cloud_min_h,cloud_width,cloud_height
    
    def most_common_type(self):
        types_count = Counter(self.types)
        #s_types=np.sort(self.types
        #return s_types[len(s_types)/2]
        return types_count.most_common(1)[0][0]
    
    def top_variation(self):
        """get the std from tops"""

        print(self.tops)
        
        return np.std(self.tops)

    def top_diffsum(self):

        return np.average(np.abs(np.diff(self.tops)))

    def length(self):
        
        return (max(self.times) - min(self.times)) * 10.0
    
    def time_length(self):
        
        return max(self.times) - min(self.times)
    
    def fill_factor(self):
        
        return self.n_profiles()*30.0/self.time_length()

    def return_values(self,name,particle_types=[]):

        values=[]
        for f in self.features:
            #if f.valid==False:
            #    continue
            if particle_types:
                cc_mask = np.isin(f.classifications, particle_types)
            else:
                cc_mask = np.full(f.classifications, True)
            var = f.measurements[name]['var'][cc_mask]
            mask = f.measurements[name]['mask'][cc_mask]
            values += var[~mask].tolist()

        values = np.array(values).ravel().astype(np.float64)
        assert np.all(np.isfinite(values)), values

        return values

    def return_values_separation(self,name,spacing):

        print("separation average", name, spacing)
        if not(name in self.features[0].measurements.keys()):
            print('name not in features')
            return np.array([])

        values=[]
        print('no_features', len(self.features))
        for f in self.features:


            #if f.valid==False:
            #    continue
            #print('precip top',f.precipitation_top, f.precipitation_top-spacing)
            #print('cc_profile', f.classifications )

            if f.precipitation_top!=-1 and f.precipitation_top-spacing>0:
                mask = f.measurements[name]['mask'][f.precipitation_top-spacing]
                values += f.measurements[name]['var'][f.precipitation_top-spacing][~mask].tolist()

        values=np.array(values).ravel()
        return values


    def combined_average(self,name,ref_name,particle_types=[]):

        if not(name in self.features[0].measurements.keys()) or not(ref_name in self.features[0].measurements.keys()):
            return 0,0,0,0

        values=[0]
        for f in self.features:

            #if f.valid==False:
            #    continue
            if not(ref_name in f.measurements.keys()) or len(f.measurements[name])!=len(f.measurements[ref_name]):
                continue

            for v in range(len(f.measurements[name])):


                if f.measurements[name][v]!=0 and f.measurements[ref_name][v]!=0:
                    value=f.measurements[name][v]
                else:
                    value=0 

                if len(particle_types)>0:
                    if f.classifications[v] in particle_types:
                        values.append(value)
                else:
                    values.append(value)

        values=np.array(values)

        if len(values)>0 and ~np.all(values==0):
            avg=np.average(values[values!=0])
            med=np.median(values[values!=0])
            std=np.std(values[values!=0])
            n=len(values[values!=0])
        else:
            avg=0.0
            med=0.0
            std=0.0
            n=0

        return avg,med,std,n


    def average(self,name,particle_types=[]):

        values=[]
        print('assemble', name, h.ts_to_dt(self.features[0].time))
        for f in self.features:
            #if f.valid==False:
            #    continue
            if particle_types:
                cc_mask = np.isin(f.classifications, particle_types)
            else:
                cc_mask = np.full(f.classifications, True)

            #print('classification ', f.classifications)
            #print('cc_mask ', cc_mask)
            if name in f.measurements.keys() and not type(f.measurements[name]) == list \
                    and any(cc_mask):
                var = f.measurements[name]['var'][cc_mask]
                mask = f.measurements[name]['mask'][cc_mask].astype(bool)
                #print('var', type(var), var.dtype, var)
                #print('mask', type(mask), mask.dtype, mask)
                if not np.all(mask) and np.all(~np.isnan(var)):
                    values += var[~mask].tolist()

        #print('average', name, values[:10], values[-10:])
        values = np.array(values).ravel().astype(np.float64)
        assert np.all(np.isfinite(values)), values
        if len(values)>0:
            avg=np.average(values)
            med=np.median(values)
            std=np.std(values)
            n=len(values)
        else:
            avg=0.0
            med=0.0
            std=0.0
            n=0
        print('assembled ', name, values[:10], values[-10:], ' avg ', avg)
        
        return avg,med,std,n
   
    def pp90(self,name,particle_types=[]):

        values=[]
        for f in self.features:
            #if f.valid==False:
            #    continue
            if particle_types:
                cc_mask = np.isin(f.classifications, particle_types)
            else:
                cc_mask = np.full(f.classifications, True)
            var = f.measurements[name]['var'][cc_mask]
            mask = f.measurements[name]['mask'][cc_mask]
            values += var[~mask].tolist()

        values = np.array(values).ravel().astype(np.float64)
        assert np.all(np.isfinite(values)), values

        if len(values)>1 and ~np.all(values==0):
            sorted_vals=np.sort(values[values!=0])
            ln=len(sorted_vals)

            return sorted_vals[int(0.1*ln)],sorted_vals[int(0.9*ln)]

        else:

            return 0,0

 
    def n_values(self,name,particle_types=[]):
        
        n_values=0
        for f in self.features:
            for v in range(len(f.measurements[name])):
                value=f.measurements[name][v]
                if len(particle_types)>0:
                    if (f.classifications[v] in particle_types) and value!=0:
                        n_values+=1
                else:
                    if value!=0:
                        n_values+=1
        
        return n_values

    def temperature_range(self):
        """loop over features, take the base and top temperature
        and return min(base), med(top), max(top)
        """ 
        top_temps=[]
        base_temps=[]

        for f in self.features:
            top_temps.append(f.measurements["T"]['var'][-1])
            base_temps.append(f.measurements["T"]['var'][0])
        
        return np.max(base_temps),np.median(top_temps),np.min(top_temps)

    def velocities_radar(self):

        v_top=[]
        v_top.append(0)

        width_top=[]
        width_top.append(0)

        Z_top=[]
        Z_top.append(0)

        for f in self.features:
            if f.valid==False:
                continue
            ll_base=-1
            if len(f.liquid_layer_base)>0:
                ll_base=f.liquid_layer_base[0]
            elif f.type=="pure liquid":
                ll_base=f.base

            if not f.measurements["Z"]["mask"][-1]:
                Z_top.append(f.measurements["Z"]['var'][-1])

            #for i in range(len(f.measurements["v"])):
            vr=f.measurements["v"]
            wt=f.measurements["width"]
            rg_valid = vr["rg"] >= ll_base
            v_top += flatten(fill_with(vr, 0)[rg_valid].tolist())
            width_top += flatten(fill_with(wt, 0)[rg_valid].tolist())

        v_top=np.array(v_top)
        v_top=v_top[v_top!=0]
        v_top_gt0=v_top[v_top>0.0]

        width_top=np.array(width_top)
        width_top=width_top[width_top!=0]

        Z_top=np.array(Z_top)
        Z_top=Z_top[Z_top!=0]

        if len(v_top)>0:
            v_mean=np.average(v_top)
            v_std=np.std(v_top)
            v_n=len(v_top)
        else:
            v_mean=0
            v_std=0
            v_n=0

        if len(v_top_gt0)>0:
            v_mean_gt0=np.average(v_top_gt0)
            v_std_gt0=np.std(v_top_gt0)
            v_n_gt0=len(v_top_gt0)
        else:
            v_mean_gt0=0
            v_std_gt0=0
            v_n_gt0=0

        if len(width_top)>0:
            width_avg=np.average(width_top)
        else:
            width_avg=0

        if len(Z_top)>1:
            z_top=np.median(Z_top)
        else:
            z_top=0

        return v_mean,v_std,v_n,v_mean_gt0,v_std_gt0,v_n_gt0,width_avg,v_top,z_top

    def velocities_liquid_radar(self, where):
        """
        retrun the velocities of the radar at the specified position of the liquid layer
        ``top``, ``whole``, ``base``
        """

        v_base=[]
        locations_of_vel = []

        for f in self.features:
            
            if f.valid==False or "v" not in f.measurements.keys() \
                or len(f.measurements['v']['rg'].shape) == 0:
                print('no measurements in this feature, valid? ', f.valid)
                continue

            ll_idx = None
            if len(f.liquid_layer_base)>0:
                #print('liquid layer bases',f.liquid_layer_base)
                #select highest liquid layer
                if where == 'base':
                    ll_idx = h.argnearest(f.measurements['v']['rg'], f.liquid_layer_base[-1])
                elif where == 'top-90':
                    ll_idx = h.argnearest(f.measurements['v']['rg'], f.liquid_layer_top[-1])-3
                elif where == 'top':
                    ll_idx = h.argnearest(f.measurements['v']['rg'], f.liquid_layer_top[-1])
                elif where == 'whole':
                    i_base = h.argnearest(f.measurements['v']['rg'], f.liquid_layer_base[-1])
                    i_top = h.argnearest(f.measurements['v']['rg'], f.liquid_layer_top[-1])
                    ll_idx = slice(i_base, i_top+1)

            elif f.type=="pure_liquid":
                if where == 'base':
                    ll_idx = h.argnearest(f.measurements['v']['rg'], f.base_range)
                elif where == 'top-90':
                    ll_idx = h.argnearest(f.measurements['v']['rg'], f.top_range)-3
                elif where == 'top':
                    ll_idx = h.argnearest(f.measurements['v']['rg'], f.top_range)
                elif where == 'whole':
                    i_base = h.argnearest(f.measurements['v']['rg'], f.base_range)
                    i_top = h.argnearest(f.measurements['v']['rg'], f.top_range)
                    ll_idx = slice(i_base, i_top+1)

            #print('velocities', f.measurements.keys())
            #print('ll_idx', ll_idx, f.type)
            #print('found indices ', ll_idx, f.liquid_layer_base, f.liquid_layer_top)
            if ll_idx is not None and len(f.measurements['v']['rg'].shape) > 0:
                v = f.measurements["v"]
                #v_lidar['mask'] = np.logical_or(v_lidar['mask'], a_lidar['var']<a_thr)
                v['mask'] = np.logical_or(v['mask'], v['var']<-990)

                if not np.all(v['mask'][ll_idx]):
                    v_base.append(v['var'][ll_idx][~v['mask'][ll_idx]].tolist())
                    locations_of_vel.append((v['ts'], v['rg'][ll_idx]))

        v_base=np.array(h.flatten(v_base))
        v_base=v_base[v_base!=0]

        if len(v_base)>0:
            v_mean=np.average(v_base)
            v_std=np.std(v_base)
            v_n=len(v_base)
        else:
            v_mean=0
            v_std=0
            v_n=0

        return v_mean,v_std,v_n,v_base,locations_of_vel

    def velocities_fixed_height(self, sample_height):
        """read the doppler lidar velocities at a fixed height of maximum backscatter
        [Suggestion by referee2]
        
        """

        a_thr=0
        #a_thr=8e4

        v_base=[]
        v_base.append(0)
        locations_of_vel = []

        for f in self.features:

            # 2020-10-12: try without the validity check
            # if f.valid==False or "v_lidar" not in f.measurements.keys():
            #     print('no measurements in this feature, valid? ', f.valid)
            #     continue

            
            if 'v_lidar' not in f.measurements:
                print('v_lidar missing at ', h.ts_to_dt(f.time))

            #print('here ', ll_base)
            
            if 'v_lidar' in f.measurements and sample_height >= 0 and len(f.measurements["v_lidar"]) > 0:
                v_lidar = f.measurements["v_lidar"]
                a_lidar = f.measurements["a_lidar"]

                v_lidar['mask'] = np.logical_or(v_lidar['mask'], a_lidar['var']<a_thr)
                a_lidar['mask'] = np.logical_or(a_lidar['mask'], a_lidar['var']<a_thr)

                #find base of liquid layer
                #print(v_lidar['rg'], v_lidar['rg'].shape)
                #print(v_lidar['var'], v_lidar['var'].shape)
                if v_lidar['var'].shape[1] > 1 and v_lidar['var'].shape[0] > 1:
                    mx_ind = h.argnearest(v_lidar['rg'], sample_height)
                else:
                    mx_ind = 0
                #print(ll_base, mx_ind)
                #print(len(a_lidar['var']), a_lidar['var'].shape)
                
                for it in range(a_lidar['var'].shape[0]):
                    #index of bsc max above liquid base
                    idx = np.argmax(fill_with(a_lidar, -99)[it,mx_ind:])
                    #print(it, idx, mx_ind+idx)
                    if not v_lidar['mask'][it, mx_ind+idx]:
                        v_base.append(v_lidar['var'][it, mx_ind+idx])
                        #print('v', v_lidar['var'][it, mx_ind+idx])
                        #print('a',a_lidar['var'][it, mx_ind+idx])
                        # sometimes rg is only a float
                        if isinstance(v_lidar['rg'], np.ndarray):
                            rg = v_lidar['rg'][mx_ind+idx]
                        else:
                            rg = v_lidar['rg']
                        locations_of_vel.append((v_lidar['ts'][it], rg, v_lidar['var'][it, mx_ind+idx], f.valid))
                    #print(v_lidar['ts'][it], v_lidar['rg'][mx_ind+idx])
                    #if v_lidar['var'].shape[1] > 1:
                        
                    #else:
                        #locations_of_vel.append((v_lidar['ts'][it], v_lidar['rg']))

                
        v_base=np.array(v_base)
        v_base=v_base[v_base != 0]
        v_base=v_base[v_base != None]
        print('v_base', v_base)

        if len(v_base)>0:
            v_mean=np.mean(v_base)
            v_std=np.std(v_base)
            v_n=len(v_base)
        else:
            v_mean=0
            v_std=0
            v_n=0

        return v_mean,v_std,v_n,v_base,locations_of_vel



    def velocities(self):
        """refactored
        reads out the doppler lidar velocities of each feature

        Returns:
            v_mean,v_std,v_n,v_base,locations_of_vel


        with locations_of_vel = v_lidar['ts'][it], v_lidar['rg'][mx_ind+idx], v_lidar['var'][it, mx_ind+idx], f.valid
        
        """
        a_thr=0
        #a_thr=8e4

        v_base=[]
        v_base.append(0)
        locations_of_vel = []

        for f in self.features:

            # 2020-10-12: try without the validity check
            # if f.valid==False or "v_lidar" not in f.measurements.keys():
            #     print('no measurements in this feature, valid? ', f.valid)
            #     continue

            ll_base=-1
            if len(f.liquid_layer_base)>0:
                #select highest liquid layer
                ll_base=f.liquid_layer_base[-1]
            elif f.type=="pure_liquid":
                ll_base=f.base_range
            #print('velocities', f.measurements.keys())
            #print('ll_base', ll_base, f.type)
            
            if 'v_lidar' not in f.measurements:
                print('v_lidar missing at ', h.ts_to_dt(f.time))

            #print('here ', ll_base)
            
            if 'v_lidar' in f.measurements and ll_base >= 0 and len(f.measurements["v_lidar"]) > 0:
                v_lidar = f.measurements["v_lidar"]
                a_lidar = f.measurements["a_lidar"]

                v_lidar['mask'] = np.logical_or(v_lidar['mask'], a_lidar['var']<a_thr)
                a_lidar['mask'] = np.logical_or(a_lidar['mask'], a_lidar['var']<a_thr)

                #find base of liquid layer
                #print(v_lidar['rg'], v_lidar['rg'].shape)
                #print(v_lidar['var'], v_lidar['var'].shape)
                if v_lidar['var'].shape[1] > 1 and v_lidar['var'].shape[0] > 1:
                    mx_ind = h.argnearest(v_lidar['rg'], ll_base)
                else:
                    mx_ind = 0
                #print(ll_base, mx_ind)
                #print(len(a_lidar['var']), a_lidar['var'].shape)
                
                for it in range(a_lidar['var'].shape[0]):
                    #index of bsc max above liquid base
                    idx = np.argmax(fill_with(a_lidar, -99)[it,mx_ind:])
                    #print(it, idx, mx_ind+idx)
                    if not v_lidar['mask'][it, mx_ind+idx]:
                        v_base.append(v_lidar['var'][it, mx_ind+idx])
                        #print('v', v_lidar['var'][it, mx_ind+idx])
                        #print('a',a_lidar['var'][it, mx_ind+idx])
                        # sometimes rg is only a float
                        if isinstance(v_lidar['rg'], np.ndarray):
                            rg = v_lidar['rg'][mx_ind+idx]
                        else:
                            rg = v_lidar['rg']
                        locations_of_vel.append((v_lidar['ts'][it], rg, v_lidar['var'][it, mx_ind+idx], f.valid))
                    #print(v_lidar['ts'][it], v_lidar['rg'][mx_ind+idx])
                    #if v_lidar['var'].shape[1] > 1:
                        
                    #else:
                        #locations_of_vel.append((v_lidar['ts'][it], v_lidar['rg']))

                
        v_base=np.array(v_base)
        v_base=v_base[v_base != 0]
        v_base=v_base[v_base != None]
        print('v_base', v_base)

        if len(v_base)>0:
            v_mean=np.mean(v_base)
            v_std=np.std(v_base)
            v_n=len(v_base)
        else:
            v_mean=0
            v_std=0
            v_n=0

        return v_mean,v_std,v_n,v_base,locations_of_vel

    def no_node_hist(self):
        """histogram over number of nodes for full cloud

        Returns:
            histogram of the node numbers
        """

        no_nodes=np.array([])

        for f in self.features:
            if 'pT_no' in f.measurements:
                var = f.measurements["pT_no"]['var'].ravel()
                no_nodes = np.append(no_nodes, var, axis=0)
        
        hist, bins = np.histogram(no_nodes, bins=[0,1,3,5,7,9,11,13,15,17,19,21,23])

        return hist.tolist()

    def no_node_hist_above_cb(self):
        """histogram over number of nodes for full cloud


        Returns:
            histogram of the node numbers
        """

        no_nodes=np.array([])

        for f in self.features:
            print(f.liquid_layer_base)
            if 'pT_no' in f.measurements:
                rg = f.measurements["pT_no"]['rg']
                print('pT range', rg)
                rg_lt_base =  np.where(rg > f.liquid_layer_base[0])[0]
                i_base = rg_lt_base[0]
                print(i_base)
                #print(f.measurements["pT_no"]['var'].shape)
                #print(f.measurements["pT_no"]['var'][:,i_base:].shape)
                var = f.measurements["pT_no"]['var'][:,i_base:].ravel()
                no_nodes = np.append(no_nodes, var, axis=0)
        
        hist, bins = np.histogram(no_nodes, bins=[0,1,3,5,7,9,11,13,15,17,19,21,23])

        return hist.tolist()


    def no_node_hist_ice_liq(self):
        """histogram over number of nodes for full cloud


        Returns:
            histogram of the node numbers
        """

        no_nodes=np.array([])

        for f in self.features:
            print(f.classifications)
            print(f.ranges)
            if 'pT_no' in f.measurements:
                rg = f.measurements["pT_no"]['rg']
                print('pT range', rg)
                lowest_ice = np.where(np.isin(f.classifications, [1,4,5]))[0]
                print('lowest ice ', lowest_ice)
                if len(lowest_ice) > 0:
                    rg_lowest = f.ranges[lowest_ice][0]
                    rg_lt_base =  np.where(rg > rg_lowest)[0]
                    i_base = rg_lt_base[0]
                    print(i_base)
                    var = f.measurements["pT_no"]['var'][:,i_base:].ravel()
                    no_nodes = np.append(no_nodes, var, axis=0)
        
        hist, bins = np.histogram(no_nodes, bins=[0,1,3,5,7,9,11,13,15,17,19,21,23])
        return hist.tolist()


    def horizontal_wind(self, cth, h_range):

        N_features=len(self.features)
        f=self.features[N_features/2]

        wp_hasl=f.measurements["wp_hasl"]
        wp_vel=f.measurements["wp_vel"]

        cth_index = np.argmin(np.abs(wp_hasl-cth))
        
        top = cth_index+h_range
        bottom = cth_index-h_range

        if top>=len(wp_vel):
            top=len(wp_vel)-1

        if bottom<0:
            bottom=0

        #Advection profile around cloud top
        vel_profile=np.array(wp_vel[bottom:top])
        hasl_profile=np.array(wp_hasl[bottom:top])

        clean_index=(~np.isnan(vel_profile))*(vel_profile>0.0)*(vel_profile<100.0)
        
        vel_profile=vel_profile[clean_index]
        hasl_profile=hasl_profile[clean_index]

        if len(vel_profile)>0:
            
            std=np.std(vel_profile)
            avg=np.average(vel_profile)
            mx=np.max(vel_profile)
            mn=np.min(vel_profile)
            dvdh=np.average(np.diff(vel_profile)/np.diff(hasl_profile))

        else:
           
            std,avg,mx,mn,dvdh = 0,0,0,0,0

        return avg, std, mx, mn, dvdh

            
    def ilr(self, spacing):

        iwc_top=[]
        lwc_top=[]

        for f in self.features:
            if f.precipitation_top!=-1 and f.precipitation_top-spacing>0:
                iwc_top.append(f.measurements["IWC"]['var'][f.precipitation_top-spacing])
                if np.sum(f.measurements["LWC"]['var'])>0:
                    #print('lwc', f.measurements["LWC"]['var'], f.measurements["LWC"]['mask'])
                    lwc_top.append(np.average(f.measurements["LWC"]['var'][f.measurements["LWC"]['var']!=0]))
                else:
                    lwc_top.append(0)

        iwc_top=np.array(iwc_top)
        lwc_top=np.array(lwc_top)


        if len(lwc_top!=0)>0:
            ilcr=iwc_top[lwc_top!=0]/lwc_top[lwc_top!=0]
        else:
            ilcr=[0]

        print('at ilr; iwc, lwc, ilcr average', 
              np.average(iwc_top), np.average(lwc_top), 
              np.average(ilcr))

        #n_ice=len(ilcr>0.9)/float(len(ilcr))
        #n_liq=len(ilcr<0.1)/float(len(ilcr))
        #n_mix=(len(ilcr)-n_ice-n_liq)/float(len(ilcr)


        bins = np.logspace(-5, 0, 20)
        if len(ilcr)>0:
            histogr=np.histogram(ilcr,bins=bins)[0]
        else:
            histogr=np.repeat(0,len(bins))

        return np.average(ilcr),np.median(ilcr),len(ilcr),list(histogr)

    def separation_location(self, spacing):

        times=[]
        ranges=[]
        for f in self.features:
            if f.precipitation_top!=-1 and f.precipitation_top-spacing>0:
                ranges.append(f.ranges[f.precipitation_top-spacing])
                times.append(f.time)

        return times, ranges


    def separation_average(self,name,spacing):

        print("separation average", name, spacing)
        if not(name in self.features[0].measurements.keys()):
            print('name not in features')
            return 0,0,0,0

        values=[]
        print('no_features', len(self.features))
        for f in self.features:


            #if f.valid==False:
            #    continue
            #print('precip top',f.precipitation_top, f.precipitation_top-spacing)
            #print('cc_profile', f.classifications )

            if f.precipitation_top!=-1 and f.precipitation_top-spacing>0:
                mask = f.measurements[name]['mask'][f.precipitation_top-spacing].astype(bool)
                #print(type(mask), mask.dtype, mask.data)
                values += f.measurements[name]['var'][f.precipitation_top-spacing][~mask].tolist()


        values=np.array(values).ravel()

        #print(values)
        if len(values)>0 and ~np.all(values==0):
            values_avg=np.average(values[values!=0])
            values_med=np.median(values[values!=0])
            values_std=np.std(values[values!=0])
            n=len(values[values!=0])
        else:
            values_avg=0
            values_med=0
            values_std=0
            n=0

        return values_avg,values_med,values_std,n


    def separation_values(self,name,spacing):

        print("separation values", name, spacing)
        if not(name in self.features[0].measurements.keys()):
            print('name not in features')
            return 0,0,0,0

        values=[]
        print('no_features', len(self.features))
        for f in self.features:


            #if f.valid==False:
            #    continue
            print('precip top',f.precipitation_top, f.precipitation_top-spacing)
            print('cc_profile', f.classifications )

            if f.precipitation_top!=-1 and f.precipitation_top-spacing>0:
                mask = f.measurements[name]['mask'][f.precipitation_top-spacing]
                values += f.measurements[name]['var'][f.precipitation_top-spacing][~mask].tolist()


        values=np.array(values).ravel()

        return values


    def liquid_layer_variation(self):

        cb=[]
        for f in self.features:
            if f.liquid_layer_base!=[]:
                cb.append(f.liquid_layer_base[0])
            else:
                cb.append(f.base)

        return np.std(cb)

    def cloud_top_thickness(self):
        """thickness of the liquid layer
        
        if a liquid base is detected?

        ! change mr:
        omit non liquid thickness

        Returns:
            np.average(thickness), np.median(thickness), np.std(thickness), thickness_with_time


        where (f.time, f.top_range-f.liquid_layer_base[0], f.liquid_layer_base[0], f.top_range, flag)

        flag 0: feature.liquid_layer_base[0]
        
        flag 1: feature.base_range [omitted as base_range is frequently the base of the ice]
        """

        thickness_with_time = []
        for f in self.features:

            print(f.time, f.top_range, f.liquid_layer_base, f.base_range)
            # 2020-10-12 f.base_range and liquid_layer_base are not equal???
            #

            if f.liquid_layer_base!=[]:
                thickness_with_time.append((f.time, f.top_range-f.liquid_layer_base[0], f.liquid_layer_base[0], f.top_range, 0))
            else:
                pass
                #
                #thickness_with_time.append((f.time, f.top_range-f.base_range, f.base_range, f.top_range, 1))

        print('cloud top thickness: no tops ', len(self.tops), ' no thickness ',  len(thickness_with_time))
        thickness = [e[1] for e in thickness_with_time]
        return np.average(thickness), np.median(thickness), np.std(thickness), thickness_with_time


    def cloud_top_avg(self, frac=0.5):
        """calculated the mean height of the mid of the liquid layer
        (based on cloud_top_thickness())
        
        Returns:
            np.average(thickness), np.median(thickness), np.std(thickness), thickness_with_time

        """

        mid_with_time = []
        for f in self.features:

            print(f.time, f.top_range, f.liquid_layer_base, f.base_range)
            # 2020-10-12 f.base_range and liquid_layer_base are not equal???
            #

            if f.liquid_layer_base!=[]:
                mid = f.liquid_layer_base[0] + frac*(f.top_range-f.liquid_layer_base[0])
                mid_with_time.append((f.time, mid))

        mids = [e[1] for e in mid_with_time]
        return np.average(mids), mid_with_time



    def average_paths(self):
    
        lwp=[]
        lwp_s=[]
        iwp=[]

        #dh=self.features[0].ranges[1]-self.features[0].ranges[0]
        #dh = self.features[0].measurements["LWC"]['rg'][1] - self.features[0].measurements["LWC"]['rg'][0]
        dh = self.features[0].dh
        for f in self.features:
            lwp.append(np.sum(get_only_valid(f.measurements["LWC"]))*dh)
            lwp_s.append(np.sum(get_only_valid(f.measurements["LWC_S"]))*dh)
            iwp.append(np.sum(get_only_valid(f.measurements["IWC"]))*dh)

        #ATTENTION!

        lwp=np.array(lwp)
        lwp_s=np.array(lwp_s)
        iwp=np.array(iwp)
        print('paths lwp, lwp_s', lwp, lwp_s)

        if lwp!=[] and ~np.all(lwp==0) and len(lwp[lwp_s!=0])>0:
            lwp_average=np.average(lwp[lwp_s!=0])
        else:
            lwp_average=0

        if lwp_s!=[] and ~np.all(lwp_s==0):
            lwp_s_average=np.average(lwp_s[lwp_s!=0])
        else:
            lwp_s_average=0

        if iwp!=[] and ~np.all(iwp==0):
            iwp_average=np.average(iwp[iwp!=0])
            iwp_std=np.std(iwp[iwp!=0])
        else:
            iwp_average=0
            iwp_std=0

        print('path averages ', lwp_average, lwp_s_average, iwp_average)

        return lwp_average, lwp_s_average, iwp_average, iwp_std

    def pressure_range(self):
        """loop over features, take the base and top pressure
        and return max(base), min(top)
        """ 
        top_pressure=[]
        base_pressure=[]

        for f in self.features:
            top_pressure.append(f.measurements["p"]['var'][-1])
            base_pressure.append(f.measurements["p"]['var'][0])

        return np.max(base_pressure),np.min(top_pressure)

    def print_values(self, name1, name2):
        for f in self.features:
            for v in range(len(f.measurements[name1])):
                if f.measurements[name1][v]!=0 and f.measurements[name2][v]!=0:
                    print(f.measurements[name1][v], f.measurements[name2][v])

    def n_melting(self):

        n_melt=0

        for f in self.features:
            if f.has_melting:
                n_melt+=1

        return n_melt

    def n_drizzle(self):

        n_driz=0

        for f in self.features:
            if f.has_drizzle:
                n_driz+=1

        return n_driz

    def correct_LDR(self,snr_co=-20,ldr_limit=-33):
        # isn't that done during reading?
        pass
        #for f in self.features:

        #    LDR_min=snr_co-f.measurements["SNR"]
        #    f.measurements["LDR"][f.measurements["LDR"]<ldr_limit]=0
        #    f.measurements["LDR"][f.measurements["LDR"]<LDR_min]=0
    
    def FA(self, spacing=3):

        values=[]
        times=[]

        with open("FA/FA.txt","a+") as of:
            for f in self.features:

                target_index=f.precipitation_top-spacing
                if target_index<0:
                    target_index=0
                    
                Zval=f.measurements["Z"][target_index]
                Ival=f.measurements["IWC"][target_index]
                time=f.time
                temp=f.measurements["T"][-1]
                of.write(str(time)+" "+str(temp)+" "+str(Zval)+" "+str(Ival)+"\n")

            of.write("\n")

    def vel_and_dir(self):
        """calculate velocity and direction from uwind and vwind
        and return the median values"""
        vels=[]
        dirs=[]
        for i in range(len(self.features)):

            u=self.features[i].measurements["uwind"]['var'][-1]
            v=self.features[i].measurements["vwind"]['var'][-1]

            vec=u-1.0j*v
            dir=np.angle(vec,deg=True)
            dir=np.mod(dir+270,360)
            vel=np.abs(vec)

            vels.append(vel)
            dirs.append(dir)

        return np.median(vels), np.median(dirs)


    def get_class_file_history(self):
       
        file_hist_all_features = [f.classification_file_history for f in self.features]
        assert len(file_hist_all_features[0]) == 1, "too many file_hist elements in list"
        processing_date = [datetime.datetime.strptime(file_hist[0][:20], "%d %b %Y %H:%M:%S")\
                for file_hist in file_hist_all_features]

        assert len(set(processing_date)) == 1, 'more than one processing date found for this day'
        #05 Apr 2019 18:11:59

        return processing_date[0].strftime("%Y%m%d_%H%M%S")
