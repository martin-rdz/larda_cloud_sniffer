#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')

import numpy as np
from scipy import interpolate
import scipy.special as sps
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap
#from IO_Functions import *
import pylab
import datetime
#import VDC_ParticleTable
import argparse
import copy
#from rayleigh_bg import *
#import dream

#Matplotlib Configuration
mpl.rcParams['xtick.major.size']=12.0 
mpl.rcParams['xtick.minor.size']=12.0
mpl.rcParams['ytick.major.size']=12.0
mpl.rcParams['ytick.minor.size']=6.0
mpl.rcParams['ytick.major.width']=1.5
mpl.rcParams['xtick.major.width']=1.5
mpl.rcParams['axes.linewidth']=2.0
mpl.rcParams['font.sans-serif'] = 'arial'
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.size'] = 24.0
mpl.rcParams['legend.fontsize'] = 18.0

mpl.rcParams['mathtext.fontset']='custom'
mpl.rcParams['mathtext.rm']='arial'
mpl.rcParams['mathtext.sf']='arial'

#parser = argparse.ArgumentParser(description='Process some integers.')
#parser.add_argument('--station', metavar='station_str', type=str, help='coudnet station for which the algorithm should be run')
#args = parser.parse_args()

ldr_colors=(
#(0.6 , 0.6 , 0.6),
#(0.6 , 0.6 , 0.6),
(0.6 , 0.6 , 0.6),

(0.0 , 0.0 , 0.5),
(0.0 , 0.0 , 0.5),
(0.0 , 0.0 , 0.5),
(0.0 , 0.0 , 0.5),

(0.0 , 0.7 , 1.0),
(0.0 , 0.7 , 1.0),
(0.0 , 0.7 , 1.0),

(0.0 , 0.9 , 0.0),
(0.0 , 0.9 , 0.0),
(0.0 , 0.9 , 0.0),
(0.0 , 0.9 , 0.0),

(1.0 , 0.8 , 0.0),
(1.0 , 0.8 , 0.0),
(1.0 , 0.8 , 0.0),
(1.0 , 0.8 , 0.0),

(1.0 , 0.0 , 0.0),
(1.0 , 0.0 , 0.0),
(1.0 , 0.0 , 0.0),
(1.0 , 0.0 , 0.0),

(0.8 , 0.8 , 0.8),
(0.8 , 0.8 , 0.8),
(0.8 , 0.8 , 0.8),
(0.8 , 0.8 , 0.8),
(0.8 , 0.8 , 0.8),
(0.8 , 0.8 , 0.8),
(0.8 , 0.8 , 0.8),
(0.8 , 0.8 , 0.8),
(0.8 , 0.8 , 0.8),
(0.8 , 0.8 , 0.8),
(0.8 , 0.8 , 0.8),
(0.8 , 0.8 , 0.8),

(1.0 , 1.0 , 1.0),
(1.0 , 1.0 , 1.0),
(1.0 , 1.0 , 1.0),
(1.0 , 1.0 , 1.0),
(1.0 , 1.0 , 1.0),
)

ldrscale=np.array([-36,-33,-32,-30,-29,-24,-21,-12,0])
ldrscale=1.0-ldrscale/(-36.0)
s=ldrscale
s[0]=0.0
s[8]=1.0

ldr_colors_scaled = {
	 'red':   ((s[0], 0.6, 0.6),
                   (s[1], 0.6, 0.0),
                   (s[2], 0.0, 0.0),
                   (s[3], 0.0, 0.0),#(s[3], 0.3, 0.3),
                   (s[4], 0.0, 0.0),
                   (s[5], 0.0, 0.0),
                   (s[6], 1.0, 1.0),
                   (s[7], 1.0, 0.6),
                   (s[8], 0.6, 0.6)),

         'green': ((s[0], 0.6, 0.6),
                   (s[1], 0.6, 0.0),
                   (s[2], 0.0, 0.0),
                   (s[3], 0.0, 0.0),#(s[3], 0.3, 0.3),
                   (s[4], 0.5, 0.5),
                   (s[5], 1.0, 1.0),
                   (s[6], 1.0, 1.0),
                   (s[7], 0.0, 0.6),
                   (s[8], 0.6, 0.6)),

         'blue':  ((s[0], 0.6, 0.6),
                   (s[1], 0.6, 0.0),
                   (s[2], 0.0, 0.0),
                   (s[3], 0.0, 0.0),#(s[3], 1.0, 1.0),
                   (s[4], 1.0, 1.0),
                   (s[5], 0.0, 0.0),
                   (s[6], 0.0, 0.0),
                   (s[7], 0.0, 0.6),
                   (s[8], 0.6, 0.6)),
        }

def magnus_i(T):
    return 61120.0*np.exp(22.46*T/(272.62+T))

def magnus_w(T):
    return 61120.0*np.exp(17.62*T/(243.12+T))

def wv(T):
    return -1.494e-04+5.025e-03*np.exp(T/1.436e+01)

def hogan_alpha(Z,T_C):

    return 10**(0.000447*Z*T_C+0.0683*Z-0.0171*T_C-3.11)
    
def hogan_m(Z,T_C):

    return 10**(0.000242*Z*T_C+0.0699*Z-0.0186*T_C-1.63) / 1000.0
    
def m_i_correction_factor(ZDR,T):

    Z_red=10*np.log10((10**(ZDR/10.0)+1.0)/2.0)

    return 10**(-0.000242*Z_red*T-0.0699*Z_red)

def IWC_Correction(IWC_raw,T_C):
   
    IWC_new=copy.deepcopy(IWC_raw)

    cf=-1

    for i in range(len(IWC_raw)):
        if T_C[i]<-22:
            cf=m_i_correction_factor(3,T_C[i]) 
        elif -16<T_C[i]<-12:
            cf=m_i_correction_factor(6,T_C[i])
        elif -6<T_C[i]<-2:
            cf=m_i_correction_factor(4,T_C[i])
        else:
            cf=1

        IWC_new[i]=IWC_new[i]*cf
    
    return IWC_new


ldr_colormap = ListedColormap(ldr_colors, "ldr_colormap")
#ldr_colormap = LinearSegmentedColormap("ldr_colormap",ldr_colors_scaled)
plt.register_cmap(cmap=ldr_colormap)

#Constants
R_air=287.058 #J/(kg*K)
g=9.81 #kg/(m*s^2)
eta_air=1.55e-5 #kg/(m*s)
delta_0=8.0
C_0=0.35

def plot_histogram(savename, data, plot_range, n_bins, x_label, y_label, normed=True):
    
    plt.hist(x=data, bins=n_bins, range=plot_range, normed=normed)
    
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    
    plt.savefig(savename, bbox_inches='tight')
    plt.clf()

def plot_color_scatter(savename,x,y,z,x_label,y_label,z_label,xlog=False, ylog=False, zlog=False, xlim=(np.NAN,np.NAN), ylim=(np.NAN,np.NAN), zlim=(np.NAN,np.NAN), cmap_name="spectral", dates=[], dpi=150, x_delta=2.5, meanplot=False):

    x=np.array(x)
    y=np.array(y)
    z=np.array(z)

    x_meanplot=np.arange(-40,10,x_delta)
    y_meanplot=[]
    y_stdplot=[]
  
    for xi in x_meanplot:
    
        index=(x>xi-x_delta/2.0)*(x<xi+x_delta/2.0)
        y_values_in_range=y[index]
        
        mean=np.NAN
        std=np.NAN

        if ylog:

            if len(y_values_in_range)>0:
                
                valid_index=(~np.isnan(y_values_in_range))*(~np.isinf(y_values_in_range))*(y_values_in_range>0)
            
                mean=np.mean(np.log10(y_values_in_range[valid_index]))
                std=np.std(np.log10(y_values_in_range[valid_index]))
                
            y_meanplot.append(10**mean)
            y_stdplot.append(10**std)
            
        else:
            
            if len(y_values_in_range)>0:
                
                valid_index=(~np.isnan(y_values_in_range))*(~np.isinf(y_values_in_range))
            
                mean=np.mean(y_values_in_range[valid_index])
                std=np.std(y_values_in_range[valid_index])
                
            y_meanplot.append(mean)
            y_stdplot.append(std)

    y_meanplot=np.array(y_meanplot)
    y_stdplot=np.array(y_stdplot)

    if ylog:
        y_error_low=(y_meanplot-y_meanplot/y_stdplot)
        y_error_high=(y_meanplot*y_stdplot-y_meanplot)
    else:
        y_error_low=y_stdplot
        y_error_high=y_stdplot

    fig, ax = plt.subplots(facecolor='white')
    
    if cmap_name=="ldr_colormap" and len([z==0])>0:
        s=ax.scatter(x[z!=0], y[z!=0] , c=z[z!=0], cmap=mpl.cm.get_cmap(cmap_name), s=175, marker='o', lw=0.0, alpha=0.8, vmin=zlim[0], vmax=zlim[1])
        s1=ax.scatter(x[z==0], y[z==0], s=175, marker='o', lw=0.5, alpha=0.8, facecolor='none')
        #alpha=0.5
    else:
        s=ax.scatter(x, y , c=z, cmap=mpl.cm.get_cmap(cmap_name), s=75, marker='o', lw=0.5, alpha=0.75, vmin=zlim[0], vmax=zlim[1])
        
        #e=ax.errorbar(x,y,xerr=None, yerr=0.5*y, mew=0, marker=None, fmt=None, lw=0.5, zorder=0)
    
    if meanplot:
        s_mean=ax.scatter(x_meanplot, y_meanplot, s=150, marker='s', zorder=100, facecolor='white', lw=2.0, alpha=0.75)
        s_std=ax.errorbar(x_meanplot, y_meanplot, xerr=None, yerr=(y_error_low,y_error_high), mew=0, marker=None, fmt=None, lw=2.0, ecolor='black', zorder=99)
    
    if xlog:
        plt.xscale('log')
    if ylog:
        plt.yscale('log')
        
    plt.xlim(xlim[0],xlim[1])
    plt.ylim(ylim[0],ylim[1])
    
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    
    if dates!=[]:
        for i in range(len(x)):
            if ~np.isnan(x[i]) and ~np.isnan(y[i]):
                ax.annotate(dates[i].strftime("%y%m%d%H"),xy=(x[i],y[i]),fontsize=2)
    
    cbar = plt.colorbar(mappable=s, ax=ax)
    if cmap_name=="ldr_colormap":
        cbar.ax.set_yticklabels(['-32','-28','-24','-20','-16','-12','-8','  ','low\nsig.'])
    cbar.set_label(z_label)

    plt.savefig(savename, bbox_inches='tight', dpi=dpi)
    plt.clf()

    return x_meanplot,y_meanplot,y_stdplot

def plot2D(savename, x, y, x_label,y_label, xerr=[], yerr=[] ,label="",xlog=False, ylog=False, xlim=(np.NAN,np.NAN), ylim=(np.NAN,np.NAN), x_delta=1.0, fmt='.'):
  

    if xerr!=[] and yerr!=[]:
        plt.errorbar(x,y,xerr=xerr, yerr=yerr, label=label, fmt=fmt, ecolor='g')
    else:
        plt.plot(x,y,fmt,label=label)
    
    if xlog:
        plt.xscale('log')
    if ylog:
        plt.yscale('log')
        
    plt.xlim(xlim[0],xlim[1])      
    plt.ylim(ylim[0],ylim[1])
        
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    
    plt.savefig(savename, bbox_inches='tight')
    
    plt.clf()

def multiple_plot2D(savename, x, y, xerr, yerr, x_label,y_label,labels=[],xlog=False, ylog=False, xlim=(np.NAN,np.NAN), ylim=(np.NAN,np.NAN)):

    for i in range(len(y)):
        plt.errorbar(x[i],y[i],xerr=xerr[i], yerr=yerr[i], label=labels[i])
      
    if xlog:
        plt.xscale('log')
    if ylog:
        plt.yscale('log')
        
    plt.xlim(xlim[0],xlim[1])      
    plt.ylim(ylim[0],ylim[1])
        
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    
    plt.legend(loc=3)
    
    plt.savefig(savename, bbox_inches='tight')
    
    plt.clf()

def plot_2D_histogram(savename, x, y, Nx, Ny, rangex, rangey, x_label, y_label, norm="", display_y=(np.NAN,np.NAN)):

    fig, ax = plt.subplots(1, figsize=(12, 8))
    hist2d=np.histogram2d(x,y, (Nx,Ny), (rangex,rangey))[0]

    if display_y==(np.NAN,np.NAN):
        display_y=rangey

    for i in range(len(hist2d)):
        hist2d[i,np.argmax(hist2d[i,:])]=0

    if norm=='column':
        for n in range(np.shape(hist2d)[0]):
            hist2d[n,:]/=np.sum(hist2d[n,:])
    elif norm=='row':
        for n in range(np.shape(hist2d)[1]):
            hist2d[:,n]/=np.sum(hist2d[:,n])

    hist2d=np.rot90(np.log10(hist2d),1)

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xlim(rangex[0],rangex[1])
    plt.ylim(display_y[0],display_y[1])

    img=plt.imshow(hist2d,extent=(rangex[0],rangex[1],rangey[1],rangey[0]),aspect='auto', origin='lowerleft', interpolation='nearest', cmap=mpl.cm.get_cmap("spectral"))
    cbar = fig.colorbar(img)

    cbar.set_label("Relative Frequency (log)")

    plt.savefig(savename, bbox_inches='tight')

    plt.clf()

    
def plot_picture(savename, data, x_range, y_range, z_range, x_label, y_label, z_label):

    fig, ax = plt.subplots(1, figsize=(12, 4))
    pcmesh = ax.pcolormesh(np.arange(x_range[0],x_range[1]), 
          np.arange(y_range[0],y_range[1]), 
          np.transpose(data),
          vmin=-3, vmax=3)
    cbar = fig.colorbar(pcmesh)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    cbar.ax.set_ylabel(z_label)
    fig.savefig(savename)
    plt.close(fig)

def estimate_N(Z,v_t,p_Pa,T_K,p_type):

    k=3
    particle_props=p_table[p_type]
    c=particle_props[3]
    d=particle_props[4]
    D_max,Re=v_D_inverter(p_Pa,T_K,v_t,particle_props)
    D0=D_max/(k+2*d)
    norm_Z = 0.174/0.93 * (6/(np.pi*917.0))**2 * c**2 * D0**(2*d) * sps.gamma(k+2*d)/sps.gamma(k) * (100**d/1000.0)**2
    
    return (10**(Z/10.0)/(norm_Z/0.001**6)), Re

def terminal_fall_velocity (particle_props, Diameter_m, p_Pa, T_K):

    #Input
    p=p_Pa #Pa
    T=T_K #K
    D=Diameter_m #m
    
    #Paremtrization Selection        
    alpha=particle_props[3]
    beta=particle_props[4]
    gamma=particle_props[5]
    sigma=particle_props[6]
    
    #parameters
    rho_air=p/(R_air*T)
    m=alpha * (D*100)**beta / 1000.0 #conversion cgs-SI
    A=gamma * (D*100)**sigma/ 10000.0 #conversion cgs-SI
    A_r=A/(np.pi*(D/2)**2)
    
    #Modified Best-Number
    Xm=rho_air * 8.0 * m * g / (eta_air**2* np.pi * A_r**0.5) 
    #Reynolds-Number
    Re=0.25 * delta_0**2 *  ((1.0 + 4.0 * Xm**0.5 / (delta_0**2 *C_0**0.5))**0.5 - 1.0)**2
    #Terminal Fall Velocity
    v_t= eta_air * Re / (rho_air * D)
    
    return v_t, m, A, Re

def v_D_inverter(p_Pa,T_K,v,particle_props):

    v_ar=[]
    D_ar=[]
    R_ar=[]
    
    for d_max in np.arange(100e-6,5000e-6,100e-6):
        D_ar.append(d_max)
        v_t,m,A,Re=terminal_fall_velocity(particle_props,d_max,p_Pa,T_K)
        v_ar.append(v_t)
        R_ar.append(Re)
    
    D_v=interpolate.interp1d(v_ar,D_ar,bounds_error=False, fill_value=0)
    R_v=interpolate.interp1d(v_ar,R_ar,bounds_error=False, fill_value=0)
    
    return D_v(v), R_v(v)

exit()
p_table=VDC_ParticleTable.get_particle_table()

station_name=args.station

print('cloud_collections/cloud_collection_'+station_name+'_all.csv')
datafile=load_csv('cloud_collections/cloud_collection_'+station_name+'_all.csv')

folder=station_name+'/'

station_lat={}
station_lon={}
station_alt={}
dust_path={}

station_lat["LEIPZIG"]=51.340333
station_lon["LEIPZIG"]=12.37475
station_alt["LEIPZIG"]=117.0
dust_path["LEIPZIG"]="dream/dream-v8b-dust_profiles-leipzig.nc"

station_lat["ICON_HDCP2"]=51.340333
station_lon["ICON_HDCP2"]=12.37475
station_alt["ICON_HDCP2"]=117.0
dust_path["ICON_HDCP2"]=None

station_lat["LINDENBERG"]=52.209444
station_lon["LINDENBERG"]=14.128333
station_alt["LINDENBERG"]=103.0
dust_path["LINDENBERG"]="dream/dream-v8b-dust_profiles-lindenberg.nc"

station_lat["MACEHEAD"]=53.333
station_lon["MACEHEAD"]=-9.9
station_alt["MACEHEAD"]=5.0
dust_path["MACEHEAD"]="dream/dream-v8b-dust_profiles-macehead.nc"

station_lat["POTENZA"]=40.60
station_lon["POTENZA"]=15.72
station_alt["POTENZA"]=760.0
dust_path["POTENZA"]="dream/dream-v8b-dust_profiles-potenza.nc"

station_lat["CABAUW"]=51.964444
station_lon["CABAUW"]=4.898333
station_alt["CABAUW"]=1.0
dust_path["CABAUW"]="dream/dream-v8b-dust_profiles-cabauw.nc"

station_lat["KRAUTHAUSEN"]=50.872194
station_lon["KRAUTHAUSEN"]=6.423852
station_alt["KRAUTHAUSEN"]=230
dust_path["KRAUTHAUSEN"]=None

station_lat["BARBADOS"]=-59.5355639
station_lon["BARBADOS"]=13.1901325
station_alt["BARBADOS"]=100
dust_path["BARBADOS"]=None

station_lat["CHILBOLTON"]=51.144593
station_lon["CHILBOLTON"]=-1.438779
station_alt["CHILBOLTON"]=17
dust_path["CHILBOLTON"]=None

station_lat["PUNTAARENAS"]=-70
station_lon["PUNTAARENAS"]=-50
station_alt["PUNTAARENAS"]=10
dust_path["PUNTAARENAS"]=None

lat=station_lat[station_name]
lon=station_lon[station_name]
alt=station_alt[station_name]
        
season={}
season["(1) spring (MAM)"]=[3,4,5]
season["(2) summer (JJA)"]=[6,7,8]
season["(3) autumn (SON)"]=[9,10,11]
season["(4) winter (DJF)"]=[12,1,2]
season["total_year"]=range(0,13)

true_liquid=0
true_mixed=0
false_liquid=0
false_mixed=0
total_SSC=0
mixed_listed=0
liquid_listed=0

LDR_scale_min=-34
LDR_scale_max=2

n_ILCR_t=np.zeros(10)

#dust_model=dream.dream(dust_path[station_name])

labels=[]

x_temp=[]
y_temp=[]
xerr_temp=[]
yerr_temp=[]

x_hgts=[]
y_hgts=[]
xerr_hgts=[]
yerr_hgts=[]

ice_positive_list=load_dict("ice_positive_list.dat")
ice_negative_list=load_dict("ice_negative_list.dat")

#target_year=2015

for season_name in sort(season.keys()):    

    ui=[]
    CTT=[]
    CTH=[]
    CBT=[]
    CTP=[]
    CBP=[]
    VEL=[]
    DIR=[]
    #BLO=[]
    Z=[]
    SNR=[]
    SNR_pp90=[]
    LDR=[]
    LDR_STD=[]
    LDR_MAX=[]
    LDR_HIST=[]
    Z_HIST=[]
    v_HIST=[]
    vl_HIST=[]
    vl_HIST_liquid=[]
    vl_HIST_single_liquid=[]
    vl_HIST_mixed=[]
    vl_HIST_single_mixed=[]
    vl_HIST_CTT=[]
    vl_HIST_CTH=[]
    vl_HIST_VEL=[]
    vl_HIST_DIR=[]
    n_LDR=[]
    IWC=[]
    n_IWC=[]
    ILCR=[]
    n_ILCR=[]
    n_ILCR_l=[]
    n_ILCR_m=[]
    n_ILCR_i=[]
    IWC_TOP_AVG=[]
    LWC=[]
    n_LWC=[]
    LWP=[]
    LWP_S=[]
    IWP=[]
    IWP_STD=[]
    v=[]
    n_v=[]
    v_lidar=[]
    std_lidar=[]
    n_lidar=[]
    gt0_std=[]
    gt0_mean=[]
    width=[]
    n_width=[]
    ct_width=[]
    alpha_hogan=[]
    beta=[]
    n_beta=[]
    delta=[]
    n_delta=[]
    n_Profiles=[]
    n_MixedPhase=[]
    n_Melting=[]
    n_Liquid=[]
    n_Drizzle=[]
    begin_datetime=[]
    N=[]
    Re=[]
    ice_positive=[]
    drizzle_positive=[]
    thickness=[]
    dust_conc=[]
    obs_time=[]
    Z_top=[]
    
    for i in range(len(datafile)):
    
        identifier=int(datafile[i]["A_Unique_Identifier"])
        current_datetime=datetime.datetime.utcfromtimestamp(float(datafile[i]["Begin_Date_Unix"]))
        begin=float(datafile[i]["Begin_Date_Unix"])
        end=float(datafile[i]["End_Date_Unix"])
        ctt=float(datafile[i]["CTT_MED"])
        cth=float(datafile[i]["CTH"])
        cbh=float(datafile[i]["CBH"])
        n_profiles=float(datafile[i]["N_Profiles"])
        n_mixed=float(datafile[i]["N_MixedPhase"])
        n_liquid=float(datafile[i]["N_Liquid"])
        z=float(datafile[i]["Z_MED"])
        N_delta=float(datafile[i]["delta_N"])
        delta_h=float(datafile[i]["Cloud_Thickness_AVG"])
        std_h=float(datafile[i]["CTH_STD"])
        cth_diffsum=float(datafile[i]["CTH_DIFFSUM"])
        llh_std=float(datafile[i]["LLH_STD"])
        Cloud_Run_Settings = datafile[i]["Cloud_Run"]

        year=current_datetime.year
        month=current_datetime.month
        day=current_datetime.day
        hour=current_datetime.hour
        
        #d_conc, d_ext = dust_model.get_dust(begin,
        print(end-begin, n_profiles, cth, delta_h, std_h, ctt)
        if Cloud_Run_Settings=="mixed-phase" and end-begin>900 and n_profiles/((end-begin)/30.0)>0.75 and cth>500 and delta_h<500.0 and std_h<300 and ctt>0: # and month in season[season_name] and current_datetime>datetime.datetime(2011,12,18):# and year==target_year and ctt<273.15:
            #print identifier
            print("identifier", identifier)

            ui.append(identifier)
            
            begin_datetime.append(current_datetime)
            obs_time.append(end-begin)
            
            CTT.append(float(datafile[i]["CTT_MED"]))
            CBT.append(float(datafile[i]["CBT"]))
            
            CTH.append(float(datafile[i]["CTH"]))

            CTP.append(float(datafile[i]["CTP"]))
            CBP.append(float(datafile[i]["CBP"]))
            
            VEL.append(float(datafile[i]["VEL"]))
            DIR.append(float(datafile[i]["DIR"]))

            #BLO.append(float(datafile[i]["CTH"])-500.0+2500*np.sin(2*np.pi/24.0*month))
    
            n_Profiles.append(n_profiles)
    
            Z.append(float(datafile[i]["Z_TOP_MED"]))
            SNR.append(float(datafile[i]["SNR_AVG"]))
            SNR_pp90.append(float(datafile[i]["SNR_90pp"]))
    
            LDR.append(float(datafile[i]["LDR_MED"]))
            LDR_STD.append(float(datafile[i]["LDR_STD"]))
            n_LDR.append(int(datafile[i]["LDR_N"]))
            #LDR_MAX.append(float(datafile[i]["LDR_MAX_MED"]))

            LDR_HIST.append(eval(datafile[i]["LDR_values"]))
            Z_HIST.append(eval(datafile[i]["Z_values"]))
            v_HIST.append(eval(datafile[i]["v_values"]))
            vl_HIST.append(eval(datafile[i]["v_lidar_histo"]))
            if identifier==201207120200456138:
            #if identifier==201108260046145389:
            #if identifier==201202210338445296:
                vl_HIST_single_mixed=np.array(eval(datafile[i]["v_lidar_histo"]))
                save_dict([vl_HIST_single_mixed],"HIST_L/vl_HIST_single_mixed")
                #print vl_HIST_single_mixed

            #if identifier==201207172249142926:
            if identifier==201208010257143519:
                vl_HIST_single_liquid=np.array(eval(datafile[i]["v_lidar_histo"]))
                save_dict([vl_HIST_single_liquid],"HIST_L/vl_HIST_single_liquid")
                #print vl_HIST_single_liquid
    
            IWC.append(float(datafile[i]["IWC_MED"]))
            n_IWC.append(int(datafile[i]["IWC_TOP_N"]))
            IWC_TOP_AVG.append(float(datafile[i]["IWC_TOP_MED"]))
            
            ILCR.append(float(datafile[i]["ILCR_MED"]))
            n_ILCR.append(float(datafile[i]["ILCR_N"]))
            
            ILCR_hist=np.array(eval(datafile[i]["ILCR_values"]))
            #if np.sum(ILCR_hist)>0:
            #    ILCR_hist/=np.sum(ILCR_hist)
            n_ILCR_t+=ILCR_hist
            n_ILCR_l.append(ILCR_hist[0])
            n_ILCR_m.append(np.sum(ILCR_hist[1:9]))
            n_ILCR_i.append(ILCR_hist[9])

            LWC.append(float(datafile[i]["LWC_MED"]))
            n_LWC.append(int(datafile[i]["LWC_N"]))
            
            LWP.append(float(datafile[i]["PATH_LWP_AVG"]))
            LWP_S.append(float(datafile[i]["PATH_LWP_S_AVG"]))
            IWP.append(float(datafile[i]["PATH_IWP_AVG"]))
            IWP_STD.append(float(datafile[i]["PATH_IWP_STD"]))
            
            v.append(float(datafile[i]["v_TOP_MED"]))
            n_v.append(int(datafile[i]["v_TOP_N"]))
            
            width.append(float(datafile[i]["width_MED"]))
            n_width.append(int(datafile[i]["width_N"]))

            n_MixedPhase.append(int(datafile[i]["N_MixedPhase"]))
            n_Liquid.append(int(datafile[i]["N_Liquid"]))
            n_Melting.append(int(datafile[i]["N_Melting"]))
            
            n_Drizzle.append(int(datafile[i]["N_Drizzle"]))
            
            alpha_hogan.append(float(datafile[i]["alpha_Hogan_AVG"]))
            
            beta.append(float(datafile[i]["beta_MED"]))
            n_beta.append(int(datafile[i]["beta_N"]))
            delta.append(float(datafile[i]["delta_AVG"]))
            n_delta.append(int(datafile[i]["delta_N"]))
            
            thickness.append(float(datafile[i]["Cloud_Thickness_MED"]))
            
            v_lidar.append(float(datafile[i]["v_lidar_AVG"]))
            std_lidar.append(float(datafile[i]["v_lidar_STD"]))
            n_lidar.append(int(datafile[i]["v_lidar_N"]))
            gt0_std.append(float(datafile[i]["v_lgt0_AVG"]))    
            gt0_mean.append(float(datafile[i]["v_lgt0_STD"]))

            ct_width.append(float(datafile[i]["v_radar_WIDTH"]))
            Z_top.append(float(datafile[i]["Z_top"]))
 
            if n_IWC[-1]>0:
                if (n_LDR[-1]>2 and LDR[-1]>-27) or (n_delta[-1]>2 and delta[-1]>0.3) or n_Melting[-1]>3:
                    ice_positive.append(1.0)
                elif n_Melting[-1]<3 and n_Drizzle[-1]>0:
                    ice_positive.append(0)
                elif n_Drizzle[-1]>0:
                    ice_positive.append(0)
                else:
                    ice_positive.append(-1)
            else:
               ice_positive.append(-1)
            
            if np.sum(eval(datafile[i]["v_lidar_histo"]))>100:# and v_lidar[-1]>-0.1 and v_lidar[-1]<0.1:
                if n_IWC[-1]>0:#/float(n_Profiles[-1])>0.1:
                    #print identifier
                    vl_HIST_mixed.append(eval(datafile[i]["v_lidar_histo"]))
                    vl_HIST_liquid.append(np.histogram([0],60,(-3.0,3.0))[0])
                    vl_HIST_CTT.append(CTT[-1])
                    vl_HIST_CTH.append(CTH[-1])
                    vl_HIST_VEL.append(VEL[-1])
                    vl_HIST_DIR.append(DIR[-1])
                else:
                    if n_Drizzle[-1]==0:
                        vl_HIST_liquid.append(eval(datafile[i]["v_lidar_histo"]))
                        vl_HIST_mixed.append(np.histogram([0],60,(-3.0,3.0))[0])
                        vl_HIST_CTT.append(CTT[-1])
                        vl_HIST_CTH.append(CTH[-1])
                        vl_HIST_VEL.append(VEL[-1])
                        vl_HIST_DIR.append(DIR[-1])
            
            
            if season_name=="total_year":
                if ctt<273.0 and ctt>267.0:
                
                    total_SSC+=1
                            
                    if identifier in ice_positive_list:
                        mixed_listed+=1
                    if identifier not in ice_positive_list:
                        liquid_listed+=1

                    if (ice_positive[-1]>0 and (identifier in ice_positive_list) ):
                        true_mixed+=1
                    elif (ice_positive[-1]<=0 and (identifier not in ice_positive_list) ):
                        true_liquid+=1
                    elif (ice_positive[-1]>0 and (identifier not in ice_positive_list) ):
                        false_mixed+=1
                    elif (ice_positive[-1]<=0 and (identifier in ice_positive_list) ):
                        false_liquid+=1
                        #print identifier
                    else:
                        print("ERROR")
                        exit()
            
            #if ctt<273.0 and ctt>267.0 and n_IWC[-1]>0 and (identifier not in ice_positive_list):
            #if ctt<273.0 and ctt>267.0 and ice_positive[-1]<=0:
            #    #print ice_positive[-1]
            #    n_Melting[-1]=0
            #    n_LDR[-1]=0
            #    n_IWC[-1]=0          
            
            #dust_conc.append(d_conc)

            p_type=14
            if v[-1]<0 and CTT[-1]>0:
               
                #if LDR[-1]==0:
                        
                ctt=CTT[-1]-273.15
    
                if ctt>-40 and ctt<-22:
                    p_type=4
                if ctt>-22 and ctt<-18:
                    p_type=15
                if ctt>-18 and ctt<-15:
                    p_type=8
                if ctt>-15 and ctt<-12:
                    p_type=12
                if ctt>-12 and ctt<-8:
                    p_type=17
                if ctt>-8 and ctt<0:
                    p_type=5

        
        
                #elif LDR[-1]>-22 and LDR[-1]<-15:
                #    p_type=6
                #elif LDR[-1]>-26 and LDR[-1]<-22:
                #    p_type=5
                #else:
                #    p_type=6
                
                N_est,Re_est=estimate_N(Z[-1],-v[-1],(CTP[-1]+CBP[-1])/2.0,(CTT[-1]+CBT[-1])/2.0,p_type)
                N.append(N_est)
                Re.append(Re_est)
            else:
                N.append(0)
                Re.append(0)

    ui=np.array(ui)
    CTT=np.array(CTT)
    CBT=np.array(CBT)
    CTH=np.array(CTH)
    CTP=np.array(CTP)
    CBP=np.array(CBP)
    VEL=np.array(VEL)
    DIR=np.array(DIR)
    Z=np.array(Z)
    SNR=np.array(SNR)
    SNR_pp90=np.array(SNR_pp90)
    LDR=np.array(LDR)
    LDR_STD=np.array(LDR_STD)
    n_LDR=np.array(n_LDR)
    LDR_MAX=np.array(LDR_MAX)
    IWC=np.array(IWC)
    IWC_TOP_AVG=np.array(IWC_TOP_AVG)
    n_IWC=np.array(n_IWC)
    ILCR=np.array(ILCR)
    n_ILCR=np.array(n_ILCR)
    LWC=np.array(LWC)
    n_LWC=np.array(n_LWC)
    LWP=np.array(LWP)
    LWP_S=np.array(LWP_S)
    IWP=np.array(IWP)
    IWP_STD=np.array(IWP_STD)
    v=np.array(v)
    n_v=np.array(n_v)
    width=np.array(width)
    n_width=np.array(n_width)
    n_Profiles=np.array(n_Profiles)
    n_MixedPhase=np.array(n_MixedPhase)
    n_Melting=np.array(n_Melting)
    n_Drizzle=np.array(n_Drizzle)
    alpha_hogan=np.array(alpha_hogan)
    beta=np.array(beta)
    n_beta=np.array(n_beta)
    delta=np.array(delta)
    n_delta=np.array(n_delta)
    n_Liquid=np.array(n_Liquid)
    N=np.array(N)
    Re=np.array(Re)
    ice_positive=np.array(ice_positive)
    thickness=np.array(thickness)
    dust_conc=np.array(dust_conc)
    v_lidar=np.array(v_lidar)
    std_lidar=np.array(std_lidar)
    n_lidar=np.array(n_lidar)
    obs_time=np.array(obs_time)
    gt0_mean=np.array(gt0_mean)
    gt0_std=np.array(gt0_std) 
    ct_width=np.array(ct_width)
    n_ILCR_l=np.array(n_ILCR_l)
    n_ILCR_m=np.array(n_ILCR_m)
    n_ILCR_i=np.array(n_ILCR_i)
    Z_top=np.array(Z_top)

    dust_conc_LOG=np.log10(dust_conc)
    dust_thr=(dust_conc>2e-9)
    LWC_LOG=np.log10(LWC)
    LWP_LOG=np.log10(LWP)  
    
    Z[Z>=0]=np.NAN
    LDR[n_LDR<30]=0
    LDR_STD[n_LDR<30]=0
    LDR_invalidation=(CTT<263.0)*(LDR>-20)
    LDR[LDR_invalidation]=0
    IWC[n_IWC<5]=0
    v[n_v<10]=0
    IWC_log=np.log10(IWC)
    delta[delta==0]=1
    #beta[n_delta>0]=0
    beta[n_beta<100]=0
    beta_log=np.log10(beta)
    std_lidar[n_lidar<30]=0
    v_lidar[n_lidar<30]=0
    #IWC_TOP_AVG=IWC_Correction(IWC_TOP_AVG,CTT-273.15)
    wvm=wv(CTT-273.15)*thickness
    SNR_pp90[SNR_pp90==0]=np.NAN
    SNR_pp90[SNR_pp90>22]=np.NAN
    SNR[SNR==0]=np.NAN
    SNR[SNR>26]=np.NAN
  
    #search for ice missclassifications
    #index=(LDR==0)*(n_Melting==0)*(CTT>263.0)
    #index=(ice_positive==0)*(CTT>260.0)
    #index[CTT<260]=0
    #index=(n_IWC<0.1*n_Profiles)
    #n_Liquid[index]+=n_MixedPhase[index]
    #n_MixedPhase[index]=0
    #IWC[index>0]=0
    #IWP[index>0]=0
    #ILCR[index>0]=0
   
    #LDR[:]=-32.5

    SNR_min=10**(-23.0/10)
    hgts=np.arange(100,10000,100)
    Z_min = 10*np.log10( SNR_min * 0.00254362123253 / 5000.0**2 * (hgts)**2)
    #Z_min_Leipzig = 10*np.log10( SNR_min * 0.00254362123253 / 5000.0**2 * (hgts-117.0)**2)
    #Z_min_Potenza = 10*np.log10( SNR_min * 0.00254362123253 / 5000.0**2 * (hgts-760.0)**2)
    #Z_min_MaceHead = 10*np.log10( SNR_min * 0.00254362123253 / 5000.0**2 * (hgts-5.0)**2)
    #Z_min=np.concatenate((Z_min_Leipzig,Z_min_Potenza,Z_min_MaceHead))
    #hgts=np.concatenate((hgts,hgts,hgts))
    
    #particle-depol
    mol_bs=molecular_backscatter(512e-9,(CTP-CBP)/2.0,(CTT-CBT)/2.0)
    delta=delta*(mol_bs/beta+1)
        
    #molecular background-correction
    beta = alpha_hogan/20.0
    #beta=beta-molecular_backscatter(1064e-9,(CTP-CBP)/2.0,(CTT-CBT)/2.0)
    lidar_snr=beta/(molecular_backscatter(2022e-9,(CTP-CBP)/2.0,(CTT-CBT)/2.0))
    lidar_snr_thr=0
    
    index=IWC>1e-99
    #LDR[:] = -20
    plot_color_scatter(folder+'ldr_iwcta_ctt_'+season_name+'.png',CTT[index]-273.15,IWC_TOP_AVG[index],LDR[index], x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Ice Water Content [$\mathsf{kg}\,\mathsf{m}^{-3}$]", z_label="LDR [dB]", ylog=True, xlim=(-40,0) ,ylim=(1e-9,1e-4), zlim=(LDR_scale_min,LDR_scale_max), cmap_name="ldr_colormap")
    #plot_color_scatter(folder+'ldr_iwc_lwc_ratio_direct_'+season_name+'.png',CTT-273.15,ILCR,LDR, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Mean ILCR Ratio", z_label="Linear Depol. Ratio [dB]", ylog=True, xlim=(-40,0) ,ylim=(1e-5,1e-0), zlim=(LDR_scale_min,LDR_scale_max), cmap_name="ldr_colormap")
    plot_color_scatter(folder+'ldr_flux_ctt_'+season_name+'.png',CTT-273.15,(IWC_TOP_AVG*(-1.0)*v),LDR, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Ice Flux [$\mathsf{kg}\,\mathsf{s}^{-1}\,\mathsf{m}^{-2}$]", z_label="LDR [dB]", ylog=True, xlim=(-40,0) ,ylim=(1e-9,1e-4), zlim=(LDR_scale_min,LDR_scale_max), cmap_name="ldr_colormap")
    
    if season_name=="total_year":
        
        print( "Ice classification above -5deg:")
        print("true_liquid:",true_liquid,", true_mixed:",true_mixed,", false_liquid",false_liquid,", false_mixed:",false_mixed)
        print("Liquid listed:", liquid_listed, ", Mixed listed:",mixed_listed)


        #exit()

        LDR_complete=[]
        Z_complete=[]
        v_complete=[]
        CTT_complete=[]
        CTT_Z_complete=[]
        CTT_v_complete=[]
        
        for n in range(len(CTT)):
            if LDR[n]!=0:
                for ldr in range(0,len(LDR_HIST[n])):
                    for k in range(LDR_HIST[n][ldr]):
                        LDR_complete.append(-35+ldr/2.0)
                        CTT_complete.append(CTT[n])
        CTT_complete=np.array(CTT_complete)
        LDR_complete=np.array(LDR_complete)

        for n in range(len(CTT)):
            for z in range(0,len(Z_HIST[n])):
                for k in range(Z_HIST[n][z]):
                    Z_complete.append(z-70)
                    CTT_Z_complete.append(CTT[n])
        CTT_Z_complete=np.array(CTT_Z_complete)
        Z_complete=np.array(Z_complete)

        for n in range(len(CTT)):
            for vi in range(0,len(v_HIST[n])):
                for k in range(v_HIST[n][vi]):
                    v_complete.append(-1.0*(-1.5+vi*1.5/60.0))
                    CTT_v_complete.append(CTT[n])
        CTT_v_complete=np.array(CTT_v_complete)
        v_complete=np.array(v_complete)
        
        x=np.arange(-2.9,3.1,0.1)
        
        output_directory="HIST_L/"
        if station_name=="LEIPZIG":
            save_dict(vl_HIST,output_directory+"vl_HIST")
            save_dict(vl_HIST_liquid,output_directory+"vl_HIST_liquid")
            save_dict(vl_HIST_mixed,output_directory+"vl_HIST_mixed")
            save_dict(vl_HIST_CTT,output_directory+"vl_HIST_CTT")
            save_dict(vl_HIST_CTH,output_directory+"vl_HIST_CTH")
            save_dict(vl_HIST_DIR,output_directory+"vl_HIST_DIR")
            save_dict(vl_HIST_VEL,output_directory+"vl_HIST_VEL")
        
        print(len(vl_HIST_liquid), len(vl_HIST_mixed))
        
        vl_HIST=np.sum(vl_HIST,axis=0)/float(np.sum(vl_HIST))
        vl_HIST_liquid=np.sum(vl_HIST_liquid,axis=0)/float(np.sum(vl_HIST_liquid))
        vl_HIST_mixed=np.sum(vl_HIST_mixed,axis=0)/float(np.sum(vl_HIST_mixed))
        vl_HIST_single_mixed=np.array(vl_HIST_single_mixed)/float(np.sum(vl_HIST_single_mixed))
        vl_HIST_single_liquid=np.array(vl_HIST_single_liquid)/float(np.sum(vl_HIST_single_liquid))
                
        plot_2D_histogram(folder+'ldr_'+season_name+'.png',CTT_complete-273.15,LDR_complete,Nx=15,Ny=60,rangex=(-40,0),rangey=(-35,-5),x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Linear Depol. Ratio [dB]", norm="column")
        plot_2D_histogram(folder+'z_'+season_name+'.png',CTT_Z_complete-273.15,Z_complete,Nx=15,Ny=90,rangex=(-40,0),rangey=(-70,20),x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Radar Reflectivity Factor [dBZ]", norm="column", display_y=(-60,0))
        plot_2D_histogram(folder+'v_'+season_name+'.png',CTT_v_complete-273.15,v_complete,Nx=15,Ny=59,rangex=(-40,0),rangey=(0,1.5),x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Doppler Velocity [$\mathsf{m}\,\mathsf{s}^{-1}$]", norm="column")
        plot_color_scatter(folder+'ldr_iwc_ctt_'+season_name+'.png',CTT-273.15,IWC,LDR, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Ice Water Content [$\mathsf{kg}\,\mathsf{m}^{-3}$]", z_label="Linear Depol. Ratio [dB]", ylog=True, xlim=(-40,0) ,ylim=(1e-9,1e-4), zlim=(LDR_scale_min,LDR_scale_max), cmap_name="ldr_colormap")
        plot_color_scatter(folder+'ice_iwc_ctt_'+season_name+'.png',CTT-273.15,IWC,ice_positive, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Ice Water Content [$\mathsf{kg}\,\mathsf{m}^{-3}$]", z_label="Definitive Ice Detection [1=Yes,0=No]", ylog=True, xlim=(-40,0) ,ylim=(1e-9,1e-4), zlim=(-1,1.25))
        plot_color_scatter(folder+'melting_iwc_ctt_'+season_name+'.png',CTT-273.15,IWC_TOP_AVG,1.0*(n_Melting>0), x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Ice Water Content [$\mathsf{kg}\,\mathsf{m}^{-3}$]", z_label="Melting Ratio", ylog=True, xlim=(-40,0) ,ylim=(1e-9,1e-4), zlim=(-0.25,1.25))
        plot_color_scatter(folder+'drizzle_iwc_ctt_'+season_name+'.png',CTT-273.15,IWC,1.0*(n_Drizzle>0), x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Ice Water Content [$\mathsf{kg}\,\mathsf{m}^{-3}$]", z_label="Drizzle Ratio", ylog=True, xlim=(-40,0) ,ylim=(1e-9,1e-4), zlim=(-0.25,1.25))
        #plot_color_scatter(folder+'ldrmax_iwc_ctt_'+season_name+'.png',CTT-273.15,IWC,LDR_MAX, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Ice Water Content [$\mathsf{kg}\,\mathsf{m}^{-3}$]", z_label="Max. Linear Depol. Ratio [dB]", ylog=True, xlim=(-40,0) ,ylim=(1e-9,1e-4), zlim=(LDR_scale_min,LDR_scale_max), cmap_name="ldr_colormap")
        #plot_color_scatter(folder+'ldr_alphahogan_ctt_'+season_name+'.png',CTT-273.15,alpha_hogan*1e6,LDR, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Alpha [1/Mm]", z_label="Linear Depol. Ratio [dB]", ylog=True, xlim=(-40,0) ,ylim=(1e-1,1e3), zlim=(LDR_scale_min,LDR_scale_max), cmap_name="ldr_colormap")
        #plot_color_scatter(folder+'ldrSTD_iwc_ctt_'+season_name+'.png',CTT-273.15,IWC,LDR_STD, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Ice Water Content [$\mathsf{kg}\,\mathsf{m}^{-3}$]", z_label="LDR STD [dB]", ylog=True, xlim=(-40,0) ,ylim=(1e-9,1e-4), zlim=(0,10))

        #plot_color_scatter(folder+'ldr_N_ctt_'+season_name+'.png',CTT[index]-273.15,N[index],LDR[index], x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Particle Number Conc. [1/m^3]", z_label="Linear Depol. Ratio [dB]", ylog=True, xlim=(-40,0) ,ylim=(0.1,1e7), zlim=(LDR_scale_min,LDR_scale_max), cmap_name="ldr_colormap")
        #plot2D(folder+'v_ctt_'+season_name+'.png',CTT[LDR!=0]-273.15,v[LDR!=0], x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Particle Fall Velocity [$\mathsf{m}\,\mathsf{s}^{-1}$]", ylog=False, xlim=(-40,0) ,ylim=(-2,0),xerr=1.5,yerr=np.abs(v[LDR!=0]/np.sqrt(n_v[LDR!=0]+1)))
        #plot2D(folder+'v_ctp_'+season_name+'.png',CTP[LDR!=0]/100.0,v[LDR!=0], x_label="Cloud Top Pressure [hPa]", y_label="Particle Fall Velocity [$\mathsf{m}\,\mathsf{s}^{-1}$]", ylog=False, xlim=(200,1000) ,ylim=(-2,0),xerr=1.5,yerr=np.abs(v[LDR!=0]/np.sqrt(n_v[LDR!=0]+1)))

        index=IWC>1e-99
	
        valid_dt=[]
        for d in range(len(begin_datetime)):
            if index[d]==True:          
    	        valid_dt.append(begin_datetime[d])
    	        
    	#for i in range(len(index)):
    	#    if index[i]==True and CTT[i]<273 and CTT[i]>260:
    	#        print thickness[i],IWC_TOP_AVG[i]

        v[v==0]=np.NAN
        plot_color_scatter(folder+'ldr_v_ctt_cor_'+season_name+'.png',CTT[index]-273.15,-1.0*v[index]*((CTP[index]*292.0)/(CTT[index]*102300.0))**(0.7),LDR[index], x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Doppler Velocity [$\mathsf{m}\,\mathsf{s}^{-1}$]", z_label="LDR [dB]", ylog=False, xlim=(-40,0) ,ylim=(0,0.8), zlim=(LDR_scale_min,LDR_scale_max), cmap_name="ldr_colormap")#,dates=valid_dt)
        plot_color_scatter(folder+'ldr_v_ctt_'+season_name+'.png',CTT[index]-273.15,-1.0*v[index],LDR[index], x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Doppler Velocity [$\mathsf{m}\,\mathsf{s}^{-1}$]", z_label="LDR [dB]", ylog=False, xlim=(-40,0) ,ylim=(0,1.2), zlim=(LDR_scale_min,LDR_scale_max), cmap_name="ldr_colormap")#,dates=valid_dt)
        plot_color_scatter(folder+'thickness_v_ctt_'+season_name+'.png',CTT[index]-273.15,-1.0*v[index]*((CTP[index]*292.0)/(CTT[index]*102300.0))**0.7,thickness[index], x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Doppler Velocity [$\mathsf{m}\,\mathsf{s}^{-1}$]", z_label="Thickness [m]", ylog=False, xlim=(-40,0) ,ylim=(0,0.8), zlim=(0,600))
        plot_color_scatter(folder+'ldr_v_ctp_'+season_name+'.png',CTT[index]-273.15,-1.0*v[index]*((CTP[index]*292.0)/(CTT[index]*102300.0))**(0.7),CTP[index]/100.0, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Doppler Velocity [$\mathsf{m}\,\mathsf{s}^{-1}$]", z_label="Top Pressure [hPa]", ylog=False, xlim=(-40,0) ,ylim=(0,1.0), zlim=(300,1000))
        LDR_indexed=LDR[index]
        LDR_indexed[LDR_indexed==0]=np.NAN
        plot_color_scatter(folder+'LDR_ctt_'+season_name+'.png',CTT[index]-273.15,LDR_indexed, SNR[index], x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Linear Depol. Ratio [dB]", z_label="Radar SNR [dB]", ylog=False, xlim=(-40,0) ,ylim=(LDR_scale_min,LDR_scale_max), zlim=(-20,30))
        plot_color_scatter(folder+'ldr_iwcta_ctt_'+season_name+'.svg',CTT[index]-273.15,IWC_TOP_AVG[index],LDR[index], x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Ice Water Content [$\mathsf{kg}\,\mathsf{m}^{-3}$]", z_label="LDR [dB]", ylog=True, xlim=(-40,0) ,ylim=(1e-9,1e-5), zlim=(LDR_scale_min,LDR_scale_max), cmap_name="ldr_colormap")#, dates=valid_dt)
        plot_color_scatter(folder+'cth_iwcta_ctt_'+season_name+'.png',CTT[index]-273.15,IWC_TOP_AVG[index],CTH[index], x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Ice Water Content [$\mathsf{kg}\,\mathsf{m}^{-3}$]", z_label="Cloud Top Height [m]", ylog=True, xlim=(-40,0) ,ylim=(1e-9,1e-4), zlim=(1000,10000))#, dates=valid_dt)
        plot_color_scatter(folder+'ldr_iwcta_cth_'+season_name+'.png',CTH[index],IWC_TOP_AVG[index],LDR[index], x_label="Cloud Top Height [m]", y_label="Ice Water Content [$\mathsf{kg}\,\mathsf{m}^{-3}$]", z_label="LDR [dB]", ylog=True, xlim=(10000,0) ,ylim=(1e-9,1e-4), zlim=(LDR_scale_min,LDR_scale_max), cmap_name="ldr_colormap", meanplot=False)#, dates=valid_dt)
        plot_color_scatter(folder+'Ztop_iwcta_cth_'+season_name+'.png',CTH[index],IWC_TOP_AVG[index],Z_top[index], x_label="Cloud Top Height [m]", y_label="Ice Water Content [$\mathsf{kg}\,\mathsf{m}^{-3}$]", z_label="Z_top [dBZ]", ylog=True, xlim=(10000,0) ,ylim=(1e-9,1e-4), zlim=(-50,-10), meanplot=False)#, dates=valid_dt)
        plot_color_scatter(folder+'Ztop_iwcta_ctt_'+season_name+'.png',CTT[index]-273.15,IWC_TOP_AVG[index],Z_top[index], x_label="Cloud Top Temperatrue [degC]", y_label="Ice Water Content [$\mathsf{kg}\,\mathsf{m}^{-3}$]", z_label="Z_top [dBZ]", ylog=True, xlim=(-40,0) ,ylim=(1e-9,1e-4), zlim=(-50,-10), meanplot=False)#, dates=valid_dt)
        #plot_color_scatter(folder+'dust_N_ctt_'+season_name+'.png',CTT-273.15,N,dust_conc_LOG, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", ylog=True, y_label="N", z_label="Dust Threshold", xlim=(-40,0) ,ylim=(1e-1,1e5), zlim=(-14,-7))
        #plot_color_scatter(folder+'dust_iwcta_ctt_'+season_name+'.png',CTT[index]-273.15,IWC_TOP_AVG[index],dust_thr[index], x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="IWCTA", z_label="Dust Threshold", ylog=True, xlim=(-40,0) ,ylim=(1e-9,1e-4), zlim=(-1.2,1.2))        
        plot_color_scatter(folder+'v_iwcta_ctt_'+season_name+'.png',CTT[index]-273.15,IWC_TOP_AVG[index],-v[index], x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Ice Water Content [$\mathsf{kg}\,\mathsf{m}^{-3}$]", z_label="Doppler velocity [$\mathsf{m}\,\mathsf{s}^{-1}$]", ylog=True, xlim=(-40,0) ,ylim=(1e-9,1e-4), zlim=(0,1.5))
        #plot_color_scatter(folder+'lwp_iwcta_ctt_'+season_name+'.png',CTT[index]-273.15,IWC_TOP_AVG[index],LWP_LOG[index], x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Ice Water Content [$\mathsf{kg}\,\mathsf{m}^{-3}$]", z_label="LWP [log(kg)]", ylog=True, xlim=(-40,0) ,ylim=(1e-9,1e-4), zlim=(-3,-0))
        #plot_color_scatter(folder+'vwidth_iwcta_ctt_'+season_name+'.png',CTT[index]-273.15,IWC_TOP_AVG[index],width[index], x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Ice Water Content [$\mathsf{kg}\,\mathsf{m}^{-3}$]", z_label="velocity width [$\mathsf{m}\,\mathsf{s}^{-1}$]", ylog=True, xlim=(-40,0) ,ylim=(1e-9,1e-4), zlim=(0.0,0.4))
        #plot_color_scatter(folder+'ctp_iwcta_ctt_'+season_name+'.png',CTT[index]-273.15,IWC_TOP_AVG[index],CTP[index]/100.0, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Ice Water Content [$\mathsf{kg}\,\mathsf{m}^{-3}$]", z_label="CTP [hPa]", ylog=True, xlim=(-40,0) ,ylim=(1e-9,1e-4), zlim=(300,900))
        plot_color_scatter(folder+'lifetime_ctt_'+season_name+'.png',CTT[index]-273.15,(LWP[index]/IWC_TOP_AVG[index]*(-1.0)*v[index]),LDR[index], x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Static Lifetime Index [s]", z_label="LDR [dB]", ylog=True, xlim=(-40,0) ,ylim=(1e3,1e7), zlim=(LDR_scale_min,LDR_scale_max), cmap_name="ldr_colormap")        
        plot_color_scatter(folder+'lifetime_obstime_ctt_'+season_name+'.png',CTT[index]-273.15,(LWP[index]+wvm[index])/(IWC_TOP_AVG[index]*(-1.0)*v[index]),obs_time[index], x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Lifetime Index [s]", z_label="Observed Time [s]", ylog=True, xlim=(-40,0) ,ylim=(1e4,1e9), zlim=(900,30000))
        plot_color_scatter(folder+'ldr_flux_ctt_'+season_name+'.png',CTT[index]-273.15,(IWC_TOP_AVG[index]*(-1.0)*v[index]),LDR[index], x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Ice Flux [$\mathsf{kg}\,\mathsf{s}^{-1}\,\mathsf{m}^{-2}$]", z_label="LDR [dB]", ylog=True, xlim=(-40,0) ,ylim=(1e-9,1e-4), zlim=(LDR_scale_min,LDR_scale_max), cmap_name="ldr_colormap")
        plot_color_scatter(folder+'thickness_iwcta_ctt_'+season_name+'.png',CTT[index]-273.15,IWC_TOP_AVG[index],thickness[index], x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Ice Water Content [$\mathsf{kg}\,\mathsf{m}^{-3}$]", z_label="Cloud Thickness [m]", ylog=True, xlim=(-40,0) ,ylim=(1e-9,1e-4), zlim=(0,500))#, dates=valid_dt)
        plot_color_scatter(folder+'thickness_lwp_ctt_'+season_name+'.png',CTT-273.15,LWP,thickness, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Liquid Water Path [$\mathsf{kg}\,\mathsf{m}^{-2}$]", z_label="Cloud Thickness [m]", ylog=True, xlim=(-40,10) ,ylim=(1e-3,1e-0), zlim=(0,500))
        plot_color_scatter(folder+'Ztop_lwp_ctt_'+season_name+'.png',CTT-273.15,LWP,Z_top, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Liquid Water Path [$\mathsf{kg}\,\mathsf{m}^{-2}$]", z_label="Z_top [dBZ]", ylog=True, xlim=(-40,10) ,ylim=(1e-3,1e-0), zlim=(-50,-10))
        plot_color_scatter(folder+'ldr_iwc_lwc_ratio_'+season_name+'.png',CTT[index]-273.15,IWC_TOP_AVG[index]/(2.0*LWC[index]),LDR[index], x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Mean ILCR Ratio", z_label="Linear Depol. Ratio [dB]", ylog=True, xlim=(-40,0) ,ylim=(1e-7,1e1), zlim=(LDR_scale_min,LDR_scale_max), cmap_name="ldr_colormap")
        #plot_color_scatter(folder+'lwp_v_ctt_'+season_name+'.png',CTT[index]-273.15,v[index],LWP_LOG[index], x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Fall Velocity [$\mathsf{m}\,\mathsf{s}^{-1}]", z_label="LWP [kg/m^2]", ylog=False, xlim=(-40,0) ,ylim=(-2,0), zlim=(-3,0))
        plot_color_scatter(folder+'drizzle_lwp_ctt_'+season_name+'.png',CTT-273.15,LWP,1.0*(n_Drizzle>0), x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Liquid Water Path [$\mathsf{kg}\,\mathsf{m}^{-2}$]", z_label="Drizzle Ratio", ylog=True, xlim=(-40,20) ,ylim=(1e-4,1e-1), zlim=(-0.25,1.25))#, dates=begin_datetime)
        

        """
        index=[]
        for i in range(len(begin_datetime)):
            if begin_datetime[i].date()==datetime.date(2014,4,1):
                print "hit"
                index.append(1)
            else:
                index.append(0)
        index=np.array(index)
        """
        
        N_log=np.log10(N)
        #N_log[index==0]=1e99
        index=(N_log>-1e99)
        plot_color_scatter(folder+'N_iwc_ctt_'+season_name+'.png',CTT-273.15,IWC_TOP_AVG,N_log, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Ice Water Content [$\mathsf{kg}\,\mathsf{m}^{-3}$]", z_label="Ice Part. Number Con. [1/m^3]", ylog=True, xlim=(-40,0) ,ylim=(1e-9,1e-4), zlim=(2,9))
        plot_color_scatter(folder+'N_ilr_ctt_'+season_name+'.png',CTT-273.15,ILCR,N_log, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Ice Liquid Ratio (direct)", z_label="Ice Part. Number Con. [1/m^3]", ylog=True, xlim=(-40,0), ylim=(1e-5,1e0), zlim=(2,9))
        plot_color_scatter(folder+'ldr_N_ctt_'+season_name+'.png',CTT-273.15,N, LDR, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="N", z_label="LDR", xlim=(-40,0), ylim=(1e-1,1e6), zlim=(LDR_scale_min, LDR_scale_max), ylog=True, cmap_name="ldr_colormap")
        plot_color_scatter(folder+'lwp_N_ctt_'+season_name+'.png',CTT-273.15,N, LWP, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="N", z_label="LWP", xlim=(-40,0), ylim=(1e-1,1e5), zlim=(1e-4,1e-1), ylog=True)
        plot_color_scatter(folder+'ldr_iwp_lwp_ratio_'+season_name+'.png',CTT-273.15,IWP/LWP,LDR, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Mean ILPR Ratio", z_label="Linear Depol. Ratio [dB]", ylog=True, xlim=(-40,0) ,ylim=(1e-6,1e1), zlim=(LDR_scale_min,LDR_scale_max), cmap_name="ldr_colormap")
        plot_color_scatter(folder+'ldr_lwp_ctt_'+season_name+'.png',CTT-273.15,LWP,LDR,x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]",y_label="Liquid Water Path [$\mathsf{kg}\,\mathsf{m}^{-2}$]", z_label="LDR [dB]", ylog=True, xlim=(-40,10), ylim=(1e-4,1e0), zlim=(LDR_scale_min,LDR_scale_max), cmap_name="ldr_colormap")
        plot_color_scatter(folder+'ldr_lwps_ctt_'+season_name+'.png',CTT-273.15,LWP_S,LDR,x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]",y_label="Liquid Water Path [$\mathsf{kg}\,\mathsf{m}^{-2}$]", z_label="LDR [dB]", ylog=True, xlim=(-40,10), ylim=(1e-4,1e0), zlim=(LDR_scale_min,LDR_scale_max), cmap_name="ldr_colormap")
        ctt_index=(CTT<(273.15-12.0))
        plot_color_scatter(folder+'ctt_iwcta_lwp_'+season_name+'.png',IWC_TOP_AVG[ctt_index],LWP[ctt_index],CTT[ctt_index]-273.15,x_label="IWCTA",y_label="LWP", z_label="CTT", xlog=True, ylog=True, xlim=(1e-8,1e-5), ylim=(1e-4,1e0), zlim=(-40,0))
        plot_color_scatter(folder+'ldr_iwp_ctt_'+season_name+'.png',CTT-273.15,IWP,LDR,x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]",y_label="Ice Water Path [kg/m^2]", z_label="Linear Depol. Ratio [dB]", ylog=True, xlim=(-40,0),    ylim=(1e-8,1e-0), zlim=(LDR_scale_min,LDR_scale_max), cmap_name="ldr_colormap")
        #std_lidar[CTH<4000]=0   
        #plot_color_scatter(folder+'vlstd_lwp_ctt_'+season_name+'.png',CTT[std_lidar>0]-273.15,LWP_S[std_lidar>0],std_lidar[std_lidar>0],x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]",y_label="Liquid Water Path [$\mathsf{kg}\,\mathsf{m}^{-2}$]", z_label="Vertical Velocity STD [$\mathsf{m}\,\mathsf{s}^{-1}$]", ylog=True, xlim=(-40,0), ylim=(1e-3,1e-1), zlim=(0,0.6), meanplot=False)
        plot_color_scatter(folder+'vrwidth_lwp_ctt_'+season_name+'.png',CTT[ct_width>0]-273.15,LWP_S[ct_width>0],ct_width[ct_width>0],x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]",y_label="Liquid Water Path [$\mathsf{kg}\,\mathsf{m}^{-2}$]", z_label="Vertical Velocity STD [$\mathsf{m}\,\mathsf{s}^{-1}$]", ylog=True, xlim=(-40,0), ylim=(1e-3,1e-1), zlim=(0,0.4), meanplot=False)
        plot_color_scatter(folder+'vlstd_iwcta_ctt_'+season_name+'.png',CTT[std_lidar>0]-273.15,IWC_TOP_AVG[std_lidar>0],std_lidar[std_lidar>0],x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]",y_label="IWC at Cloud Top [$\mathsf{kg}\,\mathsf{m}^{-3}$]", z_label="Vertical Velocity STD [$\mathsf{m}\,\mathsf{s}^{-1}$]", ylog=True, xlim=(-40,0), ylim=(1e-9,1e-4), zlim=(0,2.0), meanplot=False)
        #plot_color_scatter(folder+'vlstd_flux_ctt_'+season_name+'.png',CTT[std_lidar>0]-273.15,IWC_TOP_AVG[std_lidar>0]*-v[std_lidar>0],std_lidar[std_lidar>0],x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]",y_label="Ice Mass Flux [$\mathsf{kg}\,\mathsf{m}^{-2}\,}\mathsf{s}^{-1}$]", z_label="Vertical Velocity STD [$\mathsf{m}\,\mathsf{s}^{-1}$]", ylog=True, xlim=(-40,0), ylim=(1e-9,1e-4), zlim=(0,0.6), meanplot=False)
        plot_color_scatter(folder+'vrwidth_iwcta_ctt_'+season_name+'.png',CTT[ct_width>0]-273.15,IWC_TOP_AVG[ct_width>0],ct_width[ct_width>0],x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]",y_label="IWC at Cloud Top [$\mathsf{kg}\,\mathsf{m}^{-3}$]", z_label="Cloud Top Spectrum Width [$\mathsf{m}\,\mathsf{s}^{-1}$]", ylog=True, xlim=(-40,0), ylim=(1e-8,1e-5), zlim=(0,0.4), meanplot=False)

        #Vertical Velocity Analysis
        print( "std_lidar gt0_std gt0_mean IWC_TOP_AVG LWP LWP_S CTT CTH ui v F")
        for i in range(len(std_lidar)):
            if std_lidar[i]!=0 and CTT[i]-273.15<-12.0 and IWC_TOP_AVG[i]>0:
                print(std_lidar[i],gt0_std[i],gt0_mean[i],IWC_TOP_AVG[i], LWP[i], LWP_S[i],CTT[i]-273.15,CTH[i],ui[i], v[i], IWC_TOP_AVG[i]*(-1.0)*v[i])
                #plt.plot(vl_HIST[i])
                #plt.savefig("vl_test"+str(std_lidar[i])+"_"+str(np.log10(IWC_TOP_AVG[i]))+".png")
                #plt.clf()

        print("Z_LDR_ctt")
        print(plot_color_scatter(folder+'Z_LDR_ctt_'+season_name+'.png',CTT-273.15,Z,LDR, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Radar Reflectivity Factor [dBZ]", z_label="Linear Depol. Ratio [dB]", xlim=(-40,0) ,ylim=(-60,0), zlim=(LDR_scale_min,LDR_scale_max), cmap_name="ldr_colormap"))#, dates=begin_datetime))
        plot_color_scatter(folder+'ldr_iwc_lwc_ratio_direct_'+season_name+'.svg',CTT-273.15,ILCR,LDR, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="ice-to-liquid mass ratio", z_label="Linear Depol. Ratio [dB]", ylog=True, xlim=(-40,0) ,ylim=(1e-5,1e-0), zlim=(LDR_scale_min,LDR_scale_max), cmap_name="ldr_colormap")
        #plot_color_scatter(folder+'dust_iwc_lwc_ratio_direct_'+season_name+'.png',CTT-273.15,ILCR,dust_thr, dates=begin_datetime, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Mean ILCR Ratio (direct)", z_label="Linear Depol. Ratio [dB]", ylog=True, xlim=(-40,0) ,ylim=(1e-6,1e0), zlim=(-1.2,1.2))

        #CTD=CTP/(R_air*CTT)	
        #plot_color_scatter(folder+'ldr_iwc_lwc_ratio_direct_ctd_'+season_name+'.png',LDR,ILCR,CTD, x_label="Mean LDR [dB]", y_label="Mean ILCR Ratio (direct)", z_label="Cloud Top Density [$\mathsf{kg}\,\mathsf{m}^{-3}$]", ylog=True, xlim=(-35,-10) ,ylim=(1e-7,1e1), zlim=(0.3,1))
        #plot_color_scatter(folder+'ldr_iwc_lwc_ratio_'+season_name+'.png',CTT-273.15,IWC_TOP_AVG/(IWC_TOP_AVG+LWC),LDR, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Mean ILCR Ratio", z_label="Linear Depol. Ratio [dB]", ylog=True, xlim=(-40,0) ,ylim=(1e-5,1e1), zlim=(LDR_scale_min,LDR_scale_max), cmap_name="ldr_colormap")
        plot_color_scatter(folder+'ldr_lwc_ctt_'+season_name+'.png',CTT-273.15,LWC,LDR,x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]",y_label="Liquid Water Content [$\mathsf{kg}\,\mathsf{m}^{-3}$]", z_label="Linear Depol. Ratio [dB]", ylog=True, xlim=(-40,0), ylim=(1e-5,1e-3), zlim=(LDR_scale_min,LDR_scale_max), cmap_name="ldr_colormap", meanplot=True)
        #plot_color_scatter(folder+'delta_iwc_ctt_'+season_name+'.png',CTT-273.15,IWC,delta, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Ice Water Content [$\mathsf{kg}\,\mathsf{m}^{-3}$]", z_label="Lidar Depol. []", ylog=True, xlim=(-40,0) ,ylim=(1e-9,1e-4), zlim=(0,0.6))
        plot_color_scatter(folder+'LSNR_Z_ctt_'+season_name+'.png',CTT-273.15,Z,np.log10(lidar_snr), x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Radar Reflectivity Factor [dBZ]", z_label="Lidar SNR", xlim=(-40,0) ,ylim=(-60,0), zlim=(-2,1))
        plot_color_scatter(folder+'IWC_LSNR_ctt_'+season_name+'.png',CTT-273.15,np.log10(lidar_snr),np.log10(IWC_TOP_AVG), x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Lidar P/M SNR at 1064nm (log)", z_label="Ice Water Content [$\mathsf{kg}\,\mathsf{m}^{-3}$]", xlim=(-40,0) ,ylim=(-3,2), zlim=(-9,-4))
        #plot_color_scatter(folder+'Z_beta_ctt_'+season_name+'.png',CTT-273.15,Z,beta_log+6, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Radar Reflectivity Factor [dBZ]", z_label="Lidar Beta [log(1/(kmSr))]", xlim=(-40,0) ,ylim=(-60,0), zlim=(-1,3))
        #plot_color_scatter(folder+'LDR_Zmin_cth_'+season_name+'.png',np.concatenate((CTH-alt,hgts)),np.concatenate((Z,Z_min)),np.concatenate((LDR,np.repeat(-16.0,len(hgts)))), x_label="Cloud Top Range [m]", y_label="Radar Reflectivity Factor [dBZ]", z_label="Linear Depol. Ratio [dB]]", xlim=(10000,0) ,ylim=(-70,0), zlim=(LDR_scale_min,LDR_scale_max), meanplot=False, cmap_name="ldr_colormap")
        plot_color_scatter(folder+'LDR_Zmin_cth_'+season_name+'.png',CTH-alt,Z,LDR, x_label="Cloud Top Range [m]", y_label="Radar Reflectivity Factor [dBZ]", z_label="Linear Depol. Ratio [dB]]", xlim=(10000,0) ,ylim=(-70,0), zlim=(LDR_scale_min,LDR_scale_max), meanplot=False, cmap_name="ldr_colormap")
        plot_color_scatter(folder+'Z_ctt_cth_'+season_name+'.png',CTH-alt,CTT-273.15,Z, x_label="Cloud Top Range [m]", y_label="Cloud Top Temperature [deg]", z_label="Z [ dBZ]]", xlim=(10000,0) ,ylim=(-40,0), zlim=(-70,0), meanplot=False)
        plot_color_scatter(folder+'SNR_iwc_ctt_'+season_name+'.png',CTT-273.15,IWC,SNR, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Ice Water Content [$\mathsf{kg}\,\mathsf{m}^{-3}$]]", z_label="SNR [dB]", xlim=(-40,0) ,zlim=(-20,30), ylog=True, ylim=(1e-9,1e-4)) 
        print(len(SNR_pp90),len(CTT),len(LDR))
        plot_color_scatter(folder+'SNR_LDRthr_ctt_'+season_name+'.png',CTT-273.15,-22.0-SNR,LDR, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="LDR Detection Threshold [dB]", z_label="Linear Depol. Ratio [dB]", xlim=(-40,0) ,ylim=(0,-50), zlim=(LDR_scale_min,LDR_scale_max), cmap_name="ldr_colormap") 
        plot_color_scatter(folder+'SNR_LDR_ctt_'+season_name+'.png',CTT-273.15,SNR,LDR, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="SNR in ice virga [dB]", z_label="Linear Depol. Ratio [dB]", xlim=(-40,0) ,ylim=(-23,20), zlim=(LDR_scale_min,LDR_scale_max), cmap_name="ldr_colormap")
        #SNR_gt_m10=SNR_pp90[CTT-273.15<-10]
        #for i in range(len(SNR_gt_m10)):
        #    if ~np.isnan(SNR_gt_m10[i]):
        #        print SNR_gt_m10[i]
        #plot_color_scatter(folder+'SNR_ice_ctt_'+season_name+'.png',CTT-273.15,-22.0-SNR_pp90,ice_positive, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="LDR Detection Threshold [dB]", z_label="Ice Detection Level", xlim=(-40,0) ,ylim=(0,-50), zlim=(-0.2,1.2)) 
        #plot_color_scatter(folder+'beta_delta_ctt_'+season_name+'.png',CTT-273.15,beta*1e6,delta, x_label="Cloud Top Height [m]", y_label="Lidar beta [1/(MmSr)]", z_label="Lidar Depol. []", xlim=(-40,0) , ylog=True, ylim=(1e-1,1e3), zlim=(0,1))
        #plot_color_scatter(folder+'beta_delta_cmp_'+season_name+'.png',(CTP+CBP)/2.0/100.0,beta*1e6,delta, x_label="Mean Pressure [hPa]", y_label="Lidar beta [1/(MmSr)]", z_label="Lidar Depol. []", xlim=(300,1000) , ylog=True, ylim=(1e-1,1e3), zlim=(0,1))
        #plot_color_scatter(folder+'ldr_vs_delta_cth.png',x=LDR,y=delta, z=CTT-273.15, xlim=(-35,-15) , ylim=(0,0.6), zlim=(-40,0), x_label="Linear Depol. Ratio [dB]", y_label="delta", z_label="CTT [$\mathsf{{^\circ}C}$]")
        #plot_color_scatter(folder+'ldr_vs_delta_melting.png',x=LDR,y=delta, z=1.0*(n_Melting/n_Profiles>0.1), xlim=(-35,-15) , ylim=(0,0.6), zlim=(0,1.2), x_label="Linear Depol. Ratio [dB]", y_label="delta", z_label="Melting Layer Present (0=False, 1=True)")
        #plot_color_scatter(folder+'ldr_ctt.png',x=(CTT-273.15),y=LDR, z=1.0*(n_Drizzle>0), xlim=(-40,0) , ylim=(LDR_scale_min,LDR_scale_max), zlim=(-0.25,1.25), x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="LDR [dB]", z_label="Drizzle Detected [1=Yes,0=No]")

        x1=10**np.arange(-7,-3,0.01)
        x2=5*10**np.arange(-7,-3,0.01)
        x3=0.2*10**np.arange(-7,-3,0.01)
        xx=np.concatenate((x1,x2,x3))
        fxx=np.concatenate((x1,x1,x1))/20.0
        #plot2D(folder+'alpha_beta.png',x=np.concatenate((alpha_hogan,xx))*1e6,y=20*np.concatenate((beta,fxx))*1e6, xlim=(1e-1,1e3) , ylim=(1e-1,1e3), x_label="Alpha (Hogan) [1/Mm]", y_label="Beta*20sr (PollyXT) [1/Mm]", xlog=True, ylog=True)
        #plot_color_scatter(folder+'alpha_beta_ctt.png',x=alpha_hogan*1e6,y=20*beta*1e6, z=CTT-273.15, xlim=(1e-1,1e3) , ylim=(1e-1,1e3), zlim=(-40,0), x_label="Alpha (Hogan) [1/Mm]", y_label="Beta*20sr (PollyXT) [1/Mm]", z_label="CTH [m]", xlog=True, ylog=True)
        #plot_color_scatter(folder+'alpha_beta_Re.png',x=alpha_hogan[Re<100]*1e6,y=25*beta[Re<100]*1e6, z=Re[Re<100], xlim=(1e-1,1e3) , ylim=(1e-1,1e3), zlim=(0,200), x_label="Alpha (Hogan) [1/Mm]", y_label="Beta*25sr (PollyXT) [1/Mm]", z_label="Mean Reynolds Number []", xlog=True, ylog=True)
            
        plot_histogram(folder+'histo_ctt_'+season_name+'.png', CTT-273.15, (-50,20),18, x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Number of Cases", normed=False)
        plot_histogram(folder+'histo_ctt_MP_'+season_name+'.png', CTT[n_IWC>0]-273.15, (-50,20),18, x_label="Cloud Top Temperature (Mixed-Phase Clouds) [$\mathsf{{^\circ}C}$]", y_label="Number of Cases", normed=False)
        plot_histogram(folder+'histo_ctt_LP_'+season_name+'.png', CTT[n_IWC==0]-273.15, (-50,20),18, x_label="Cloud Top Temperature (Liquid Clouds) [$\mathsf{{^\circ}C}$]", y_label="Number of Cases", normed=False)
        plot_histogram(folder+'histo_cth_'+season_name+'.png', CTH, (0,11000),22, x_label="Cloud Top Height [m]", y_label="PDF [1/m]", normed=False)

    
        ldr_40_22=[]
        ldr_22_10=[]
        ldr_10_3=[]
        ldr_3_0=[]
        
        for i in range(len(CTT)):
            
            ctt=CTT[i]-273.15
            
            if ctt>-40 and ctt<-22:
                ldr_40_22.append(LDR[i])
            if ctt>-22 and ctt<-10:
                ldr_22_10.append(LDR[i])
            if ctt>-10 and ctt<-3:
                ldr_10_3.append(LDR[i])
            if ctt>-3 and ctt<0:
                ldr_3_0.append(LDR[i])
                
        ldr_40_22=np.array(ldr_40_22)
        ldr_22_10=np.array(ldr_22_10)
        ldr_10_3=np.array(ldr_10_3)
        ldr_3_0=np.array(ldr_3_0)
        
        steps=20
        ldr_range=(-40,0)
        plot_histogram(folder+'histo_ldr_40_22.png', ldr_40_22[ldr_40_22!=0], ldr_range,steps, x_label="Linear Depol. Ratio [dB]", y_label="PDF [1/dB]")
        plot_histogram(folder+'histo_ldr_22_10.png', ldr_22_10[ldr_22_10!=0], ldr_range,steps, x_label="Linear Depol. Ratio [dB]", y_label="PDF [1/dB]")
        plot_histogram(folder+'histo_ldr_10_3.png', ldr_10_3[ldr_10_3!=0], ldr_range,steps, x_label="Linear Depol. Ratio [dB]", y_label="PDF [1/dB]")
        plot_histogram(folder+'histo_ldr_3_0.png', ldr_3_0[ldr_3_0!=0], ldr_range,steps, x_label="Linear Depol. Ratio [dB]", y_label="PDF [1/dB]")


    
    mp_temps=[]
    liquid_temps=[]
    mp_hgts=[]
    liquid_hgts=[]

    """
    for i in range(len(CTT)):
        if n_MixedPhase[i]>0 and IWC[i]>0:
            #for n in range(n_MixedPhase[i]):
            mp_temps.append(CTT[i])
            mp_hgts.append(CTH[i])
        elif n_MixedPhase[i]==0 and n_Liquid[i]>0:
            #for n in range(n_Liquid[i]):
            liquid_temps.append(CTT[i])
            liquid_hgts.append(CTH[i])
    """



    #delta_limit=0.2
    #n_Liquid[delta<delta_limit]+=n_MixedPhase[delta<delta_limit]
    #n_MixedPhase[delta<delta_limit]=0
    
    #lsnr_limit=lidar_snr_thr
    #n_Liquid[lidar_snr<lsnr_limit]+=n_MixedPhase[lidar_snr<lsnr_limit]
    #n_MixedPhase[lidar_snr<lsnr_limit]=0
    #n_IWC[lidar_snr<lsnr_limit]=0

    #IWC_limit=2e-7
    #n_Liquid[IWC<IWC_limit]+=n_MixedPhase[IWC<IWC_limit]
    #n_MixedPhase[IWC<IWC_limit]=0
  
    #ILR_limit=1e-5
    #ILR=ILCR
    #n_Liquid[ILR<ILR_limit]+=n_MixedPhase[ILR<ILR_limit]
    #n_MixedPhase[ILR<ILR_limit]=0
    
    #n_Liquid[n_IWC<2*n_Profiles]+=n_MixedPhase[n_IWC<2*n_Profiles]
    #n_MixedPhase[n_IWC<2*n_Profiles]=0
  
    #alpha_limit=15e-6
    #n_Liquid[alpha_hogan<alpha_limit]+=n_MixedPhase[alpha_hogan<alpha_limit]
    #n_MixedPhase[alpha_hogan<alpha_limit]=0

    
    for i in range(len(CTT)):
    
        if n_IWC[i]/n_Profiles[i]>0.1:
            for n in range(n_IWC[i]):
                mp_temps.append(CTT[i])
                mp_hgts.append(CTH[i])
        else:
           n_Liquid[i]=n_Profiles[i]

        if n_Liquid[i]>0:
            for n in range(n_Liquid[i]):
                liquid_temps.append(CTT[i])
                liquid_hgts.append(CTH[i])    
    
    """             
    for i in range(len(CTT)):
    
        if n_IWC[i]>0:
            mp_temps.append(CTT[i])
            mp_hgts.append(CTH[i])
        else:
            liquid_temps.append(CTT[i])
            liquid_hgts.append(CTH[i])
    """
    """
    for i in range(len(CTT)):

        if lidar_snr[i]>lsnr_limit:
            mp_temps.append(CTT[i])
            mp_hgts.append(CTH[i])
        else:
            liquid_temps.append(CTT[i])
            liquid_hgts.append(CTH[i])
    """

    n_bins=10
    low_T=233.15
    high_T=283.15
    mp_histo_temp=np.histogram(mp_temps, n_bins, (low_T,high_T))
    liquid_histo_temp=np.histogram(liquid_temps, n_bins, (low_T,high_T))
    x_axis_temp=mp_histo_temp[1][0:-1]-273.15+(mp_histo_temp[1][1]-mp_histo_temp[1][0])/2.0
    
    n_bins=10
    low_H=0
    high_H=10000
    mp_histo_hgts=np.histogram(mp_hgts, n_bins, (low_H,high_H))
    liquid_histo_hgts=np.histogram(liquid_hgts, n_bins, (low_H,high_H))
    x_axis_hgts=mp_histo_hgts[1][0:-1]+(mp_histo_hgts[1][1]-mp_histo_hgts[1][0])/2.0
    
    total_temp=liquid_histo_temp[0]+mp_histo_temp[0]
    mpc_fraction_temp=np.divide(mp_histo_temp[0].astype(float),total_temp)

    total_hgts=liquid_histo_hgts[0]+mp_histo_hgts[0]
    mpc_fraction_hgts=np.divide(mp_histo_hgts[0].astype(float),total_hgts)
    
    labels.append(season_name)
    
    xerr_temp.append(2.0)
    yerr_temp.append(1.0/np.sqrt(liquid_histo_temp[0]+np.sqrt(mp_histo_temp[0])))
    x_temp.append(x_axis_temp)
    y_temp.append(mpc_fraction_temp)

    xerr_hgts.append(100.0)
    yerr_hgts.append(1.0/np.sqrt(liquid_histo_hgts[0]+np.sqrt(mp_histo_hgts[0])))
    x_hgts.append(x_axis_hgts)
    y_hgts.append(mpc_fraction_hgts)

    if season_name=="total_year":
        print("High-Temp. Misclassification Rate:", 100.0*mp_histo_temp[0][8]/(mp_histo_temp[0][8]+liquid_histo_temp[0][8]),"%")

multiple_plot2D(folder+'mpc_fraction_temp.png',x=x_temp,y=y_temp,xerr=xerr_temp,yerr=yerr_temp, labels=labels, xlim=(-40,10) , ylim=(0,1), x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Mixed Phase Fraction")
plot2D(folder+'mpc_fraction_temp_total_year.png',x=x_temp[4],y=y_temp[4], xlim=(-40,10) , ylim=(0,1), x_label="Cloud Top Temperature [$\mathsf{{^\circ}C}$]", y_label="Mixed Phase Fraction", fmt='-')
multiple_plot2D(folder+'mpc_fraction_hgts.png',x=x_hgts,y=y_hgts,xerr=xerr_hgts,yerr=yerr_hgts, labels=labels, xlim=(10000,0) , ylim=(0,1), x_label="Cloud Top Height [m]", y_label="Mixed Phase Fraction")

for i in range(len(n_ILCR_t)):
    print( n_ILCR_t[i])
print("Low-Temp. Misclassification Rate:", 100.0*len(CTT[CTT<233.15])/float(len(CTT)),"%")
print("Length of dataset at", station_name,":",len(CTT))
