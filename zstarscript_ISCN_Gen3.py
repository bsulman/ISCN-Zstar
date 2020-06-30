# -*- coding: utf-8 -*-
"""
Script for Z*. Processing ISCN dataset
direct exponential fit to pctC

Created on Fri Jun 19 10:43:46 2015

@author: happysk8er

Edited by Benjamin Sulman Sept 2016
"""
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pylab
from numba import jit

import pandas as pd

import time

@jit
def expfunc(x, K, I):
    '''
    Z* function. not forcing through Csurf.
    @params: K, I
        z*      = -1/K
        csurf   = I
    '''
    return I*np.exp(K*x)

def linfunc(x, a, b):
    '''
    Z* function. not forcing throught Csurf. pass in log trans pctC
    z*      = -1/b
    csurf   = np.exp(a)
    '''
    return a + b*x

@jit
def piecewise_func(z,inv_zstar,csurf,zmax):
    '''
    Z* function. not forcing through Csurf.
    inv_zstar is 1/z*
    Constant below z=zmax
    Cmin = Csurf*exp(-zmax/z*)
    @params: inv_zstar, csurf, zmax
    '''
    out=np.zeros_like(z,dtype=float)
    out[z<=zmax]=csurf*np.exp(-z[z<=zmax]*inv_zstar)
    out[z>zmax]=csurf*np.exp(-zmax*inv_zstar)
    return out


failcodes={-1:'No mineral soil',-2:'Less than 3 layers',-3:'No available layers',-999:'Optimization failed'}

#@jit
def zstarfunc(depth, pctC, Cvalues, fit_type = 'exp', org_C_cutoff=20.0):
    '''
    Pass in observations of each profile, fit func, return yhat (if plott=False), zstar, stat
    parameters:
        fit_type      : String.  Current possibilities: 'exp','lin','piecewise'
                        by default, fit exponential func ('exp')
    output:
        fitted value: 1-d vec
        failure code:
            -1   : no mineral soil
            -2   : layer number < 3, inside zstarfunc
            -3   : no avalable layer, in raw data
            -999 : optimization failed
    '''
    from scipy.optimize import curve_fit
    # print 'fit profile ', profid

    # define failure code
    nomin    = -1
    toofewly = -2
    optifd   = -999


    if min(pctC) >= org_C_cutoff: # no minearl soil layer
        return nomin
    Csurf         = Cvalues[pctC<org_C_cutoff][0] # defined by not used. override with fitted value
    Zsurf         = depth[pctC<org_C_cutoff][0]
    depth_mineral = depth - Zsurf # depth vec starts from Zsurf
    Cmin          = np.nanmin(Cvalues)
    Zmin          = np.nanmax(depth_mineral[Cvalues==Cmin]) # max or min?
    Cdeep         = np.nanmean(Cvalues[depth_mineral>=Zmin])
    if fit_type != 'piecewise':
        idx           = np.logical_and(depth_mineral >= 0, depth_mineral <= Zmin)
    else:      # piecewise fit should calculate its own Zmin
        idx           = (depth_mineral >= 0)
    idx           = np.logical_and(idx, Cvalues>0)
    fitdepth      = depth_mineral[idx]
    fitC          = Cvalues[idx]
    nlayer        = fitdepth.shape[0]


    if nlayer < 3:
        return toofewly

    try:
        if fit_type=='exp':
            popt, pcov = curve_fit(expfunc, fitdepth, fitC, maxfev=500,
                                   p0=(-0.01,fitC[0]))
        elif fit_type == 'lin':
            popt, pcov = curve_fit(linfunc, fitdepth, np.log(fitC), maxfev=500)

        elif fit_type == 'piecewise':
            popt, pcov = curve_fit(piecewise_func, fitdepth, fitC, maxfev=500,
                                   p0=(0.01,fitC[0],max(fitdepth)*0.75))

        else:
            raise ValueError('fit_type %s not implemented'%fit_type)


    except RuntimeError:
        return optifd

    if fit_type == 'exp':
        Csurf       = popt[1]
        zstar       = -1./popt[0]
        yhat        = expfunc(fitdepth, *popt)
        prop = [zstar,Csurf,Cmin,Zmin,Cdeep]

    elif fit_type == 'lin':
        Csurf       = np.exp(popt[0])
        zstar       = -1./popt[1]
        yhat        = np.exp(linfunc(fitdepth, *popt))
        prop = [zstar,Csurf,Cmin,Zmin,Cdeep]

    elif fit_type == 'piecewise':
        Csurf = popt[1]
        zstar = 1./popt[0]
        zmax  = popt[2]
        yhat  = piecewise_func(fitdepth,*popt)
        Cdeep = Csurf*np.exp(-min(zmax,fitdepth.max())/zstar)
        prop = [zstar,Csurf,Cmin,zmax,Cdeep]

    else:
        raise ValueError('fit_type %s not implemented'%fit_type)

    # plt.plot(fitdepth, fitC);plt.plot(yhat, fitdepth)
    z_r2        = np.corrcoef(fitC, yhat)[0,1]**2  #mysm.cal_R2(fitC, yhat)
    z_rmse      = np.sqrt(sum((yhat-fitC)**2))     #mysm.cal_RMSE(fitC, yhat)
    z_pcterr    = np.mean((yhat-fitC)/fitC)        #mysm.cal_pctERR(fitC, yhat)

    return {'prop':prop, 'stat':[z_r2, z_rmse, z_pcterr, nlayer], 'fitting':[yhat,fitC,fitdepth]}

def mainrun(layerdata, profdata, uniqprofname, plott=False, ppf=9,maxprofiles=100000, fit_type = 'exp', do_C_stock=False, org_C_cutoff=20.0):
    '''
    fit Z* function to ISCN dataset
    parameters:
        plott: whether to plot (True) or not (default, False)
        ppf: plots per figure, default 9
        failure code:
            -1   : no mineral soil
            -2   : layer number < 3, inside zstarfunc
            -3   : no avalable layer, in raw data
            -999 : optimization failed
    '''
    norawly     = -3
    out_fitting = []
    out_prop    = []
    out_stat    = []
    failed      = []
    out_method  = []
    out_name    = []

    n_OC=0
    n_Ctot=0
    n_stock=0

    tstart=time.clock()

    failures={-1:'No mineral soil',-2:'Less than 3 layers (in zstarfunc)',-3:'No available layers',-999:'Optimization failed'}

    if fit_type not in ['exp','lin','piecewise']:
        raise ValueError('fit_type %s not implemented.  Use exp, lin, or piecewise'%fit_type)

    for profid, profname in enumerate(uniqprofname[:maxprofiles]):
        fitstart=time.clock()

        layer_bot    = np.array(layerdata.loc[profname]['layer_bot (cm)'],ndmin=1)
        layer_top    = np.array(layerdata.loc[profname]['layer_top (cm)'],ndmin=1)
        rawdepth     = 0.5*(layer_bot+layer_top)
        rawpctCtot  = np.array(layerdata.loc[profname]['c_tot (percent)'],ndmin=1)
        rawpctCoc   = np.array(layerdata.loc[profname]['oc (percent)'],ndmin=1)

        layer_thickness = layer_bot-layer_top
        rawCstock   = np.array(layerdata.loc[profname]['soc (g cm-2)'],ndmin=1)/layer_thickness #gC/m3

        totalprofs=len(uniqprofname[:maxprofiles])

        if len(rawpctCoc[~np.isnan(rawpctCoc)])>2:
            rawpctC = rawpctCoc
            C_measure='OC'
        else:
            rawpctC = rawpctCtot
            C_measure='Ctot'

        if do_C_stock:
            rawCvalues=rawCstock
            notNaNs = ~np.isnan(rawdepth) & ~np.isnan(rawCvalues) & ~np.isnan(rawpctC)
        else:
            notNaNs = ~np.isnan(rawdepth) & ~np.isnan(rawpctC)
            rawCvalues=rawpctC

        depth = rawdepth[notNaNs]; Cvalues = rawCvalues[notNaNs]; pctC = rawpctC[notNaNs]

        if depth.shape[0] > 0:
            res = zstarfunc(depth, pctC, Cvalues, fit_type=fit_type, org_C_cutoff=org_C_cutoff)
            fitend=time.clock()
            if (not isinstance(res, int)) and (not isinstance(res, list)):
                # print 'Fit number %d, profile name %s [%s]. Took %1.1g seconds.'%(profid,profname,C_measure,fitend-fitstart)
                if plott:
                    out_fitting.append(res['fitting'])
                    fign = (profid+1)%ppf
                    if profid%ppf == 0:
                        fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(12,10))
                    ax = fig.axes[fign-1]
                    ax.invert_yaxis()
                    ax.scatter(Cvalues, depth-depth[Cvalues<20][0], label='obs', c='g')
                    ax.plot(res['fitting'][0], res['fitting'][2], label='fitted')
                    if do_C_stock:
                        ax.set_xlabel('SOC (gC m${-3}$) ' + str(profid) + ':(' + profname + ')')
                    elif C_measure=='OC':
                        ax.set_xlabel('OC(%) ' + str(profid) + ':(' + profname + ')')
                    else:
                        ax.set_xlabel('Ctot(%) ' + str(profid) + ':(' + profname + ')')
                    ax.set_ylabel('depth (cm)')
                    if fit_type == 'piecewise':
                        ax.text(0.7,0.1,'z* = %1.1f\nZmax = %1.1f\nr$^2$ = %1.2f'%(res['prop'][1],res['prop'][3],res['stat'][0]),transform=ax.transAxes)
                    else:
                        ax.text(0.7,0.1,'z* = %1.1f\nr$^2$ = %1.2f'%(res['prop'][1],res['stat'][0]),transform=ax.transAxes)
                    plt.legend()
                    #pylab.text(.8, .8, df.loc[profname]['ecoregion'],fontsize=8)
                    if fign == 0:
                        print ('save plot, profid ', profid)
                        plt.tight_layout()
                        if fit_type == 'lin':
                            fig.savefig(pathh + 'figures_ben_gen3/lmfit/fig_%s.png'%(profid))
                        elif fit_type == 'exp':
                            fig.savefig(pathh + 'figures_ben_gen3/expfit/fig_%s.png'%(profid))
                        elif fit_type == 'piecewise':
                            fig.savefig(pathh + 'figures_ben_gen3/piecewisefit/fig_%s.png'%(profid))
                        plt.close()
                out_prop.append(res['prop'])
                out_stat.append(res['stat'])
                if do_C_stock:
                    # Sometimes the top layer did not have data but lower layers did
                    # Fix so it uses layers that had data used for the calculation
                    out_method.append(layerdata.loc[profname].dropna(subset=['soc (g cm-2)'])['soc_carbon_flag'].iloc[0])
                else:
                    out_method.append(C_measure)
                out_name.append(profname)
                if C_measure=='OC':
                    n_OC=n_OC+1
                else:
                    n_Ctot=n_Ctot+1
            if isinstance(res, int):
                failed.append([profname,res,failures[res]])
                # print 'Fit number %d, profile name %s [FAILED: %s]'%(profid,profname,failures[res])
        else:
            failed.append([profname,norawly,failures[norawly]])
            # print 'Fit number %d, profile name %s [FAILED: %s]'%(profid,profname,failures[norawly])
        if profid%100==0 and profid>0:
            tt=time.clock()
            print ('\n')
            print ('*** Summary after %d of %d total profiles:'%(profid,totalprofs))
            print ('*** %d OC, %d Ctot, %d failed. '%(n_OC,n_Ctot,len(failed)))
            time_per_profile=(tt-tstart)/profid
            time_left=time_per_profile*(totalprofs-profid)
            print ('Total time=%1.1f s. Mean time per profile=%1.1g s. Est time to finish: %1.1f minutes'%(tt-tstart,time_per_profile,time_left/60))
            print ('\n')
    if len(out_prop)>0:
        prop_array=np.asarray(out_prop)
        stat_array=np.asarray(out_stat)
    else:
        prop_array=np.zeros((1,5))+np.nan
        stat_array=np.zeros((1,5))+np.nan
    out_df=pd.DataFrame({
        'zstar':prop_array[:,0],
        'Csurf':prop_array[:,1],
        'Cmin' :prop_array[:,2],
        'Zmin' :prop_array[:,3],
        'Cdeep':prop_array[:,4],
        'r2'   :stat_array[:,0],
        'rmse' :stat_array[:,1],
        'pcterr':stat_array[:,2],
        'npoints':stat_array[:,3],
        'Cmeasure':out_method,
        },index=out_name)
    failed_array=np.asarray(failed)
    if len(failed_array)>0:
        return out_df, pd.DataFrame({'failcode':failed_array[:,1],'reason':failed_array[:,2]},index=failed_array[:,0])
    else:
        return out_df,pd.DataFrame({'failcode':failed_array})

#%%
if __name__ == "__main__":

    pathh       = '.'
    #pathh = 'C:\\Users\\happysk8er\\Google Drive\\manuscripts\\C14_synthesis\\JenMeetingApr10_Zstar\\'
    # layerfn     = '/Users/bsulman/Documents/Zstar project/ISCN-data/ISCNLayerData_LATEST.csv'
    # proffn      = '/Users/bsulman/Documents/Zstar project/ISCN-data/ISCNProfileData_LATEST.csv'
    proffn  = '../ISCN-data/ISCN_ALL-DATA_PROFILE_1-1.csv'
    profdata    = pd.read_csv(proffn,encoding='iso-8859-1',index_col='profile_name')

    dataset_all=profdata['dataset_name_soc']
    profdata=profdata[(dataset_all=='ISCN SOC stock computation')|(dataset_all=='ISCN No SOC stock computation')]

    for do_stock in [True,False]:
        fit_type='piecewise'

        for num in [1,2,3,4]:
            layerfn = '../ISCN-data/ISCN_ALL_DATA_LAYER_C%d_1-1.csv'%num
            layerdata   = pd.read_csv(layerfn,encoding='iso-8859-1')
            layerdata.set_index('profile_name',inplace=True)
            dataset_all=layerdata['dataset_name_soc']
            layerdata=layerdata[(dataset_all=='ISCN SOC stock computation')|(dataset_all=='ISCN No SOC stock computation')]
            uniqprofname = layerdata.index.unique()

            print ('Finished reading datasets for segment C%d'%num)
            if do_stock:
                print ("Using C stocks")
            else:
                print ("Using percent C")

            maxruns = 1000000

            # run exp fit
            out_df, failed_exp = mainrun(layerdata, profdata, uniqprofname, plott=False, fit_type=fit_type, maxprofiles=maxruns, do_C_stock=do_stock)
            # np.savez('out_exp_C1',out_stat_exp=out_stat_exp, out_prop_exp=out_prop_exp, failed_exp=failed_exp,out_method=out_method)

            if do_stock:
                out_df.to_csv('../Zstar_params_output/C%d_Cstock_%s.csv'%(num,fit_type),index_label='profile_name')
                failed_exp.to_csv('../Zstar_params_output/C%d_Cstock_%s_failed.csv'%(num,fit_type),index_label='profile_name')
            else:
                out_df.to_csv('../Zstar_params_output/C%d_pctC_%s.csv'%(num,fit_type),index_label='profile_name')
                failed_exp.to_csv('../Zstar_params_output/C%d_pctC_%s_failed.csv'%(num,fit_type),index_label='profile_name')


    #
    # # run linear fit
    # out_fitting_lm, out_prop_lm, \
    #     out_stat_lm, failed_lm = mainrun(layerdata, profdata, uniqprofname, plott=True, fit_type = 'lin', maxprofiles=maxruns)
    # out_stat_lm = np.array(out_stat_lm)
    # out_prop_lm = np.array(out_prop_lm)
    # np.savez('out_lm',out_stat_lm=out_stat_lm, out_prop_lm=out_prop_lm, failed_lm=failed_lm)

    # run piecewise fit
    # out_fitting_piecewise, out_prop_piecewise, \
    #     out_stat_piecewise, failed_piecewise = mainrun(layerdata, profdata, uniqprofname, plott=True, fit_type='piecewise', maxprofiles=maxruns)
    # out_stat_piecewise = np.array(out_stat_piecewise)
    # out_prop_piecewise = np.array(out_prop_piecewise)
    # np.savez('out_piecewise',out_stat_piecewise=out_stat_piecewise, out_prop_piecewise=out_prop_piecewise, failed_piecewise=failed_piecewise)
