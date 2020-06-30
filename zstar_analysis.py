import pandas,numpy
from pylab import *

Cmeasure='Cstock'
Zmax=30
minr2=0.90
min_npoints=5
only_oc=False

soil_order_mappings = {
'alf':'Alfisol',
'and':'Andisol',
'el':'Gelisol',
'ent':'Entisol',
'ept':'Inceptisol',
'ert':'Vertisol',
'id':'Aridisol',
'ist':'Histosol',
'od':'Spodosol',
'oll':'Mollisol',
'ox':'Oxisol',
'ult':'Ultisol',
}

def get_soil_order(soil_name):
    if not isinstance(soil_name,str):
        return 'Unknown'
    elif len(soil_name.split())>1:
        # Multiple words in soil taxon. Try to get order, starting from the end
        for word in soil_name.split()[::-1]:
            order=get_soil_order(word.strip(' ,.;'))
            if order != 'Unknown':
                return order
        return 'Unknown'
    elif soil_name.capitalize() in soil_order_mappings.values():
        return soil_name.capitalize()
    elif soil_name.lower().endswith('s'):
        return get_soil_order(soil_name[:-1])
    else:
        name=soil_name.lower()
        if name[-2:] in ['el','id','od','ox']:
            return soil_order_mappings.get(name[-2:],'Unknown')
        else:
            return soil_order_mappings.get(name[-3:],'Unknown')


NLCD_mappings={
    31:'Barren Land',
    82:'Cultivated Crops',
    41:'Deciduous Forest',
    24:'Developed, High Intensity',
    22:'Developed, Low Intensity',
    23:'Developed, Medium Intensity',
    21:'Developed, Open Space',
    95:'Emergent Herbanceous Wetlands',
    42:'Evergreen Forest',
    71:'Grassland/Herbaceous',
    43:'Mixed Forest',
    11:'Open Water',
    81:'Pasture/Hay',
    12:'Perennial Ice/Snow',
    52:'Shrub/Scrub',
    90:'Woody Wetlands',
    32:'Quarries/Strip Mines/Gravel Pits',
    33:'Transitional',
    61:'Orchards/Vineyards/Other',
    83:'Small Grains',
    84:'Fallow',
    85:'Urban/Recreational Grasses',
    91:'Woody Wetlands',
    92:'Emergent Herbaceous Wetlands',
    51:'Shrubland',
    0:'Unknown',
}



def plot_profile(layerdata,profiledata,profid=None,fig_text=True,clearaxis=True,Cmeasure=Cmeasure,*args,**kwargs):
    if clearaxis:
        cla()
    if profid is None:
        while True:
            profid=profiledata.index[numpy.random.randint(0,len(profiledata))]
            if profid in profiledata.index:
                break
        print ('Plotting profile '+profid)
    profile=layerdata[layerdata.index==profid]
    top=profile['layer_top (cm)']
    bottom=profile['layer_bot (cm)']
    middle=(bottom+top)/2.0
    ctot=profile['c_tot (percent)']
    corg=profile['oc (percent)']
    cstock=profile['soc (g cm-2)']/(bottom-top)
    if corg.isnull().all():
        cpct=ctot
    else:
        cpct=corg

    if Cmeasure is 'Cstock':
        c=cstock
    else:
        c=cpct

    dots=plot(c,-middle,'o',**kwargs)
    kwargs.pop('c',None)
    kwargs.pop('marker',None)
    z=linspace(0,profiledata['Zmin'][profid],100)
    start_depth = middle[cpct<20].iloc[0]
    plot(profiledata['Csurf'][profid]*exp(-z/profiledata['zstar'][profid]),-z-start_depth,'-',c=dots[0].get_color(),**kwargs)
    plot([profiledata['Cdeep'][profid],profiledata['Cdeep'][profid]],[-profiledata['Zmin'][profid]-start_depth,-bottom.max()-start_depth],'--',c=dots[0].get_color(),**kwargs)
    textstr='''Profile {profile}:
               Soil order: {order}
               Latitude: {latitude:1.1f}
               Country: {country}
               zstar = {zstar:1.2f}, Zmin = {Zmin:1.2f}
               Csurf = {Csurf:1.2e}, Cdeep = {Cdeep:1.2e}'''.format(
                                    profile=profid,
                                    # profiledata.landuse[profid],
                                    order=get_soil_order(profiledata.soil_taxon[profid]),
                                    zstar=profiledata.zstar[profid],
                                    Csurf=profiledata.Csurf[profid],Cdeep=profiledata.Cdeep[profid],Zmin=profiledata.Zmin[profid],
                                    latitude=lat[profid],country=profiledata['country (country)'][profid])
    if fig_text:

        text(0.8,0.2,textstr,transform=gca().transAxes,horizontalalignment='right')
    else:
        print (textstr)

    ylabel('Depth (cm)')
    if Cmeasure is "Cstock":
        xlabel('C content (gC cm$^{-3}$)')
    else:
        xlabel('C content (%)')
    return (middle.values,c.values,dots[0])



def profile_totalcarbon(Csurf,Cdeep,Zstar,Zmin,Zmax):
    '''Integrate the C in the profile based on profile stats:
       Csurf and Cdeep in gC cm-3
       Zstar, Zmin, Zmax in cm
       Returns C stock in kgC m-2 to depth Zmax'''
    # To do this with %C we need C stock units
    if Cmeasure != 'Cstock':
        raise RuntimeError('Profiles can only be calculated when using C stocks')
    if isinstance(Zmin,float):
        if Zmax<Zmin:
            totalcarbon= -Csurf*Zstar*(exp(-Zmax/Zstar)-1)
        else:
            totalcarbon= -Csurf*Zstar*(exp(-Zmin/Zstar)-1) + Cdeep*(Zmax-Zmin)
    else:
        totalcarbon = where(Zmin>Zmax,-Csurf*Zstar*(exp(-Zmax/Zstar)-1),-Csurf*Zstar*(exp(-Zmin/Zstar)-1) + Cdeep*(Zmax-Zmin))
        if isinstance(Zmin,pandas.Series):
            totalcarbon=pandas.Series(totalcarbon,index=Zmin.index)
    return totalcarbon*1e4*1e-3


# Read in profile and z* data

zstardata=pandas.concat([pandas.read_csv('../Zstar_params_output/C%d_Cstock_piecewise.csv'%(num),index_col='profile_name',low_memory=False) for num in [1,2,3,4]])
zstardata_percent=pandas.concat([pandas.read_csv('../Zstar_params_output/C%d_pctC_piecewise.csv'%num,index_col='profile_name',low_memory=False) for num in [1,2,3,4]])
profiledata=pandas.read_csv('../ISCN-data/profiledata_gen3_with_climate_and_nlcd.csv',encoding='latin-1',index_col='profile_name',low_memory=False)
profiledata.mask(profiledata==-9999.0,inplace=True)

failed_Cstock=pandas.concat([pandas.read_csv('../Zstar_params_output/C%d_Cstock_piecewise_failed.csv'%(num),index_col='profile_name',low_memory=False) for num in [1,2,3,4]])
failed_pctC=pandas.concat([pandas.read_csv('../Zstar_params_output/C%d_pctC_piecewise_failed.csv'%(num),index_col='profile_name',low_memory=False) for num in [1,2,3,4]])

# Add number of layers to this figure?
fig=figure('Fit stats',figsize=(12,5.5));clf()
fig,axes=subplots(2,5,num='Fit stats');
ax=axes[0,0]
ax.hist(zstardata_percent.r2,arange(0,1.05,0.05))
ax.set_title('R$^2$')
ax.set_ylabel('Soil C concentration\nNumber of profiles')
ax.set_xlabel('R$^2$')
labx=0.02
laby=1.05
ax.text(labx,laby,'(a)',transform=ax.transAxes)
ax=axes[1,0]
ax.hist(zstardata.r2,arange(0,1.05,0.05))
# ax.set_title('Correlation coefficient')
ax.set_ylabel('Soil C density\nNumber of profiles')
ax.set_xlabel('R$^2$')
ax.text(labx,laby,'(b)',transform=ax.transAxes)

ax=axes[0,1]
ax.hist(zstardata_percent.zstar,arange(0,200,10))
# title('C concentration')
ax.set_ylabel('Number of profiles')
ax.set_xlabel('Z* (cm)')
ax.set_title('Z*')
ax.text(labx,laby,'(c)',transform=ax.transAxes)
ax=axes[1,1]
ax.hist(zstardata.zstar,arange(0,200,10))
# title('C density')
ax.set_ylabel('Number of profiles')
ax.set_xlabel('Z* (cm)')
ax.text(labx,laby,'(d)',transform=ax.transAxes)

ax=axes[0,2]
ax.hist(zstardata_percent.Csurf,linspace(0,25,30))
# title('C concentration')
ax.set_ylabel('Number of profiles')
ax.set_xlabel('$C_{s}$ (%)')
ax.set_title('$C_{s}$')
ax.text(labx,laby,'(e)',transform=ax.transAxes)
ax=axes[1,2]
ax.hist(zstardata.Csurf,linspace(0,0.3,30))
# title('C density')
ax.set_ylabel('Number of profiles')
ax.set_xlabel('$C_{s}$ (gC cm$^{-3}$)')
ax.text(labx,laby,'(f)',transform=ax.transAxes)

ax=axes[0,3]
ax.hist(zstardata_percent.Zmin,linspace(0,300,30))
# title('C concentration')
ax.set_ylabel('Number of profiles')
ax.set_xlabel('$Z_{min}$ (cm)')
ax.set_title('$Z_{min}$')
ax.text(labx,laby,'(g)',transform=ax.transAxes)
ax=axes[1,3]
ax.hist(zstardata.Zmin,linspace(0,300,30))
# title('C density')
ax.set_ylabel('Number of profiles')
ax.set_xlabel('$Z_{min}$ (cm)')
ax.text(labx,laby,'(h)',transform=ax.transAxes)

ax=axes[0,4]
ax.hist(zstardata_percent.Cdeep,linspace(0,4.5,30))
# title('C concentration')
ax.set_ylabel('Number of profiles')
ax.set_xlabel('$C_{deep}$ (%)')
ax.set_title('$C_{deep}$')
ax.text(labx,laby,'(i)',transform=ax.transAxes)
ax=axes[1,4]
ax.hist(zstardata.Cdeep,linspace(0,0.07,30))
# title('C density')
ax.set_ylabel('Number of profiles')
ax.set_xlabel('$C_{deep}$ (gC cm$^{-3}$)')
ax.text(labx,laby,'(j)',transform=ax.transAxes)


# tight_layout()

# Merge into one dataset with profile properties and calculated stats
alldata=pandas.merge(zstardata,profiledata,how='inner',left_on='profile_name',right_on='profile_name')
alldata_pctC=pandas.merge(zstardata_percent,profiledata,how='inner',left_on='profile_name',right_on='profile_name')

max_zstar=200
min_zstar=0
lon=alldata['long (dec. deg)']
lat=alldata['lat (dec. deg)']
good_zstar=abs(lon<360)&abs(lat<=90)&(alldata['zstar']>min_zstar)
if only_oc:
    good_zstar=good_zstar&~alldata.Cmeasure.str.contains('c_tot')

alldata_good=alldata[good_zstar].copy()

# Make a subset of just profiles with Z* fit r2 greater than a minimum value

alldata_goodr2=alldata_good[(alldata_good['r2']>minr2)&(alldata_good['npoints']>=min_npoints)].copy()

nprofs=len(profiledata)
nolayers_pctC=len(failed_pctC[(failed_pctC.failcode==-1)|(failed_pctC.failcode==-3)])  
nolayers_Cstock=len(failed_Cstock[(failed_Cstock.failcode==-1)|(failed_Cstock.failcode==-3)])  
lt3layers_pctC=len(failed_pctC[failed_pctC.failcode==-2])
lt3layers_Cstock=len(failed_Cstock[failed_Cstock.failcode==-2])
lt5layers_pctC=lt3layers_pctC+len(zstardata_percent[zstardata_percent.npoints<5])
lt5layers_Cstock=lt3layers_Cstock+len(zstardata[zstardata.npoints<5])
badfit_pctC=len(zstardata_percent[zstardata_percent.zstar<0])+len(failed_pctC[failed_pctC.failcode==-999])
badfit_Cstock=len(zstardata[zstardata.zstar<0])+len(failed_Cstock[failed_Cstock.failcode==-999])
noloc_pctC=sum((alldata_pctC['long (dec. deg)']>=360)|(alldata_pctC['lat (dec. deg)']>90))
noloc_Cstock=sum((lon>=360)|(lat>90))
lowr2_pctC=len(alldata_pctC[(alldata_pctC['r2']<0.9)&(alldata_pctC['npoints']>=5)&(alldata_pctC['zstar']>0)])
lowr2_Cstock=len(alldata[(alldata['r2']<0.9)&(alldata['npoints']>=5)&(alldata['zstar']>0)])
r2_75_pctC=len(alldata_pctC[(alldata_pctC['r2']<0.75)&(alldata_pctC['npoints']>=5)&(alldata_pctC['zstar']>0)])
r2_75_Cstock=len(alldata[(alldata['r2']<0.75)&(alldata['npoints']>=5)&(alldata['zstar']>0)])
ctot_pctC=len(alldata_pctC[(alldata_pctC['r2']>=0.9)&(alldata_pctC['npoints']>=5)&(alldata_pctC['zstar']>0)&(alldata_pctC.Cmeasure.str.contains('Ctot'))])
ctot_Cstock=len(alldata[(alldata['r2']>=0.9)&(alldata['npoints']>=5)&(alldata['zstar']>0)&(alldata.Cmeasure.str.contains('c_tot'))])

reject_table=pandas.DataFrame(columns=['Total','At least one layer','At least 3 layers','At least 5 layers','Successful fit','R2>0.75','R2>0.9','Only OC'],index=pandas.MultiIndex.from_product([['Pct C','C stock'],['Remaining','Rejected','Rejected (% of total)']]),dtype=int,data=-999)
reject_table['Total'].loc[:,'Rejected']=0
reject_table['Total'].loc[:,'Remaining']=nprofs
reject_table['At least one layer'].loc['Pct C','Rejected']=nolayers_pctC
reject_table['At least 3 layers'].loc['Pct C','Rejected']=lt3layers_pctC
reject_table['At least 5 layers'].loc['Pct C','Rejected']=lt5layers_pctC-lt3layers_pctC
reject_table['Successful fit'].loc['Pct C','Rejected']=badfit_pctC
reject_table['R2>0.9'].loc['Pct C','Rejected']=lowr2_pctC-r2_75_pctC
reject_table['R2>0.75'].loc['Pct C','Rejected']=r2_75_pctC
reject_table['Only OC'].loc['Pct C','Rejected']=ctot_pctC

reject_table['At least one layer'].loc['C stock','Rejected']=nolayers_Cstock
reject_table['At least 3 layers'].loc['C stock','Rejected']=lt3layers_Cstock
reject_table['At least 5 layers'].loc['C stock','Rejected']=lt5layers_Cstock-lt3layers_Cstock
reject_table['Successful fit'].loc['C stock','Rejected']=badfit_Cstock
reject_table['R2>0.9'].loc['C stock','Rejected']=lowr2_Cstock-r2_75_Cstock
reject_table['R2>0.75'].loc['C stock','Rejected']=r2_75_Cstock
reject_table['Only OC'].loc['C stock','Rejected']=ctot_Cstock

reject_table['At least one layer'].loc['Pct C','Remaining']=nprofs-nolayers_pctC
reject_table['At least 3 layers'].loc['Pct C','Remaining']=nprofs-nolayers_pctC-lt3layers_pctC
reject_table['At least 5 layers'].loc['Pct C','Remaining']=lt5layers=nprofs-nolayers_pctC-lt5layers_pctC
reject_table['Successful fit'].loc['Pct C','Remaining']=nprofs-nolayers_pctC-lt5layers_pctC-badfit_pctC
reject_table['R2>0.9'].loc['Pct C','Remaining']=nprofs-nolayers_pctC-lt5layers_pctC-badfit_pctC-lowr2_pctC
reject_table['R2>0.75'].loc['Pct C','Remaining']=nprofs-nolayers_pctC-lt5layers_pctC-badfit_pctC-r2_75_pctC
reject_table['Only OC'].loc['Pct C','Remaining']=nprofs-nolayers_pctC-lt5layers_pctC-badfit_pctC-lowr2_pctC-ctot_pctC

reject_table['At least one layer'].loc['C stock','Remaining']=nprofs-nolayers_Cstock
reject_table['At least 3 layers'].loc['C stock','Remaining']=nprofs-nolayers_Cstock-lt3layers_Cstock
reject_table['At least 5 layers'].loc['C stock','Remaining']=lt5layers=nprofs-nolayers_Cstock-lt5layers_Cstock
reject_table['Successful fit'].loc['C stock','Remaining']=nprofs-nolayers_Cstock-lt5layers_Cstock-badfit_Cstock
reject_table['R2>0.9'].loc['C stock','Remaining']=nprofs-nolayers_Cstock-lt5layers_Cstock-badfit_Cstock-lowr2_Cstock
reject_table['R2>0.75'].loc['C stock','Remaining']=nprofs-nolayers_Cstock-lt5layers_Cstock-badfit_Cstock-r2_75_Cstock
reject_table['Only OC'].loc['C stock','Remaining']=nprofs-nolayers_Cstock-lt5layers_Cstock-badfit_Cstock-lowr2_Cstock-ctot_Cstock

reject_table.loc[('Pct C','Rejected (% of total)')]=reject_table.loc[('Pct C','Rejected')]/reject_table.loc[('Pct C','Remaining')]['Total']*100
reject_table.loc[('C stock','Rejected (% of total)')]=reject_table.loc[('C stock','Rejected')]/reject_table.loc[('C stock','Remaining')]['Total']*100


print(reject_table.round(1))


zstar=alldata_goodr2['zstar'].copy()
# Z* sometimes gets very large values for profiles that don't decline much with depth, so limit it to a maximum value
zstar[zstar>max_zstar]=max_zstar

in_US=alldata_goodr2['country (country)']=='United States'
data_US=alldata_goodr2[in_US]

# Read in data of individual soil layers from four files and combine into one dataset
layerdata_list=[]
for num in [1,2,3,4]:
    layerfn = '../ISCN-data/ISCN_ALL_DATA_LAYER_C%d_1-1.csv'%num
    layerdata   = pandas.read_csv(layerfn,encoding='iso-8859-1',low_memory=False)
    layerdata.set_index('profile_name',inplace=True)
    dataset_all=layerdata['dataset_name_soc']
    layerdata_list.append(layerdata[(dataset_all=='ISCN SOC stock computation')|(dataset_all=='ISCN No SOC stock computation')])

layerdata=pandas.concat(layerdata_list)

# Make a map of where profiles are
soilorders=alldata_goodr2['soil_taxon'].apply(get_soil_order)
lat_good=alldata_goodr2['lat (dec. deg)']
lon_good=alldata_goodr2['long (dec. deg)']

alldata_goodr2['soilorder']=soilorders

import cartopy

# Colors for soil orders
cols = {
'Alfisol':'C0',
'Andisol':'C1',
'Gelisol':'C2',
'Entisol':'C3',
'Inceptisol':'C4',
'Vertisol':'C5',
'Aridisol':'C6',
'Histosol':'C7',
'Spodosol':'C8',
'Mollisol':'C9',
'Oxisol':'#F0E68C',
'Ultisol':'#4B0082',
'Unknown':'k'
}

zstar_vmax=150
zmin_vmax=250
if Cmeasure is "Cstock":
    Cunits='gC cm$^{-3}$'
    C_vmax=0.3
    Cdeep_vmax=0.03
    Cstocks=profile_totalcarbon(alldata_goodr2.Csurf,alldata_goodr2.Cdeep,alldata_goodr2.zstar,alldata_goodr2.Zmin,Zmax)
    alldata_goodr2['Cstocks_zstar_%dcm'%Zmax]=Cstocks
    Cstocks_vmax=15.0
else:
    Cunits='%'
    C_vmax=20
    Cdeep_vmax=4.0

figure('Profiles map global',figsize=(11.8,6.5));clf()
set_cmap('YlOrBr')
ax=subplot(111,projection=cartopy.crs.PlateCarree())
ax.coastlines(lw=0.5,color='gray')

orders_sorted=sorted(soilorders.unique())
orders_sorted.remove('Unknown') 
orders_sorted.append('Unknown') 
for order in orders_sorted:
    xx=soilorders==order
    ax.plot(lon_good[xx],lat_good[xx],'o',label=order,ms=1.5,alpha=0.7,mec='None',c=cols[order])

legend(markerscale=3,bbox_to_anchor=(-0.01,1.0),loc='upper right')
ax.set_title('Soil orders')
# 
# ax=subplot(322,projection=cartopy.crs.PlateCarree())
# ax.coastlines(lw=0.5,color='gray')
# h=ax.scatter(lon_good,lat_good,label=order,s=3,alpha=0.7,edgecolor='None',c=alldata_goodr2['zstar'],vmax=zstar_vmax)
# ax.set_title('Z*')
# cb=colorbar(h,extend='max');cb.set_label('Z* (cm)')
# 
# ax=subplot(324,projection=cartopy.crs.PlateCarree())
# ax.coastlines(lw=0.5,color='gray')
# h=ax.scatter(lon_good,lat_good,label=order,s=3,alpha=0.7,edgecolor='None',c=alldata_goodr2['Csurf'],vmax=C_vmax)
# ax.set_title('Csurf')
# cb=colorbar(h,extend='max');cb.set_label('Csurf (%s)'%Cunits)
# 
# ax=subplot(326,projection=cartopy.crs.PlateCarree())
# ax.coastlines(lw=0.5,color='gray')
# h=ax.scatter(lon_good,lat_good,label=order,s=3,alpha=0.7,edgecolor='None',c=alldata_goodr2['Cdeep'],vmax=Cdeep_vmax)
# ax.set_title('Cdeep')
# cb=colorbar(h,extend='max');cb.set_label('Cdeep (%s)'%Cunits)
# 
# if Cmeasure=='Cstock':
#     ax=subplot(323,projection=cartopy.crs.PlateCarree())
#     ax.coastlines(lw=0.5,color='gray')
#     h=ax.scatter(lon_good,lat_good,label=order,s=3,alpha=0.7,edgecolor='None',c=Cstocks,vmax=Cstocks_vmax)
#     ax.set_title('C stocks to %d cm'%Zmax)
#     cb=colorbar(h,extend='max');cb.set_label('C stock (%s)'%'kgC m$^{-2}$')
# 
# ax=subplot(325,projection=cartopy.crs.PlateCarree())
# ax.coastlines(lw=0.5,color='gray')
# h=ax.scatter(lon_good,lat_good,label=order,s=3,alpha=0.7,edgecolor='None',c=alldata_goodr2['zstar'],vmax=zmin_vmax)
# ax.set_title('Zmin')
# cb=colorbar(h,extend='max');cb.set_label('Zmin (cm)')

# subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.8,hspace=0.1,wspace=0.1)

# Same maps, but zoomed to US

f,axs=subplots(3,2,num='Profiles map US',figsize=(11.8,6.5),clear=True,subplot_kw=dict(projection=cartopy.crs.PlateCarree()))
ax=axs[0,0]
ax.add_feature(cartopy.feature.BORDERS,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.COASTLINE,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.STATES,lw=0.25,facecolor='none',edgecolor='gray')

for order in unique(soilorders):
    xx=soilorders==order
    ax.plot(lon_good[xx],lat_good[xx],'o',label=order,ms=1.5,alpha=0.7,mec='None',c=cols[order])

ax.legend(ncol=1,loc=(-.2,0),markerscale=3,fontsize='small')
ax.set_title('Soil orders')
ax.set_xlim(-129,-65)
ax.set_ylim(23,52)
labx=0.02
laby=0.90
ax.text(labx,laby,'(a)',transform=ax.transAxes)


# ax=subplot(322,projection=cartopy.crs.PlateCarree())
ax=axs[0,1]
ax.add_feature(cartopy.feature.BORDERS,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.COASTLINE,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.STATES,lw=0.25,facecolor=[0.7,0.7,0.7],edgecolor='gray')
h=ax.scatter(lon_good,lat_good,label=order,s=3,alpha=0.7,edgecolor='None',c=alldata_goodr2['zstar'],vmax=zstar_vmax,zorder=10)
ax.set_title('Z*')
cb=colorbar(h,extend='max',ax=ax);cb.set_label('Z* (cm)')
ax.set_xlim(-129,-65)
ax.set_ylim(23,52)
ax.text(labx,laby,'(b)',transform=ax.transAxes)

# ax=subplot(324,projection=cartopy.crs.PlateCarree())
ax=axs[1,1]
ax.add_feature(cartopy.feature.BORDERS,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.COASTLINE,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.STATES,lw=0.25,facecolor=[0.7,0.7,0.7],edgecolor='gray')
h=ax.scatter(lon_good,lat_good,label=order,s=3,alpha=0.7,edgecolor='None',c=alldata_goodr2['Csurf'],vmax=C_vmax,zorder=10)
ax.set_title('$C_{surf}$')
cb=colorbar(h,extend='max',ax=ax);cb.set_label('$C_{surf}$ (%s)'%Cunits)
ax.set_xlim(-129,-65)
ax.set_ylim(23,52)
ax.text(labx,laby,'(d)',transform=ax.transAxes)

# ax=subplot(326,projection=cartopy.crs.PlateCarree())
ax=axs[2,1]
ax.add_feature(cartopy.feature.BORDERS,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.COASTLINE,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.STATES,lw=0.25,facecolor=[0.7,0.7,0.7],edgecolor='gray')
h=ax.scatter(lon_good,lat_good,label=order,s=3,alpha=0.7,edgecolor='None',c=alldata_goodr2['Cdeep'],vmax=Cdeep_vmax,zorder=10)
ax.set_title('$C_{deep}$')
cb=colorbar(h,extend='max',ax=ax);cb.set_label('$C_{deep}$ (%s)'%Cunits)
ax.set_xlim(-129,-65)
ax.set_ylim(23,52)
ax.text(labx,laby,'(f)',transform=ax.transAxes)

if Cmeasure=='Cstock':
    # ax=subplot(323,projection=cartopy.crs.PlateCarree())
    ax=axs[1,0]
    ax.add_feature(cartopy.feature.BORDERS,lw=0.25,facecolor='none',edgecolor='gray')
    ax.add_feature(cartopy.feature.COASTLINE,lw=0.25,facecolor='none',edgecolor='gray')
    ax.add_feature(cartopy.feature.STATES,lw=0.25,facecolor=[0.7,0.7,0.7],edgecolor='gray')
    h=ax.scatter(lon_good,lat_good,label=order,s=3,alpha=0.7,edgecolor='None',c=Cstocks,vmax=Cstocks_vmax,zorder=10)
    ax.set_title('Soil C stocks to %d cm'%Zmax)
    cb=colorbar(h,extend='max',ax=ax);cb.set_label('Soil C stock (%s)'%'kgC m$^{-2}$')
    ax.set_xlim(-129,-65)
    ax.set_ylim(23,52)
    ax.text(labx,laby,'(c)',transform=ax.transAxes)

# ax=subplot(325,projection=cartopy.crs.PlateCarree())
ax=axs[2,0]
ax.add_feature(cartopy.feature.BORDERS,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.COASTLINE,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.STATES,lw=0.25,facecolor=[0.7,0.7,0.7],edgecolor='gray')
h=ax.scatter(lon_good,lat_good,label=order,s=3,alpha=0.7,edgecolor='None',c=alldata_goodr2['Zmin'],vmax=zmin_vmax,zorder=10)
ax.set_title('$Z_{min}$')
cb=colorbar(h,extend='max',ax=ax);cb.set_label('$Z_{min}$ (cm)')
ax.set_xlim(-129,-65)
ax.set_ylim(23,52)
ax.text(labx,laby,'(e)',transform=ax.transAxes)


# subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.8,hspace=0.1,wspace=0.1)


has_Ap=layerdata.hzn_desgn.str.upper().str.contains('AP').groupby(layerdata.index).any()
max_meas_depth=layerdata['layer_bot (cm)'].groupby(layerdata.index).max()
meanclay=layerdata['clay_tot_psa (percent)'].groupby(layerdata.index).mean()

alldata_goodr2['has_Ap']=has_Ap[alldata_goodr2.index]
alldata_goodr2['max_meas_depth']=max_meas_depth[alldata_goodr2.index]
alldata_goodr2['meanclay']=meanclay[alldata_goodr2.index]

f,axs=subplots(5,2,num='Ap effect maps',figsize=(9,9),clear=True,subplot_kw=dict(projection=cartopy.crs.PlateCarree()))
# ax=subplot(521,projection=cartopy.crs.PlateCarree())
ax=axs[0,0]
ax.add_feature(cartopy.feature.BORDERS,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.COASTLINE,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.STATES,lw=0.25,facecolor='none',edgecolor='gray')
h=ax.scatter(lon_good,lat_good,label=order,s=3,alpha=0.7,edgecolor='None',c=alldata_goodr2['zstar'].mask(has_Ap),vmin=0,vmax=200)
ax.set_title('Z* (no Ap)')
cb=colorbar(h,extend='max',ax=ax);cb.set_label('Z* (cm)')
ax.set_xlim(-129,-65)
ax.set_ylim(23,52)

# ax=subplot(522,projection=cartopy.crs.PlateCarree())
ax=axs[0,1]
ax.add_feature(cartopy.feature.BORDERS,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.COASTLINE,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.STATES,lw=0.25,facecolor='none',edgecolor='gray')
h=ax.scatter(lon_good,lat_good,label=order,s=3,alpha=0.7,edgecolor='None',c=alldata_goodr2['zstar'].mask(~has_Ap),vmin=0,vmax=200)
ax.set_title('Z* (has Ap)')
cb=colorbar(h,extend='max',ax=ax);cb.set_label('Z* (cm)')
ax.set_xlim(-129,-65)
ax.set_ylim(23,52)

# ax=subplot(523,projection=cartopy.crs.PlateCarree())
ax=axs[1,0]
ax.add_feature(cartopy.feature.BORDERS,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.COASTLINE,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.STATES,lw=0.25,facecolor='none',edgecolor='gray')
h=ax.scatter(lon_good,lat_good,label=order,s=3,alpha=0.7,edgecolor='None',c=alldata_goodr2['Zmin'].mask(has_Ap),vmin=0,vmax=250)
ax.set_title('Zmin (no Ap)')
cb=colorbar(h,extend='max',ax=ax);cb.set_label('Zmin (cm)')
ax.set_xlim(-129,-65)
ax.set_ylim(23,52)

# ax=subplot(524,projection=cartopy.crs.PlateCarree())
ax=axs[1,1]
ax.add_feature(cartopy.feature.BORDERS,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.COASTLINE,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.STATES,lw=0.25,facecolor='none',edgecolor='gray')
h=ax.scatter(lon_good,lat_good,label=order,s=3,alpha=0.7,edgecolor='None',c=alldata_goodr2['Zmin'].mask(~has_Ap),vmin=0,vmax=250)
ax.set_title('Zmin (has Ap)')
cb=colorbar(h,extend='max',ax=ax);cb.set_label('Zmin (cm)')
ax.set_xlim(-129,-65)
ax.set_ylim(23,52)

# ax=subplot(525,projection=cartopy.crs.PlateCarree())
ax=axs[2,0]
ax.add_feature(cartopy.feature.BORDERS,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.COASTLINE,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.STATES,lw=0.25,facecolor='none',edgecolor='gray')
h=ax.scatter(lon_good,lat_good,label=order,s=3,alpha=0.7,edgecolor='None',c=alldata_goodr2['Csurf'].mask(has_Ap),vmin=0,vmax=C_vmax)
ax.set_title('Csurf (no Ap)')
cb=colorbar(h,extend='max',ax=ax);cb.set_label('Csurf (%s)'%Cunits)
ax.set_xlim(-129,-65)
ax.set_ylim(23,52)

# ax=subplot(526,projection=cartopy.crs.PlateCarree())
ax=axs[2,1]
ax.add_feature(cartopy.feature.BORDERS,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.COASTLINE,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.STATES,lw=0.25,facecolor='none',edgecolor='gray')
h=ax.scatter(lon_good,lat_good,label=order,s=3,alpha=0.7,edgecolor='None',c=alldata_goodr2['Csurf'].mask(~has_Ap),vmin=0,vmax=C_vmax)
ax.set_title('Csurf (has Ap)')
cb=colorbar(h,extend='max',ax=ax);cb.set_label('Csurf (%s)'%Cunits)
ax.set_xlim(-129,-65)
ax.set_ylim(23,52)

# ax=subplot(527,projection=cartopy.crs.PlateCarree())
ax=axs[3,0]
ax.add_feature(cartopy.feature.BORDERS,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.COASTLINE,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.STATES,lw=0.25,facecolor='none',edgecolor='gray')
h=ax.scatter(lon_good,lat_good,label=order,s=3,alpha=0.7,edgecolor='None',c=alldata_goodr2['Cdeep'].mask(has_Ap),vmax=Cdeep_vmax)
ax.set_title('Cdeep (no Ap)')
cb=colorbar(h,extend='max',ax=ax);cb.set_label('Cdeep (%s)'%Cunits)
ax.set_xlim(-129,-65)
ax.set_ylim(23,52)

# ax=subplot(528,projection=cartopy.crs.PlateCarree())
ax=axs[3,1]
ax.add_feature(cartopy.feature.BORDERS,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.COASTLINE,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.STATES,lw=0.25,facecolor='none',edgecolor='gray')
h=ax.scatter(lon_good,lat_good,label=order,s=3,alpha=0.7,edgecolor='None',c=alldata_goodr2['Cdeep'].mask(~has_Ap),vmax=Cdeep_vmax)
ax.set_title('Cdeep (has Ap)')
cb=colorbar(h,extend='max',ax=ax);cb.set_label('Cdeep (%s)'%Cunits)
ax.set_xlim(-129,-65)
ax.set_ylim(23,52)

if Cmeasure=='Cstock':
    # ax=subplot(529,projection=cartopy.crs.PlateCarree())
    ax=axs[4,0]
    ax.add_feature(cartopy.feature.BORDERS,lw=0.25,facecolor='none',edgecolor='gray')
    ax.add_feature(cartopy.feature.COASTLINE,lw=0.25,facecolor='none',edgecolor='gray')
    ax.add_feature(cartopy.feature.STATES,lw=0.25,facecolor='none',edgecolor='gray')
    h=ax.scatter(lon_good,lat_good,label=order,s=3,alpha=0.7,edgecolor='None',c=Cstocks.mask(has_Ap),vmax=Cstocks_vmax)
    ax.set_title('C stocks to %d cm (no Ap)'%Zmax)
    cb=colorbar(h,extend='max',ax=ax);cb.set_label('C stock (%s)'%'kg C m$^{-2}$')
    ax.set_xlim(-129,-65)
    ax.set_ylim(23,52)

    # ax=subplot(5,2,10,projection=cartopy.crs.PlateCarree())
    ax=axs[4,1]
    ax.add_feature(cartopy.feature.BORDERS,lw=0.25,facecolor='none',edgecolor='gray')
    ax.add_feature(cartopy.feature.COASTLINE,lw=0.25,facecolor='none',edgecolor='gray')
    ax.add_feature(cartopy.feature.STATES,lw=0.25,facecolor='none',edgecolor='gray')
    h=ax.scatter(lon_good,lat_good,label=order,s=3,alpha=0.7,edgecolor='None',c=Cstocks.mask(~has_Ap),vmax=Cstocks_vmax)
    ax.set_title('C stocks to %d cm (has Ap)'%Zmax)
    cb=colorbar(h,extend='max',ax=ax);cb.set_label('C stock (%s)'%'kg C m$^{-2}$')
    ax.set_xlim(-129,-65)
    ax.set_ylim(23,52)

# subplots_adjust(hspace=0.2,wspace=0.1,bottom=0.07,top=0.95,left=0.05,right=0.97)


f,axs=subplots(3,2,num='Ap horizon histograms',figsize=(8,8),clear=True)
labx=0.02
laby=0.93
bins=linspace(0,150,30)
sca(axs[0,0])
plot(bins[:-1]+bins[1]/2,histogram(alldata_goodr2['zstar'][has_Ap],density=True,bins=bins)[0]*bins[1],'r-',label='Has Ap')
plot(bins[:-1]+bins[1]/2,histogram(alldata_goodr2['zstar'][~has_Ap],density=True,bins=bins)[0]*bins[1],'b-',label='No Ap')
plot(bins[:-1]+bins[1]/2,histogram(alldata_goodr2['zstar'][has_Ap&(soilorders=='Mollisol')],density=True,bins=bins)[0]*bins[1],'r--',label='Has Ap, Mollisol')
plot(bins[:-1]+bins[1]/2,histogram(alldata_goodr2['zstar'][(~has_Ap)&(soilorders=='Mollisol')],density=True,bins=bins)[0]*bins[1],'b--',label='No Ap, Mollisol')
legend(fontsize='small')
title('Z*')
xlabel('Z* values (cm)')
ylabel('Distribution')
text(labx,laby,'(a)',transform=gca().transAxes)

bins=linspace(0,C_vmax,30)
sca(axs[0,1])
plot(bins[:-1]+bins[1]/2,histogram(alldata_goodr2['Csurf'][has_Ap],density=True,bins=bins)[0]*bins[1],'r-',label='Has Ap')
plot(bins[:-1]+bins[1]/2,histogram(alldata_goodr2['Csurf'][~has_Ap],density=True,bins=bins)[0]*bins[1],'b-',label='No Ap')
plot(bins[:-1]+bins[1]/2,histogram(alldata_goodr2['Csurf'][has_Ap&(soilorders=='Mollisol')],density=True,bins=bins)[0]*bins[1],'r--')
plot(bins[:-1]+bins[1]/2,histogram(alldata_goodr2['Csurf'][(~has_Ap)&(soilorders=='Mollisol')],density=True,bins=bins)[0]*bins[1],'b--')
# legend()
title('$C_{surf}$')
xlabel('$C_{surf}$ values (%s)'%Cunits)
ylabel('Distribution')
text(labx,laby,'(b)',transform=gca().transAxes)

bins=linspace(0,Cdeep_vmax,30)
# subplot(324)
sca(axs[1,1])
plot(bins[:-1]+bins[1]/2,histogram(alldata_goodr2['Cdeep'][has_Ap],density=True,bins=bins)[0]*bins[1],'r-',label='Has Ap')
plot(bins[:-1]+bins[1]/2,histogram(alldata_goodr2['Cdeep'][~has_Ap],density=True,bins=bins)[0]*bins[1],'b-',label='No Ap')
plot(bins[:-1]+bins[1]/2,histogram(alldata_goodr2['Cdeep'][has_Ap&(soilorders=='Mollisol')],density=True,bins=bins)[0]*bins[1],'r--')
plot(bins[:-1]+bins[1]/2,histogram(alldata_goodr2['Cdeep'][(~has_Ap)&(soilorders=='Mollisol')],density=True,bins=bins)[0]*bins[1],'b--')
# legend()
title('$C_{deep}$')
xlabel('$C_{deep}$ values (%s)'%Cunits)
ylabel('Distribution')
text(labx,laby,'(d)',transform=gca().transAxes)

bins=linspace(0,250,30)
# subplot(323)
sca(axs[1,0])
plot(bins[:-1]+bins[1]/2,histogram(alldata_goodr2['Zmin'][has_Ap],density=True,bins=bins)[0]*bins[1],'r-',label='Has Ap')
plot(bins[:-1]+bins[1]/2,histogram(alldata_goodr2['Zmin'][~has_Ap],density=True,bins=bins)[0]*bins[1],'b-',label='No Ap')
plot(bins[:-1]+bins[1]/2,histogram(alldata_goodr2['Zmin'][has_Ap&(soilorders=='Mollisol')],density=True,bins=bins)[0]*bins[1],'r--',label='Has Ap, Mollisol')
plot(bins[:-1]+bins[1]/2,histogram(alldata_goodr2['Zmin'][(~has_Ap)&(soilorders=='Mollisol')],density=True,bins=bins)[0]*bins[1],'b--',label='No Ap, Mollisol')
title('$Z_{min}$')
xlabel('$Z_{min}$ values (cm)')
ylabel('Distribution')
text(labx,laby,'(c)',transform=gca().transAxes)

if Cmeasure=='Cstock':
    bins=linspace(0,Cstocks_vmax,30)
    # subplot(325)
    sca(axs[2,0])
    plot(bins[:-1]+bins[1]/2,histogram(Cstocks[has_Ap],density=True,bins=bins)[0]*bins[1],'r-',label='Has Ap')
    plot(bins[:-1]+bins[1]/2,histogram(Cstocks[~has_Ap],density=True,bins=bins)[0]*bins[1],'b-',label='No Ap')
    plot(bins[:-1]+bins[1]/2,histogram(Cstocks[has_Ap&(soilorders=='Mollisol')],density=True,bins=bins)[0]*bins[1],'r--')
    plot(bins[:-1]+bins[1]/2,histogram(Cstocks[(~has_Ap)&(soilorders=='Mollisol')],density=True,bins=bins)[0]*bins[1],'b--')
    # legend()
    title('Soil C stocks to %d cm'%Zmax)
    xlabel('Soil C stock values (%s)'%'kgC m$^{-2}$')
    ylabel('Distribution')
    text(labx,laby,'(e)',transform=gca().transAxes)

axs[2,1].set_visible(False)
# tight_layout()

# Try to sort into nearby profiles

Cstock_df=pandas.DataFrame({'Cstock':Cstocks,'lat':lat_good,'lon':lon_good,'has_Ap':has_Ap[Cstocks.index],'soilorder':soilorders})

grid_spacing=0.5 # Degrees
lon_bins=arange(-180,180,grid_spacing)
lat_bins=arange(-90,90,grid_spacing)
lat_cut=pandas.cut(lat_good,lat_bins)
lon_cut=pandas.cut(lon_good,lon_bins)

# Categorize into groups for each 2-dimensional bin.
Cstocks_grouped=Cstock_df.groupby([lon_cut,lat_cut,'soilorder'])

# Function to calculate mean C stock in a group for either with or without Ap horizon
def Ap_mean(group,hasAp,col='Cstock'):
    vals=group[col]
    has_Ap=group['has_Ap']
    return vals[has_Ap==hasAp].mean()

def Ap_mean_stats(group,hasAp,col='Cstock'):
    vals=alldata_goodr2[col][group.index]
    has_Ap=group['has_Ap']
    return vals[has_Ap==hasAp].mean()

Ap_comparison=pandas.DataFrame(
    {'Cstock_Ap':Cstocks_grouped.apply(Ap_mean,hasAp=True),'Cstock_noAp':Cstocks_grouped.apply(Ap_mean,hasAp=False)}
    )

Ap_comparison['zstar_Ap']=Cstocks_grouped.apply(Ap_mean_stats,hasAp=True,col='zstar')
Ap_comparison['Zmin_Ap']=Cstocks_grouped.apply(Ap_mean_stats,hasAp=True,col='Zmin')
Ap_comparison['Csurf_Ap']=Cstocks_grouped.apply(Ap_mean_stats,hasAp=True,col='Csurf')
Ap_comparison['Cdeep_Ap']=Cstocks_grouped.apply(Ap_mean_stats,hasAp=True,col='Cdeep')
Ap_comparison['zstar_noAp']=Cstocks_grouped.apply(Ap_mean_stats,hasAp=False,col='zstar')
Ap_comparison['Zmin_noAp']=Cstocks_grouped.apply(Ap_mean_stats,hasAp=False,col='Zmin')
Ap_comparison['Csurf_noAp']=Cstocks_grouped.apply(Ap_mean_stats,hasAp=False,col='Csurf')
Ap_comparison['Cdeep_noAp']=Cstocks_grouped.apply(Ap_mean_stats,hasAp=False,col='Cdeep')

Ap_comp_reindexed=Ap_comparison.reset_index().rename(columns={'long (dec. deg)':'lon_interval','lat (dec. deg)':'lat_interval'}).dropna()
Ap_comp_reindexed['lat']=Ap_comp_reindexed.lat_interval.apply(lambda x: x.mid).astype(float)
Ap_comp_reindexed['lon']=Ap_comp_reindexed.lon_interval.apply(lambda x: x.mid).astype(float)

Ap_comp_nosoilorder=Ap_comp_reindexed.groupby(['lon','lat']).mean().reset_index()

f,axs=subplots(3,1,num='Ap direct comparison',figsize=(6.8,6.9),clear=True,subplot_kw=dict(projection=cartopy.crs.PlateCarree()))
# ax=subplot(311,projection=cartopy.crs.PlateCarree())
ax=axs[0];sca(ax)
ax.add_feature(cartopy.feature.BORDERS,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.COASTLINE,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.STATES,lw=0.25,facecolor='none',edgecolor='gray')
h=ax.scatter(Ap_comp_nosoilorder.lon,Ap_comp_nosoilorder.lat,s=3,alpha=0.7,edgecolor='None',c=Ap_comp_nosoilorder['Cstock_Ap'],vmax=Cstocks_vmax)
ax.set_title('C stocks to %d cm (has Ap)'%Zmax)
cb=colorbar(h,extend='max');cb.set_label('C stock (%s)'%'kg C m$^{-2}$')
xlim(-129,-65)
ylim(23,52)

# ax=subplot(312,projection=cartopy.crs.PlateCarree())
ax=axs[1];sca(ax)
ax.add_feature(cartopy.feature.BORDERS,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.COASTLINE,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.STATES,lw=0.25,facecolor='none',edgecolor='gray')
h=ax.scatter(Ap_comp_nosoilorder.lon,Ap_comp_nosoilorder.lat,s=3,alpha=0.7,edgecolor='None',c=Ap_comp_nosoilorder['Cstock_noAp'],vmax=Cstocks_vmax)
ax.set_title('C stocks to %d cm (no Ap)'%Zmax)
cb=colorbar(h,extend='max');cb.set_label('C stock (%s)'%'kg C m$^{-2}$')
xlim(-129,-65)
ylim(23,52)

# ax=subplot(313,projection=cartopy.crs.PlateCarree())
ax=axs[2];sca(ax)
ax.add_feature(cartopy.feature.BORDERS,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.COASTLINE,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.STATES,lw=0.25,facecolor='none',edgecolor='gray')
h=ax.scatter(Ap_comp_nosoilorder.lon,Ap_comp_nosoilorder.lat,s=3,alpha=0.7,edgecolor='None',c=Ap_comp_nosoilorder['Cstock_Ap']-Ap_comp_nosoilorder['Cstock_noAp'],vmax=10,vmin=-10,cmap=get_cmap('RdBu_r'))
ax.set_title('C stocks to %d cm (Ap minus no-Ap)'%Zmax)
cb=colorbar(h,extend='both');cb.set_label('C stock (%s)'%'kg C m$^{-2}$')
xlim(-129,-65)
ylim(23,52)

f=figure('Ap direct comparison histograms',figsize=(11.8,5.5),clear=True)
gs=f.add_gridspec(2,2)

# ax=subplot(221,projection=cartopy.crs.PlateCarree())
ax=f.add_subplot(gs[0,0],projection=cartopy.crs.PlateCarree())
ax.add_feature(cartopy.feature.BORDERS,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.COASTLINE,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.STATES,lw=0.25,facecolor=[0.75,0.75,0.75],edgecolor='gray')
h=ax.scatter(Ap_comp_nosoilorder.lon,Ap_comp_nosoilorder.lat,s=4,alpha=1.0,edgecolor='None',
        c=(Ap_comp_nosoilorder['Cstock_Ap']-Ap_comp_nosoilorder['Cstock_noAp']),
        vmax=7,vmin=-7,cmap=get_cmap('RdBu_r'),zorder=10)
ax.set_title('Soil C stocks to %d cm (Ap minus no-Ap)'%Zmax)
cb=colorbar(h,extend='both',ax=ax);cb.set_label('Soil C stock (%s)'%'kg C m$^{-2}$')
ax.set_xlim(-129,-65)
ax.set_ylim(23,52)
ax.text(0.03,1.03,'(a)',transform=ax.transAxes)

# ax=subplot(222,projection=cartopy.crs.PlateCarree())
ax=f.add_subplot(gs[0,1],projection=cartopy.crs.PlateCarree())
ax.add_feature(cartopy.feature.BORDERS,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.COASTLINE,lw=0.25,facecolor='none',edgecolor='gray')
ax.add_feature(cartopy.feature.STATES,lw=0.25,facecolor=[0.75,0.75,0.75],edgecolor='gray')
h=ax.scatter(Ap_comp_nosoilorder.lon,Ap_comp_nosoilorder.lat,s=4,alpha=1.0,edgecolor='None',
        c=(Ap_comp_nosoilorder['Cstock_Ap']-Ap_comp_nosoilorder['Cstock_noAp'])/Ap_comp_nosoilorder['Cstock_noAp']*100,
        vmax=75,vmin=-75,cmap=get_cmap('RdBu_r'),zorder=10)
ax.set_title('Soil C stocks to %d cm (%% difference Ap minus no-Ap)'%Zmax)
cb=colorbar(h,extend='both',ax=ax);cb.set_label('% difference')
ax.set_xlim(-129,-65)
ax.set_ylim(23,52)
ax.text(0.03,1.03,'(b)',transform=ax.transAxes)

bins=linspace(-Cstocks_vmax,Cstocks_vmax,30)
binwidth=bins[1]-bins[0]
# subplot(223)
ax=f.add_subplot(gs[1,0]);sca(ax)
# h=plot(bins[:-1]+binwidth/2,histogram(Ap_comp_reindexed.Cstock_Ap-Ap_comp_reindexed.Cstock_noAp,density=False,bins=bins)[0]*binwidth,'-')
hist(Ap_comp_reindexed.Cstock_Ap-Ap_comp_reindexed.Cstock_noAp,density=False,bins=bins)
# legend()
plot([0,0],[0,200],'k:',lw=0.75)
median_diff=(Ap_comp_reindexed.Cstock_Ap-Ap_comp_reindexed.Cstock_noAp).median()
plot([median_diff,median_diff],[0,200],'r:',lw=1.5)
text(median_diff-1,25,'Median difference',rotation=90,color='r',va='bottom')
title('Difference between Ap and no Ap')
xlabel('Soil C stock to %d cm (%s)'%(Zmax,'kg C m$^{-2}$'))
ylabel('Number of cells')
text(-Cstocks_vmax*0.7,100,'Less soil C with Ap',ha='center')
text(Cstocks_vmax*0.7,100,'More soil C with Ap',ha='center')
text(0.03,1.03,'(c)',transform=gca().transAxes)

bins=linspace(-100,100,30)
binwidth=bins[1]-bins[0]
# subplot(224)
ax=f.add_subplot(gs[1,1]);sca(ax)
# h=plot(bins[:-1]+binwidth/2,histogram(Ap_comp_reindexed.Cstock_Ap-Ap_comp_reindexed.Cstock_noAp,density=False,bins=bins)[0]*binwidth,'-')
hist((Ap_comp_reindexed.Cstock_Ap-Ap_comp_reindexed.Cstock_noAp)/Ap_comp_reindexed.Cstock_noAp*100,density=False,bins=bins)
# legend()
plot([0,0],[0,100],'k:',lw=0.75)
median_diff=((Ap_comp_reindexed.Cstock_Ap-Ap_comp_reindexed.Cstock_noAp)/Ap_comp_reindexed.Cstock_noAp*100).median()
plot([median_diff,median_diff],[0,100],'r:',lw=1.5)
# text(median_diff-1,25,'Median difference',rotation=90,color='r',va='bottom')
title('Percent difference between Ap and no Ap')
xlabel('Percent difference (%)')
ylabel('Number of cells')
text(-75,75,'Less soil C with Ap',ha='center')
text(75,75,'More soil C with Ap',ha='center')
text(0.03,1.03,'(d)',transform=gca().transAxes)

# tight_layout()
# subplots_adjust(top=0.95,left=0.15,right=0.95,hspace=0.5,bottom=0.08)  

# Scatter plot of C loss by SOC content and soil order
fig,ax=subplots(num='C_diff_scatter',clear=True)
# ax=fig.add_subplot()
sca(ax)
for order in unique(Ap_comp_reindexed.soilorder):
    xx=Ap_comp_reindexed.soilorder==order
    ax.plot(Ap_comp_reindexed.Cstock_noAp[xx],Ap_comp_reindexed.Cstock_Ap[xx],'o',label=order,ms=3.0,alpha=0.7,mec='None',c=cols[order])
# reg=scipy.stats.linregress(Ap_comp_reindexed.Cstock_noAp,Ap_comp_reindexed.Cstock_Ap) 
import statsmodels.formula.api as smf
reg=smf.ols('Cstock_Ap ~ Cstock_noAp',Ap_comp_reindexed).fit()
reg_rlm=smf.rlm('Cstock_Ap ~ Cstock_noAp',Ap_comp_reindexed).fit()
reg_orders=smf.ols('Cstock_Ap ~ Cstock_noAp+C(soilorder,Sum(omit="Unknown"))',Ap_comp_reindexed).fit()
reg_orders2=smf.ols('Cstock_Ap ~ Cstock_noAp*C(soilorder,Sum(omit="Unknown"))',Ap_comp_reindexed).fit()
x=linspace(0.1,29,10)
ax.plot(x,reg.params['Intercept']+x*reg.params['Cstock_noAp'],'k--',label='Regression')
fmt='C(soilorder, Sum(omit="Unknown"))[S.%s]'
rmse=sqrt(((reg.model.predict(reg.params)-Ap_comp_reindexed['Cstock_Ap'])**2).mean()) 
# for order in unique(Ap_comp_reindexed.soilorder):
#     if fmt%order in reg_orders.params:
#         if True or reg_orders.pvalues[fmt%order]<0.05:
#             ax.plot(x,reg_orders.params['Intercept']+reg_orders.params[fmt%order]+x*reg_orders.params['Cstock_noAp'],'--',c=cols[order])
# plot(x,reg2.params['Intercept']+x*reg2.params['noAp'],'k.-',label='Regression')
# fmt='C(soilorder, Sum(omit="Unknown"))[S.%s]'
# fmt2='Cstock_noAp:C(soilorder, Sum(omit="Unknown"))[S.%s]'
# for order in unique(Ap_comp_reindexed.soilorder):
#     if fmt%order in reg_orders.params:
#         if reg_orders2.pvalues[fmt%order]<0.05 or reg_orders2.pvalues[fmt2%order]<0.05:
#             ax.plot(x,reg_orders2.params['Intercept']+reg_orders2.params[fmt%order]+x*(reg_orders2.params['Cstock_noAp']+reg_orders2.params[fmt2%order]),'--',c=cols[order])

ax.plot(x,x,'k:',label='1:1 line')
legend(markerscale=2.0,ncol=3)
ax.set_aspect('equal')

xlabel('Soil C stock without Ap (kg C m$^{-2}$)')
ylabel('Soil C stock with Ap (kg C m$^{-2}$)')
title('Soil C stock differences')
ylim(-1,31)
xlim(-1,31)

# Bar plots of land cover comparison
obs_date=pandas.DatetimeIndex(pandas.to_datetime(alldata_goodr2['observation_date (YYYY-MM-DD)']))

alldata_goodr2['landuse']=''
alldata_goodr2['landuse'][obs_date.year<=1995]=alldata_goodr2['NLCD_92'].map(NLCD_mappings)
alldata_goodr2['landuse'][(obs_date.year>1995)&(obs_date.year<=2003)]=alldata_goodr2['NLCD_01'].map(NLCD_mappings)
alldata_goodr2['landuse'][(obs_date.year>2003)&(obs_date.year<=2009)]=alldata_goodr2['NLCD_06'].map(NLCD_mappings)
alldata_goodr2['landuse'][(obs_date.year>2009)]=alldata_goodr2['NLCD_11'].map(NLCD_mappings)

lu_to_use = [
    'Cultivated Crops',
    'Pasture/Hay',
    'Grassland/Herbaceous',
    'Shrub/Scrub',
    'Deciduous Forest',
    'Mixed Forest',
    'Evergreen Forest',
    ]
    

def bar_plot(m,se,c,labels=None,colors=None):
    bar(arange(len(m)),m,yerr=se,ecolor='k',color=colors)
    if labels is None:
        xticks(arange(len(m)),m.index,rotation=90)
    else:
        xticks(arange(len(m)),labels,rotation=90)
    ylims=ylim()
    # dy=(ylims[1]-ylims[0])/10
    # for ii in arange(len(c)):
        # text(ii+0.5,m[ii]+dy,'n=%d'%c[ii],horizontalalignment='center',rotation=90,verticalalignment='bottom')
    # ylim(ylims[0],ylims[1]+dy)

def groups_bar(group,includeonly=None,labels=None,colors=None):
    m=group.median()
    s=group.std()
    c=group.count()
    se=group.sem()
    if includeonly is not None:
        m=m.loc[includeonly]
        s=s.loc[includeonly]
        c=c.loc[includeonly]
        se=se.loc[includeonly]
    bar_plot(m,se,c,labels,colors)
    
f=figure('Land use bar plot',clear=True)
# ax=f.add_axes([0.15,0.35,0.35,0.55])
axs=f.subplots(ncols=2)
sca(axs[0])
colors=[]
for name in lu_to_use:
    if "Forest" in name:
        colors=colors+['green']
    else:
        colors=colors+['brown']
groups_bar(alldata_goodr2[obs_date>datetime.datetime(1990,1,1)].Csurf.groupby(alldata_goodr2.landuse),lu_to_use,colors=colors)
ylabel('Csurf')
title('Csurf')

sca(axs[1])
# ax=f.add_axes([0.6,0.35,0.35,0.55])
groups_bar((alldata_goodr2[obs_date>datetime.datetime(1990,1,1)]['zstar']).groupby(alldata_goodr2.landuse),lu_to_use,colors=colors)
ylabel('Z*')
title('Z*')

def plot_cultiv_history(group,categories=lu_to_use):
    m=group.mean()[:,True][categories]
    se=group.sem()[:,True][categories]
    c=group.count()[:,True][categories]
    b_Ap=bar(arange(len(m)),m,yerr=se,ecolor='k',width=0.3,color='red',label='has Ap')
    # ylims=ylim()
    # dy=(ylims[1]-ylims[0])/10
    # for ii in arange(len(c)):
    #     text(ii+0.15,m[ii]+dy,'n=%d'%c[ii],horizontalalignment='center',verticalalignment='bottom',rotation=90)

    m=group.mean()[:,False][categories]
    se=group.sem()[:,False][categories]
    c=group.count()[:,False][categories]
    b_notAp=bar(arange(len(m))+0.4,m,yerr=se,ecolor='k',width=0.3,color='blue',label='no Ap')
    # ylims=ylim()
    # dy=(ylims[1]-ylims[0])/10
    # for ii in arange(len(c)):
    #     text(ii+0.15+.4,m[ii]+dy,'n=%d'%c[ii],horizontalalignment='center',verticalalignment='bottom',rotation=90)

    xticks(arange(len(m))+0.5*(0.4+0.3),categories,rotation=90)

    # ylim(ylims[0],ylims[1]+dy)
    

f,axs=subplots(2,3,num='Ap effect by land use and soil order',clear=True);
forest_delim=3.75
# subplot(231)
sca(axs[0,0])
dg=alldata_goodr2[obs_date>datetime.datetime(1990,1,1)].zstar.groupby([alldata_goodr2.landuse,alldata_goodr2.has_Ap])
plot_cultiv_history(dg)
ylabel('Z* (cm)')
title('Z*')
ylims=ylim()
plot([forest_delim,forest_delim],[0,ylims[1]],'k--')
text(forest_delim-0.25,ylims[1]*0.99,'Non-forest',fontsize='small',ha='right',va='top')
text(forest_delim+0.25,ylims[1]*0.99,'Forest',fontsize='small',ha='left',va='top')
ylim(ylims)
legend(loc='upper right',bbox_to_anchor=(0.95,0.95))
text(0.03,0.93,'(a)',transform=gca().transAxes)
# subplot(232)
sca(axs[0,1])
dg=alldata_goodr2[obs_date>datetime.datetime(1990,1,1)].Csurf.groupby([alldata_goodr2.landuse,alldata_goodr2.has_Ap])
plot_cultiv_history(dg)
ylims=ylim()
plot([forest_delim,forest_delim],[0,ylims[1]],'k--')
text(forest_delim-0.25,ylims[1]*0.99,'Non-forest',fontsize='small',ha='right',va='top')
text(forest_delim+0.25,ylims[1]*0.99,'Forest',fontsize='small',ha='left',va='top')
ylim(ylims)
ylabel('C$_{s}$ (%s)'%Cunits)
title('C$_{s}$')
text(0.03,0.93,'(b)',transform=gca().transAxes)
# subplot(233)
sca(axs[0,2])
dg=Cstocks.groupby([alldata_goodr2[obs_date>datetime.datetime(1990,1,1)].landuse,alldata_goodr2.has_Ap])
plot_cultiv_history(dg)
ylims=ylim()
plot([forest_delim,forest_delim],[0,ylims[1]],'k--')
text(forest_delim-0.25,ylims[1]*0.99,'Non-forest',fontsize='small',ha='right',va='top')
text(forest_delim+0.25,ylims[1]*0.99,'Forest',fontsize='small',ha='left',va='top')
ylim(ylims)
ylabel('Soil C stocks to 30 cm (kgC m$^{-2}$)')
title('Soil C stock')
text(0.03,0.93,'(c)',transform=gca().transAxes)
# tight_layout()

# figure('Ap effect by soil order');clf()
# subplot(234)
sca(axs[1,0])
dg=alldata_goodr2[obs_date>datetime.datetime(1990,1,1)].zstar.groupby([alldata_goodr2.soilorder,alldata_goodr2.has_Ap])
soilorders_to_use=[]
for so in unique(soilorders):
    if (so,True) in dg.groups and (so,False) in dg.groups and len(dg.groups[so,True])>=10 and len(dg.groups[so,False])>=10:
        soilorders_to_use.append(so)
if 'Unknown' in soilorders_to_use:
    soilorders_to_use.remove('Unknown')
    soilorders_to_use.append('Unknown')
plot_cultiv_history(dg,soilorders_to_use)
ylabel('Z* (cm)')
title('Z*')
ylims=ylim()

ylim(ylims)
legend(loc='upper center')
text(0.03,0.93,'(d)',transform=gca().transAxes)
# subplot(235)
sca(axs[1,1])
dg=alldata_goodr2[obs_date>datetime.datetime(1990,1,1)].Csurf.groupby([alldata_goodr2.soilorder,alldata_goodr2.has_Ap])
plot_cultiv_history(dg,soilorders_to_use)
ylims=ylim()

ylim(ylims)
ylabel('C$_{s}$ (%s)'%Cunits)
title('C$_{s}$')
text(0.03,0.93,'(e)',transform=gca().transAxes)
# subplot(236)
sca(axs[1,2])
dg=Cstocks.groupby([alldata_goodr2[obs_date>datetime.datetime(1990,1,1)].soilorder,alldata_goodr2.has_Ap])
plot_cultiv_history(dg,soilorders_to_use)
ylims=ylim()

ylim(ylims)
ylabel('Soil C stocks to 30 cm (kgC m$^{-2}$)')
title('Soil C stock')
text(0.03,0.93,'(f)',transform=gca().transAxes)
# tight_layout()

def get_group_stats(num):
    lat_int=Ap_comp_reindexed.lat_interval.iloc[num]
    lon_int=Ap_comp_reindexed.lon_interval.iloc[num]
    soilorder=Ap_comp_reindexed.soilorder.iloc[num]
    group_index=(lon_int,lat_int,soilorder)
    group=Cstocks_grouped.get_group(group_index).copy()
    group['zstar']=alldata_goodr2['zstar'][group.index]
    group['Zmin']=alldata_goodr2['Zmin'][group.index]
    group['Csurf']=alldata_goodr2['Csurf'][group.index]
    group['Cdeep']=alldata_goodr2['Cdeep'][group.index]
    return group_index,group


Ap_comp_all=Ap_comparison.reset_index().rename(columns={'long (dec. deg)':'lon_interval','lat (dec. deg)':'lat_interval'})
cellnumbers=pandas.Series(index=alldata_goodr2.index) 
cell_lats=pandas.Series(index=alldata_goodr2.index) 
cell_lons=pandas.Series(index=alldata_goodr2.index) 
for num in range(len(Ap_comp_all)):
    lat_int=Ap_comp_all.lat_interval.iloc[num]
    lon_int=Ap_comp_all.lon_interval.iloc[num]
    soilorder=Ap_comp_all.soilorder.iloc[num]
    group_index=(lon_int,lat_int,soilorder)
    group=Cstocks_grouped.get_group(group_index).copy()
    ind=group.index
    cellnumber=Ap_comp_all.index[num]
    cellnumbers[ind]=cellnumber
    
alldata_goodr2['cell_number']=cellnumbers
alldata_goodr2['lon_cell_bin']=lon_cut
alldata_goodr2['lat_cell_bin']=lat_cut

def plot_profiles_grouped(num=None,clearaxis=True):
    if num is None:
        num=numpy.random.randint(len(Ap_comp_reindexed))
    if clearaxis:
        cla()
    group_index,group=get_group_stats(num)
    print(num,group_index)
    styles=['o','+','x','^','*','>']
    alldepths=[]
    alldots=[]
    labels=[]
    for profnum,prof in enumerate(group[group['has_Ap']].index):
        depths,Cvalues,dots=plot_profile(layerdata,alldata_goodr2,profid=prof,fig_text=False,clearaxis=False,c='r',lw=0.5,ms=4.0,marker=styles[profnum%len(styles)])
        alldepths.append(depths)
        alldots.append(dots)
        labels.append('{profname}: Z*={zstar:1.1f}, Zmin={Zmin:1.1f} Csurf={Csurf:1.1e} Cdeep={Cdeep:1.1e}'\
                .format(profname=prof,order=soilorders[prof],zstar=alldata_goodr2.zstar[prof],
                Csurf=alldata_goodr2.Csurf[prof],Cdeep=alldata_goodr2.Cdeep[prof],Zmin=alldata_goodr2.Zmin[prof],
                latitude=lat[prof]))
    for profnum,prof in enumerate(group[~group['has_Ap']].index):
        depths,Cvalues,dots=plot_profile(layerdata,alldata_goodr2,profid=prof,fig_text=False,clearaxis=False,c='b',lw=0.5,ms=4.0,marker=styles[profnum%len(styles)])
        alldepths.append(depths)
        alldots.append(dots)
        labels.append('{profname}: Z*={zstar:1.1f}, Zmin={Zmin:1.1f} Csurf={Csurf:1.1e} Cdeep={Cdeep:1.1e}'\
                .format(profname=prof,order=soilorders[prof],zstar=alldata_goodr2.zstar[prof],
                Csurf=alldata_goodr2.Csurf[prof],Cdeep=alldata_goodr2.Cdeep[prof],Zmin=alldata_goodr2.Zmin[prof],
                latitude=lat[prof]))

    alldepths=concatenate(alldepths)
    means_Ap=group.groupby('has_Ap').mean().loc[True]
    means_noAp=group.groupby('has_Ap').mean().loc[False]
    start_depth=alldepths.min()
    bottom=alldepths.max()

    z=linspace(0,means_Ap['Zmin'],100)
    h_Ap_mean=plot(means_Ap['Csurf']*exp(-z/means_Ap['zstar']),-z-start_depth,'-',c='r',lw=2.0)[0]
    plot([means_Ap['Cdeep'],means_Ap['Cdeep']],[-means_Ap['Zmin']-start_depth,-bottom-start_depth],'--',c='r',lw=2.0)
    alldots.append(h_Ap_mean)
    labels.append('Mean Ap stats: Z*={zstar:1.1f}, Zmin={Zmin:1.1f} Csurf={Csurf:1.1e} Cdeep={Cdeep:1.1e}'\
            .format(zstar=means_Ap['zstar'],
            Csurf=means_Ap['Csurf'],Cdeep=means_Ap['Cdeep'],Zmin=means_Ap['Zmin']))

    z=linspace(0,means_noAp['Zmin'],100)
    h_noAp_mean=plot(means_noAp['Csurf']*exp(-z/means_noAp['zstar']),-z-start_depth,'-',c='b',lw=2.0)[0]
    plot([means_noAp['Cdeep'],means_noAp['Cdeep']],[-means_noAp['Zmin']-start_depth,-bottom-start_depth],'--',c='b',lw=2.0)
    alldots.append(h_noAp_mean)
    labels.append('Mean no Ap stats: Z*={zstar:1.1f}, Zmin={Zmin:1.1f} Csurf={Csurf:1.1e} Cdeep={Cdeep:1.1e}'\
            .format(zstar=means_noAp['zstar'],
            Csurf=means_noAp['Csurf'],Cdeep=means_noAp['Cdeep'],Zmin=means_noAp['Zmin']))

    legend(handles=alldots,labels=labels,fontsize='small')
    title('Latitude: {lat:1.1f}, Soil order: {order}'.format(lat=lat[prof],order=soilorders[prof]))
    

def make_label(profid,profiledata=alldata_goodr2):
    textstr='Profile ID: {profile}\nSoil order: {order}'
    if profiledata.has_Ap[profid]:
        textstr=textstr + ' (has Ap)'
    else:
        textstr=textstr + ' (no Ap)'
    return (textstr+'''
Z* = {zstar:1.1f} cm
$Z_{{min}}$ = {Zmin:1.1f} cm
$C_s$ = {Csurf:1.2g} gC cm$^{{-3}}$
$C_{{deep}}$ = {Cdeep:1.2g} gC cm$^{{-3}}$
Soil C stock to 30 cm = {Cstock:1.1f} kgC m$^{{-2}}$''').format(
                                        profile=profid,
                                        # profiledata.landuse[profid],
                                        order=get_soil_order(profiledata.soil_taxon[profid]),
                                        zstar=profiledata.zstar[profid],
                                        Csurf=profiledata.Csurf[profid],Cdeep=profiledata.Cdeep[profid],Zmin=profiledata.Zmin[profid],Cstock=Cstocks[profid],
                                        latitude=lat[profid],country=profiledata['country (country)'][profid])
        
# figure('Profile examples');clf()
fig,axes=subplots(2,2,num='Profile examples',squeeze=False,clear=True)
sca(axes[0,0])
plot_profile(layerdata,alldata_goodr2,profid='09N0987',fig_text=False)
text(0.3,0.05,make_label('09N0987'),transform=gca().transAxes)
text(0.03,0.93,'(a)',transform=gca().transAxes)
ylim(-201,0)
xlim(-0.001,0.075)
xlabel('Soil C density (gC cm$^{-3}$)')
sca(axes[0,1])
plot_profile(layerdata,alldata_goodr2,profid='60135',fig_text=False)
text(0.3,0.05,make_label('60135'),transform=gca().transAxes)
text(0.03,0.93,'(b)',transform=gca().transAxes)
ylim(-201,0)
xlim(-0.001,0.075)
xlabel('Soil C density (gC cm$^{-3}$)')

sca(axes[1,0])
plot_profile(layerdata,alldata_goodr2,profid='79P0163',fig_text=False)
text(0.3,0.05,make_label('79P0163'),transform=gca().transAxes)
text(0.03,0.93,'(c)',transform=gca().transAxes)
ylim(-201,0)
xlim(-0.001,0.075)
xlabel('Soil C density (gC cm$^{-3}$)')

sca(axes[1,1])
plot_profile(layerdata,alldata_goodr2,profid='uiuc196200221',fig_text=False)
text(0.3,0.05,make_label('uiuc196200221'),transform=gca().transAxes)
text(0.03,0.93,'(d)',transform=gca().transAxes)
ylim(-201,0)
xlim(-0.001,0.075)
xlabel('Soil C density (gC cm$^{-3}$)')

# tight_layout()
figlist=['Profile examples','Fit stats','Profiles map global','Profiles map US','Ap horizon histograms','Ap effect by land use and soil order',
        'Ap direct comparison histograms','C_diff_scatter']
def save_all_figs(dirname,fmt='pdf',savedata=True,fignumbers=False):
    for fname in get_figlabels():
        print(fname)
        if fignumbers and fname in figlist:
            figure(fname).savefig(dirname+'/'+'Figure %d %s'%(figlist.index(fname)+1,fname)+'.'+fmt)
        else:
            figure(fname).savefig(dirname + '/' + fname + '.' + fmt)
    if savedata:
        alldata_goodr2.to_csv(dirname+'/alldata_goodr2.csv')

show()
