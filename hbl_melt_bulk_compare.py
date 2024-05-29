#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   hbl_melt_bulk_compare.py
@Time    :   2023/12/07 20:53:43
@Author  :   Ke Gao 
@Version :   1.0
@Contact :   earthgaoke@gmail.com
@Desc    :   Comparison of crystals and melts produced by simulations of whole-rock compositions ( both from georoc volcanic rocks and our compiled whole-rock-hornblende pairing dataset) as starting compositions.
'''

#%% here put the import lib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# from matplotlib.patches import Rectangle
import seaborn as sns

# bulk_melt = pd.read_csv('hbl_bulk_georoc_vol.csv')
# bulk_melt = pd.read_csv('bulk_georoc_plu_all.csv')
# bulk_melt = pd.read_csv('bulk_georoc_vol_all.csv')

# simulated from georoc volcanic rocks
bulk_melt1 = pd.read_csv('bulk_georoc_plu_all.csv')
bulk_melt2 = pd.read_csv('bulk_georoc_vol_all.csv')

# simulated from our compiled whole-rock-hornblende pairing dataset
# bulk_melt1 = pd.read_csv('hbl_bulk_georoc_plu.csv')
# bulk_melt2 = pd.read_csv('hbl_bulk_georoc_vol.csv')
bulk_melt = pd.concat([bulk_melt1,bulk_melt2])
# bulk_melt.to_csv('georoc_hbl_result.csv', index=False)
bulk_elements = ['SiO2', 'TiO2', 'Al2O3', 'FeO', 'MgO', 'CaO', 'K2O']
# bulk_elements = ['FeO', 'SiO2']
# %% simulate crystallization
bulk_melt = bulk_melt.loc[(bulk_melt['melt_vol'] > 5) & (bulk_melt['melt_vol'] < 95) & (bulk_melt['hbl_SiO2'].notna()) & (bulk_melt['bulk_SiO2'] > 40)]
bulk_melt = bulk_melt.sample(n=10000)

for b in bulk_elements:
    fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(5, 4))

    sns.kdeplot(x = bulk_melt['bulk_' + b],y = bulk_melt['melt_' + b], fill=True, color = '#7E2FFF', ax = axs)
    axs.scatter(bulk_melt['bulk_' + b], bulk_melt['melt_' + b], color = 'gray', alpha = .1, s = .2)
    axs.plot(np.arange(0,100),np.arange(0,100), c='gray')
    axs.set_title(r'$Simulated$ ' + r'$plutonic$')
    # axs.set_xlabel('Perple_X system ' + b)
    # axs.set_ylabel('Perple_X pure melt ' + b)
    
    if b == 'SiO2':
        axs.set_xlim(40, 80)
        axs.set_ylim(40, 80) 
        axs.set_xlabel(r'$Bulk$ ' + r'$SiO_2$')
        axs.set_ylabel(r'$Melt$ ' + r'$SiO_2$')


    elif b == 'TiO2':
        axs.set_xlim(0, 4)
        axs.set_ylim(0, 4)
        axs.set_xlabel(r'$Bulk$ ' + r'$TiO_2$')
        axs.set_ylabel(r'$Melt$ ' + r'$TiO_2$')
        
    elif b == 'Al2O3':
        axs.set_xlim(7, 25)
        axs.set_ylim(7, 25) 
        axs.set_xlabel(r'$Bulk$ ' + r'$Al_2O_3$')
        axs.set_ylabel(r'$Melt$ ' + r'$Al_2O_3$')
        
    elif b == 'MgO':
        axs.set_xlim(0, 16)
        axs.set_ylim(0, 16)
        axs.set_xlabel(r'$Bulk$ ' + r'$MgO$')
        axs.set_ylabel(r'$Melt$ ' + r'$MgO$')
        
    elif b == 'CaO':
        axs.set_xlim(0, 15)
        axs.set_ylim(0, 15)    
        axs.set_xlabel(r'$Bulk$ ' + r'$CaO$')
        axs.set_ylabel(r'$Melt$ ' + r'$CaO$')
        
    elif b == 'K2O':
        axs.set_xlim(0, 8)
        axs.set_ylim(0, 8)   
        axs.set_xlabel(r'$Bulk$ ' + r'$K_2O$')
        axs.set_ylabel(r'$Melt$ ' + r'$K_2O$')

    elif b == 'FeO':
        axs.set_xlim(0, 23)
        axs.set_ylim(0, 23)
        axs.set_xlabel(r'$Bulk$ ' + r'$FeO_t$')
        axs.set_ylabel(r'$Melt$ ' + r'$FeO_t$')

    plt.savefig('./fig/' + b + '_sim_plu_cry_georoc_all_5-95vol' + '.svg', format="svg", transparent=True)
    plt.show()
# %% Simulating volcanics formed by melt extraction
bulk_melt = bulk_melt.loc[(bulk_melt['melt_vol'] > 30) & (bulk_melt['melt_vol'] < 50) & (bulk_melt['hbl_SiO2'].notna()) & (bulk_melt['bulk_SiO2'] > 40)]
bulk_melt = bulk_melt.sample(n=10000)

# add some crystals to the melt according to the separation efficiency of crystals and melts
crystal_percent = abs(bulk_melt['melt_vol'] - 40) / 100
crystal_percent_1 = np.array([np.random.normal(e, 0.1) for e in crystal_percent])
crystal_percent_1 = abs(crystal_percent_1)
crystal_percent_2 = np.array([np.random.normal(e, 0.1) for e in crystal_percent])
crystal_percent_2 = abs(crystal_percent_2)

# build a DF to store melts (Mixed with a few crystals)
melt_melt = pd.DataFrame()
for e in bulk_elements:
    melt_melt['melt_1_'+e] = bulk_melt['melt_'+e] * (1 - crystal_percent_1) + bulk_melt['solid_'+e] * crystal_percent_1
    melt_melt['melt_2_'+e] = bulk_melt['melt_'+e] * (1 - crystal_percent_2) + bulk_melt['solid_'+e] * crystal_percent_2

    fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(5, 4))

    sns.kdeplot(x = melt_melt['melt_1_' + e],y = melt_melt['melt_2_' + e], fill=True, color = '#f8ac8c', ax = axs, levels=20)
    axs.scatter(melt_melt['melt_1_' + e], melt_melt['melt_2_' + e], color = 'gray', alpha = .1, s = .1)
    axs.plot(np.arange(0,100),np.arange(0,100), c='gray')
    axs.set_title(r'$Simulated$ ' + r'$volcanic$')
    # axs.set_xlabel('Perple_X system ' + b)
    # axs.set_ylabel('Perple_X pure melt ' + b)
    
    if e == 'SiO2':
        axs.set_xlim(40, 80)
        axs.set_ylim(40, 80) 
        axs.set_xlabel(r'$Bulk$ ' + r'$SiO_2$')
        axs.set_ylabel(r'$Melt$ ' + r'$SiO_2$')


    elif e == 'TiO2':
        axs.set_xlim(0, 2.5)
        axs.set_ylim(0, 2.5)
        axs.set_xlabel(r'$Bulk$ ' + r'$TiO_2$')
        axs.set_ylabel(r'$Melt$ ' + r'$TiO_2$')
        
    elif e == 'Al2O3':
        axs.set_xlim(8, 22)
        axs.set_ylim(8, 22) 
        axs.set_xlabel(r'$Bulk$ ' + r'$Al_2O_3$')
        axs.set_ylabel(r'$Melt$ ' + r'$Al_2O_3$')
        
    elif e == 'MgO':
        axs.set_xlim(0, 11)
        axs.set_ylim(0, 11)
        axs.set_xlabel(r'$Bulk$ ' + r'$MgO$')
        axs.set_ylabel(r'$Melt$ ' + r'$MgO$')
        
    elif e == 'CaO':
        axs.set_xlim(0, 11)
        axs.set_ylim(0, 11)    
        axs.set_xlabel(r'$Bulk$ ' + r'$CaO$')
        axs.set_ylabel(r'$Melt$ ' + r'$CaO$')
        
    elif e == 'K2O':
        axs.set_xlim(0, 6)
        axs.set_ylim(0, 6)   
        axs.set_xlabel(r'$Bulk$ ' + r'$K_2O$')
        axs.set_ylabel(r'$Melt$ ' + r'$K_2O$')

    elif e == 'FeO':
        axs.set_xlim(0, 18)
        axs.set_ylim(0, 18)
        axs.set_xlabel(r'$Bulk$ ' + r'$FeO_t$')
        axs.set_ylabel(r'$Melt$ ' + r'$FeO_t$')

    # plt.savefig('./fig/' + b + '_bulk_melt_georoc_all' + '.svg', format="svg", transparent=True)
    plt.savefig('./fig/' + e + '_sim_vol_melt_extaract_georoc_all' + '.svg', format="svg", transparent=True)
    plt.show()
# %% Simulating plutons formed by melt extraction
# melt = pd.read_csv('melt_extraction_simulated.csv')
# melt = melt.loc[(melt['melt_vol'] >= 30) & (melt['melt_vol'] <= 50) & (melt['hbl_SiO2'].notna()) 
#                      & (melt['bulk_SiO2'] >= 40) & (melt['bulk_SiO2'] <= 80)]
# # add ~30% melts to crystals
# melt_pst = np.random.normal(20, 3, size=len(melt)) / 100
# crystal_pst = 1 - melt_pst
# # build a DF to store melts (Mixed with a few crystals)
# crystal_melt = pd.DataFrame()
# for e in bulk_elements:
#     crystal_melt['melt_'+e] = melt['melt_'+e]
#     crystal_melt['residue_'+e] = melt['melt_'+e] * melt_pst + melt['solid_'+e] * crystal_pst

#     fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(5, 4))
#     sns.kdeplot(x = crystal_melt['residue_' + e], y = crystal_melt['melt_' + e], fill=True, color = '#7E2FFF', ax = axs, levels=20)
#     axs.scatter(crystal_melt['residue_' + e], y = crystal_melt['melt_' + e], color = 'gray', alpha = .1, s = .1)
#     axs.plot(np.arange(0,100),np.arange(0,100), c='gray')
#     axs.set_title(r'$Simulated$ ' + r'$plutonic$')
#     # axs.set_xlabel('Perple_X system ' + b)
#     # axs.set_ylabel('Perple_X pure melt ' + b)
    
#     if e == 'SiO2':
#         axs.set_xlim(40, 80)
#         axs.set_ylim(40, 80) 
#         axs.set_xlabel(r'$Bulk$ ' + r'$SiO_2$')
#         axs.set_ylabel(r'$Melt$ ' + r'$SiO_2$')


#     elif e == 'TiO2':
#         axs.set_xlim(0, 4)
#         axs.set_ylim(0, 4)
#         axs.set_xlabel(r'$Bulk$ ' + r'$TiO_2$')
#         axs.set_ylabel(r'$Melt$ ' + r'$TiO_2$')
        
#     elif e == 'Al2O3':
#         axs.set_xlim(7, 25)
#         axs.set_ylim(7, 25) 
#         axs.set_xlabel(r'$Bulk$ ' + r'$Al_2O_3$')
#         axs.set_ylabel(r'$Melt$ ' + r'$Al_2O_3$')
        
#     elif e == 'MgO':
#         axs.set_xlim(0, 16)
#         axs.set_ylim(0, 16)
#         axs.set_xlabel(r'$Bulk$ ' + r'$MgO$')
#         axs.set_ylabel(r'$Melt$ ' + r'$MgO$')
        
#     elif e == 'CaO':
#         axs.set_xlim(0, 15)
#         axs.set_ylim(0, 15)    
#         axs.set_xlabel(r'$Bulk$ ' + r'$CaO$')
#         axs.set_ylabel(r'$Melt$ ' + r'$CaO$')
        
#     elif e == 'K2O':
#         axs.set_xlim(0, 8)
#         axs.set_ylim(0, 8)   
#         axs.set_xlabel(r'$Bulk$ ' + r'$K_2O$')
#         axs.set_ylabel(r'$Melt$ ' + r'$K_2O$')

#     elif e == 'FeO':
#         axs.set_xlim(0, 23)
#         axs.set_ylim(0, 23)
#         axs.set_xlabel(r'$Bulk$ ' + r'$FeO_t$')
#         axs.set_ylabel(r'$Melt$ ' + r'$FeO_t$')

#     # plt.savefig('./fig/' + b + '_melt_georoc_all' + '.svg', format="svg", transparent=True)
#     plt.savefig('./fig/' + e + '_sim_plu_melt_extract_georoc_all' + '.svg', format="svg", transparent=True)
#     plt.show()
# %% Simulating plutons formed by melt extraction
bulk_melt = bulk_melt.loc[(bulk_melt['melt_vol'] > 30) & (bulk_melt['melt_vol'] < 50) & (bulk_melt['hbl_SiO2'].notna()) & (bulk_melt['bulk_SiO2'] > 40)]
bulk_melt = bulk_melt.sample(n=10000)

# add ~30% melts to crystals
melt_pst = np.random.normal(40, 3, size=len(bulk_melt)) / 100
crystal_pst = 1 - melt_pst
# build a DF to store melts (Mixed with a few crystals)
crystal_melt = pd.DataFrame()
for e in bulk_elements:
    crystal_melt['melt_'+e] = bulk_melt['melt_'+e]
    crystal_melt['residue_'+e] = bulk_melt['melt_'+e] * melt_pst + bulk_melt['solid_'+e] * crystal_pst

    fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(5, 4))
    sns.kdeplot(x = crystal_melt['residue_' + e],y = crystal_melt['melt_' + e], fill=True, color = '#7E2FFF', ax = axs, levels=20)
    axs.scatter(crystal_melt['residue_' + e], y = crystal_melt['melt_' + e], color = 'gray', alpha = .1, s = .1)
    axs.plot(np.arange(0,100),np.arange(0,100), c='gray')
    axs.set_title(r'$Simulated$ ' + r'$plutonic$')
    # axs.set_xlabel('Perple_X system ' + b)
    # axs.set_ylabel('Perple_X pure melt ' + b)
    
    if e == 'SiO2':
        axs.set_xlim(40, 80)
        axs.set_ylim(40, 80) 
        axs.set_xlabel(r'$Bulk$ ' + r'$SiO_2$')
        axs.set_ylabel(r'$Melt$ ' + r'$SiO_2$')


    elif e == 'TiO2':
        axs.set_xlim(0, 4)
        axs.set_ylim(0, 4)
        axs.set_xlabel(r'$Bulk$ ' + r'$TiO_2$')
        axs.set_ylabel(r'$Melt$ ' + r'$TiO_2$')
        
    elif e == 'Al2O3':
        axs.set_xlim(7, 25)
        axs.set_ylim(7, 25) 
        axs.set_xlabel(r'$Bulk$ ' + r'$Al_2O_3$')
        axs.set_ylabel(r'$Melt$ ' + r'$Al_2O_3$')
        
    elif e == 'MgO':
        axs.set_xlim(0, 16)
        axs.set_ylim(0, 16)
        axs.set_xlabel(r'$Bulk$ ' + r'$MgO$')
        axs.set_ylabel(r'$Melt$ ' + r'$MgO$')
        
    elif e == 'CaO':
        axs.set_xlim(0, 15)
        axs.set_ylim(0, 15)    
        axs.set_xlabel(r'$Bulk$ ' + r'$CaO$')
        axs.set_ylabel(r'$Melt$ ' + r'$CaO$')
        
    elif e == 'K2O':
        axs.set_xlim(0, 8)
        axs.set_ylim(0, 8)   
        axs.set_xlabel(r'$Bulk$ ' + r'$K_2O$')
        axs.set_ylabel(r'$Melt$ ' + r'$K_2O$')

    elif e == 'FeO':
        axs.set_xlim(0, 23)
        axs.set_ylim(0, 23)
        axs.set_xlabel(r'$Bulk$ ' + r'$FeO_t$')
        axs.set_ylabel(r'$Melt$ ' + r'$FeO_t$')

    # plt.savefig('./fig/' + b + '_bulk_melt_georoc_all' + '.svg', format="svg", transparent=True)
    plt.savefig('./fig/' + e + '_sim_plu_melt_extract_georoc_all' + '.svg', format="svg", transparent=True)
    plt.show()
# %% Simulating reequliburate
bulk_melt = bulk_melt.loc[(bulk_melt['melt_vol'] > 0) & (bulk_melt['melt_vol'] <= 5) & (bulk_melt['hbl_SiO2'].notna()) & (bulk_melt['bulk_SiO2'] > 40)]
resampling_weights = (bulk_melt['bulk_SiO2'] - 40) / 40
bulk_melt = bulk_melt.sample(n=10000, weights=resampling_weights, replace=True)
# bulk_melt.to_csv('bulk_melt_resample_si_0-5_good.csv') # Save the samples that look good
bulk_melt = pd.read_csv('bulk_melt_resample_si_0-5_good.csv')

for b in bulk_elements:
    fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(5, 4))

    sns.kdeplot(x = bulk_melt['bulk_' + b],y = bulk_melt['melt_' + b], fill=True, color = '#7E2FFF', ax = axs)
    axs.scatter(bulk_melt['bulk_' + b], bulk_melt['melt_' + b], color = 'gray', alpha = .1, s = .1)
    axs.plot(np.arange(0,100),np.arange(0,100), c='gray')
    axs.set_title(r'$Simulated$ ' + r'$plutonic$')
    # axs.set_xlabel('Perple_X system ' + b)
    # axs.set_ylabel('Perple_X pure melt ' + b)
    
    if b == 'SiO2':
        axs.set_xlim(40, 80)
        axs.set_ylim(40, 80) 
        axs.set_xlabel(r'$Bulk$ ' + r'$SiO_2$')
        axs.set_ylabel(r'$Melt$ ' + r'$SiO_2$')


    elif b == 'TiO2':
        axs.set_xlim(0, 4)
        axs.set_ylim(0, 2)
        axs.set_xlabel(r'$Bulk$ ' + r'$TiO_2$')
        axs.set_ylabel(r'$Melt$ ' + r'$TiO_2$')
        
    elif b == 'Al2O3':
        axs.set_xlim(7, 25)
        axs.set_ylim(7, 25) 
        axs.set_xlabel(r'$Bulk$ ' + r'$Al_2O_3$')
        axs.set_ylabel(r'$Melt$ ' + r'$Al_2O_3$')
        
    elif b == 'MgO':
        axs.set_xlim(0, 16)
        axs.set_ylim(0, 6)
        axs.set_xlabel(r'$Bulk$ ' + r'$MgO$')
        axs.set_ylabel(r'$Melt$ ' + r'$MgO$')
        
    elif b == 'CaO':
        axs.set_xlim(0, 17.5)
        axs.set_ylim(0, 6)    
        axs.set_xlabel(r'$Bulk$ ' + r'$CaO$')
        axs.set_ylabel(r'$Melt$ ' + r'$CaO$')
        
    elif b == 'K2O':
        axs.set_xlim(0, 6)
        axs.set_ylim(0, 6)   
        axs.set_xlabel(r'$Bulk$ ' + r'$K_2O$')
        axs.set_ylabel(r'$Melt$ ' + r'$K_2O$')

    elif b == 'FeO':
        axs.set_xlim(0, 17.5)
        axs.set_ylim(0, 10)
        axs.set_xlabel(r'$Bulk$ ' + r'$FeO_t$')
        axs.set_ylabel(r'$Melt$ ' + r'$FeO_t$')

    # plt.savefig('./fig/' + b + '_bulk_melt_georoc_all' + '.svg', format="svg", transparent=True)
    plt.savefig('./fig/' + b + '_sim_plu_reeq' + '.svg', format="svg", transparent=True)
    plt.show()
# %% no melt extraction volcanic
bulk_melt = bulk_melt.loc[(bulk_melt['melt_vol'] > 80) & (bulk_melt['melt_vol'] <= 200) & (bulk_melt['hbl_SiO2'].notna()) & (bulk_melt['bulk_SiO2'] > 40)]

resampling_weights = (bulk_melt['bulk_SiO2'] - 40) / 40
bulk_melt = bulk_melt.sample(n=10000, weights=resampling_weights, replace=True)

for b in bulk_elements:
    fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(5, 4))

    sns.kdeplot(x = bulk_melt['bulk_' + b],y = bulk_melt['melt_' + b], fill=True, color = '#7E2FFF', ax = axs)
    axs.scatter(bulk_melt['bulk_' + b], bulk_melt['melt_' + b], color = 'gray', alpha = .1, s = .1)
    axs.plot(np.arange(0,100),np.arange(0,100), c='gray')
    axs.set_title(r'$Simulated$ ' + r'$plutonic$')
    # axs.set_xlabel('Perple_X system ' + b)
    # axs.set_ylabel('Perple_X pure melt ' + b)
    
    if b == 'SiO2':
        axs.set_xlim(40, 80)
        axs.set_ylim(40, 80) 
        axs.set_xlabel(r'$Bulk$ ' + r'$SiO_2$')
        axs.set_ylabel(r'$Melt$ ' + r'$SiO_2$')


    elif b == 'TiO2':
        axs.set_xlim(0, 4)
        axs.set_ylim(0, 4)
        axs.set_xlabel(r'$Bulk$ ' + r'$TiO_2$')
        axs.set_ylabel(r'$Melt$ ' + r'$TiO_2$')
        
    elif b == 'Al2O3':
        axs.set_xlim(7, 25)
        axs.set_ylim(7, 25) 
        axs.set_xlabel(r'$Bulk$ ' + r'$Al_2O_3$')
        axs.set_ylabel(r'$Melt$ ' + r'$Al_2O_3$')
        
    elif b == 'MgO':
        axs.set_xlim(0, 16)
        axs.set_ylim(0, 16)
        axs.set_xlabel(r'$Bulk$ ' + r'$MgO$')
        axs.set_ylabel(r'$Melt$ ' + r'$MgO$')
        
    elif b == 'CaO':
        axs.set_xlim(0, 17.5)
        axs.set_ylim(0, 10)    
        axs.set_xlabel(r'$Bulk$ ' + r'$CaO$')
        axs.set_ylabel(r'$Melt$ ' + r'$CaO$')
        
    elif b == 'K2O':
        axs.set_xlim(0, 6)
        axs.set_ylim(0, 6)   
        axs.set_xlabel(r'$Bulk$ ' + r'$K_2O$')
        axs.set_ylabel(r'$Melt$ ' + r'$K_2O$')

    elif b == 'FeO':
        axs.set_xlim(0, 17.5)
        axs.set_ylim(0, 15)
        axs.set_xlabel(r'$Bulk$ ' + r'$FeO_t$')
        axs.set_ylabel(r'$Melt$ ' + r'$FeO_t$')

    # plt.savefig('./fig/' + b + '_bulk_melt_georoc_all' + '.svg', format="svg", transparent=True)
    # plt.savefig('./fig/' + b + '_sim_plu_reeq' + '.svg', format="svg", transparent=True)
    plt.show()
# %%
