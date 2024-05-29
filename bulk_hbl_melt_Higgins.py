#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   bulk_hbl_melt_Higgins.py
@Time    :   2022/12/07 11:18:09
@Author  :   Ke Gao 
@Version :   1.0
@Contact :   earthgaoke@gmail.com
@Desc    :   None
'''

#%% import
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

bulk_melt = pd.read_csv('vol_plu_all_higgins2021.csv')
bulk_meltCM = bulk_melt[bulk_melt['TECTONIC SETTING'] == 'CONVERGENT MARGIN']
vol = bulk_meltCM[bulk_meltCM['ROCK TYPE'] == 'VOL']
plu = bulk_meltCM[bulk_meltCM['ROCK TYPE'] == 'PLU'] 

bulk_elements = ['SIO2(WT%)', 'TIO2(WT%)', 'AL2O3(WT%)', 'FEOT(WT%)', 'MGO(WT%)', 'CAO(WT%)', 'K2O(WT%)', 'NA2O(WT%)']
melt_elements = ['SiO2_pred', 'Al2O3_pred', 'CaO_pred', 'Na2O_pred', 'K2O_pred', 'FeO_pred', 'MgO_pred', 'TiO2_pred']
#%% 
for b in bulk_elements:
    for m in melt_elements:
        if b.lower()[:-6] in m.lower():
            fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(10, 4))

            sns.kdeplot(x = vol[b],y = vol[m], fill = True, color = '#E60000',
            # scatter_kws={'alpha':0.4, 's':4,}, fit_reg = False, x_jitter = 0.5, y_jitter = 0.2, 
                        ax = axs[0])
            axs[0].scatter(vol[b], vol[m], color = 'gray', alpha = .3, s = 1)
            axs[0].plot(np.arange(0,100),np.arange(0,100), c='gray')
            axs[0].set_title(r'$Volcanic$')
#             axs[0].set_xlabel('Bulk ' + b)
#             axs[0].set_ylabel('Melt ' + b)
            
            sns.kdeplot(x = plu[b],y =  plu[m], fill = True, color = '#003153',
            # scatter_kws={'alpha':0.4, 's':4,}, fit_reg = False, x_jitter = 0.5, y_jitter = 0.2, 
                        ax = axs[1])
            axs[1].scatter(plu[b], plu[m], color = 'gray', alpha = .3, s = 1)
            axs[1].plot(np.arange(0,100),np.arange(0,100), c='gray')
            axs[1].set_title(r'$Plutonic$')
            axs[1].set_xlabel('Bulk ' + b)
            axs[1].set_ylabel('Melt ' + b)
            
            if b == 'SIO2(WT%)':
                axs[0].set_xlim(45, 80)
                axs[0].set_ylim(45, 80) 
                axs[0].set_xlabel(r'$Bulk$ ' + r'$SiO_2$')
                axs[0].set_ylabel(r'$Melt$ ' + r'$SiO_2$')

                axs[1].set_xlim(40, 80)
                axs[1].set_ylim(50, 80)
                axs[1].set_xlabel(r'$Bulk$ ' + r'$SiO_2$')
                axs[1].set_ylabel(r'$Melt$ ' + r'$SiO_2$')


                
            elif b == 'TIO2(WT%)':
                axs[0].set_xlim(0, 2)
                axs[0].set_ylim(0, 2)
                axs[0].set_xlabel(r'$Bulk$ ' + r'$TiO_2$')
                axs[0].set_ylabel(r'$Melt$ ' + r'$TiO_2$')
                
                axs[1].set_xlim(0, 2)
                axs[1].set_ylim(0, 2)
                axs[1].set_xlabel(r'$Bulk$ ' + r'$TiO_2$')
                axs[1].set_ylabel(r'$Melt$ ' + r'$TiO_2$')
                
            elif b == 'AL2O3(WT%)':
                axs[0].set_xlim(10, 23)
                axs[0].set_ylim(10, 22) 
                axs[0].set_xlabel(r'$Bulk$ ' + r'$Al_2O_3$')
                axs[0].set_ylabel(r'$Melt$ ' + r'$Al_2O_3$')
                
                axs[1].set_xlim(10, 23)
                axs[1].set_ylim(10, 22)
                axs[1].set_xlabel(r'$Bulk$ ' + r'$Al_2O_3$')
                axs[1].set_ylabel(r'$Melt$ ' + r'$Al_2O_3$')
                
            elif b == 'MGO(WT%)':
                axs[0].set_xlim(0, 12)
                axs[0].set_ylim(0, 7)
                axs[0].set_xlabel(r'$Bulk$ ' + r'$MgO$')
                axs[0].set_ylabel(r'$Melt$ ' + r'$MgO$')
                
                axs[1].set_xlim(0, 15)
                axs[1].set_ylim(0, 6)
                axs[1].set_xlabel(r'$Bulk$ ' + r'$MgO$')
                axs[1].set_ylabel(r'$Melt$ ' + r'$MgO$')
                
            elif b == 'CAO(WT%)':
                axs[0].set_xlim(0, 15)
                axs[0].set_ylim(0, 12.5)    
                axs[0].set_xlabel(r'$Bulk$ ' + r'$CaO$')
                axs[0].set_ylabel(r'$Melt$ ' + r'$CaO$')
                
                axs[1].set_xlim(0, 17.5)
                axs[1].set_ylim(0, 15)
                axs[1].set_xlabel(r'$Bulk$ ' + r'$CaO$')
                axs[1].set_ylabel(r'$Melt$ ' + r'$CaO$')
                
            elif b == 'K2O(WT%)':
                axs[0].set_xlim(0, 6)
                axs[0].set_ylim(0, 6)   
                axs[0].set_xlabel(r'$Bulk$ ' + r'$K_2O$')
                axs[0].set_ylabel(r'$Melt$ ' + r'$K_2O$')
                
                axs[1].set_xlim(0, 6)
                axs[1].set_ylim(0, 6)
                axs[1].set_xlabel(r'$Bulk$ ' + r'$K_2O$')
                axs[1].set_ylabel(r'$Melt$ ' + r'$K_2O$')
            elif b == 'FEOT(WT%)':
                axs[0].set_xlim(0, 12)
                axs[0].set_ylim(0, 12)
                axs[0].set_xlabel(r'$Bulk$ ' + r'$FeO_t$')
                axs[0].set_ylabel(r'$Melt$ ' + r'$FeO_t$')
                
                axs[1].set_xlim(0, 17.5)
                axs[1].set_ylim(0, 10)
                axs[1].set_xlabel(r'$Bulk$ ' + r'$FeO_t$')
                axs[1].set_ylabel(r'$Melt$ ' + r'$FeO_t$')
            plt.savefig('./figures/figures_geology/' + m + '_bulk_melt' + '_Higgins.pdf', transparent=True)
            plt.show()

#%%
melt_bulk_dict = dict(zip(melt_elements, bulk_elements))
for e in melt_elements:
    fig, axs = plt.subplots()
    if e == 'SiO2_pred':
        elm_min = 50
        elm_max = 80
        axs.set_xlim(elm_min, elm_max)
        axs.set_xlabel('$SiO_2$')
    elif e == 'FeO_pred':
        elm_min = 0
        elm_max = 11
        axs.set_xlim(elm_min, elm_max)
        axs.set_xlabel('$FeO_t$')
    elif e == 'MgO_pred':
        elm_min = 0
        elm_max = 5
        axs.set_xlim(elm_min, elm_max)
        axs.set_xlabel('$MgO$')
    elif e == 'CaO_pred':
        elm_min = 0
        elm_max = 12
        axs.set_xlim(elm_min, elm_max)
        axs.set_xlabel('$CaO$')
    elif e == 'Al2O3_pred':
        elm_min = 12
        elm_max = 21
        axs.set_xlim(elm_min, elm_max)
        axs.set_xlabel('$Al_2O_3$')
    elif e == 'K2O_pred':
        elm_min = 0
        elm_max = 5        
        axs.set_xlim(elm_min, elm_max)
        axs.set_xlabel('$K_2O$')
    elif e == 'TiO2_pred':
        elm_min = 0
        elm_max = 1.5        
        axs.set_xlim(elm_min, elm_max)
        axs.set_xlabel('$TiO_2$')
    elif e == 'Na2O_pred':
        elm_min = 2
        elm_max = 6        
        axs.set_xlim(elm_min, elm_max)
        axs.set_xlabel('$Na_2O$')
    sns.histplot(vol[e][(vol[e]<=elm_max) & (vol[e]>=elm_min)], color = '#E60000', bins = np.arange(elm_min,elm_max,(elm_max-elm_min)/40), stat = 'probability', common_norm=False)
    sns.histplot(plu[e][(plu[e]<=elm_max) & (plu[e]>=elm_min)], color = '#003153', bins = np.arange(elm_min,elm_max,(elm_max-elm_min)/40), stat = 'probability', common_norm=False)
    sns.kdeplot(vol[melt_bulk_dict[e]][(vol[melt_bulk_dict[e]]<=elm_max) & (vol[melt_bulk_dict[e]]>=elm_min)], color = '#db4646', alpha=0.5)
    sns.kdeplot(plu[melt_bulk_dict[e]][(plu[melt_bulk_dict[e]]<=elm_max) & (plu[melt_bulk_dict[e]]>=elm_min)], color = '#0090ab', alpha=0.5)
    axs.set_ylabel(r'$Probability$')
    axs.legend(labels=['Volcanic bulk',"Plutonic bulk", "Volcanic melt", "Plutonic melt"])
    plt.savefig('./figures/figures_geology/melts_hist_' + e + '_Higgins.pdf', format="pdf", transparent=True)
    plt.show()            
# %%
