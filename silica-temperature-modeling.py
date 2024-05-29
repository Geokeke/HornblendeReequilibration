#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   silica-temperature-modeling.py
@Time    :   2024/05/01 14:13:55
@Author  :   Ke Gao 
@Version :   1.0
@Contact :   earthgaoke@gmail.com
@Desc    :   Comparison of the chemical behavior of simulated volcanic and plutonic rocks in cooling crystallization, melt extraction, and reequilibration modes.
'''

# %% import lib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
# %%
nbins = 20
SiO2min = 40
SiO2max = 80
binedges = np.linspace(SiO2min, SiO2max, num=(nbins + 1))

bulk_melt1 = pd.read_csv('bulk_georoc_plu_all.csv')
bulk_melt2 = pd.read_csv('bulk_georoc_vol_all.csv')
bulk_melt = pd.concat([bulk_melt1,bulk_melt2])
# %% make crystallization
cry = bulk_melt.loc[(bulk_melt['melt_vol'] >= 5) & (bulk_melt['melt_vol'] <= 95) & (bulk_melt['hbl_SiO2'].notna()) 
                        & (bulk_melt['bulk_SiO2'] >= 40) & (bulk_melt['bulk_SiO2'] <= 80)]

resampling_weights = (cry['bulk_SiO2'] - 40) / 40
cry = cry.sample(n=10000, weights=resampling_weights, replace=True)
cry_si_t = cry[['bulk_SiO2', 'T(K)']]
cry_si_t.loc[:, 'T(K)'] = cry_si_t['T(K)'] - 273.15
cry_si_t.loc[:,'silica_bins'] = pd.cut(cry_si_t['bulk_SiO2'], binedges)
cry_avg = cry_si_t.groupby(['silica_bins'], observed=True).mean()
cry_sem = cry_si_t.groupby(['silica_bins'], observed=True).sem()

# %% make reeq
reeq = bulk_melt.loc[(bulk_melt['melt_vol'] > 0) & (bulk_melt['melt_vol'] <= 5) & (bulk_melt['hbl_SiO2'].notna()) 
                     & (bulk_melt['bulk_SiO2'] >= 40) & (bulk_melt['bulk_SiO2'] <= 80)]

resampling_weights = (reeq['bulk_SiO2'] - 40) / 40
reeq = reeq.sample(n=10000, weights=resampling_weights, replace=True)
# reeq = pd.read_csv('bulk_melt_resample_si_0-5_good.csv')
si_t = reeq[['bulk_SiO2', 'T(K)']]
si_t.loc[:, 'T(K)'] = si_t['T(K)'] - 273.15
si_t.loc[:,'silica_bins'] = pd.cut(si_t['bulk_SiO2'], binedges)
avg = si_t.groupby(['silica_bins'], observed=True).mean()
sem = si_t.groupby(['silica_bins'], observed=True).sem()
# %% melt extraction
melt = bulk_melt.loc[(bulk_melt['melt_vol'] >= 30) & (bulk_melt['melt_vol'] <= 50) & (bulk_melt['hbl_SiO2'].notna()) 
                        & (bulk_melt['bulk_SiO2'] >= 40) & (bulk_melt['bulk_SiO2'] <= 80)]

resampling_weights = (melt['bulk_SiO2'] - 40) / 40
melt = melt.sample(n=10000, weights=resampling_weights, replace=True)
# calculate the volcanic part
melt_vol_si_t = melt[['melt_SiO2', 'T(K)']]
melt_vol_si_t.loc[:, 'T(K)'] = melt_vol_si_t['T(K)'] - 273.15
melt_vol_si_t.loc[:,'silica_bins'] = pd.cut(melt_vol_si_t['melt_SiO2'], binedges)

melt_vol_avg = melt_vol_si_t.groupby(['silica_bins'], observed=True).mean()
melt_vol_sem = melt_vol_si_t.groupby(['silica_bins'], observed=True).sem()
# add ~30% melts to crystals
melt_pst = np.random.normal(30, 3, size=len(melt)) / 100
crystal_pst = 1 - melt_pst
# build a DF to store melts (Mixed with a few crystals)
bulk_elements = ['SiO2', 'TiO2', 'Al2O3', 'FeO', 'MgO', 'CaO', 'K2O']
for e in bulk_elements:
    melt['residue_'+e] = melt['melt_'+e] * melt_pst + melt['solid_'+e] * crystal_pst
melt_plu_si_t = melt[['melt_SiO2', 'residue_SiO2', 'T(K)']]
melt_plu_si_t.loc[:, 'T(K)'] = melt_plu_si_t['T(K)'] - 273.15
# melt_plu_si_t.loc[:,'silica_bins'] = pd.cut(melt_plu_si_t['residue_SiO2'], binedges)

melt_plu_si_t.loc[:,'silica_bins'] = pd.cut(melt_plu_si_t['melt_SiO2'], binedges)
melt_plu_avg = melt_plu_si_t.groupby(['silica_bins'], observed=True).mean()
melt_plu_sem = melt_plu_si_t.groupby(['silica_bins'], observed=True).sem()
# %% no melt extraction volcanic (with hbl)
# no_melt_w_hbl = bulk_melt.loc[(bulk_melt['melt_vol'] >= 30) & (bulk_melt['melt_vol'] <= 60) & (bulk_melt['hbl_SiO2'].notna()) & (bulk_melt['bulk_SiO2'] >= 40) & (bulk_melt['bulk_SiO2'] <= 80)]
# resampling_weights = (no_melt_w_hbl['bulk_SiO2'] - 40) / 40
# no_melt_w_hbl = no_melt_w_hbl.sample(n=10000, weights=resampling_weights, replace=True)
no_melt_w_hbl = pd.read_csv('no_melt_extraction_vol_36-64.csv')
# resampling_weights = (no_melt_w_hbl['bulk_SiO2'] - 40) / 40
# no_melt_w_hbl = no_melt_w_hbl.sample(n=10000, weights=resampling_weights, replace=True)
# no_melt_w_hbl = no_melt_w_hbl.sample(n=10000, replace=True)
no_melt_w_hbl_si_t = no_melt_w_hbl[['bulk_SiO2', 'T(K)']]
no_melt_w_hbl_si_t.loc[:, 'T(K)'] = no_melt_w_hbl_si_t['T(K)'] - 273.15
no_melt_w_hbl_si_t.loc[:,'silica_bins'] = pd.cut(no_melt_w_hbl_si_t['bulk_SiO2'], binedges)
no_melt_w_hbl_avg = no_melt_w_hbl_si_t.groupby(['silica_bins'], observed=True).median()
no_melt_w_hbl_sem = no_melt_w_hbl_si_t.groupby(['silica_bins'], observed=True).std()
# %% no melt extraction volcanic (without hbl)
# no_me_vol1 = pd.read_csv('hbl_bulk_georoc_vol_all_highT.csv')
# no_me_vol2 = pd.read_csv('hbl_bulk_georoc_plu_all_highT.csv')
# no_me_vol = pd.concat([no_me_vol1,no_me_vol2])

# no_me_vol = no_me_vol.loc[(no_me_vol['melt_vol'] >= 90) & (no_me_vol['melt_vol'] <= 98.5) & (no_me_vol['bulk_SiO2'] >= 40) & (no_me_vol['bulk_SiO2'] <= 80)]

# no_me_vol = no_me_vol.sample(n=10000, replace=True)
# no_me_vol_si_t = no_me_vol[['bulk_SiO2', 'T(K)']]
# no_me_vol_si_t.loc[:, 'T(K)'] = no_me_vol_si_t['T(K)'] - 273.15
# no_me_vol_si_t.loc[:,'silica_bins'] = pd.cut(no_me_vol_si_t['bulk_SiO2'], binedges)
# no_me_vol_avg = no_me_vol_si_t.groupby(['silica_bins'], observed=True).median()
# no_me_vol_sem = no_me_vol_si_t.groupby(['silica_bins'], observed=True).sem()

# %% read natural samples
natural_samples = pd.read_csv('vol_plu_all_higgins2021_cleaned.csv')
natural_vol = natural_samples[natural_samples['ROCK TYPE'] == 'VOL']
natural_plu = natural_samples[natural_samples['ROCK TYPE'] == 'PLU']

natural_mean = pd.read_csv('resampled_bined_vol_plu_mean.csv')
natural_err = pd.read_csv('resampled_bined_vol_plu_err.csv')
# %% plot
fig, ax = plt.subplots(figsize=(8, 12))
ax2 = ax.twinx()

# sns.scatterplot(data=si_t, x='bulk_SiO2', y='T(K)',
#                 s=1, alpha=0.1, c='blue', legend=False, ax=ax)

# sns.scatterplot(data=melt_plu_si_t, x='residue_SiO2', y='T(K)',
#                 s=1, alpha=0.1, c='black', legend=False, ax=ax)

# sns.scatterplot(data=melt_vol_si_t, x='melt_SiO2', y='T(K)',
#                 s=1, alpha=0.1, c='red', legend=False, ax=ax)

# sns.scatterplot(data=no_melt_w_hbl_si_t, x='bulk_SiO2', y='T(K)',
#                 s=1, alpha=0.1, c='pink', legend=False, ax=ax)


# ax.errorbar(x=avg['bulk_SiO2'], y=avg['T(K)'], yerr=sem['T(K)']*2*3.5,
#             fmt="o", capsize=4, markersize=3, c='blue', label='Reequilibrium')

# ax.errorbar(x=melt_plu_avg['residue_SiO2'], y=melt_plu_avg['T(K)'], yerr=melt_plu_sem['T(K)']*2*3.5,
#             fmt="o", capsize=4, markersize=3, c='black', label='Melt extraction residue')

# ax.errorbar(x=melt_vol_avg['melt_SiO2'], y=melt_vol_avg['T(K)'], yerr=melt_vol_sem['T(K)']*2*3.5,
#             fmt="o", capsize=4, markersize=3, c='red', label='Melt extraction melt')

# ax.errorbar(x=no_melt_w_hbl_avg['bulk_SiO2'], y=no_melt_w_hbl_avg['T(K)'], yerr=no_melt_w_hbl_sem['T(K)'],
#             fmt="o", capsize=4, markersize=3, c='pink', label='No melt extraction melt (with Hbl)')

sns.scatterplot(data=natural_samples, x='SIO2(WT%)', y='T', hue='ROCK TYPE', s=4, alpha=0.3, legend='auto', palette='RdBu', ax=ax2, zorder=1)
# ax2.scatter(x=natural_vol['SIO2(WT%)'], y=natural_vol['T'], s=2, alpha=0.03, label='Volcanic', color='red', zorder=1)
# ax2.scatter(x=natural_plu['SIO2(WT%)'], y=natural_plu['T'], s=2, alpha=0.03, label='Plutonic', color='blue', zorder=1)
for row in range(melt_plu_avg.shape[0]):
    sio2_list = [melt_plu_avg.iloc[row, 0], melt_plu_avg.iloc[row, 1]]
    temp_list = [melt_plu_avg.iloc[row, 2], melt_plu_avg.iloc[row, 2]]
    ax.plot(sio2_list, temp_list, c='gray')
ax.fill_between(melt_plu_avg['residue_SiO2'], melt_plu_avg['T(K)']-melt_plu_sem['T(K)']*2*3.5, melt_plu_avg['T(K)']+melt_plu_sem['T(K)']*2*3.5, color = 'gray', alpha=.2, label='Plutonic - melt extraction residue', linewidth=0, zorder=2)
ax.fill_between(cry_avg['bulk_SiO2'], cry_avg['T(K)']-cry_sem['T(K)']*2, cry_avg['T(K)']+cry_sem['T(K)']*2, color = '#134857', alpha=.3, label='Plutonic - crystallization only', linewidth=0, zorder=2)
ax.fill_between(avg['bulk_SiO2'], avg['T(K)']-sem['T(K)']*2*3.5, avg['T(K)']+sem['T(K)']*2*3.5, alpha=.4, label='Plutonic - reequilibration', linewidth=0, zorder=2)
ax.fill_between(melt_vol_avg['melt_SiO2'], melt_vol_avg['T(K)']-melt_vol_sem['T(K)']*2*3.5, melt_vol_avg['T(K)']+melt_vol_sem['T(K)']*2*3.5, color = 'orange', alpha=.2, label='Volcanic - extracted melt', linewidth=0, zorder=2)
ax.fill_between(no_melt_w_hbl_avg['bulk_SiO2'], no_melt_w_hbl_avg['T(K)']-no_melt_w_hbl_sem['T(K)'], no_melt_w_hbl_avg['T(K)']+no_melt_w_hbl_sem['T(K)'], color = 'red', alpha=.3, label='Volcanic - no melt extraction', linewidth=0, zorder=2)


ax2.errorbar(
            x=natural_mean[natural_mean['ROCK TYPE'] == 'VOL']['SIO2(WT%)'],
            y=natural_mean[natural_mean['ROCK TYPE'] == 'VOL']['T'],
            yerr=natural_err[natural_err['ROCK TYPE'] == 'VOL']['T'],
            fmt="o",
            capsize=4,
            label='Volcanic',
            markersize=3,
            color='#E60000',
            zorder=3
        )
ax2.errorbar(
            x=natural_mean[natural_mean['ROCK TYPE'] == 'PLU']['SIO2(WT%)'],
            y=natural_mean[natural_mean['ROCK TYPE'] == 'PLU']['T'],
            yerr=natural_err[natural_err['ROCK TYPE'] == 'PLU']['T'],
            fmt="o",
            capsize=4,
            label='Plutonic',
            markersize=3,
            color='#003153',
            zorder=3
        )
# ax2.plot(
#     natural_mean[natural_mean['ROCK TYPE'] == 'VOL']['SIO2(WT%)'],
#     natural_mean[natural_mean['ROCK TYPE'] == 'VOL']['T'],
#     label='Natural volcanic',
#     markersize=3,
#     color='#E60000',
#     alpha = 0.5
# )
# ax2.plot(
#     natural_mean[natural_mean['ROCK TYPE'] == 'PLU']['SIO2(WT%)'],
#     natural_mean[natural_mean['ROCK TYPE'] == 'PLU']['T'],
#     label='Natural plutonic',
#     markersize=3,
#     color='#003153',
#     alpha=0.5
# )

ax.legend(title='Simulated', loc='upper left', bbox_to_anchor=(0.0, 1.0))
ax2.legend(title='Natural', loc='upper right', bbox_to_anchor=(1, 1.0),)

ax.set_xlabel('$SiO_2\;[wt.\%]$')
ax.set_ylabel('Simulated magma temperature [$^\circ$C]')
ax2.set_ylabel('Amphibole temperature [$^\circ$C]')

plt.xlim(40, 80)
ax.set_ylim(650, 1100)
ax2.set_ylim(750, 1200)
# plt.savefig('fig/pluton_silica_temperature_compare_natural_w_modeling.svg')
plt.savefig('fig/pluton_silica_temperature_compare_natural_w_modeling_another_extraction.svg')

plt.show
# %%