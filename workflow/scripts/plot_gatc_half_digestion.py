import logging 
import sys
import pandas as pd
import numpy as np
import sqlite3
import seaborn as sns
import matplotlib.pyplot as plt
import os
import gc

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
streamHandler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter("[%(asctime)s] %(levelname)s:%(name)s:#%(lineno)d:%(message)s")
streamHandler.setFormatter(formatter)
logger.addHandler(streamHandler)
fileHandler = logging.FileHandler(snakemake.log[0])
fileHandler.setFormatter(formatter)
logger.addHandler(fileHandler)

gatc_db = snakemake.params['data_db']
filter_size = int(snakemake.params['filter_size'])
sample_sheet = snakemake.params['sample_sheet']
prefix = snakemake.params['prefix']

if not os.path.exists(gatc_db):
    raise ValueError(f'{gatc_db} does not exist')
if not os.path.exists(sample_sheet):
    raise ValueError(f'{sample_sheet} does not exist')
if not os.path.exists('results/figures'):
    os.makedirs('results/figures')
    logger.info(f'Created results/figures')
if not os.path.exists('results/tables'):
    os.makedirs('results/tables')
    logger.info(f'Created results/tables')

gatc_con = sqlite3.connect(gatc_db)


samples = []
with open(sample_sheet, 'r') as f:
    for line in f:
        samples.append(line.strip().split('\t')[-1])

tmp_data = []
tmp_filter_data = []
tmp_occ_data = []
tmp_occ_filter_data = []

for sample in samples:
    dpni_pd = pd.read_sql_query('select * from "{}"'.format(sample), gatc_con)
    pos_dist = dpni_pd['Apos'].diff()
    pos_dist[pos_dist < 0] = np.nan
    dpni_pd['Left_Flank'] = pos_dist
    dpni_pd['Right_Flank'] = pos_dist.shift(-1)
    # remove CpG overlap
    dpni_pd = dpni_pd.loc[dpni_pd['is_overlap'] == 0] 
    # stats for the filtered locations
    tmp = dpni_pd.loc[(dpni_pd['Left_Flank'] >= filter_size)]['DpnI_Left'].values  
    tmp_occ = dpni_pd.loc[(dpni_pd['Left_Flank'] >= filter_size)]['Occ_G'].values
    tmp = tmp[~np.isnan(tmp)]
    tmp_occ = tmp_occ[~np.isnan(tmp_occ)]
    if len(tmp) > 0:
        filter_1l = np.quantile(tmp, 0.25)
        filter_2l = np.quantile(tmp, 0.5)
        filter_3l = np.quantile(tmp, 0.75)
    else:
        filter_1l = np.nan
        filter_2l = np.nan
        filter_3l = np.nan
    if len(tmp_occ) > 0:
        filter_1l_occ = np.quantile(tmp_occ, 0.25)
        filter_2l_occ = np.quantile(tmp_occ, 0.5)
        filter_3l_occ = np.quantile(tmp_occ, 0.75)
    else:
        filter_1l_occ = np.nan
        filter_2l_occ = np.nan
        filter_3l_occ = np.nan
    tmp_occ_g = tmp_occ
    del tmp
    tmp = dpni_pd.loc[(dpni_pd['Right_Flank'] >= filter_size)]['DpnI_Right'].values
    tmp = tmp[~np.isnan(tmp)]
    tmp_occ = dpni_pd.loc[(dpni_pd['Right_Flank'] >= filter_size)]['Occ_C'].values
    if len(tmp) > 0:
        filter_1r = np.quantile(tmp, 0.25)
        filter_2r = np.quantile(tmp, 0.5)
        filter_3r = np.quantile(tmp, 0.75)
    else:
        filter_1r = np.nan
        filter_2r = np.nan
        filter_3r = np.nan
    if len(tmp_occ) > 0:
        filter_1r_occ = np.quantile(tmp_occ, 0.25)
        filter_2r_occ = np.quantile(tmp_occ, 0.5)
        filter_3r_occ = np.quantile(tmp_occ, 0.75)
    else:
        filter_1r_occ = np.nan
        filter_2r_occ = np.nan
        filter_3r_occ = np.nan
    tmp_occ_c = tmp_occ
    del tmp
    tmp_filter_data.append([sample, filter_1l, filter_2l, filter_3l, filter_1r, filter_2r, filter_3r])
    tmp_occ_filter_data.append([sample, filter_1l_occ, filter_2l_occ, filter_3l_occ, filter_1r_occ, filter_2r_occ, filter_3r_occ])

    # draw Occ_G and Occ_C plots after filtering
    max_occ = max(tmp_occ_g.max(), tmp_occ_c.max())
    if max_occ > 60:
        max_occ = 60
    g = sns.histplot(data={'Occ-G': tmp_occ_g, 'Occ-C': tmp_occ_c}, 
                     bins=range(0, max_occ+1, 1), element='step', fill=False)
    g.set(xlim=(0, max_occ))
    g.set(xlabel='Occupancy', ylabel='Number of GATC Sites')
    g.set_title(f'{sample} Occupancy with Flank >= {filter_size} bp')
    plt.tight_layout()
    plt.savefig(f'results/figures/{prefix}_{sample}_Occupancy_Filtered_GATC.png', bbox_inches='tight',
                transparent=False, dpi=300, facecolor='white')
    plt.close()

    del tmp_occ_g
    del tmp_occ_c
    gc.collect()

    for i in range(1, 601):
        tmp = dpni_pd.loc[(dpni_pd['Left_Flank'] == i)]['DpnI_Left'].values
        tmp = tmp[~np.isnan(tmp)]
        if len(tmp) > 0:
            q1l = np.quantile(tmp, 0.25)
            q2l = np.quantile(tmp, 0.5)
            q3l = np.quantile(tmp, 0.75)
        else:
            q1l = np.nan
            q2l = np.nan
            q3l = np.nan

        tmp = dpni_pd.loc[(dpni_pd['Right_Flank'] == i)]['DpnI_Right'].values
        tmp = tmp[~np.isnan(tmp)]        
        if len(tmp) > 0:
            q1r = np.quantile(tmp, 0.25)
            q2r = np.quantile(tmp, 0.5)
            q3r = np.quantile(tmp, 0.75)
        else:
            q1r = np.nan
            q2r = np.nan
            q3r = np.nan
        tmp_data.append([sample, i, q1l, q2l, q3l, q1r, q2r, q3r])
        del tmp
    del dpni_pd
    gc.collect()

tmp_data = pd.DataFrame(tmp_data, columns=['sample', 'FlankLength', 'q1l', 'q2l', 'q3l', 'q1r', 'q2r', 'q3r'])

tmp_data.to_csv('results/tables/{}_flank_quantile.csv'.format(prefix), index=False)

tmp_filter_data = pd.DataFrame(tmp_filter_data, columns=['sample', 'filter_q1l', 'filter_q2l', 'filter_q3l', 'filter_q1r', 'filter_q2r', 'filter_q3r'])

tmp_occ_filter_data = pd.DataFrame(tmp_occ_filter_data, columns=['sample', 'filter_q1l_occG', 'filter_q2l_occG', 'filter_q3l_occG', 
                                                                 'filter_q1r_occC', 'filter_q2r_occC', 'filter_q3r_occC'])

tmp_filter_data.to_csv('results/tables/{}_filter_quantile.csv'.format(prefix), index=False)
tmp_occ_filter_data.to_csv('results/tables/{}_filter_occ_quantile.csv'.format(prefix), index=False)
logger.info(f'Finished reading {prefix} for half site digestion statistics')

# plot the figures

# set color pellete to Tab10
sns.set_palette(sns.color_palette("tab10"))
for sample in samples:
    tmp = tmp_data.loc[tmp_data['sample'] == sample]
    tmpm = tmp['q2l'].values    
    plt.plot(tmp['FlankLength'], tmp['q2l'], label='{}'.format(sample), lw=1, zorder=1)
plt.axvline(x=filter_size, color='grey', ls='--', zorder=7, lw=1)
plt.xlabel('Left Flanking Distance (bp)')
plt.ylabel('Left DpnI Digestion (%)')
plt.xlim(1, 600)
plt.ylim(-1, 100)
plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1), frameon=False)
plt.tight_layout()
plt.savefig(f'results/figures/{prefix}_DpnI_Digestion_LeftHalf.png', bbox_inches='tight',
            transparent=False, dpi=300, facecolor='white')
plt.close()

for sample in samples:
    tmp = tmp_data.loc[tmp_data['sample'] == sample]
    tmpm = tmp['q2r'].values    
    plt.plot(tmp['FlankLength'], tmp['q2r'], label='{}'.format(sample), lw=1, zorder=1)
plt.axvline(x=filter_size, color='grey', ls='--', zorder=7, lw=1)
plt.xlabel('Right Flanking Distance (bp)')
plt.ylabel('Right DpnI Digestion (%)')
plt.xlim(1, 600)
plt.ylim(-1, 100)
plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1), frameon=False)
plt.tight_layout()
plt.savefig(f'results/figures/{prefix}_DpnI_Digestion_RightHalf.png', bbox_inches='tight',
            transparent=False, dpi=300, facecolor='white')
plt.close()
logger.info(f'Finished plotting {prefix} for half site digestion statistics')

gatc_con.close()