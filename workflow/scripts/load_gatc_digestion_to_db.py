import os
import pickle
import logging 
import sys
import pandas as pd
import copy
import numpy as np
import sqlite3
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
dirn = snakemake.params['bamdir']
chrom_fn = snakemake.params['chrom_size']

if not os.path.exists(dirn):
    raise ValueError(f'{dirn} does not exist')
if not os.path.exists(chrom_fn):
    raise ValueError(f'{chrom_fn} does not exist')
if not os.path.exists(gatc_db):
    raise ValueError(f'{gatc_db} does not exist')

sample_sheet = snakemake.params['sample_sheet']

chrom_size_pd = pd.read_csv(chrom_fn, sep='\t', header=None, 
                            names=['Chrom', 'Ref', 'Size'])
chrom_sizes = {}
for _, row in chrom_size_pd.iterrows():
    chrom_sizes[row['Ref']] = (row['Size'], row['Chrom'])


gatc_con = sqlite3.connect(gatc_db)
gatc_pos = pd.read_sql('select * from gatc_pos', gatc_con)

if sample_sheet == 'None':
    sample_func = lambda x: x
else:
    sample_pd = pd.read_csv(sample_sheet, sep='\t', header=None, 
                            names=['Prefix', 'Sample'])
    sample_dict = {}
    for _, row in sample_pd.iterrows():
        sample_dict[row['Prefix']] = row['Sample']
    sample_func = lambda x: sample_dict[x] if x in sample_dict else None


for fn in sorted(os.listdir(dirn)):
    if not fn.endswith('occ.pickle'):
        continue
    prefix = fn.split('.')[0]
    sample = sample_func(prefix)
    if sample is None:
        logger.warn(f'{prefix} not in sample sheet')
        continue

    # load occupancy data
    occ_fn = os.path.join(dirn, fn)
    with open(occ_fn, 'rb') as f:
        occ_dict = pickle.load(f)
    # load left site data
    left_fn = os.path.join(dirn, prefix + '.gatc.left.pickle')
    right_fn = os.path.join(dirn, prefix + '.gatc.right.pickle')
    with open(os.path.join(dirn, prefix + '.gatc.left.pickle'), 'rb') as f:
        left_dict = pickle.load(f)
    # load right site data
    with open(os.path.join(dirn, prefix + '.gatc.right.pickle'), 'rb') as f:
        right_dict = pickle.load(f)
    tmp_dpni_pd = []
    for ref, (size, chrom) in chrom_sizes.items():
        tmp_pd = copy.deepcopy(gatc_pos[gatc_pos['chrom'] == chrom])
        tmp_pd['Gpos'] = tmp_pd['Apos'] - 1
        tmp_pd['Tpos'] = tmp_pd['Apos'] + 1
        tmp_pd['Cpos'] = tmp_pd['Apos'] + 2

        for nt in 'GA':
            occ_nt = tmp_pd['{}pos'.format(nt)]
            occ = occ_dict[ref][occ_nt]
            cut = left_dict[ref][occ_nt]
            tmp_pd['Occ_{}'.format(nt)] = occ
            tmp_pd['Cut_{}'.format(nt)] = cut 
        for nt in 'TC':
            occ_nt = tmp_pd['{}pos'.format(nt)]
            occ = occ_dict[ref][occ_nt]
            cut = right_dict[ref][occ_nt]
            tmp_pd['Occ_{}'.format(nt)] = occ
            tmp_pd['Cut_{}'.format(nt)] = cut 
        left_cut = tmp_pd['Cut_G'] + tmp_pd['Cut_A']
        right_cut = tmp_pd['Cut_T'] + tmp_pd['Cut_C']
        left_occ = tmp_pd['Occ_G']
        right_occ = tmp_pd['Occ_C']
        zeroleft = (left_occ == 0)
        zeroright = (right_occ == 0)         
        left = left_cut / left_occ * 100
        left[zeroleft] = np.nan 
        right = right_cut / right_occ * 100
        right[zeroright] = np.nan
        tmp_pd['DpnI_Left'] = left
        tmp_pd['DpnI_Right'] = right
        valboth = ~np.logical_or(zeroleft, zeroright)
        tmp_pd['DpnI_Both'] = np.nan 
        tmp_pd.loc[valboth, 'DpnI_Both'] = (left[valboth] + right[valboth])/2
        valleft = np.logical_and(~zeroleft, ~valboth)
        valright = np.logical_and(~zeroright, ~valboth)
        tmp_pd.loc[valleft, 'DpnI_Both']  = left[valleft]
        tmp_pd.loc[valright, 'DpnI_Both']  = right[valright]
        for nt in 'TCG':
            tmp_pd = tmp_pd.drop(['{}pos'.format(nt)], axis=1)
        tmp_dpni_pd.append(tmp_pd)
    dpni_pd = pd.concat(tmp_dpni_pd)
    dpni_pd.to_sql(sample, gatc_con, if_exists='append', index=False)
    logger.info('{} Finished processing {}'.format(prefix, sample))
   
    del occ_dict
    del left_dict
    del right_dict
    gc.collect()
    os.remove(occ_fn)
    os.remove(left_fn)
    os.remove(right_fn)
gatc_con.close()