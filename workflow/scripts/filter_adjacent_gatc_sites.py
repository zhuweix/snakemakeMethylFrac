import logging 
import sys
import pandas as pd
import gc
import sqlite3
import numpy as np
import os

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
streamHandler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter("[%(asctime)s] %(levelname)s:%(name)s:#%(lineno)d:%(message)s")
streamHandler.setFormatter(formatter)
logger.addHandler(streamHandler)
fileHandler = logging.FileHandler(snakemake.log[0])
fileHandler.setFormatter(formatter)
logger.addHandler(fileHandler)


gatc_db = snakemake.params['gatc_db']
flank_limit = int(snakemake.params['filter_size'])
sample_sheet = snakemake.params['sample_sheet']
filter_gatc = snakemake.params['filter_overlap']
if filter_gatc not in [True, False]:
    raise ValueError(f'filter_gatc must be True or False, not {filter_gatc}')

if not os.path.exists(gatc_db):
    raise ValueError(f'{gatc_db} does not exist')
if not os.path.exists(sample_sheet):
    raise ValueError(f'{sample_sheet} does not exist')

samples = []
with open(sample_sheet, 'r') as f:
    for line in f:
        samples.append(line.strip().split('\t')[-1])

gatc_con = sqlite3.connect(gatc_db)

for sample in samples:
    # test is sample in the database
    if not gatc_con.execute('select count(*) from sqlite_master where type="table" and name="{}"'.format(sample)).fetchone()[0]:
        logger.warn(f'{sample} not in database')
        continue
    dpni_pd = pd.read_sql_query('select * from "{}"'.format(sample), gatc_con)

    left_dist = dpni_pd['Apos'].diff()
    left_dist[left_dist < 0] = np.nan
    right_dist = left_dist.shift(-1)  
    dpni_pd['Left_Flank'] = left_dist
    dpni_pd['Right_Flank'] = right_dist
    if filter_gatc:
        no_cpg = (dpni_pd['is_overlap'] == 0)
        dpni_pd = dpni_pd.loc[no_cpg]  
    
    tmp_dpn = dpni_pd['DpnI_Both'].values
    tmp_occ = (dpni_pd['Occ_G'].values + dpni_pd['Occ_C'].values)/2
    mask_left = (dpni_pd['Left_Flank'] < flank_limit)
    mask_right = (dpni_pd['Right_Flank'] < flank_limit)
    mask_both = ~(mask_left & mask_right)
    tmp_dpn[mask_left] = dpni_pd.loc[mask_left]['DpnI_Right'].values
    tmp_occ[mask_left] = dpni_pd.loc[mask_left]['Occ_C'].values
    tmp_dpn[mask_right] = dpni_pd.loc[mask_right]['DpnI_Left'].values
    tmp_occ[mask_right] = dpni_pd.loc[mask_right]['Occ_G'].values
    tmp_dpn = tmp_dpn[mask_both]
    tmp_occ = tmp_occ[mask_both]
    dpni_pd = dpni_pd.loc[mask_both]
    tmp_digest = pd.DataFrame({
        'chrom': dpni_pd['chrom'].values,
        'Apos': dpni_pd['Apos'].values,        
        'DpnI': tmp_dpn,
        'Occ': tmp_occ,
    })
    tmp_digest.to_sql('Filter_{}'.format(sample), gatc_con, if_exists='replace', index=False)
    logger.info(f'Finished {sample}')
    del tmp_digest
    del dpni_pd
    gc.collect()

gatc_con.close()
