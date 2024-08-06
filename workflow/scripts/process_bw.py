import logging 
import sys
import gc
import os
import pandas as pd
import sqlite3
import numpy as np
import pyBigWig

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
sample_sheet = snakemake.params['sample_sheet']
ref_size_fn = snakemake.params['ref_size_sheet']
bw_prefix = snakemake.params['bw_prefix']

if not os.path.exists(gatc_db):
    raise ValueError(f'{gatc_db} does not exist')
if not os.path.exists(sample_sheet):
    raise ValueError(f'{sample_sheet} does not exist')
if not os.path.exists(ref_size_fn):
    raise ValueError(f'{ref_size_fn} does not exist')

filter_samples = []
chrom_sizes = {}

with open(sample_sheet, 'r') as f:
    for line in f:
        sample = line.strip().split('\t')[1]
        filter_samples.append(f"Filter_{sample}")

with open(ref_size_fn) as f:
    for line in f:
        line = line.strip().split('\t')
        chrom_sizes[line[0]] = int(line[-1])

gatc_con = sqlite3.connect(gatc_db)

for sample in filter_samples:
    bw_fn = f"{bw_prefix}.{sample}.bw"

    with pyBigWig.open(bw_fn, 'w') as bw_p:
        header = [(c, s) for (c, s) in chrom_sizes.items()]
        bw_p.addHeader(header)
        tmp_dpni = pd.read_sql_query(f'select chrom, Apos, DpnI from "{sample}"', gatc_con)
        tmp_dpni.dropna(inplace=True)
        for chrom, size in header:
            tmp_pd = tmp_dpni.loc[tmp_dpni['chrom'] == chrom]
            pos = tmp_pd['Apos'].values
            dpni = tmp_pd['DpnI'].values
            chroms = np.repeat(chrom, len(pos))
            bw_p.addEntries(chroms, pos, ends=pos + 1, values=dpni)
        logger.info('%s done', sample)
gatc_con.close()