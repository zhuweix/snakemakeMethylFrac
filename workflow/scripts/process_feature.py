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



anno_db = snakemake.input[1]

gatc_db = snakemake.params['gatc_db']
flank = int(snakemake.params['flank_size'])
sample_sheet = snakemake.params['sample_sheet']
ref_size_fn = snakemake.params['ref_size_sheet']
org = snakemake.params['org']

if not os.path.exists(gatc_db):
    raise ValueError(f'{gatc_db} does not exist')
if not os.path.exists(anno_db):
    raise ValueError(f'{anno_db} does not exist')
if not os.path.exists(sample_sheet):
    raise ValueError(f'{sample_sheet} does not exist')
if not os.path.exists(ref_size_fn):
    raise ValueError(f'{ref_size_fn} does not exist')
if not org in ['MCF7', 'MCF10']:
    raise ValueError(f'{org} is not supported, only MCF7 and MCF10 are supported')

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
anno_con = sqlite3.connect(anno_db)
if org == 'MCF7':
    gene_pd = pd.read_sql_query('SELECT * from "Gene_Best_Transcript_GSE201262"', anno_con)
elif org == 'MCF10':
    gene_pd = pd.read_sql_query('SELECT * from "Gene_Best_Transcript_GSE237066"', anno_con)
other_feature_pd = pd.read_sql_query('SELECT * from "other_features"', anno_con)
anno_con.close()

tmp_gene_pd = []
tmp_other_pd = []

for chrom in chrom_sizes:   
    for sample in filter_samples:
        # test is sample in the database
        if not gatc_con.execute(f'select count(*) from sqlite_master where type="table" and name="{sample}"').fetchone()[0]:
            logger.warn(f'{sample} not in database')
            continue
        chrom_loc = (gene_pd['Chrom'] == chrom)
        tmp_dpni = pd.read_sql_query('SELECT * from "{}" where chrom="{}"'.format(sample, chrom), gatc_con)   
        sname = sample.replace('Filter_', '')
        tmp_pd = gene_pd.loc[chrom_loc]
        for ent in tmp_pd.itertuples():
            start = ent.Start
            end = ent.End
            gid = ent.ID
            tmp = tmp_dpni.loc[(tmp_dpni['Apos'] >= start -flank) & (tmp_dpni['Apos'] < end + flank)] 
            tmp_genebody_loc = (tmp['Apos'] >= start) & (tmp['Apos'] < end)
            tmp_ratio = tmp['DpnI'].to_numpy()
            tmp_occ = tmp['Occ'].to_numpy()
            tmp_apos = tmp['Apos'].to_numpy()
            # change to relative position to start site
            if ent.Strand == "+":
                tmp_apos = tmp_apos - start
            elif ent.Strand == "-":
                tmp_apos = end - 1 - tmp_apos
            tmp_chrom = chrom
            tmp_gene_pd.append(pd.DataFrame(
                {'sample': sname, 'feature': 'genebody', 'DpnI': tmp_ratio[tmp_genebody_loc], 
                 'Occ':tmp_occ[tmp_genebody_loc], 'geneid': gid,
                 'Apos': tmp_apos[tmp_genebody_loc]}
            ))
            if ent.Strand == "+":
                tmp_gene = tmp.loc[(tmp['Apos'] >= start) & (tmp['Apos'] < start + flank)]
                tmp_ratio = tmp_gene['DpnI'].to_numpy()
                tmp_occ = tmp_gene['Occ'].to_numpy()  
                tmp_apos = tmp_gene['Apos'].to_numpy() - start       
                tmp_gene_pd.append(pd.DataFrame(
                    {'sample': sname, 'feature': f'gene_{flank}bp', 'DpnI': tmp_ratio, 'Occ':tmp_occ, 'geneid': gid,
                     'Apos': tmp_apos}
                ))            
                tmp_promoter = tmp.loc[(tmp['Apos'] >= start - flank) & 
                                            (tmp['Apos'] < start )] 
                tmp_ratio = tmp_promoter['DpnI'].to_numpy()
                tmp_occ = tmp_promoter['Occ'].to_numpy()
                tmp_apos = tmp_promoter['Apos'].to_numpy() - start
                tmp_gene_pd.append(pd.DataFrame(
                    {'sample': sname, 'feature': f'promoter_{flank}bp', 'DpnI': tmp_ratio, 'Occ':tmp_occ,
                     'geneid': gid, 'Apos': tmp_apos}
                ))     
            elif ent.Strand == "-":
                tmp_gene = tmp.loc[(tmp['Apos'] >= end -flank) & (tmp['Apos'] < end)] 
                tmp_ratio = tmp_gene['DpnI'].to_numpy()
                tmp_occ = tmp_gene['Occ'].to_numpy()
                tmp_apos = end - 1 - tmp_gene['Apos']
                tmp_gene_pd.append(pd.DataFrame(
                    {'sample': sname, 'feature': f'gene_{flank}bp', 'DpnI': tmp_ratio, 'Occ':tmp_occ, 
                     'geneid': gid, 'Apos': tmp_apos}
                ))            
                tmp_promoter = tmp.loc[(tmp['Apos'] >= end) & 
                                            (tmp['Apos'] < end + flank )] 
                tmp_ratio = tmp_promoter['DpnI'].to_numpy()
                tmp_occ = tmp_promoter['Occ'].to_numpy()
                tmp_apos = end - 1 - tmp_promoter['Apos']
                tmp_gene_pd.append(pd.DataFrame(
                    {'sample': sname, 'feature': f'promoter_{flank}bp', 'DpnI': tmp_ratio, 'Occ':tmp_occ, 
                     'geneid': gid, 'Apos': tmp_apos}
                ))

        chrom_loc = (other_feature_pd['Chrom'] == chrom)
        tmp_pd = other_feature_pd.loc[chrom_loc]  
        for ent in tmp_pd.itertuples():
            start = ent.Start
            end = ent.End
            tp = ent.Type
            feat_id = ent.FeatureIndex
            tmp = tmp_dpni.loc[(tmp_dpni['Apos'] >= start) & (tmp_dpni['Apos'] < end)] 
            tmp_ratio = tmp['DpnI'].to_numpy()
            tmp_occ = tmp['Occ'].to_numpy()
            tmp_apos = tmp['Apos'].to_numpy()
            if ent.Strand == "+":
                tmp_apos = tmp_apos - start
            else:
                tmp_apos = end - 1 - tmp_apos
            tmp_other_pd.append(pd.DataFrame(
                {'sample': sname, 'feature': tp, 'DpnI': tmp_ratio, 'Occ':tmp_occ, 'FeatureIndex': feat_id,
                 'Apos': tmp_apos}
            ))                   
    logger.info(f'{chrom} done')
    gc.collect()
tmp_gene_pd = pd.concat(tmp_gene_pd)
tmp_gene_pd.to_sql('gene_raw_site', gatc_con, if_exists='replace', index=False)
tmp_other_pd = pd.concat(tmp_other_pd)
tmp_other_pd.to_sql('other_feature_raw_site', gatc_con, if_exists='replace', index=False)

del tmp_gene_pd
del tmp_other_pd
gc.collect()
gatc_con.close()




