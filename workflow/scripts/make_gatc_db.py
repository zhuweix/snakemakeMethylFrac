import logging 
import sys
import sqlite3 
import re
import pandas as pd
from Bio import SeqIO


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
streamHandler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter("[%(asctime)s] %(levelname)s:%(name)s:#%(lineno)d:%(message)s")
streamHandler.setFormatter(formatter)
logger.addHandler(streamHandler)
fileHandler = logging.FileHandler(snakemake.log[0])
fileHandler.setFormatter(formatter)
logger.addHandler(fileHandler)

ref_size_fn = snakemake.input['ref_size_sheet']
ref_fn = snakemake.input['ref_file']
db_fn = snakemake.output[0]


db_con = sqlite3.connect(db_fn)

chrom_sizes = {}
chrom_names = {}


ref_size_pd = pd.read_csv(ref_size_fn, sep='\t', header=None, 
                          names=['Chrom', 'RefName', 'Size'])

ref_size_pd.to_sql('chrom_size_table', db_con, if_exists='replace', index=False)

for _, row in ref_size_pd.iterrows():
    chrom_sizes[row['Chrom']] = row['Size']
    chrom_names[row['RefName']] = row['Chrom']

ref_seq = SeqIO.parse(ref_fn, 'fasta')

seq_pd = []
for record in ref_seq:
    if not record.id in chrom_names:
        continue 
    seq = str(record.seq).upper()
    chrom = chrom_names[record.id]
    chrom_size = chrom_sizes[chrom]
    for match in re.finditer(r'GATC', seq):
        start = match.start()
        is_overlap = 0
        if start > 0 and seq[start - 1] == 'C':
            is_overlap += 1
        if start + 4 < chrom_size and seq[start + 4] == 'G':
            is_overlap += 1
        seq_pd.append((chrom, start + 1, is_overlap))
    logger.info(f'Finished processing {chrom}')
seq_pd = pd.DataFrame(seq_pd, columns=['chrom', 'Apos', 'is_overlap'])
seq_pd.to_sql('gatc_pos', db_con, if_exists='replace', index=False)

db_con.close()