import os
import logging 
import sys
import pandas as pd
import gc
import process_raw

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
streamHandler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter("[%(asctime)s] %(levelname)s:%(name)s:#%(lineno)d:%(message)s")
streamHandler.setFormatter(formatter)
logger.addHandler(streamHandler)
fileHandler = logging.FileHandler(snakemake.log[0])
fileHandler.setFormatter(formatter)
logger.addHandler(fileHandler)


dirn = snakemake.params['bamdir']
if not os.path.exists(dirn):
    raise ValueError(f'{dirn} does not exist')
chrom_fn = snakemake.params['chrom_size']
if not os.path.exists(chrom_fn):
    raise ValueError(f'{chrom_fn} does not exist')

chrom_size_pd = pd.read_csv(chrom_fn, sep='\t', header=None, 
                            names=['Chrom', 'Ref', 'Size'])
chrom_sizes = {}
for _, row in chrom_size_pd.iterrows():
    chrom_sizes[row['Ref']] = row['Size']

for fn in sorted(os.listdir(dirn)):
    if fn.endswith('.count.forward.1kb.gz'):
        prefix = fn.split('.')[0]        
        process_raw.bed_end_to_pickle(fn=os.path.join(dirn, fn), 
                                      chrom_size=chrom_sizes,
                                      out_fn=os.path.join(dirn, "{}.gatc.right.pickle".format(prefix)))
        logger.info(f'Finished processing {prefix} Forward')
        os.remove(os.path.join(dirn, fn))        
        gc.collect()
    elif fn.endswith('.count.reverse.1kb.gz'):
        process_raw.bed_end_to_pickle(fn=os.path.join(dirn, fn), 
                                      chrom_size=chrom_sizes,
                                      out_fn=os.path.join(dirn, "{}.gatc.left.pickle".format(prefix)))
        logger.info(f'Finished processing {prefix} Reverse')  
        os.remove(os.path.join(dirn, fn))        
        gc.collect()        
    else:
        continue

for fn in sorted(os.listdir(dirn)):
    if fn.endswith('.1kb.bedgraph.gz'):
        prefix = fn.split('.')[0]        
        process_raw.bg_to_pickle(fn=os.path.join(dirn, fn), 
                                       chrom_size=chrom_sizes,
                                       out_fn=os.path.join(dirn, "{}.occ.pickle".format(prefix)))
        os.remove(os.path.join(dirn, fn))
        logger.info(f'Finished processing {prefix} Occupancy')    
        gc.collect()        
    else:
        continue