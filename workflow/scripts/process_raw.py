import numpy as np
import pickle
import gzip

def count_bg(fn: str, chrom_size: dict) -> dict:
    tmp = {}
    for chrom_name, size in chrom_size.items():
        tmp[chrom_name] = np.zeros(size, dtype=np.int32)
    if fn.endswith('.gz'):
        with gzip.open(fn, 'rt') as filep:
            for line in filep:
                ent = line.strip().split('\t')
                chrom_name = ent[0]
                if chrom_name not in chrom_size:
                    continue
                start, end, count = map(int, ent[1:4])
                tmp[chrom_name][start: end] = count
    else:
        with open(fn) as filep:
            for line in filep:
                ent = line.strip().split('\t')
                chrom_name = ent[0]
                if chrom_name not in chrom_size:
                    continue
                start, end, count = map(int, ent[1:4])
                tmp[chrom_name][start: end] = count
    return tmp

def count_bed_end_gzip(fn: str, chrom_size: dict) -> dict:
    tmp = {}
    for chrom_name, size in chrom_size.items():
        tmp[chrom_name] = np.zeros(size, dtype=np.int32)
    with gzip.open(fn, 'rt') as filep:
        for line in filep:
            ent = line.strip().split('\t')
            chrom_name = ent[0]
            if chrom_name not in chrom_size:
                continue
            pos, count = map(int, ent[1:3])
            tmp[chrom_name][pos] = count
    return tmp

def bg_to_pickle(fn: str, chrom_size: dict, out_fn: str):
    tmp = count_bg(fn, chrom_size)
    with open(out_fn, 'wb') as filep:
        pickle.dump(tmp, filep)
    del tmp

def bed_end_to_pickle(fn: str, chrom_size: dict, out_fn: str):
    tmp = count_bed_end_gzip(fn, chrom_size)
    with open(out_fn, 'wb') as filep:
        pickle.dump(tmp, filep)
    del tmp


