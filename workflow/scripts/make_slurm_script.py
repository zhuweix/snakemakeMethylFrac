import os 
import logging 
import sys

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
streamHandler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter("[%(asctime)s] %(levelname)s:%(name)s:#%(lineno)d:%(message)s")
streamHandler.setFormatter(formatter)
logger.addHandler(streamHandler)
fileHandler = logging.FileHandler(snakemake.log[0])
fileHandler.setFormatter(formatter)
logger.addHandler(fileHandler)

input_dir = snakemake.input['input_dir']

suffix1 = snakemake.params['suffix1']
suffix2 = snakemake.params['suffix2']
jobname = snakemake.params['jobname']
output_dir = snakemake.params['output_dir']
ref_index = snakemake.params['ref_index']
jobtime = snakemake.params['jobtime']
maxinsert = snakemake.params['maxinsert']
log_dir = 'logs/'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    logger.info(f'Created {output_dir}')

if not os.path.exists('results/slurm'):
    os.makedirs('results/slurm')
    logger.info(f'Created results/slurm')

fn_list = []
i = 0

for fn in sorted(os.listdir(input_dir)):
    if not fn.endswith(suffix1):
        continue
    i += 1
    header = f'''#!/bin/bash
#SBATCH --job-name={jobname}-{i}
#SBATCH --cpus-per-task=48
#SBATCH --mem=72g
#SBATCH --gres=lscratch:340
#SBATCH --output={os.path.join(log_dir, f"{jobname}_{i}.out")} 
#SBATCH --error={os.path.join(log_dir, f"{jobname}_{i}.err")} 
#SBATCH --time={jobtime}

module load samtools
module load bowtie/2
module load bedtools

set -e

input_dir={input_dir}
output_dir={output_dir}
mkdir -p $output_dir

export TMPDIR=/lscratch/$SLURM_JOB_ID

ref_genome={ref_index}

bowtie2_params="-q --no-mixed --no-discordant --no-unal --no-contain --no-overlap -X {maxinsert}"
'''

    base = fn.split(suffix1)[0]
    fn2 = base + suffix2
    fn = os.path.join(input_dir, fn)
    fn2 = os.path.join(input_dir, fn2)
    if not os.path.exists(fn2):
        raise ValueError(f'{fn2} does not exist')
    
    content = [header, 
               f'''
base="{base}"
fn="{fn}"
fn2="{fn2}"
bowtie2  --threads $SLURM_CPUS_PER_TASK -x $ref_genome \
    -1 $fn -2 $fn2 \
    $bowtie2_params  -S $TMPDIR/$base.sam --un-conc-gz $output_dir/$base.unmapped.R0%.fq.gz

samtools view -h -f 3 -F 260 $TMPDIR/$base.sam | samtools sort -@ $SLURM_CPUS_PER_TASK -T $TMPDIR   -o $output_dir/$base.sort.1kb.bam
samtools index -@ $SLURM_CPUS_PER_TASK $output_dir/$base.sort.1kb.bam

samtools flagstat -@ $SLURM_CPUS_PER_TASK $TMPDIR/$base.sam > $output_dir/$base.flagstat.txt

genomeCoverageBed  -ibam $output_dir/$base.sort.1kb.bam -pc -bg | gzip > $output_dir/$base.1kb.bedgraph.gz

samtools view -h -f 3 -F 288 $output_dir/$base.sort.1kb.bam | \
    bedtools genomecov -ibam stdin -dz -5 | \
    gzip > $output_dir/$base.count.reverse.1kb.gz 

samtools view -h -f 35 -F 256 $output_dir/$base.sort.1kb.bam | \
    bedtools genomecov -ibam stdin -dz -5 | \
    gzip > $output_dir/$base.count.forward.1kb.gz                  
               ''']
    
    outfn = f'results/slurm/{jobname}_{i}.sh'
    fn_list.append(outfn)
    with open(outfn, 'w') as filep:
        filep.write('\n'.join(content))

with open(snakemake.output[0], 'w') as filep:
    filep.write('\n'.join(fn_list))

logger.info(f'Finished writing scripts in {snakemake.output[0]}')