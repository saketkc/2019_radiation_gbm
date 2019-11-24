GTF = '/home/cmb-panasas2/skchoudh/genomes/hg38/annotation/Homo_sapiens.GRCh38.96.gtf'
STAR_INDEX = '/home/cmb-panasas2/skchoudh/genomes/hg38/star_annotated_ribopod_v96'
HTSEQ_STRANDED='no'
RAWDATA_DIR ='/home/cmb-06/as/skchoudh/dna/Choudhary_Burns_et_al_2019_radiation_GEO_submission/RNA_seq'
OUT_DIR = '/staging/as/skchoudh/rna/Feb_02_2017_Radiation_GBM_merged_fastq'

from collections import defaultdict
from riboraptor.helpers import mkdir_p
from itertools import chain, combinations
from os.path import join
import glob
import os
import re

ADAPTER = 'default'

mkdir_p(join(OUT_DIR, 'slurm-logs'))
workdir: OUT_DIR


SAMPLES = glob.glob('{}**/*.fastq.gz'.format(RAWDATA_DIR), recursive=False)
SAMPLE_LANE_R1 = []
SAMPLE_LANE_R2 = []
SAMPLE_NAMES = []
for sample in SAMPLES:
    sample = sample.replace('{}/'.format(RAWDATA_DIR),'')
    if 'R1_001' in sample:
      sample_name = re.split(r'_L\d\d\d_', sample)[0]
      lane_name = re.search(r'L\d\d\d', sample).group()
      SAMPLE_LANE_R1.append((sample_name, lane_name))
    elif 'R2_001' in sample:
      sample_name = re.split(r'_L\d\d\d_', sample)[0]
      lane_name = re.search(r'L\d\d\d', sample).group()
      SAMPLE_LANE_R2.append((sample_name, lane_name))
    else:
      print('Neother R1 nor R2')

    SAMPLE_NAMES.append(sample_name.replace('_R1_001','').replace('_R2_001', ''))

SAMPLE_LANE_R1 = sorted(SAMPLE_LANE_R1, key=lambda tup: tup[0])
SAMPLE_LANE_R2 = sorted(SAMPLE_LANE_R2, key=lambda tup: tup[0])

SAMPLE_NAMES_R1, LANE_NAMES_R1 = zip(*SAMPLE_LANE_R1)
SAMPLE_NAMES_R2, LANE_NAMES_R2 = zip(*SAMPLE_LANE_R2)

SAMPLEWISE_LANES_R1 = defaultdict(list)
for sample_name, lane in SAMPLE_LANE_R1:
    SAMPLEWISE_LANES_R1[sample_name].append(lane)

SAMPLEWISE_LANES_R2 = defaultdict(list)
for sample_name, lane in SAMPLE_LANE_R2:
    SAMPLEWISE_LANES_R2[sample_name].append(lane)


def merge_fastq_input_R1(wildcards):
    return [os.path.join(RAWDATA_DIR, '') + '{}_{}_R1_001'.format(wildcards.sample_name, lane)+ '.fastq.gz' for lane in SAMPLEWISE_LANES_R1[wildcards.sample_name] ]

def merge_fastq_input_R2(wildcards):
    return [os.path.join(RAWDATA_DIR, '') + '{}_{}_R2_001'.format(wildcards.sample_name, lane)+ '.fastq.gz' for lane in SAMPLEWISE_LANES_R2[wildcards.sample_name] ]

print(SAMPLEWISE_LANES_R1)
print(SAMPLEWISE_LANES_R2)

rule all:
  input:
        expand('merged_fastq/{sample_name}_R1.fastq.gz', sample_name=SAMPLE_NAMES_R1),
        expand('merged_fastq/{sample_name}_R2.fastq.gz', sample_name=SAMPLE_NAMES_R2),
        expand('md5/{sample_name}_R1.txt', sample_name=SAMPLE_NAMES_R1),
        expand('md5/{sample_name}_R2.txt', sample_name=SAMPLE_NAMES_R2),
        expand('bams_sortedbyname/{sample_name}.sortedByName.bam', sample_name=SAMPLE_NAMES),
        expand('HTSeq-counts/byCDS/{sample_name}.CDS.counts.tsv', sample_name=SAMPLE_NAMES),
        expand('HTSeq-counts/byExon/{sample_name}.exon.counts.tsv', sample_name=SAMPLE_NAMES),
        expand('size_metrics/{sample_name}.pdf', sample_name=SAMPLE_NAMES),

rule merge_fastqs_R1:
    input: merge_fastq_input_R1
    output: 'merged_fastq/{sample_name}_R1.fastq.gz'
    shell: 
        r'''cat {input} > {output}'''

rule merge_fastqs_R2:
    input: merge_fastq_input_R2
    output: 'merged_fastq/{sample_name}_R2.fastq.gz'
    shell: 
        r'''cat {input} > {output}'''

rule md5_1:
    input: 'merged_fastq/{sample_name}_R1.fastq.gz'
    output: 'md5/{sample_name}_R1.txt',
    shell: 
        r'''md5sum {input} > {output}'''


rule md5_2:
    input: 'merged_fastq/{sample_name}_R2.fastq.gz'
    output: 'md5/{sample_name}_R2.txt',
    shell: 
        r'''md5sum {input} > {output}'''


rule perfom_trimming:
    input:
        R1='merged_fastq/{sample_name}_R1.fastq.gz',
        R2='merged_fastq/{sample_name}_R2.fastq.gz'
    params:
        out_dir='preprocessed/{sample_name}',
        phred_cutoff=5
    output:
        'preprocessed/{sample_name}/{sample_name}_R1_val_1.fq.gz',
        'preprocessed/{sample_name}/{sample_name}_R2_val_2.fq.gz',
    shell:
        r'''
            trim_galore --paired -o {params.out_dir} -q {params.phred_cutoff} {input.R1} {input.R2}
        '''

rule map_starlong:
    input:
        R1='preprocessed/{sample_name}/{sample_name}_R1_val_1.fq.gz',
        R2='preprocessed/{sample_name}/{sample_name}_R2_val_2.fq.gz',
    output:
        bam = 'bams/{sample_name}.bam',
    params:
        prefix = 'bams/{sample_name}',
        starlogs = 'mapped/starlogs'
    threads: 16
    shell:
        r'''
        STARlong --runThreadN {threads}\
             --genomeDir {STAR_INDEX}\
             --outFileNamePrefix {params.prefix} --readFilesIn {input.R1} {input.R2}\
             --outSAMtype BAM SortedByCoordinate\
             --readFilesCommand zcat\
            --outFilterMultimapScoreRange 20\
            --outFilterScoreMinOverLread 0\
            --outFilterMatchNminOverLread 0.66\
            --outFilterMismatchNmax 100\
            --winAnchorMultimapNmax 200\
            --seedPerReadNmax 100000 --seedPerWindowNmax 100\
            && mv {params.prefix}Aligned.sortedByCoord.out.bam {output.bam} &&\
            mkdir -p {params.starlogs} && mv {params.prefix}Log.final.out\
            {params.prefix}Log.out {params.prefix}Log.progress.out {params.starlogs}\
            && samtools index {output.bam}
        '''

rule sort_by_name:
    input: 'bams/{sample_name}.bam'
    output: 'bams_sortedbyname/{sample_name}.sortedByName.bam'
    shell:
        r'''
            samtools sort -on {input} -T /tmp/ -o {output}
        '''


rule count_byCDS:
    input: 'bams_sortedbyname/{sample_name}.sortedByName.bam'
    params:
        annotation=GTF,
        phred_cutoff=5
    output: 'HTSeq-counts/byCDS/{sample_name}.CDS.counts.tsv'
    shell:
        r'''
        htseq-count --order=name --format=bam --mode=intersection-strict --stranded={HTSEQ_STRANDED} --minaqual={params.phred_cutoff} --type=CDS --idattr=gene_id {input} {params.annotation} > {output}.temp && cat {output}.temp | sed -E 's/\.[0-9]+//' > {output} && rm {output}.temp
        '''

rule count_byExon:
    input: 'bams_sortedbyname/{sample_name}.sortedByName.bam'
    params:
        annotation=GTF,
        phred_cutoff=5
    output: 'HTSeq-counts/byExon/{sample_name}.exon.counts.tsv'
    shell:
        r'''
        htseq-count --order=name --format=bam --mode=intersection-strict --stranded={HTSEQ_STRANDED} --minaqual={params.phred_cutoff} --type=exon --idattr=gene_id {input} {params.annotation} > {output}.temp && cat {output}.temp | sed -E 's/\.[0-9]+//' > {output} && rm {output}.temp
        '''

rule picard_sizemetrics:
    input: 'bams/{sample_name}.bam'
    output: 
      size_metrics = 'size_metrics/{sample_name}.size_metrics.txt',
      hist = 'size_metrics/{sample_name}.pdf'
    shell:
      r'''
      picard CollectInsertSizeMetrics I={input} O={output.size_metrics} H={output.hist}
      '''
