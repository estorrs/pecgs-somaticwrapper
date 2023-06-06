'''
from https://github.com/ding-lab/GATK4SCNA/blob/main/src/segment_to_geneLevel.py

    Hua Sun
    2021-01-29
    gatk4cnv segment cn to gene-level
    Install:
    conda install bedtools=2.30.0
    # gene-level - >=50% gene region must overlap to query region
    INPUT
    // gatk.called.igv.seg
    Sample  Chromosome      Start   End     Num_Probes      Call    Segment_Mean
    // gene.bed
    chr star end gene
    OUTPUT
    // output
    Sample   Chromosome      Start   End     Gene    Segment_Mean    Call
    --seg <file>         # GATK cnv *.called.igv.seg file
    --gene <file>        # gene region file
    --name <str>         # new sample name (default: use *.igv.seg 'Sample' name)
    --prefix <str>       # prefix for output name (default: none)
    -o,--outdir <str>    # out dir. 
    python segment_to_geneLevel.py --name new_sampleID --prefix sampleID.T --seg gatk.called.igv.seg --gene gene.bed -o outdir
'''

import os
import sys
import pandas as pd
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--seg', required=True, help='cnv file')
parser.add_argument('--name', default='', help='new sample id')
parser.add_argument('--gene', required=True, help='gene bed file')
parser.add_argument('--prefix', default='', help='prefix for output name')
parser.add_argument('-o', '--outdir', required=True, help='outdir')

args = parser.parse_args()



def main():
    Check_File(args.seg)
    Check_File(args.gene)
    Check_Dir(args.outdir)

    SegCN_to_GeneLevel(args.name, args.seg, args.gene, args.prefix, args.outdir)


'''
    Set Functions
'''

def Check_File(file):
    if not os.path.exists(file):
        message = '[Error] File not exists ...' + str(file)
        sys.exit(message)


def Check_Dir(dir):
    if not os.path.isdir(dir):
        os.mkdir(dir)


'''Make gene-level cn'''
def SegCN_to_GeneLevel(sampleName, f_seg, f_gene, prefix, outdir):
    df_seg = pd.read_csv(f_seg, sep='\t', low_memory=False)
    df_seg = df_seg[['Chromosome','Start','End','Segment_Mean','Call','Sample']]
    df_seg.to_csv(f'{outdir}/seg_region.tmp', sep='\t', header = False, index=False)
    
    cmd = f'bedtools intersect -wa -wb -F 0.50 -a {outdir}/seg_region.tmp -b {f_gene} > {outdir}/geneLevel.from_seg.tmp'
    os.system(cmd)

    df_temp = pd.read_csv(f'{outdir}/geneLevel.from_seg.tmp', sep='\t', header=None, usecols=[*range(3, 10)], low_memory=False)
    df_temp.columns = ['Segment_Mean','Call','Sample','Chromosome','Start','End','Gene']
    df_temp = df_temp[['Sample', 'Chromosome','Start','End','Gene', 'Segment_Mean', 'Call']].drop_duplicates(['Gene'])

    if sampleName != '':
        df_temp['Sample'] = sampleName
    
    # finial output
    outfile = f'{outdir}/geneLevel.from_seg.cn'
    if prefix != '':
        outfile = f'{outdir}/{prefix}.geneLevel.from_seg.cn'
    
    df_temp.to_csv(outfile, sep='\t', index=False)

    os.remove(f'{outdir}/seg_region.tmp')
    os.remove(f'{outdir}/geneLevel.from_seg.tmp')
    



if __name__ == '__main__':
    main()
