'''
    Hua Sun
    2021-02-24 Added segment merging
    2021-01-30

    1.Merage all of segment-level to one file
    2.Merage all of gene-level to one file
    

    // *.called.igv.seg
    Sample  Chromosome      Start   End     Num_Probes      Call    Segment_Mean
    
    // *.geneLevel.from_seg.cn
    Sample   Chromosome      Start   End     Gene    Segment_Mean    Call

    INPUT
    /path/dir

    OUTPUT
    `merged.geneLevel.from_seg.hg38.log2ratio.tsv`
    e.g.
    Gene Sample1 Sample2 ...

    python3 mergeMultipleFilesToOne.py ./gatk4scna

'''

import os
import sys
import pandas as pd


dir = sys.argv[1]


#----- 1.Merge segment
cmd = f'cat {dir}/*/*.T.called.igv.seg | grep -v Segment_Mean > {dir}/merged.segment_level.cn.tmp'
os.system(cmd)

df = pd.read_csv(f'{dir}/merged.segment_level.cn.tmp', sep='\t', header=None)
df.columns = ['Sample','Chromosome','Start','End','Num_Probes','Call','Segment_Mean']
df = df[['Sample','Chromosome','Start','End','Num_Probes','Segment_Mean','Call']]

df.to_csv(f'{dir}/merged.segment_level.hg38.log2ratio.tsv', sep='\t', index=False)
os.remove(f'{dir}/merged.segment_level.cn.tmp')



#----- 2.Merge gene level
cmd = f'cat {dir}/*/*.geneLevel.from_seg.cn | grep -v Segment_Mean > {dir}/merged.gene_level.cn.tmp'
os.system(cmd)

df = pd.read_csv(f'{dir}/merged.gene_level.cn.tmp', sep='\t', header=None)
df.columns = ['Sample','Chromosome','Start','End','Gene','Segment_Mean','Call']

df_spread = df.pivot_table(values='Segment_Mean', index=['Gene'], columns='Sample', aggfunc='first')
df_spread = df_spread.reset_index()

df_spread.to_csv(f'{dir}/merged.gene_level.from_seg.hg38.log2ratio.tsv', sep='\t', index=False)

df = pd.read_csv(f'{dir}/merged.gene_level.cn.tmp', sep='\t', header=None)
df.columns = ['Sample','Chromosome','Start','End','Gene','Segment_Mean','Call']

df_spread = df.pivot_table(values='Call', index=['Gene'], columns='Sample', aggfunc='first')
df_spread = df_spread.reset_index()

df_spread.to_csv(f'{dir}/merged.gene_level.from_seg.hg38.calls.tsv', sep='\t', index=False)
os.remove(f'{dir}/merged.gene_level.cn.tmp')

