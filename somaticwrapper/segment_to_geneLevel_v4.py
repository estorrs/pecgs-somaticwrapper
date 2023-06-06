'''
    Hua Sun
    Austin Southard-Smith
    2022-07-07

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

#from __future__ import generator_stop
#from importlib.metadata import entry_points
import os
#from subprocess import call
import sys
import pandas as pd
import math
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--seg', required=True, help='cnv file')
parser.add_argument('--name', default='', help='new sample id')
parser.add_argument('--gene', required=True, help='gene bed file')
parser.add_argument('--prefix', default='', help='prefix for output name')
parser.add_argument('-o', '--outdir', required=True, help='outdir')

# added in addition to austins pipeline so it works with cromwell
parser.add_argument('--step3-stderr-filepath', default='./stderr', help='Where stderr stored from step 3 is located. Default is ./stderr which is where it will be during a cromwell execution.')

args = parser.parse_args()



def main():
    Check_File(args.seg)
    Check_File(args.gene)
    Check_Dir(args.outdir)

    SegCN_to_GeneLevel(args.name, args.seg, args.gene, args.prefix, args.outdir, args.step3_stderr_filepath)


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
def SegCN_to_GeneLevel(sampleName, f_seg, f_gene, prefix, outdir, step3_stderr_filepath):
    df_seg = pd.read_csv(f_seg, sep='\t', low_memory=False)
    #print(df_seg)
    df_seg = df_seg[['Chromosome','Start','End','Segment_Mean','Call','Sample']] #this line converts the dataframe to bed format so it can be written
    df_seg.to_csv(f'{outdir}/seg_region.tmp', sep='\t', header = False, index=False)
    
    cmd = f'bedtools intersect -wo -a {outdir}/seg_region.tmp -b {f_gene} > {outdir}/geneLevel.from_seg.tmp'
    os.system(cmd)

    file_path = f'{outdir}/geneLevel.from_seg.tmp'
    file=open(file_path, 'r')
    gene_dict = {}
    gene_list = []
    #this will form a dict with the format: {gene: [[line1], [line2], [line3]]}
    for line in file:
        line_list = line.strip().split('\t')
        gene=line_list[9]
        if (gene in gene_list):
            gene_dict[gene].append(line_list)
        else:
            gene_dict[gene] = [line_list]
            gene_list.append(gene)
    file.close()

    #here you are going to essentially replicate the logic from the SimpleCopyRatioCaller.java script that is used by the CallCopyRatioSegments.java script that is run when your run the GATK4SCNA pipeline as you will need some of the values that it generates below. Good luck:
    cmd2 = f'grep "Length-weighted" {step3_stderr_filepath} > {outdir}/copy-number-calling-criteria-from-step3-log-file.txt'
    os.system(cmd2)

    cr_calling_path = f'{outdir}/copy-number-calling-criteria-from-step3-log-file.txt'
    file = open(cr_calling_path, 'r')
    file.readline()
    file.readline()
    length_weigthed_mean_for_zscore_calling = float(file.readline().strip().split(' ')[-1])
    length_weigthed_standard_deviation_for_zscore_calling = float(file.readline().strip().split(' ')[-1])
    z_score_neutral_call_upper_bound = length_weigthed_mean_for_zscore_calling + (2*length_weigthed_standard_deviation_for_zscore_calling)
    z_score_neutral_call_lower_bound = length_weigthed_mean_for_zscore_calling - (2*length_weigthed_standard_deviation_for_zscore_calling)
    # the cr neutral region is all segments with a copy ratio greater than or equal to 0.9 and less than or equal to 1.1
    # if the cr is greater than 1.1 then it is automatically called as an amplification (+)
    # if the cr is less than 0.9 then it is automatically called as a deletion (-)
    # The z-score threshold neutral call bound is only used if the segment falls outside of the bounds of 0.9 or 1.1. 
    # If the standard deviation or mean are such that the z-score threshold never falls outside of the cr neutral region bounds [0.9, 1.1],
    # then the threshold used is only the cr neutral region bounds.
    if z_score_neutral_call_lower_bound < 0.9:
        neutral_call_lower_bound = z_score_neutral_call_lower_bound
    else:
        neutral_call_lower_bound = 0.9
    if z_score_neutral_call_upper_bound > 1.1:
        neutral_call_upper_bound = z_score_neutral_call_upper_bound
    else:
        neutral_call_upper_bound = 1.1
    print("neutral_call_upper_bound: ",neutral_call_upper_bound)
    print("neutral_call_lower_bound: ",neutral_call_lower_bound)
    file.close()

    #the cr neutral region is all segments with a copy ratio greater than or equal to 0.9 and less than or equal to 1.1
    #if the cr is greater than 1.1 then it is automatically called as an amplification (+)
    #if the cr is less than 0.9 then it is automatically called as a deletion (-)
    #see https://github.com/broadinstitute/gatk/blob/2e6045a259ed2ded3e9036a5b44a1f8ba330860d/src/main/java/org/broadinstitute/hellbender/tools/copynumber/CallCopyRatioSegments.java
    #and https://github.com/broadinstitute/gatk/blob/2e6045a259ed2ded3e9036a5b44a1f8ba330860d/src/main/java/org/broadinstitute/hellbender/tools/copynumber/caller/SimpleCopyRatioCaller.java#L30
    # 0       1               2               3               4       5               6       7               8               9       10
    # chr7    54530002        54804000        2.028522        +       CPT0292720006   chr7    54542325        54571080        VSTM2A  28755
    # chr7    54530002        54804000        2.028522        +       CPT0292720006   chr7    54752250        54759974        SEC61G  7724
    # chr7    54953002        55110000        5.534218        +       CPT0292720006   chr7    55019017        55211628        EGFR    90983
    # chr7    55110002        55155000        2.035126        +       CPT0292720006   chr7    55019017        55211628        EGFR    44998
    # chr7    55155002        55237000        5.593449        +       CPT0292720006   chr7    55019017        55211628        EGFR    56626
    # chr     seg_start       seg_stop        segment_mean    call    name_from_bam   chr     gene_start      gene_stop       gene    intersect_length

    out_dict = dict()
    genes_with_more_than_2_segment_intersects = []
    file_path_to_track_genes_with_both_amplifications_and_deletions = f'{outdir}/{prefix}.amplified_and_deleted_genes.from_seg.cn'
    file = open(file_path_to_track_genes_with_both_amplifications_and_deletions, 'a')
    file.write("Length-weighted mean for z-score calling (CR space): "+str(length_weigthed_mean_for_zscore_calling)+'\n')
    file.write("Length-weighted standard deviation for z-score calling (CR space): "+str(length_weigthed_standard_deviation_for_zscore_calling)+'\n')
    file.write("Chromsome\tSegment_Start\tSegment_Stop\tSegment_Mean\tCall\tName_from_bam\tChromosome\tGene_Start\tGene_Stop\tGene_Symbol\tDistance_Overlap\n")
    for gene in gene_list:
        segment_intersects_of_gene = gene_dict[gene] #this is a list of lists [[chr     seg_start       seg_stop        copy_ratio      call    name_from_bam   chr     gene_start      gene_stop       gene    intersect_length]]
        number_of_segments_that_gene_intersects = len(segment_intersects_of_gene)
        chromosome = segment_intersects_of_gene[0][0]
        gene_start = segment_intersects_of_gene[0][7]
        gene_stop = segment_intersects_of_gene[0][8]
        gene_name = segment_intersects_of_gene[0][9]
        if number_of_segments_that_gene_intersects==1:
            segment_copy_ratio = segment_intersects_of_gene[0][3]
            segment_call = segment_intersects_of_gene[0][4]
            out_dict[gene] = [sampleName, chromosome, gene_start, gene_stop, gene_name, segment_copy_ratio, segment_call]
        else:
            copy_ratios = [] #these are log2 space 
            length_of_gene_overlaps = []
            segment_intersect_calls = []
            # max_intersect_length = int(segment_intersects_of_gene[0][10])
            # max_intersect_length_call = int(segment_intersects_of_gene[0][4])
            # max_intersect_length_index = 0
            total_length_amplified = 0
            total_length_deleted = 0
            total_length_neutral = 0
            for i in range(0,number_of_segments_that_gene_intersects):
                copy_ratios.append(2**(float(segment_intersects_of_gene[i][3]))) #this converts it to copy ratio space and adds it to a list (the weighting is the same no matter if the calculation is done in the copy ratio space or copy number space. However you must not take the weighted segment means. That does not work.)
                current_length=int(segment_intersects_of_gene[i][-1]) #this the length of the current segment intersection with the gene of interest
                length_of_gene_overlaps.append(current_length) #this adds the current intersection length to a list of all lengths
                current_call = segment_intersects_of_gene[i][4] #
                segment_intersect_calls.append(current_call)
                if current_call == "+":
                    total_length_amplified += current_length
                elif current_call == '-':
                    total_length_deleted += current_length
                else:
                    total_length_neutral += current_length
                # if current_length > max_intersect_length:
                #     if current_call != "0":
                #         max_intersect_length = current_length
                #         max_intersect_length_call = current_call
                #         max_intersect_length_index = i
            amplification_indices = [i for i, x in enumerate(segment_intersect_calls) if x == "+"]
            deletion_indices = [i for i, x in enumerate(segment_intersect_calls) if x == "-"]
            neutral_indices = [i for i, x in enumerate(segment_intersect_calls) if x == "0"]
            weighted_copy_ratio = 0
            if ("+" in segment_intersect_calls) and ("-" in segment_intersect_calls):
                total_length = total_length_amplified + total_length_deleted + total_length_neutral
                for i, copy_ratio in enumerate(copy_ratios):
                    weighted_copy_ratio += copy_ratio*(length_of_gene_overlaps[i]/total_length) #this build the weighted copy number for if there are amplifications or deletions (but not both <- exclusive or) and neutral segments
                # The below block can be uncommented to do the following:
                # For genes that do overlap a amplified seg and a deleted segment then I find if the amplification or deletion has 
                # a larger overlap for that gene and the one that does is used for the calculation along with any neutral segments 
                # that are present. This was done to prevent any genes that are amplified + deleted as being called as neutral.
                # if total_length_amplified >= total_length_deleted:
                #     adjusted_total_length = total_length_amplified + total_length_neutral
                #     if total_length_neutral != 0:
                #         for index in amplification_indices:
                #             weighted_copy_ratio += copy_ratios[index]*(length_of_gene_overlaps[index]/adjusted_total_length) #builds the weighted copy number for if there is an intersecting amplification, deletion, and neutral segments and the amplification is the largest (using only the amplified and neutral segments)
                #         for index in neutral_indices:
                #             weighted_copy_ratio += copy_ratios[index]*(length_of_gene_overlaps[index]/adjusted_total_length) #this adds the neutral contribution to the copy number calculated on line 180

                #     else:
                #         for index in amplification_indices:
                #             weighted_copy_ratio += copy_ratios[index]*(length_of_gene_overlaps[index]/adjusted_total_length) #this builds the weighted copy number for if there are only segments that overlap amplifications or deletions (no neutral) and the amplification is the largest (using only segments labeled as amplified)
                # else: #total_length_amplified < total_length_deleted:
                #     adjusted_total_length = total_length_deleted + total_length_neutral
                #     if total_length_neutral != 0:
                #         for index in deletion_indices:
                #             weighted_copy_ratio += copy_ratios[index]*(length_of_gene_overlaps[index]/adjusted_total_length) #builds the wighted copy number for if there is an intersecting amplification, deletion, and neutral segments and the deletion is the largest (using only the deleted and neutral segments)
                #         for index in neutral_indices:
                #             weighted_copy_ratio += copy_ratios[index]*(length_of_gene_overlaps[index]/adjusted_total_length) #this adds the neutral contribution to the copy number calculated on line 191
                #     else:
                #         for index in deletion_indices:
                #             weighted_copy_ratio += copy_ratios[index]*(length_of_gene_overlaps[index]/adjusted_total_length) #this builds the weighted copy number for if there are only segments that overlap amplifications or deletions (no neutral) and the deletion is the largest (using only segments labeled as deleted)
                for line in gene_dict[gene]:
                    for i in range(len(line)):
                        if i == (len(line)-1):
                            line[i] = str(line[i])+'\n'
                        else:
                            line[i] = str(line[i])+'\t'
                        #line[i] = str(line[i]) 

                    file.writelines(line)
            else:
                total_length = total_length_amplified + total_length_deleted + total_length_neutral
                for i, copy_ratio in enumerate(copy_ratios):
                    weighted_copy_ratio += copy_ratio*(length_of_gene_overlaps[i]/total_length) #this build the weighted copy number for if there are amplifications or deletions (but not both <- exclusive or) and neutral segments
            print("cutoffs and copy ratios calculated for genes that overlap multiple segments")
            print(chromosome,gene_name,"neutral_call_lower_bound:",neutral_call_lower_bound,"; weighted_copy_ratio:",weighted_copy_ratio, "; neutral_call_upper_bound:", neutral_call_upper_bound)
            if weighted_copy_ratio > neutral_call_upper_bound:
                new_call = "+" 
            elif weighted_copy_ratio < neutral_call_lower_bound:
                new_call = "-"
            else:
                new_call = "0"
            new_weighted_gene_mean = str(math.log2(weighted_copy_ratio))
            out_dict[gene] = [sampleName, chromosome, gene_start, gene_stop, gene_name, new_weighted_gene_mean, new_call]
    file.close()
    # TODO: here you will write the output to the file and drop duplicates
    # finial output
    outtmp_path = f'{outdir}/geneLevel.from_seg.trim.tmp'
    if prefix != '':
        outtmp_path = f'{outdir}/{prefix}.geneLevel.from_seg.trim.tmp'
    outtmp = open(outtmp_path, 'w')
    for gene in gene_list:
        sampleName_out, chromosome_out, gene_start_out, gene_stop_out, gene_name_out, new_weighted_gene_mean_out, new_call_out = out_dict[gene]
        #print(sampleName_out, chromosome_out, gene_start_out, gene_stop_out, gene_name_out, new_weighted_gene_mean_out, new_call_out)
        outtmp.writelines([sampleName_out+'\t', chromosome_out+'\t', gene_start_out+'\t', gene_stop_out+'\t', gene_name_out+'\t', new_weighted_gene_mean_out+'\t', new_call_out+'\n'])
    outtmp.close()
    
    # loading in the temp file and using pandas to remove any possible duplicate values (there shouldn't be)
    df_temp = pd.read_csv(f'{outdir}/{prefix}.geneLevel.from_seg.trim.tmp', sep='\t', header=None, low_memory=False)
    df_temp.columns = ['Sample', 'Chromosome','Start','End','Gene', 'Segment_Mean', 'Call']
    print(df_temp)
    df_temp = df_temp[['Sample', 'Chromosome','Start','End','Gene', 'Segment_Mean', 'Call']].drop_duplicates(subset=['Gene'], keep='first')
    if sampleName != '':
        df_temp['Sample'] = sampleName
    # finial output
    outfile = f'{outdir}/geneLevel.from_seg.cn'
    if prefix != '':
        outfile = f'{outdir}/{prefix}.geneLevel.from_seg.cn'
    df_temp.to_csv(outfile, sep='\t', index=False)

    # df_temp = pd.read_csv(f'{outdir}/geneLevel.from_seg.tmp', sep='\t', header=None, usecols=[*range(3, 10)], low_memory=False)
    # df_temp.columns = ['Segment_Mean','Call','Sample','Chromosome','Start','End','Gene']
    # df_temp = df_temp[['Sample', 'Chromosome','Start','End','Gene', 'Segment_Mean', 'Call']].drop_duplicates(subset=['Gene'], keep='first')

    # if sampleName != '':
    #     df_temp['Sample'] = sampleName
    
    # # finial output
    # outfile = f'{outdir}/geneLevel.from_seg.cn'
    # if prefix != '':
    #     outfile = f'{outdir}/{prefix}.geneLevel.from_seg.cn'
    
    # df_temp.to_csv(outfile, sep='\t', index=False)

    # os.remove(f'{outdir}/seg_region.tmp')
    # os.remove(f'{outdir}/geneLevel.from_seg.tmp')
    



if __name__ == '__main__':
    main()


