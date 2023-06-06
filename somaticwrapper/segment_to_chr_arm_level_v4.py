"""
    Hua Sun
    Austin Southard-Smith
    2022-07-07

    gatk4cnv segment cn to band-level

    Install:
    conda install bedtools=2.30.0

    # band-level - >=50% band region must overlap to query region


    INPUT
    // gatk.called.igv.seg
    Sample  Chromosome      Start   End     Num_Probes      Call    Segment_Mean

    // band.bed
    chr star end band

    OUTPUT
    // output
    Sample   Chromosome      Start   End     band    Segment_Mean    Call


    --seg <file>         # GATK cnv *.called.igv.seg file
    --band <file>        # band region file
    --name <str>         # new sample name (default: use *.igv.seg 'Sample' name)
    --prefix <str>       # prefix for output name (default: none)
    -o,--outdir <str>    # out dir. 

    python segment_to_bandLevel.py --name new_sampleID --prefix sampleID.T --seg gatk.called.igv.seg --band band.bed -o outdir

"""

#from __future__ import bandrator_stop
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
parser.add_argument('--band', required=True, help='cytoBand.txt bed file')
parser.add_argument('--prefix', default='', help='prefix for output name')
parser.add_argument('-o', '--outdir', required=True, help='outdir')

args = parser.parse_args()
def main():
    Check_File(args.seg)
    Check_File(args.band)
    Check_Dir(args.outdir)

    SegCN_to_bandLevel(args.name, args.seg, args.band, args.prefix, args.outdir)


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


'''Make band-level cn'''
def SegCN_to_bandLevel(sampleName, f_seg, f_band, prefix, outdir):
    df_seg = pd.read_csv(f_seg, sep='\t', low_memory=False)
    #print(df_seg)
    df_seg = df_seg[['Chromosome','Start','End','Segment_Mean','Call','Sample']] #this line converts the dataframe to bed format so it can be written
    df_seg.to_csv(f'{outdir}/seg_region.tmp', sep='\t', header = False, index=False)
    #f_band should be the cytoBand.txt file.
    cmd = f'bedtools intersect -wo -a {outdir}/seg_region.tmp -b {f_band} > {outdir}/bandLevel.from_seg.tmp'
    os.system(cmd)

    file_path = f'{outdir}/bandLevel.from_seg.tmp'
    file=open(file_path, 'r')
    band_dict = {}
    band_list = []
    chr_arm_dict = {}
    chr_arm_list = []
    #this will form a dict with the format: {band: [[line1], [line2], [line3]]}
    for line in file:
        line_list = line.strip().split('\t')
        if "p" in line_list[9]:
            band_no_chr = line_list[9]
            arm = band_no_chr[0] #the first character is the arm
            chromosome = line_list[0]
            chr_arm = str(chromosome)+str(arm)
            line_list[9] = chr_arm 
        elif "q" in line_list[9]:
            band_no_chr = str(line_list[9])
            arm = band_no_chr[0] #the first character is the arm
            chromosome = str(line_list[0])
            chr_arm = chromosome+str(arm) #chr1p
            band = chromosome+band_no_chr #chr1_p36.33
            line_list[9] = chr_arm 
        else:
            continue
        if chr_arm in chr_arm_list:
            chr_arm_dict[chr_arm].append(line_list)
        else:
            chr_arm_dict[chr_arm] = [line_list]
            chr_arm_list.append(chr_arm)
    file.close()

    band_dict = chr_arm_dict
    band_list = chr_arm_list
    print(chr_arm_list)
    #here you are going to essentially replicate the logic from the SimpleCopyRatioCaller.java script that is used by the CallCopyRatioSegments.java script that is run when your run the GATK4SCNA pipeline as you will need some of the values that it bandrates below. Good luck:
    cmd2 = f'grep "Length-weighted" {outdir}/../logs/gatk4cn.s3.{sampleName}.err > {outdir}/copy-number-calling-criteria-from-step3-log-file.txt'
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
# 0	    1	    2		        3				        4	    5				        6	    7	    8	    9	    10      11		
# chr1    925692  248918613       0.029639999999999996    0       CE321E1-C1S2D1_1.T      chr1    0       2300000 p36.33  gneg    1374308
# chr1    925692  248918613       0.029639999999999996    0       CE321E1-C1S2D1_1.T      chr1    2300000 5300000 p36.32  gpos25  3000000
# chr1    925692  248918613       0.029639999999999996    0       CE321E1-C1S2D1_1.T      chr1    5300000 7100000 p36.31  gneg    1800000
# chr   seg_start       seg_stop        segment_mean      call    name_from_bam           chr     band_start      band_stop       band    g_status    intersect_length
# 0		1		        2		        3		        4		5			            6	    7		        8			    9	    10		    11
# the reported segment mean is the log2(CR) (it is the log base 2 of the copy ratio) To convert it back to copy ratio space we do the following calculation:
# (2^segment_mean) <- the copy-ratio space to the copy number space.

    out_dict = dict()
    charm_out_dict = dict()
    bands_with_more_than_2_segment_intersects = []
    file_path_to_track_bands_with_both_amplifications_and_deletions = f'{outdir}/{prefix}.amplified_and_deleted_bands.from_seg.cn'
    file = open(file_path_to_track_bands_with_both_amplifications_and_deletions, 'w')
    file.write("Length-weighted mean for z-score calling (CR space): "+str(length_weigthed_mean_for_zscore_calling)+'\n')
    file.write("Length-weighted standard deviation for z-score calling (CR space): "+str(length_weigthed_standard_deviation_for_zscore_calling)+'\n')
    file.write("Chromsome\tSegment_Start\tSegment_Stop\tSegment_Mean\tCall\tName_from_bam\tChromosome\tband_Start\tband_Stop\tband_Symbol\tg_status\tDistance_Overlap\n")
    for band in band_list:
        segment_intersects_of_band = band_dict[band] #this is a list of lists [[chr     seg_start       seg_stop        segment_mean      call    name_from_bam   chr     band_start      band_stop       band    g_status intersect_length]]
        number_of_segments_that_band_intersects = len(segment_intersects_of_band)
        chromosome = segment_intersects_of_band[0][0]
        band_start = segment_intersects_of_band[0][7]
        band_stop = segment_intersects_of_band[0][8]
        band_name = segment_intersects_of_band[0][9]
        print("number_of_segments_that_band_intersects: ", number_of_segments_that_band_intersects)
        if number_of_segments_that_band_intersects==1:
            segment_mean = segment_intersects_of_band[0][3]
            segment_call = segment_intersects_of_band[0][4]
            out_dict[band] = [sampleName, chromosome, band_start, band_stop, band_name, segment_mean, segment_call]
        else:
            copy_ratios = [] #these are log2 space 
            length_of_band_overlaps = []
            segment_intersect_calls = []
            # max_intersect_length = int(segment_intersects_of_band[0][10])
            # max_intersect_length_call = int(segment_intersects_of_band[0][4])
            # max_intersect_length_index = 0
            total_length_amplified = 0
            total_length_deleted = 0
            total_length_neutral = 0
            for i in range(0,number_of_segments_that_band_intersects):
                copy_ratios.append(2**(float(segment_intersects_of_band[i][3]))) #this converts it to copy ratio space and adds it to a list (the weighting is the same no matter if the calculation is done in the copy ratio space or copy number space. However you must not take the weighted segment means. That does not work.)
                current_length=int(segment_intersects_of_band[i][-1]) #the the length of the current segment intersection with the band of interest
                length_of_band_overlaps.append(current_length) # this adds the current intersection length to a list of all lengths
                current_call = segment_intersects_of_band[i][4] #
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
                    weighted_copy_ratio += copy_ratio*(length_of_band_overlaps[i]/total_length) #this build the weighted copy number for if there are amplifications or deletions (but not both <- exclusive or) and neutral segments
                # The below block can be uncommented to do the following:
                # For bands that do overlap a amplified seg and a deleted segment then I find if the amplification or deletion has 
                # a larger overlap for that band and the one that does is used for the calculation along with any neutral segments 
                # that are present. This is done to reduce the chance that any bands that are amplified + deleted as being called as neutral.
                # if total_length_amplified >= total_length_deleted:
                #     adjusted_total_length = total_length_amplified + total_length_neutral
                #     if total_length_neutral != 0:
                #         for index in amplification_indices:
                #             weighted_copy_ratio += copy_ratios[index]*(length_of_band_overlaps[index]/adjusted_total_length) #builds the weighted copy number for if there is an intersecting amplification, deletion, and neutral segments and the amplification is the largest (using only the amplified and neutral segments)
                #         for index in neutral_indices:
                #             weighted_copy_ratio += copy_ratios[index]*(length_of_band_overlaps[index]/adjusted_total_length) #this adds the neutral contribution to the copy number calculated on line 180

                #     else:
                #         for index in amplification_indices:
                #             weighted_copy_ratio += copy_ratios[index]*(length_of_band_overlaps[index]/adjusted_total_length) #this builds the weighted copy number for if there are only segments that overlap amplifications or deletions (no neutral) and the amplification is the largest (using only segments labeled as amplified)
                # else: #total_length_amplified < total_length_deleted:
                #     adjusted_total_length = total_length_deleted + total_length_neutral
                #     if total_length_neutral != 0:
                #         for index in deletion_indices:
                #             weighted_copy_ratio += copy_ratios[index]*(length_of_band_overlaps[index]/adjusted_total_length) #builds the wighted copy number for if there is an intersecting amplification, deletion, and neutral segments and the deletion is the largest (using only the deleted and neutral segments)
                #         for index in neutral_indices:
                #             weighted_copy_ratio += copy_ratios[index]*(length_of_band_overlaps[index]/adjusted_total_length) #this adds the neutral contribution to the copy number calculated on line 191
                #     else:
                #         for index in deletion_indices:
                #             weighted_copy_ratio += copy_ratios[index]*(length_of_band_overlaps[index]/adjusted_total_length) #this builds the weighted copy number for if there are only segments that overlap amplifications or deletions (no neutral) and the deletion is the largest (using only segments labeled as deleted)
                for line in band_dict[band]:
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
                    # print(i, copy_ratio)
                    weighted_copy_ratio += copy_ratio*(length_of_band_overlaps[i]/total_length) #this build the weighted copy number for if there are amplifications or deletions (but not both <- exclusive or) and neutral segments
            print("cutoffs and copy ratios calculated for bands that overlap multiple segments")
            print(chromosome,band_name,"neutral_call_lower_bound:",neutral_call_lower_bound,"; weighted_copy_ratio:",weighted_copy_ratio, "; neutral_call_upper_bound:", neutral_call_upper_bound)
            #new_weighted_copy_ratio = math.log2(weighted_copy_ratio) #-length_weigthed_mean_for_zscore_calling
            if weighted_copy_ratio > neutral_call_upper_bound:
                new_call = "+"
            elif weighted_copy_ratio < neutral_call_lower_bound:
                new_call = "-"
            else:
                new_call = "0"
            new_weighted_band_mean = str(math.log2(weighted_copy_ratio))
            out_dict[band] = [sampleName, chromosome, band_start, band_stop, band_name, new_weighted_band_mean, new_call]
    file.close()
    # finial output
    outtmp_path = f'{outdir}/bandLevel.from_seg.trim.tmp'
    if prefix != '':
        outtmp_path = f'{outdir}/{prefix}.bandLevel.from_seg.trim.tmp'
    outtmp = open(outtmp_path, 'w')
    for band in band_list:
        sampleName_out, chromosome_out, band_start_out, band_stop_out, band_name_out, new_weight_segment_mean, new_call_out = out_dict[band]
        # print(sampleName_out, chromosome_out, band_start_out, band_stop_out, band_name_out, new_weight_segment_mean, new_call_out)
        outtmp.writelines([sampleName_out+'\t', chromosome_out+'\t', band_start_out+'\t', band_stop_out+'\t', band_name_out+'\t', new_weight_segment_mean+'\t', new_call_out+'\n'])
    outtmp.close()

    # loading in the temp file and using pandas to remove any possible duplicate values (there shouldn't be)
    df_temp = pd.read_csv(f'{outdir}/{prefix}.bandLevel.from_seg.trim.tmp', sep='\t', header=None, low_memory=False)
    df_temp.columns = ['Sample', 'Chromosome','Start','End','band', 'Segment_Mean', 'Call']
    print(df_temp)
    df_temp = df_temp[['Sample', 'Chromosome','Start','End','band', 'Segment_Mean', 'Call']].drop_duplicates(subset=['band'], keep='first')
    if sampleName != '':
        df_temp['Sample'] = sampleName
    # finial output
    outfile = f'{outdir}/bandLevel.from_seg.cn'
    if prefix != '':
        outfile = f'{outdir}/{prefix}.bandLevel.from_seg.cn'
    df_temp.to_csv(outfile, sep='\t', index=False)

    # df_temp = pd.read_csv(f'{outdir}/bandLevel.from_seg.tmp', sep='\t', header=None, usecols=[*range(3, 10)], low_memory=False)
    # df_temp.columns = ['Segment_Mean','Call','Sample','Chromosome','Start','End','band']
    # df_temp = df_temp[['Sample', 'Chromosome','Start','End','band', 'Segment_Mean', 'Call']].drop_duplicates(subset=['band'], keep='first')

    # if sampleName != '':
    #     df_temp['Sample'] = sampleName

    # # finial output
    # outfile = f'{outdir}/bandLevel.from_seg.cn'
    # if prefix != '':
    #     outfile = f'{outdir}/{prefix}.bandLevel.from_seg.cn'

    # df_temp.to_csv(outfile, sep='\t', index=False)

    # os.remove(f'{outdir}/seg_region.tmp')
    # os.remove(f'{outdir}/bandLevel.from_seg.tmp')




if __name__ == '__main__':
    main()

