import argparse
import os
import logging
import subprocess
import shutil
from pathlib import Path

import pandas as pd

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

parser = argparse.ArgumentParser()

parser.add_argument('sample', type=str,
    help='Sample id')

parser.add_argument('tumor_bam', type=str,
    help='tumor wxs bam')

parser.add_argument('normal_bam', type=str,
    help='normal WXS bam')

parser.add_argument('--run-dir', type=str, default='output',
    help='Directory that somaticwrapper is run in. Default is output/')

parser.add_argument('--log-dir', type=str, default='logs',
    help='Directory that somaticwrapper is run in. Default is logs/')

parser.add_argument('--reference', type=str,
    help='reference location')

parser.add_argument('--rescue-genes', type=str, required=True,
    help='Which genes to rescue from vaf filter. Either a .bed file with 4 columns: <chrom> <start> <stop> <gene>. Or a .txt file with one gene per line.')

parser.add_argument('--script-dir', type=str, default='/pecgs-somaticwrapper/somaticwrapper',
    help='Directory with all the somaticwrapper perl scripts')

args = parser.parse_args()


SMG_FP = 'smg.txt'


def call_somaticwrapper(run_dir, smg_fp, log_dir, step):
    pieces = [
        'perl somaticwrapper.pl',
        '--rdir', run_dir,
        '--log', log_dir,
        '--ref', args.reference,
        '--smg', smg_fp,
        '--step', str(step)
    ]
    return ' '.join(pieces)


def setup_run():
    # make run dir
    Path(args.run_dir).mkdir(parents=True, exist_ok=True)
    sample_dir = os.path.join(args.run_dir, args.sample)
    Path(sample_dir).mkdir(parents=True, exist_ok=True)
    logging.info('linking bams')
    os.symlink(args.tumor_bam, os.path.join(sample_dir, f'{args.sample}.T.bam'))
    os.symlink(f'{args.tumor_bam}.bai', os.path.join(sample_dir, f'{args.sample}.T.bam.bai'))
    os.symlink(args.normal_bam, os.path.join(sample_dir, f'{args.sample}.N.bam'))
    os.symlink(f'{args.normal_bam}.bai', os.path.join(sample_dir, f'{args.sample}.N.bam.bai'))
    # shutil.copy(args.tumor_bam, os.path.join(sample_dir, f'{args.sample}.T.bam'))
    # shutil.copy(f'{args.tumor_bam}.bai', os.path.join(sample_dir, f'{args.sample}.T.bam.bai'))
    # shutil.copy(args.normal_bam, os.path.join(sample_dir, f'{args.sample}.N.bam'))
    # shutil.copy(f'{args.normal_bam}.bai', os.path.join(sample_dir, f'{args.sample}.N.bam.bai'))
    logging.info('path exists:')
    logging.info(os.path.exists('/storage1/fs1/dinglab/Active/Projects/estorrs/pecgs_resources/somaticwrapper/softwae/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py'))


    # make smg list
    df = pd.read_csv(args.rescue_genes, sep='\t', header=None)
    if len(df.columns) == 1:
        shutil.copy(args.rescue_genes, SMG_FP)
    elif len(df.columns) == 4:
        df = df.iloc[:, -1:]
        df.to_csv(SMG_FP, sep='\t', header=None, index=None)


def run_somaticwrapper():
    logging.info('setting up run directory')
    setup_run()

    run_dir = os.path.abspath(args.run_dir)
    smg_fp = os.path.abspath(SMG_FP)
    log_dir = os.path.abspath(args.log_dir)
    Path(log_dir).mkdir(parents=True, exist_ok=True)

    logging.info(f'run dir: {run_dir}')
    logging.info(f'smg filepath: {smg_fp}')
    logging.info(f'log dir: {log_dir}')

    logging.info(f'run dir contents: {os.listdir(run_dir)}')

    os.chdir(args.script_dir)

    for i in range(1, 15, 1):
        logging.info(f'running step {i}')
        cmd = call_somaticwrapper(run_dir, smg_fp, log_dir, i)
        logging.info(f'executing command: {cmd}')
        subprocess.check_output(cmd, shell=True)

    logging.info('somaticwrapper run finished')


def main():
    run_somaticwrapper()


if __name__ == '__main__':
    main()