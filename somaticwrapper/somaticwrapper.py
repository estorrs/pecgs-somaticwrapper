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

parser.add_argument('--reference', type=str,
    help='reference location')

parser.add_argument('--rescue-genes', type=str, required=True,
    help='Which genes to rescue from vaf filter. Either a .bed file with 4 columns: <chrom> <start> <stop> <gene>. Or a .txt file with one gene per line.')

parser.add_argument('--script-dir', type=str, default='/pecgs-somaticwrapper/somaticwrapper',
    help='Directory with all the somaticwrapper perl scripts')

args = parser.parse_args()


SMG_FP = 'smg.txt'


def call_somaticwrapper(run_dir, smg_fp, log_fp, step):
    pieces = [
        'perl somaticwrapper.pl',
        '-rdir', run_dir,
        '--log', log_fp,
        '--ref', args.reference,
        '--smg', smg_fp,
        '--step', str(step)
    ]
    return ' '.join(pieces)


def setup_run():
    # make run dir
    Path(args.run_dir).mkdir(parents=True, exist_ok=True)
    os.symlink(args.tumor_bam, os.path.join(args.run_dir, f'{args.sample}.T.bam'))
    os.symlink(f'{args.tumor_bam}.bai', os.path.join(args.run_dir, f'{args.sample}.T.bam.bai'))
    os.symlink(args.normal_bam, os.path.join(args.run_dir, f'{args.sample}.N.bam'))
    os.symlink(f'{args.normal_bam}.bai', os.path.join(args.run_dir, f'{args.sample}.N.bam.bai'))

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

    logging.info(f'run dir: {run_dir}')
    logging.info(f'smg filepath: {smg_fp}')

    os.chdir(args.script_dir)

    for i in range(1, 15, 1):
        logging.info(f'running step {i}')
        log_fp = os.path.join(os.path.abspath(args.run_dir), f'step{i}.log')
        cmd = call_somaticwrapper(run_dir, smg_fp, log_fp, i)
        logging.info(f'executing command: {cmd}')
        subprocess.check_output(cmd, shell=True)

    logging.info('somaticwrapper run finished')


def main():
    run_somaticwrapper()


if __name__ == '__main__':
    main()