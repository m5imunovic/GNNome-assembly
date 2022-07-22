import argparse
import multiprocessing as mp
import os
import subprocess
from datetime import datetime
from pathlib import Path

from Bio import SeqIO
from hydra import initialize, compose
from omegaconf import OmegaConf
from typeguard import typechecked
from cmd_utils import get_read_files

import graph_dataset_dbg

#from homo_sapiens import chm13_chr_lens as chr_lens
#from homo_sapiens import get_chr_dirs, species_specific_dirs, species_reference

from assembler import assembler_factory
from bioutils.simulator import simulator_factory
from drosophila_melanogaster import release_6_plus_iso1_mt_chr_lens as chr_lens
from drosophila_melanogaster import get_chr_dirs, species_specific_dirs, species_reference


@typechecked
def validate_chrs(chrs: dict) -> bool:
    diff = set(chrs) - set(get_chr_dirs())
    if len(diff):
        raise Exception(f'Not using valid chromosome name(s): {diff}')

    return True


@typechecked
def change_description(file_path : Path, multiline : bool = True):
    new_fasta = []
    # TODO: add assertion for format
    # 'fasta' for FASTA file, 'fastq' for FASTQ file
    format = file_path.suffix.strip('.').lower()
    for record in SeqIO.parse(file_path, format=format):
        des = record.description.split(",")
        record.id = des[0][5:]
        strand = '+' if des[1] == 'forward' else '-'
        start, end = des[2][9:].split("-")
        record.description = f'strand={strand}, start={start}, end={end}'
        new_fasta.append(record)
    
    if multiline:
        SeqIO.write(new_fasta, file_path, "fasta")
    else:
        # output entire DNA sequence in a single line
        new_fasta = [SeqIO.FastaIO.as_fasta_2line(record) for record in new_fasta]
        with open(file_path, 'w') as handle:
            handle.writelines(new_fasta)


@typechecked        
def create_chr_dirs(pth : Path):
    for chr_dir in get_chr_dirs():
        chr_dir_full = pth / chr_dir
        os.makedirs(chr_dir_full, exist_ok=True)
        for subdir in ['mkdir', 'raw', 'processed', 'info', 'asm_output', 'graphia']:
            os.makedirs(chr_dir_full / subdir, exist_ok=True)


# -1. Set up the data file structure
@typechecked
def file_structure_setup(data_path: Path, ref_path: Path):
    print(f'SETUP::filesystem:: Create directories for storing data')

    os.makedirs(data_path, exist_ok=True)
    os.makedirs(ref_path, exist_ok=True)

    species_specific_dirs(ref_path)

    chromosomes_dir = ref_path / 'chromosomes'
    os.makedirs(chromosomes_dir, exist_ok=True)
    simulated_dir = data_path / 'simulated'
    os.makedirs(simulated_dir, exist_ok=True)
    create_chr_dirs(simulated_dir)

    real_dir = data_path / 'real'
    os.makedirs(real_dir, exist_ok=True)
    create_chr_dirs(real_dir)
    #subprocess.run(f'bash download_dataset.sh {data_path}', shell=True)

    experiments_dir = data_path / 'experiments'
    os.makedirs(experiments_dir, exist_ok=True)


# 0. Download reference if necessary
@typechecked
def download_reference(ref_path: Path):
    species_reference(ref_path)


@typechecked
def simulate_reads_mp(chr_seq_path: Path, chr_dist_path: Path, chr_save_path: Path, chr_len, i: int):
    print(f'\nSegment {i}: Simulating reads {chr_save_path}')
    subprocess.run(f'./vendor/seqrequester/build/bin/seqrequester simulate -genome {chr_seq_path} ' \
                    f'-genomesize {chr_len} -coverage 42.4 -distribution {chr_dist_path} > {chr_save_path}',
                    shell=True, stdout=subprocess.DEVNULL)
    change_description(chr_save_path, multiline=False)


# 1. Simulate the sequences
@typechecked
def simulate_reads(data_path: Path, ref_path: Path, chr_dict: dict, simulator: dict):
    # Dict saying how much of simulated datasets for each chromosome do we need
    # E.g., {'chr1': 4, 'chr6': 2, 'chrX': 4}

    print(f'SETUP::simulate')
    validate_chrs(chr_dict)

    simulator = simulator_factory(simulator['name'], simulator['params'])
    simulator.run(ref_path, data_path, chr_dict, chr_lens)


# 2. Generate the graphs
@typechecked
def generate_graphs(data_path: Path, chr_dict: dict, assembler: dict = {}):
    print(f'SETUP::generate')

    validate_chrs(chr_dict)

    assembler = assembler_factory(assembler['name'], assembler['params'])
    data_path = data_path.resolve()

    for chrN in chr_dict:
        if '_r' in chrN:
            continue
        chr_sim_path = data_path / 'simulated' / f'{chrN}'
        chr_raw_path = chr_sim_path / 'raw'
        chr_prc_path = chr_sim_path / 'processed'
        n_raw = len(get_read_files(chr_raw_path, ['*.fastq'], override=True))
        n_prc = len(get_read_files(chr_prc_path, ['*.fastq'], override=True))
        n_diff = n_raw - n_prc
        print(f'SETUP::generate:: Generate {n_diff} graphs for {chrN}')

        graph_dataset_dbg.AssemblyGraphDataset(chr_sim_path, nb_pos_enc=None, assembler=assembler, generate=True)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--data', type=Path, default='data', help='Path to directory with simulated and real data')
    parser.add_argument('--refs', type=Path, default='data/references', help='Path to directory with reference information')
    parser.add_argument('--out', type=str, default=None, help='Output name for figures and models')
    parser.add_argument('--overfit', action='store_true', default=False, help='Overfit on the chromosomes in the train directory')
    #parser.add_argument('--species', type=str, default='drosophila_melanogaster', help='Species to analyse')

    return parser.parse_args()


def main(config_dir: str, config_name: str):
    args = parse_args()

    initialize(config_path=config_dir, job_name='DBG Pipeline', version_base='1.1')
    # +var_name=var_value syntax appends to config
    vargs = [f'+{k}={v}' for k,v in vars(args).items()]
    # compose but still use args for now until introducing structured config
    # use hydra only for new config options
    cfg = compose(config_name=config_name, overrides=vargs)

    data_path = args.data / cfg.species
    ref_path = args.refs / cfg.species
    out = args.out

    time_start = datetime.now()
    timestamp = time_start.strftime('%Y-%b-%d-%H-%M-%S')
    if out is None:
        out = timestamp

    file_structure_setup(data_path, ref_path)
    download_reference(ref_path)
    chr_dict = { "chr2L" : 3 }
    simulator = OmegaConf.to_container(cfg.simulator)
    simulate_reads(data_path, ref_path, chr_dict=chr_dict, simulator=simulator)
    assembler = OmegaConf.to_container(cfg.assembler)
    generate_graphs(data_path, chr_dict=chr_dict, assembler=assembler)


if __name__ == '__main__':
    main("config", "config.yaml")