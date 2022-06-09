import argparse
import gzip
import multiprocessing as mp
import os
import subprocess
from datetime import datetime

from Bio import SeqIO

import graph_dataset_dbg

#from homo_sapiens import chm13_chr_lens as chr_lens
#from homo_sapiens import get_chr_dirs, species_specific_dirs, species_reference

from drosophila_melanogaster import release_6_plus_iso1_mt_chr_lens as chr_lens
from drosophila_melanogaster import get_chr_dirs, species_specific_dirs, species_reference


def validate_chrs(chrs):
    diff = set(chrs) - set(get_chr_dirs())
    if len(diff):
        raise Exception(f'Not using valid chromosome name(s): {diff}')
    return len(diff) == 0


def change_description(file_path, multiline=True):
    new_fasta = []
    for record in SeqIO.parse(file_path, file_path[-5:]): # 'fasta' for FASTA file, 'fastq' for FASTQ file
        des = record.description.split(",")
        id = des[0][5:]
        if des[1] == "forward":
            strand = '+'
        else:
            strand = '-'
        position = des[2][9:].split("-")
        start = position[0]
        end = position[1]
        record.id = id
        record.description = f'strand={strand}, start={start}, end={end}'
        new_fasta.append(record)
    if multiline:
        SeqIO.write(new_fasta, file_path, "fasta")
    else:
        # output entire DNA sequence in a single line
        new_fasta = [SeqIO.FastaIO.as_fasta_2line(record) for record in new_fasta]
        with open(file_path, 'w') as handle:
            handle.writelines(new_fasta)

        
def create_chr_dirs(pth):
    for chr_dir in get_chr_dirs():
        subprocess.run(f'mkdir {chr_dir}', shell=True, cwd=pth)
        subprocess.run(f'mkdir raw processed info asm_output graphia', shell=True, cwd=os.path.join(pth, f'{chr_dir}'))


# -1. Set up the data file structure
def file_structure_setup(data_path, ref_path):
    print(f'SETUP::filesystem:: Create directories for storing data')
    if not os.path.isdir(data_path):
        os.makedirs(data_path, exist_ok=True)

    if not os.path.isdir(ref_path):
        os.makedirs(ref_path, exist_ok=True)

    species_specific_dirs(ref_path)

    if 'chromosomes' not in os.listdir(ref_path):
        os.mkdir(os.path.join(ref_path, 'chromosomes'))
            
    if 'simulated' not in os.listdir(data_path):
        os.mkdir(os.path.join(data_path, 'simulated'))
        create_chr_dirs(os.path.join(data_path, 'simulated'))
    if 'real' not in os.listdir(data_path):
        pass
        #subprocess.run(f'bash download_dataset.sh {data_path}', shell=True)
        # os.mkdir(os.path.join(data_path, 'real'))
        # create_chr_dirs(os.path.join(data_path, 'real'))
    if 'experiments' not in os.listdir(data_path):
        os.mkdir(os.path.join(data_path, 'experiments'))


# 0. Download reference if necessary
def download_reference(ref_path):
    species_reference(ref_path)

def simulate_reads_mp(chr_seq_path, chr_dist_path, chr_save_path, chr_len, i):
    print(f'\nSegment {i}: Simulating reads {chr_save_path}')
    subprocess.run(f'./vendor/seqrequester/build/bin/seqrequester simulate -genome {chr_seq_path} ' \
                    f'-genomesize {chr_len} -coverage 42.4 -distribution {chr_dist_path} > {chr_save_path}',
                    shell=True, stdout=subprocess.DEVNULL)
    change_description(chr_save_path, multiline=False)


# 1. Simulate the sequences
def simulate_reads(data_path, ref_path, chr_dict):
    # Dict saying how much of simulated datasets for each chromosome do we need
    # E.g., {'chr1': 4, 'chr6': 2, 'chrX': 4}


    print(f'SETUP::simulate')
    validate_chrs(chr_dict)
    if 'vendor' not in os.listdir():
        os.mkdir('vendor')
    if 'seqrequester' not in os.listdir('vendor'):
        print(f'SETUP::simulate:: Download seqrequester')
        subprocess.run(f'git clone https://github.com/marbl/seqrequester.git', shell=True, cwd='vendor')
        subprocess.run(f'make', shell=True, cwd='vendor/seqrequester/src')


    

    data_path = os.path.abspath(data_path)
    chr_path = os.path.join(ref_path, 'chromosomes')
    len_path = os.path.join(ref_path, 'lengths')
    sim_path = os.path.join(data_path, 'simulated')
    simulation_data = []
    for chrN, n_need in chr_dict.items():
        if '_r' in chrN:
            continue
        chr_raw_path = os.path.join(sim_path, f'{chrN}/raw')
        n_have = len(os.listdir(chr_raw_path))
        if n_need <= n_have:
            continue
        else:
            n_diff = n_need - n_have
            print(f'SETUP::simulate:: Simulate {n_diff} datasets for {chrN}')
            # Simulate reads for chrN n_diff times
            chr_seq_path = os.path.join(chr_path, f'{chrN}.fasta')
            chr_dist_path = os.path.join(len_path, f'{chrN}.txt')
            if not os.path.exists(chr_dist_path):
                chr_dist_path = os.path.normpath(os.path.join(data_path, '../../defaults/pacbio-hifi-custom.txt'))
            if not os.path.exists(chr_dist_path):
                raise FileNotFoundError
            chr_len = chr_lens[chrN]
            for i in range(n_diff):
                idx = n_have + i
                chr_save_path = os.path.join(chr_raw_path, f'{idx}.fasta')
                simulation_data.append((chr_seq_path, chr_dist_path, chr_save_path, chr_len, i))

    # leave one processor free
    with mp.Pool(os.cpu_count() - 1) as pool:
        pool.starmap(simulate_reads_mp, simulation_data)


# 2. Generate the graphs
def generate_graphs(data_path, chr_dict):
    print(f'SETUP::generate')

    validate_chrs(chr_dict)

    if 'rust-mdbg' not in os.listdir('vendor'):
        print(f'SETUP::generate:: Download rust-mdbg')
        subprocess.run(f'git clone https://github.com/ekimb/rust-mdbg.git', shell=True, cwd='vendor')
    subprocess.run(f'cargo build --release', shell=True, cwd='vendor/rust-mdbg')


    data_path = os.path.abspath(data_path)

    for chrN in chr_dict:
        if '_r' in chrN:
            continue
        chr_sim_path = os.path.join(data_path, 'simulated', f'{chrN}')
        chr_raw_path = os.path.join(chr_sim_path, 'raw')
        chr_prc_path = os.path.join(chr_sim_path, 'processed')
        n_raw = len(os.listdir(chr_raw_path))
        n_prc = len(os.listdir(chr_prc_path))
        n_diff = n_raw - n_prc
        print(f'SETUP::generate:: Generate {n_diff} graphs for {chrN}')
        specs = {
            'threads': 8,
            'k': 31,
            'l': 8,
            'd': 0.001,
        }

        graph_dataset_dbg.AssemblyGraphDataset(chr_sim_path, nb_pos_enc=None, specs=specs, generate=True)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type=str, default='data', help='Path to directory with simulated and real data')
    parser.add_argument('--refs', type=str, default='data/references', help='Path to directory with reference information')
    parser.add_argument('--out', type=str, default=None, help='Output name for figures and models')
    parser.add_argument('--overfit', action='store_true', default=False, help='Overfit on the chromosomes in the train directory')
    parser.add_argument('--species', type=str, default='drosophila_melanogaster', help='Species to analyse')
    parser.add_argument('--assembler', type=str, default='rust-mdbg', help='Assembler to use for creating assembly graph')
    args = parser.parse_args()

    data_path = os.path.join(args.data, args.species)
    ref_path = os.path.join(args.refs, args.species)
    out = args.out
    overfit = args.overfit

    time_start = datetime.now()
    timestamp = time_start.strftime('%Y-%b-%d-%H-%M-%S')
    if out is None:
        out = timestamp


    file_structure_setup(data_path, ref_path)
    download_reference(ref_path)
    chr_dict = { "chr3R" : 1 }
    simulate_reads(data_path, ref_path, chr_dict=chr_dict)
    generate_graphs(data_path, chr_dict=chr_dict)