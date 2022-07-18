import multiprocessing as mp
import os
import subprocess
from abc import abstractmethod
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from typeguard import typechecked

from cmd_utils import compose_cmd_params, get_read_files


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



class Simulator:
    def __init__(self, cfg: dict, vendor_dir: Path):
        self.cfg = cfg
        self.simulator_root = self._install(vendor_dir)

    @typechecked
    @abstractmethod
    def _install(self, vendor_dir: Path):
        pass

    def pre_simulation_step(self, *args, **kwargs):
        pass

    def post_simulation_step(self, *args, **kwargs):
        pass

    def run(self, raw_path: Path, tmp_path: Path, save_path: Path, *args, **kwargs):
        pass


class Seqrequester(Simulator):
    @typechecked
    def _install(self, vendor_dir: Path) -> Path:
        simulator_root = vendor_dir / 'seqrequester'
        if not simulator_root.exists():
            print(f'SETUP::generate:: Download Seqrequester')
            subprocess.run(f'git clone https://github.com/marbl/seqrequester.git', shell=True, cwd='vendor')
            subprocess.run(f'make', shell=True, cwd='vendor/seqrequester/src')

        return simulator_root

    @typechecked
    def _construct_exec_cmd(self, raw_path: Path, prefix = str) -> list[str]:
        simulator_path = self.simulator_root / 'src/seqrequester'

    
    @typechecked
    def run(self, ref_path: Path, data_path: Path, chr_dict: dict, chr_lens: dict, *args, **kwargs):
        data_path = data_path.resolve()

        chr_path = ref_path / 'chromosomes'
        len_path = ref_path / 'lengths'
        sim_path = data_path / 'simulated'
        simulation_data = []
        for chrN, n_need in chr_dict.items():
            if '_r' in chrN:
                continue
            chr_raw_path = sim_path / f'{chrN}' / 'raw'
            n_have = len(get_read_files(chr_raw_path, ['.fasta'], override=True))
            if n_need <= n_have:
                continue
            else:
                n_diff = n_need - n_have
                print(f'SETUP::simulate:: Simulate {n_diff} datasets for {chrN}')
                # Simulate reads for chrN n_diff times
                chr_seq_path = chr_path / f'{chrN}.fasta'
                chr_dist_path = len_path / f'{chrN}.txt'
                if not chr_dist_path.exists():
                    chr_dist_path = (data_path / '../../defaults/pacbio-hifi-custom.txt').resolve()
                    if not chr_dist_path.exists():
                        raise FileNotFoundError

                chr_len = chr_lens[chrN]
                for i in range(n_diff):
                    idx = n_have + i
                    chr_save_path = chr_raw_path / f'{idx}.fastq'
                    simulation_data.append((chr_seq_path, chr_dist_path, chr_save_path, chr_len, i))

            # leave one processor free
            with mp.Pool(os.cpu_count() - 1) as pool:
                pool.starmap(self.simulate_reads_mp, simulation_data)

    @typechecked
    def simulate_reads_mp(self, chr_seq_path: Path, chr_dist_path: Path, chr_save_path: Path, chr_len, i: int):
        print(f'\nSegment {i}: Simulating reads {chr_save_path}')
        subprocess.run(f'{self.simulator_root}/build/bin/seqrequester simulate -genome {chr_seq_path} ' \
                        f'-genomesize {chr_len} -coverage 42.4 -distribution {chr_dist_path} > {chr_save_path}',
                        shell=True, stdout=subprocess.DEVNULL)
        change_description(chr_save_path, multiline=False)



class PbSim2(Simulator):
    @typechecked
    def _install(self, vendor_dir: Path) -> Path:
        simulator_root = vendor_dir / 'pbsim2'
        if not simulator_root.exists():
            print(f'SETUP::generate:: Download PbSim2')
            subprocess.run('git clone https://github.com/yukiteruono/pbsim2.git', shell=True, cwd=str(vendor_dir))
            subprocess.run(f'./configure', shell=True, cwd='vendor/pbsim2')
            subprocess.run(f'make', shell=True, cwd='vendor/pbsim2')

        return simulator_root

    @typechecked
    def _construct_exec_cmd(self, ref_path: Path, chr_save_path: Path, prefix: str) -> list[str]:
        simulator_path = self.simulator_root / 'src/pbsim'
        read_files = get_read_files(ref_path, pattern=['.fasta]'], override=True)
        read_params = ' '.join(f'{str(read_file)}' for read_file in read_files)
        option_params = compose_cmd_params(self.cfg)
        prefix_param = f'--prefix {prefix}'

        return [
            f'{simulator_path} {option_params} {prefix_param} {read_params}',
            f'mv {prefix}_0001.* {chr_save_path}'
        ]

    @typechecked
    def run(self, ref_path: Path, data_path: Path, chr_dict: dict, *args, **kwargs):
        data_path = data_path.resolve()

        chr_path = ref_path / 'chromosomes'
        sim_path = data_path / 'simulated'
        simulation_data = []
        for chrN, n_need in chr_dict.items():
            if '_r' in chrN:
                continue
            chr_raw_path = sim_path / f'{chrN}' / 'raw'
            n_have = len(get_read_files(chr_raw_path, pattern=['.fastq'], override=True))
            if n_need <= n_have:
                continue
            else:
                n_diff = n_need - n_have
                print(f'SETUP::simulate:: Simulate {n_diff} datasets for {chrN}')
                # Simulate reads for chrN n_diff times

                chr_seq_path = chr_path / f'{chrN}.fasta'
                for i in range(n_diff):
                    idx = n_have + i
                    chr_save_path = chr_raw_path
                    simulation_data.append((chr_seq_path, chr_save_path, str(idx), i))

            # leave one processor free
            with mp.Pool(os.cpu_count() - 1) as pool:
                pool.starmap(self.simulate_reads_mp, simulation_data)

    @typechecked
    def simulate_reads_mp(self, chr_seq_path: Path, chr_save_path: Path, prefix: str, i: int):
        print(f'\nSegment {i}: Simulating reads {chr_seq_path}')
        commands = self._construct_exec_cmd(chr_seq_path, chr_save_path, prefix)
        for cmd in commands:
            subprocess.run(cmd, shell=True)


@typechecked
def simulator_factory(simulator: str, params: Optional[dict] = None) -> Simulator:
    vendor_dir: Path = Path('vendor').resolve()
    if simulator == 'seqrequester':
        return Seqrequester(cfg=params, vendor_dir=vendor_dir)
    if simulator == 'pbsim2':
        return PbSim2(cfg=params, vendor_dir=vendor_dir)