import os
import subprocess
from pathlib import Path

from typeguard import typechecked

from assembler import Assembler, compose_cmd_params, get_read_files


class RustMDBG(Assembler):
    @typechecked
    def _install(self, vendor_dir: Path):
        assembler_root = vendor_dir / 'rust-mdbg'
        if not assembler_root.exists():
            print(f'SETUP::generate:: Download rust-mdbg')
            subprocess.run('git clone https://github.com/ekimb/rust-mdbg.git', shell=True, cwd=str(vendor_dir))
            subprocess.run('cargo build --release', shell=True, cwd=str(assembler_root))

        return assembler_root

    @typechecked
    def _construct_exec_cmd(self, raw_path: Path, prefix = str) -> list[str]:
        assembler_path = self.assembler_root / 'target/release/rust-mdbg'
        to_basespace = self.assembler_root / 'utils/complete_gfa.py'
        read_files = get_read_files(raw_path, pattern=['.fa.gz', 'fasta.gz'])
        reads_params = ' '.join([f'{str(read_file)}' for read_file in read_files])
        option_params = compose_cmd_params(self.cfg)
        prefix_param = f'--prefix {prefix}'

        return [
            f'{assembler_path} {option_params} {prefix_param} {reads_params}',
            f'python {to_basespace} {prefix}. {prefix}.gfa'
        ]

    @typechecked
    def run(self, raw_path: Path, tmp_path: Path, save_path: Path, *args, **kwargs):
        n_have = len(os.listdir(save_path))
        n_need = len(os.listdir(raw_path))
        for cnt, idx in enumerate(range(n_have, n_need)):
            prefix = f'graph_{idx}_{cnt}'
            print(f'Path to the reads: {str(raw_path)}')
            commands = self._construct_exec_cmd(raw_path, prefix)
            for cmd in commands:
                subprocess.run(cmd, shell=True, cwd=str(tmp_path))