import os
import subprocess
from abc import abstractmethod
from pathlib import Path
from typing import Optional

from typeguard import typechecked

from cmd_utils import compose_cmd_params



def setup_rust_mdbg_assembler(vendor_dir: Path) -> Path:
    assembler_root = vendor_dir / 'rust-mdbg'
    if 'rust-mdbg' not in os.listdir('vendor'):
        print(f'SETUP::generate:: Download rust-mdbg')
        subprocess.run('git clone https://github.com/ekimb/rust-mdbg.git', shell=True, cwd=str(vendor_dir))
        subprocess.run('cargo build --release', shell=True, cwd=str(assembler_root))

    return assembler_root


@typechecked
def rust_mdbg(params: dict, vendor_dir: Path) -> dict:
    assembler_root = setup_rust_mdbg_assembler(vendor_dir=vendor_dir)
    assembler_path = assembler_root / 'target/release/rust-mdbg'
    to_basespace = assembler_root / 'utils/complete_gfa.py'

    if params is None:
        params = {
            'threads': 8,
            'k': 31,
            'l': 8,
            'd': 0.001,
        }
    return {
        'asm_cmd': f'{assembler_path} {compose_cmd_params(params)}',
        'to_basespace_cmd': f'{to_basespace}'
    }


@typechecked
def setup_la_jolla_assembler(vendor_dir: Path) -> Path:
    assembler_root = vendor_dir / 'LJA'
    if not assembler_root.exists():
        print(f'SETUP::generate:: Download La Jolla Assembler')
        subprocess.run('git clone https://github.com/AntonBankevich/LJA.git', shell=True, cwd=str(vendor_dir))
        subprocess.run('cmake .', shell=True, cwd=str(assembler_root))
        subprocess.run('make -j 8', shell=True, cwd=str(assembler_root))

    return assembler_root


@typechecked
def la_jolla(params: dict, vendor_dir: Path) -> dict:
    assembler_root = setup_la_jolla_assembler(vendor_dir=vendor_dir)
    assembler_path = assembler_root / 'bin/lja'
    if params is None:
        params = {
            'threads': 8,
            'k': 51,
            'K': 2001,
        }

    return {
        'asm_cmd': f'{assembler_path} {compose_cmd_params(params)}'
    }


class Assembler:
    def __init__(self, cfg: dict, vendor_dir: Path):
        self.cfg = cfg
        self.assembler_root = self._install(vendor_dir)

    @typechecked
    @abstractmethod
    def _install(self, vendor_dir: Path):
        pass

    def pre_assembly_step(self, *args, **kwargs):
        pass

    def post_assembly_step(self, *args, **kwargs):
        pass

    def run(self, raw_path: Path, tmp_path: Path, save_path: Path, *args, **kwargs):
        pass


@typechecked
def assembler_factory(assembler: str, params: Optional[dict] = None) -> Assembler:
    vendor_dir: Path = Path('vendor').resolve()
    if assembler == 'rust-mdbg':
        from rust_mdbg import RustMDBG
        return RustMDBG(cfg=params, vendor_dir=vendor_dir)
    if assembler == 'LJA':
        from la_jolla import LaJolla
        return LaJolla(cfg=params, vendor_dir=vendor_dir)
    
