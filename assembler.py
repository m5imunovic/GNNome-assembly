import os
import subprocess
from pathlib import Path
from typing import Optional

from typeguard import typechecked


@typechecked
def compose_cmd_params(params: dict) -> str:
    append = params.pop('append', None)
    cmd_params = ' '.join([f'--{k} {v}' for k, v in params.items()])
    append_params = f' {append}' if 'append' in params and params['append'] is not None else ''
    return cmd_params + append_params


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
def assembler_factory(assembler: str, params: Optional[dict] = None) -> dict:
    vendor_dir: Path = Path('vendor').resolve()
    if assembler == 'rust-mdbg':
        return rust_mdbg(params, vendor_dir=vendor_dir)
    