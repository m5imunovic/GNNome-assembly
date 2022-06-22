import subprocess
from pathlib import Path

from typeguard import typechecked

from assembler import Assembler, compose_cmd_params, get_read_files


class LaJolla(Assembler):
    @typechecked
    def _install(vendor_dir: Path) -> Path:
        assembler_root = vendor_dir / 'LJA'
        if not assembler_root.exists():
            print(f'SETUP::generate:: Download La Jolla Assembler')
            subprocess.run('git clone https://github.com/AntonBankevich/LJA.git', shell=True, cwd=str(vendor_dir))
            subprocess.run('cmake .', shell=True, cwd=str(assembler_root))
            subprocess.run('make -j 8', shell=True, cwd=str(assembler_root))

        return assembler_root

    @typechecked
    def _construct_exec_cmd(self, reads_path: Path, output_path: Path) -> list[str]:
        read_files = get_read_files(reads_path)
        reads_cmd_params = ' '.join([f'--reads {str(read_file)}' for read_file in read_files])
        out_cmd_param = f'-o {str(output_path)}'
        option_params = compose_cmd_params(self.cfg)

        asm_executable = self.assembler_root / 'bin/lja'
        return [f'{asm_executable} {option_params} {reads_cmd_params} {out_cmd_param}']

    @typechecked
    def run(self, raw_path: Path, tmp_path: Path, save_path: Path, *args, **kwargs):
        commands = self._construct_exec_cmd(raw_path, save_path)
        for cmd in commands:
            subprocess.run(cmd, shell=True)



