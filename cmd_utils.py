from pathlib import Path
from typing import Optional
from typeguard import typechecked


@typechecked
def get_read_files(read_path: Path, pattern: Optional[list[str]] = None) -> list[Path]:
    if read_path.is_file():
        return [read_path]
    else:
        default_pattern = ['.fasta', '.fa', '.fastq', '.fq']
        default_pattern += pattern if pattern else []
        result = [p.resolve() for p in read_path.glob('**/*') if p.suffix in set(default_pattern)]
        return result


@typechecked
def compose_cmd_params(params: dict) -> str:
    short_params = ' '.join([f'-{k} {v}' for k, v in params['short'].items()]) if 'short' in params else ''
    long_params = ' '.join([f'--{k} {v}' for k, v in params['long'].items()]) if 'long' in params else ''
    append_params = f'{params["append"]}' if 'append' in params and params['append'] is not None else ''
    combined_params = ' '.join([short_params,  long_params,  append_params])
    return combined_params.strip()
