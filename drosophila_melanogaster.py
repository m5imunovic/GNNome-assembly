import os
from pathlib import Path

from Bio import Entrez
from typeguard import typechecked

Entrez.email = "entrez.user@gmail.com"

REFERENCE_NAME = 'release_6_plus_iso1_mt'

download = {
  'accessions': {
    'chrX':  'AE014298.5',
    'chr2L': 'AE014134.6',
    'chr2R': 'AE013599.5',
    'chr3L': 'AE014296.5',
    'chr3R': 'AE014297.3',
    'chr4':  'AE014135.4',
    'chrY':  'CP007106.1'
  },
  'format': 'fasta',
  'db' : 'nuccore'
}

# download_aws = "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR11906526/SRR11906526"
# https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR10238607/SRR10238607
# source: https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR11906526


release_6_plus_iso1_mt_chr_lens = {
    'chr2L': 23513712,
    'chr2R': 25286936,
    'chr3L': 28110227,
    'chr3R': 32079331,
    'chr4' : 1348131,
    'chrX' : 23542271,
    'chrY' : 3667352
}


def get_chr_dirs():
    return list(release_6_plus_iso1_mt_chr_lens.keys())


@typechecked
def species_specific_dirs(ref_path: Path):
    if REFERENCE_NAME not in os.listdir(ref_path):
        os.makedirs(ref_path / REFERENCE_NAME, exist_ok=True)


@typechecked
def species_reference(ref_path: Path):
    chr_path = ref_path / 'chromosomes'
    chr_dirs = get_chr_dirs()

    for chr_dir in chr_dirs:
        fasta_path = chr_path / f'{chr_dir}.fasta'
        if fasta_path.exists():
            print(f'Path {fasta_path} already exists. Skipping!')
            continue
        handle = Entrez.efetch(db=download['db'], id=download['accessions'][chr_dir], rettype=download['format'], retmode='text')
        lines = [line for line in handle]
        handle.close()

        with open(fasta_path, 'w') as f:
            f.writelines(lines)
