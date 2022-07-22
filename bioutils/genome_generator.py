from pathlib import Path
from typing import Dict, Optional, Sequence

import hydra
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from omegaconf import OmegaConf
from typeguard import typechecked


@typechecked
def get_random_chromosome(sequence_length: int, base_p: Optional[Sequence[float]] = None) -> str:
    bases = np.random.choice(['A', 'C', 'G', 'T'], size=sequence_length, p=base_p)
    return ''.join(bases)



@typechecked
def get_random_genome(chromosomes: Dict[str, int]) -> Dict[str, str]:

    genome = {}
    for chr_name, chr_length in chromosomes.items():
        chr_seq = get_random_chromosome(chr_length)
        genome[chr_name] = chr_seq

    return genome


@typechecked
def save_chr_to_fasta(output_path: Path, chr_name: str, chr_seq: str, description: str = '', multiline: bool = True) -> Path:
    file_path = output_path / f'{chr_name}.fasta'
    new_fasta = [SeqIO.SeqRecord(seq=Seq(chr_seq), id=chr_name, description=description)]
    if multiline:
        SeqIO.write(new_fasta, file_path, 'fasta')
    else:
        # output entire DNA sequence in a single line
        new_fasta = [SeqIO.FastaIO.as_fasta_2line(record) for record in new_fasta]
        with open(file_path, 'w') as handle:
            handle.writelines(new_fasta)

    return file_path


@typechecked
def save_genome_to_fasta(output_path: Path, genome: Dict[str, str], description: str = '', multiline: bool = True) -> None:
    if not output_path.exists():
        output_path.mkdir(parents=True)
    for chr_name, chr_seq in genome.items():
        save_chr_to_fasta(output_path, chr_name, chr_seq, description, multiline)


@hydra.main(version_base=None, config_path=".", config_name="config")
def main(cfg):
    genome = get_random_genome(dict(cfg.chromosomes))

    species_path = Path(cfg.output_path) / cfg.species
    save_genome_to_fasta(species_path / 'chromosomes', genome, multiline=False)


if __name__ == '__main__':
    main()

