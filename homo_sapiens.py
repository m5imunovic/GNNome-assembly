import os
import gzip

import requests
from Bio import SeqIO
from tqdm import tqdm


chm13_chr_lens = {
    'chr1' : 248387328,
    'chr2' : 242696752,
    'chr3' : 201105948,
    'chr4' : 193574945,
    'chr5' : 182045439,
    'chr6' : 172126628,
    'chr7' : 160567428,
    'chr8' : 146259331,
    'chr9' : 150617247,
    'chr10': 134758134,
    'chr11': 135127769,
    'chr12': 133324548,
    'chr13': 113566686,
    'chr14': 101161492,
    'chr15': 99753195,
    'chr16': 96330374,
    'chr17': 84276897,
    'chr18': 80542538,
    'chr19': 61707364,
    'chr20': 66210255,
    'chr21': 45090682,
    'chr22': 51324926,
    'chrX' : 154259566,
}

def get_chr_dirs():
    return list(chm13_chr_lens.keys())


def species_specific_dirs(ref_path):
    if 'CHM13' not in os.listdir(ref_path):
        os.mkdir(os.path.join(ref_path, 'CHM13'))


def species_reference(ref_path):
    chm_path = os.path.join(ref_path, 'CHM13')
    chr_path = os.path.join(ref_path, 'chromosomes')
    chm13_url = 'https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.1.fasta.gz'
    chm13_path = os.path.join(chm_path, 'chm13.draft_v1.1.fasta.gz')

    if len(os.listdir(chm_path)) == 0:
        # Download the CHM13 reference
        # Code for tqdm from: https://stackoverflow.com/questions/37573483/progress-bar-while-download-file-over-http-with-requests
        print(f'SETUP::download:: CHM13 not found! Downloading...')
        response = requests.get(chm13_url, stream=True)
        total_size_in_bytes= int(response.headers.get('content-length', 0))
        block_size = 1024 #1 Kibibyte
        progress_bar = tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True)

        with open(chm13_path, 'wb') as file:
            for data in response.iter_content(block_size):
                progress_bar.update(len(data))
                file.write(data)
        progress_bar.close()
        if total_size_in_bytes != 0 and progress_bar.n != total_size_in_bytes:
            print("ERROR, something went wrong")

    if len(os.listdir(chr_path)) == 0:
        # Parse the CHM13 into individual chromosomes
        print(f'SETUP::download:: Split CHM13 per chromosome')
        with gzip.open(chm13_path, 'rt') as f:
            for record in SeqIO.parse(f, 'fasta'):
                SeqIO.write(record, os.path.join(chr_path, f'{record.id}.fasta'), 'fasta')