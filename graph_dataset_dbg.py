import os
import subprocess
from pathlib import Path
from typing import Optional

import dgl
from dgl.data import DGLDataset
from assembler import Assembler


class AssemblyGraphDataset(DGLDataset):
    """
    A dataset to store the assembly graphs.

    A class that inherits from the DGLDataset and extends the
    functionality by adding additional attributes and processing the
    graph data appropriately.

    Attributes
    ----------
    root : str
        Root directory consisting of other directories where the raw
        data can be found (reads in FASTQ format), and where all the
        processing results are stored.
    tmp_dir : str
        Directory inside root where mid-results (output of the raven 
        assembler) is stored
    info_dir : str
        Directory where additional graph information is stored
    assembler_path : str
        Path to the DNA assembler
    """

    def __init__(self, root: str, assembler: Assembler, nb_pos_enc: int=10, generate: bool=False):
        """
        Parameters
        ----------
        root : str
            Root directory consisting of other directories where the raw
            data can be found (reads in FASTQ format), and where all the
            processing results are stored
        nb_pos_end : int
            Dimension of vector for storing the positional encoding of nodes
        specs : dict
            Dictionary with specifications for running assembler
        generate : bool
            Boolean which, if true, skips computing positional encoding and
            loading of the graphs.
        """
        self.root = Path(root).resolve()
        self.assembler = assembler
        required_dirs = ['raw', 'asm_output', 'processed', 'info']
        for required_dir in required_dirs:
            os.makedirs(self.root / required_dir, exist_ok=True)
        raw_dir = self.root / 'raw'
        save_dir = self.root / 'processed'
        self.tmp_dir = self.root / 'asm_output'
        self.info_dir = self.root / 'info'
        super().__init__(name='assembly_graphs', raw_dir=raw_dir, save_dir=save_dir)

        self.graph_list = []
        if not generate:
            pass
            # for file in os.listdir(self.save_dir):
            #     idx = int(file[:-4])
            #     graph = dgl.load_graphs(os.path.join(self.save_dir, file))[0][0]
            #     graph = preprocess_graph(graph, self.root, idx)
            #     if nb_pos_enc is not None:
            #         graph = add_positional_encoding(graph, nb_pos_enc) 
            #     #graph, _ = dgl.khop_in_subgraph(graph, 390, k=20) # DEBUG !!!!
            #     print(f'DGL graph idx={idx} info:\n',graph)
            #     self.graph_list.append((idx, graph))
            # self.graph_list.sort(key=lambda x: x[0])

    def has_cache(self):
        """Check if the raw data is already processed and stored."""
        return len(os.listdir(self.save_dir)) >= len(os.listdir(self.raw_dir))

    def __len__(self):
        return len(os.listdir(self.save_dir))

    def __getitem__(self, idx):
        i, graph = self.graph_list[idx]
        return i, graph

    def process(self):
        """Process the raw data and save it on the disk."""

        graphia_dir = os.path.join(self.root, 'graphia')
        if not os.path.isdir(graphia_dir):
            os.makedirs(graphia_dir, exist_ok=True)

        print(f'====> ASSEMBLY\n')

        self.assembler.run(self.raw_dir, self.tmp_dir, self.save_dir)

            
        print(f'\nAssembler generated the graph! Processing...')
        #processed_path = os.path.join(self.save_dir, f'{idx}.dgl')
        #graph, pred, succ, reads, edges, labels = graph_parser.from_csv(os.path.join(self.tmp_dir, f'{idx}_graph_1.csv'), reads_path)
        print(f'Parsed assembler output! Saving files...')

        #dgl.save_graphs(processed_path, graph)
        #pickle.dump(pred, open(f'{self.info_dir}/{idx}_pred.pkl', 'wb'))
        #pickle.dump(succ, open(f'{self.info_dir}/{idx}_succ.pkl', 'wb'))
        #pickle.dump(reads, open(f'{self.info_dir}/{idx}_reads.pkl', 'wb'))
        #pickle.dump(edges, open(f'{self.info_dir}/{idx}_edges.pkl', 'wb'))
        #pickle.dump(labels, open(f'{self.info_dir}/{idx}_labels.pkl', 'wb'))

        #graphia_path = os.path.join(graphia_dir, f'{idx}_graph.txt')
        #graph_parser.print_pairwise(graph, graphia_path)
        #print(f'Processing of graph {idx} generated from {fastq} done!\n')

