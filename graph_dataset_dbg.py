import os
import subprocess

import dgl
from dgl.data import DGLDataset


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

    def __init__(self, root, nb_pos_enc=10, specs=None, generate=False):
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
        self.root = os.path.abspath(root)
        self.specs = specs
        required_dirs = ['raw', 'asm_output', 'processed', 'info']
        for required_dir in required_dirs:
            if required_dir not in os.listdir(self.root):
                subprocess.run(f'mkdir {required_dir}', shell=True, cwd=self.root)
        raw_dir = os.path.join(self.root, 'raw')
        save_dir = os.path.join(self.root, 'processed')
        self.tmp_dir = os.path.join(self.root, 'asm_output')
        self.info_dir = os.path.join(self.root, 'info')
        self.assembler_root = os.path.abspath('vendor/rust-mdbg/')
        self.assembler_path = os.path.join(self.assembler_root, 'target/release/rust-mdbg')
        self.to_basespace = os.path.join(self.assembler_root, 'utils/complete_gfa.py')
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
        if self.specs is None:
            threads = 8
            k, l, d = 31, 7, 0.01
        else:
            threads = self.specs['threads']
            k, l, d = self.specs['k'], self.specs['l'], self.specs['d']

        graphia_dir = os.path.join(self.root, 'graphia')
        if not os.path.isdir(graphia_dir):
            os.mkdir(graphia_dir)

        print(f'====> ASSEMBLY\n')

        n_have = len(os.listdir(self.save_dir))
        n_need = len(os.listdir(self.raw_dir))
        for cnt, idx in enumerate(range(n_have, n_need)):
            fastq = f'{idx}.fasta'
            prefix = f'graph_{idx}_{cnt}'
            print(f'Step {cnt}: generating graphs for reads in {fastq}')
            reads_path = os.path.abspath(os.path.join(self.raw_dir, fastq))
            print(f'Path to the reads: {reads_path}')
            print(f'Starting assembler at: {self.assembler_path}')
            print(f'Parameters: k={k} l={l} d={d} t={threads}')
            subprocess.run(f'{self.assembler_path} -k {k} -l {l} -d {d} --threads {threads} --prefix {prefix} {reads_path}', shell=True, cwd=self.tmp_dir)
            to_basespace_path = os.path.join(self.assembler_root, 'utils/complete_gfa.py')
            subprocess.run(f'python {to_basespace_path} {prefix}. {prefix}.gfa', shell=True, cwd=self.tmp_dir)

            
            print(f'\nAssembler generated the graph! Processing...')
            processed_path = os.path.join(self.save_dir, f'{idx}.dgl')
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

