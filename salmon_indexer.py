import os
import sys


class SalmonIndexer(object):
    def __init__(self, genome, transcriptome, outdir):
        self.genome = genome
        self.transcriptome = transcriptome
        self.outdir = outdir
        self.catome = f'{self.outdir}/cat_transcriptome_genome.fasta'
        self.__check_dir()

    def __check_dir(self):
        if not os.path.exists(self.outdir):
            os.system(f'mkdir {self.outdir}')

    def create_decoys(self):
        cmd = f'grep \'^>\' {self.genome} | cut -d \" \" -f 1 > {self.outdir}/decoys.txt'
        os.system(cmd)
        os.system(f'sed -i.bak -e \'s/>//g\' {self.outdir}/decoys.txt')

    def cat_omes(self):
        cmd = f'cat {self.transcriptome} {self.genome} > {self.catome}'
        os.system(cmd)

    def index(self):
        cmd = f'salmon index -t {self.catome} -d {self.outdir}/decoys.txt -p 12 -i {self.outdir}/salmon_index'
        os.system(cmd)


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: salmon_indexer.py <genome_fasta> <transcriptome_fasta> <outdir>')
    else:
        data = SalmonIndexer(genome=sys.argv[1], transcriptome=sys.argv[2], outdir=sys.argv[3])
        data.create_decoys()
        data.cat_omes()
        data.index()
