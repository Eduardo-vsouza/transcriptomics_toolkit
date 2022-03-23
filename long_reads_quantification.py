import os
import sys


class LongQuantifier(object):
    def __init__(self, nanopore_reads, transcriptome_fasta, threads, outdir):
        self.nanoporeReads = nanopore_reads.split(",")
        self.transcriptome = transcriptome_fasta
        self.threads = threads
        self.outdir = outdir
        self.__check_dir(self.outdir)

    def __check_dir(self, folder):
        if not os.path.exists(folder):
            os.system(f'mkdir {folder}')

    def align(self):
        for file in self.nanoporeReads:
            name = file.split("/")[-1]
            print(f"Aligning the nanopore reads from {file} to the transcriptome using minimap2.\n")
            cmd = f'minimap2 -a -N 10 -t {self.threads} -x map-ont {self.transcriptome} {file} | samtools view -b | ' \
                  f'samtools sort -@ {self.threads} -o {self.outdir}/{name.replace(".fastq", "_sorted.bam")}'
            # running minimap2 on -x map-ont mode because the splice mode is not necessary here. Since we are mapping
            # reads to a transcriptome fasta file that was generated from the spliced exons of the GTf file, there is
            # no need to align using a splice-aware mode.
            os.system(cmd)

    def quantify_with_salmon(self):
        outdir = f'{self.outdir}/salmon_quantification'
        self.__check_dir(outdir)
        files = os.listdir(f'{self.outdir}/')
        for file in files:

            if file.endswith('.bam'):
                print(f"Quantifying transcripts for {file}\n")
                cmd = f'salmon quant -t {self.transcriptome} -l A --noErrorModel -a {self.outdir}/{file} -o ' \
                      f'{file[:-4]}/quant'
                ## the --noErrorModel is necessary for long reads as this model was developed for short reads.
                # long reads have much higher error and minimap2 outputs the strings without the CIGAR characters
                os.system(cmd)


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: long_reads_quantification.py <nanopore_reads_comma_separated> <transcriptome.fasta> '
              '<threads> <outdir>')
    else:
        data = LongQuantifier(nanopore_reads=sys.argv[1], transcriptome_fasta=sys.argv[2], threads=sys.argv[3],
                              outdir=sys.argv[4])
        data.align()
        # data.quantify_with_salmon()