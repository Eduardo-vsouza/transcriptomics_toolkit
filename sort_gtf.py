import sys


class GTFSorting(object):
    def __init__(self, gtf):
        """

        :param gtf: GTF file containing unordered features
        Written by: Eduardo Vieira de Souza
        """

        self.gtf = gtf

    def parse(self):
        genes = {}  # contains gene_id and the GTF line
        transcripts = {}   # contains gene_id and a dictionary with key, value: transcript_id: GTF line
        exons = {}   # contains transcript_id and the GTF line
        with open(self.gtf, 'r') as handler:
            lines = handler.readlines()
            for line in lines:
                cols = line.split("\t")
                attrs = cols[8].split(";")
                if cols[2] == 'gene':

                    for a in attrs:
                        if 'gene_id' in a:
                            gene_id = self.__get_gene_id(a)
                            if gene_id not in genes:
                                genes[gene_id] = []
                            # hierarchy[gene_id].append(line)
                            genes[gene_id].append(line)

            for line in lines:
                cols, attrs = self.__get_attrs(line)
                if cols[2] == 'transcript':
                    for a in attrs:
                        if 'gene_id' in a:
                            gene_id = self.__get_gene_id(a)
                            if gene_id not in transcripts:
                                transcripts[gene_id] = {}
                                for at in attrs:
                                    if 'transcript_id' in at:
                                        transcript_id = self.__get_gene_id(at)
                                        transcripts[gene_id][transcript_id] = line
            for line in lines:
                cols, attrs = self.__get_attrs(line)
                if cols[2] == 'exon' or cols[2] == 'CDS':
                    for a in attrs:
                        if 'gene_id' in a:
                            gene_id = self.__get_gene_id(a)
                            for at in attrs:
                                if 'transcript_id' in at:
                                    transcript_id = self.__get_gene_id(at)
                                    if gene_id in transcripts:
                                        if transcript_id not in exons:
                                            exons[transcript_id] = []
                                        exons[transcript_id].append(line)
        self.genes = genes
        self.transcripts = transcripts
        self.exons = exons

    def write(self, output):
        lines = []
        for gene in self.genes:
            lines.append(self.genes[gene][0])
            transcripts = self.transcripts[gene]
            for rna in transcripts:
                lines.append(transcripts[rna])
                exons = self.exons[rna]
                for exon in exons:
                    lines.append(exon)
        with open(output, 'w') as outfile:
            outfile.writelines(lines)

    @staticmethod
    def __get_attrs(line):
        cols = line.split("\t")
        attrs = cols[8].split(";")
        return cols, attrs

    @staticmethod
    def __get_gene_id(attribute):
        if attribute.startswith(" "):
            gene_id = attribute.split(" ")[2].replace("\"", "")
        else:
            gene_id = attribute.split(" ")[1].replace("\"", "")

        return gene_id


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('\nThis script sorts a GTF file with the following hierarchy: gene > transcript > exon/CDS. Every gene will'
              'be followed by its transcripts and each transcript by its exons and CDSs. This will exclude any other '
              'features from the GTF file.\n\nUsage: sort_gtf.py <gtf_to_fix> <output>\n')
    else:
        data = GTFSorting(gtf=sys.argv[1])
        data.parse()
        data.write(output=sys.argv[2])