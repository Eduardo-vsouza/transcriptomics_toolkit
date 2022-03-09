import sys
import numpy as np
import os


class WithinSampleComparison(object):
    def __init__(self, gtf_list):
        self.gtf_list = gtf_list.split(",")
        self.tpms = {}
        self.isoforms = {}
        self.tpmMeans = {}

        self.pctTPMs = {}
        self.genes = {}

    def parse_gtf(self):
        for gtf in self.gtf_list:
            print(f"Calculating TPMs for the isoforms in {gtf}\n")
            with open(gtf, 'r') as handler:
                lines = handler.readlines()
                for line in lines:
                    if not line.startswith("#"):
                        # if line.startswith("PCMN"):
                        cols = line.split("\t")
                        if cols[2] == 'transcript':
                            attrs = cols[8].split(";")
                            for attr in attrs:
                                if attr.startswith(" transcript_id"):
                                    transcript = attr.split(" ")[2].replace("\"", "")
                                    if 'PCMN' in transcript:
                                        original = transcript
                                    else:
                                        original = '.'.join(transcript.split(".")[:2])
                                    if original not in self.isoforms:
                                        self.isoforms[original] = []
                                    isoform = transcript
                                    if isoform not in self.isoforms[original]:
                                        self.isoforms[original].append(isoform)
                                if attr.startswith(" TPM"):
                                    tpm = float(attr.split(" ")[2].replace("\"", "")[:5])
                                    # if not tpm >= 0:
                                    #     print(tpm)
                                    # if tpm == 'nan':
                                    #     tpm = 0
                                    if isoform not in self.tpms:
                                        self.tpms[isoform] = []
                                    self.tpms[isoform].append(tpm)
                                if attr.startswith(" ref_gene_id") or attr.startswith(" reference_id"):
                                    gene = attr.split(" ")[2].replace("\"", "")
                                    if original not in self.genes:
                                        self.genes[original] = gene

        # print(self.isoforms)
        # print(self.tpms)

    def get_means(self):
        for tpm in self.tpms:
            # if 'PCMN' in tpm:
                # if 0.0 in self.tpms[tpm]:
                #     for i in self.tpms[tpm]:
                #         if i > 0:
                # print(tpm, self.tpms[tpm])
                # print(np.mean(self.tpms[tpm]))
            self.tpmMeans[tpm] = np.mean(self.tpms[tpm])

    def get_pcts(self):

        for original in self.isoforms:

            total = 0
            for isoform in self.isoforms[original]:
                total += self.tpmMeans[isoform]
            for isoform in self.isoforms[original]:
                pct = (self.tpmMeans[isoform]/total)*100
                # print(pct)
                # if isoform not in self.pctTPMs:
                #     self.pctTPMs[isoform] =
                self.pctTPMs[isoform] = pct


    def create_df(self, output):
        file = ['transcript\t', 'isoform_tpms\n']
        for original in self.isoforms:
            if original in self.genes:
                line = f'{original}_({self.genes[original]})\t'
            else:
                line = f'{original}\t'
            for isoform in self.isoforms[original]:
                # for tpm in self.tpms[isoform]:
                # if isoform in self.genes:
                #     line += f
                    # {'{isoform} ({self.genes[isoform]}):} {self.pctTPMs[isoform]},'
                # else:
                line += f'{isoform}: {self.pctTPMs[isoform]}\t'
                # line += f'{isoform}: {self.tpmMeans[isoform]},'

                # line += f'{self.genes[isoform]}'
            # line = line[:-1]
            line += '\n'
            file.append(line)
        with open(output, 'w') as outfile:
            outfile.writelines(file)

def grab_gtf_list(folder):
    gtf_list = ''
    files = os.listdir(folder)
    for file in files:
        gtf_list += f'{folder}/{file},'
    return gtf_list[:-1]

if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('\nExtracts the TPM values for each transcript after grouping them by isoforms. Assumes TPM is present in the attributes of a GTF file, just like StringTie outputs.\nusage: within_sample_isoform.py <gtf_folder> <output>')
    else:
        gtf_list = grab_gtf_list(sys.argv[1])
        data = WithinSampleComparison(gtf_list=gtf_list)
        data.parse_gtf()
        data.get_means()
        data.get_pcts()
        data.create_df(output=sys.argv[2])
