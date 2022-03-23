import sys
import os


class NanoCountQuantifier(object):
    def __init__(self, folders):
        self.folders = folders

    def count(self):
        for folder in self.folders:
            files = os.listdir(folder)
            for file in files:
                if file.endswith(".bam"):
                    if not os.path.exists(f'{folder}/{file}.bai'):
                        os.system(f'samtools index {folder}/{file}')
                    outdir = f'{folder}/nanocount'
                    if not os.path.exists(outdir):
                        os.system(f'mkdir {outdir}')
                    cmd = f'NanoCount -i {folder}/{file} -o {outdir}/{file[:-4]}.tsv'
                    os.system(cmd)


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print("usage: nanocount.py <folder_containing_alignment> <folder_2> ... <folder_x>")
    else:
        args = sys.argv[1:]
        folders = []
        for arg in args:
            folders.append(arg)
        data = NanoCountQuantifier(folders=folders)
        data.count()