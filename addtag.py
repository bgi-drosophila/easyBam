#!/usr/bin/env python3
"""usage: {} -b <bam> -o <outdir> -x <offsetx> -y <offsety> -t <threads>
"""

import os
import sys
import getopt
import pysam
import numpy as np
import pandas as pd
import time
import subprocess

class Options:
    sargs = 'hl:b:o:x:y:t:'
    largs = ['help','bam=', 'outdir=', 'offsetx=','offsety=','threads=']
    def __init__(self, cmdline):
        self.cmdline = cmdline
        self.out_file = sys.stdout
        self.keep_tag = False
        try:
            self.options, self.not_options = getopt.gnu_getopt(self.cmdline, 
                    self.sargs, self.largs)
        except getopt.GetoptError as err:
            sys.stderr.write('*** Error: {}\n'.format(err))
            sys.exit()
        if self.not_options or not self.options:
            sys.stderr.write(self.usage)
            sys.exit()
        for opt, arg in self.options:
            self.threads = 4
            if opt == '-h':
                print("unsplice.py -b <bam> -o <outdir> -x <offsetx> -y <offsety>-t <threads>")
                sys.exit()
            # elif opt in ['-l', '--lasso']:
            #     self.lasso_file = arg
            elif opt in ['-b', '--bam']:
                self.bam_file = arg
            elif opt in ['-o', '--outdir']:
                self.out_dir = arg
            elif opt in ['-x', '--offsetx']:
                self.offsetx = int(arg)
            elif opt in ['-y', '--offsety']:
                self.offsety = int(arg)
            elif opt in ['-t','--threads']:
                self.threads = int(arg)
    @property
    def usage(self):
        return __doc__.format(__file__)

# def read_bgi_as_dataframe(path: str) -> pd.DataFrame:
#     """Read a BGI read file as a pandas DataFrame.

#     Args:
#         path: Path to read file.

#     Returns:
#         Pandas Dataframe with column names `gene`, `x`, `y`, `total` and
#         additionally `spliced` and `unspliced` if splicing counts are present.
#     """
#     return pd.read_csv(
#         path,
#         sep="\t",
#         dtype={
#             "geneID": "category",  # geneID
#             "x": np.uint32,  # x
#             "y": np.uint32,  # y
#             "MIDCounts": np.uint16,  # total
#         },
#         comment="#",
#     )

if __name__ == '__main__':
    start = time.time()
    opts = Options(sys.argv[1:])
    
    save_folder = opts.out_dir
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)
    
    threads = int(opts.threads)

    ibam_file = opts.bam_file
    bamfile = os.path.basename(opts.bam_file)
    chip_id = bamfile[:-53]
    # lasso_file = opts.lasso_file
    # lasso_id = os.path.basename(lasso_file)

    bamfile_filter = os.path.join(save_folder, f"{chip_id}.filter.bam")
    bamfile_sort = os.path.join(save_folder, f"{chip_id}.filter.sorted.bam")
    bamfile_addtag = os.path.join(save_folder, f"{chip_id}.filter.sorted.addtag.bam")
   # coords_file = os.path.join(save_folder, f"{lasso_id[:-7]}.coords_file.csv")

    # df = read_bgi_as_dataframe(opts.lasso_file)
    # coords = df['x'].apply(str)+"_"+df['y'].apply(str)
    # del df
    # coords = set(list(coords))
    #coords.to_csv(coords_file)
    


    ## filter
    print("filter======")
    if not os.path.exists(bamfile_filter):
        subprocess.call(f'/share/app/samtools/1.11/bin/samtools view -h -@ {threads} -b -q 30 -F 4 -F 256 -F 1024 {ibam_file} >{bamfile_filter}', shell=True)
    time.sleep(10)
    ## sort
    print("sort======")
    if not os.path.exists(bamfile_sort):
        subprocess.call(f'/share/app/sambamba/0.7.1/sambamba sort {bamfile_filter} -t {threads} && rm {bamfile_filter}', shell=True)
    time.sleep(10)
    ## read coords
    sam = pysam.AlignmentFile(bamfile_sort, 'rb')
    out_bam = pysam.AlignmentFile(bamfile_addtag, 'wb', template=sam)
    
    
    offsetx = opts.offsetx
    offsety = opts.offsety
 
 
    bin_tag_UB_list = []
    print("addtag======")
    if not os.path.exists(bamfile_addtag):
        for read in sam.fetch():
            ### remove duplicate
            if not read.has_tag('GE'):
                continue
            x = read.get_tag('Cx') - offsetx
            y = read.get_tag('Cy') - offsety

            bin_tag = '{}_{}'.format(x, y)
            # if bin_tag not in coords:
            #     continue

            if not read.has_tag('UB'):
                UR = read.get_tag('UR')
                read.set_tag('UB',UR)
    
            # bin_tag_UB = bin_tag + read.get_tag('UB')

            # ## remove duplicate
            # if bin_tag_UB in bin_tag_UB_list:
            #         continue			
            # bin_tag_UB_list.append(bin_tag_UB)

            # if len(bin_tag_UB_list)==500:
            #     bin_tag_UB_list = bin_tag_UB_list[-300:]
            read.set_tag('CB', bin_tag)
            out_bam.write(read)
        sam.close()
        out_bam.close()
    subprocess.call(f'rm {bamfile_sort}')
    # bamfile_sort_final = os.path.join(save_folder,f"cellsorted_{chip_id}.filter.sorted.addtag.bam")
    # if not os.path.exists(bamfile_sort_final):
    #     subprocess.call(f'/share/app/samtools/1.11/bin/samtools sort -l 7 -m 100000M -t CB -O BAM -@ {threads} -o {bamfile_sort_final} bamfile_addtag && rm {bamfile_addtag}', shell=True)
    
    # time.sleep(10)
    
    #subprocess.call(f'source /home/xiangrong1/.bashrc && /hwfssz5/ST_SUPERCELLS/P21Z10200N0206/xiangrong/software/miniconda3/envs/spateo0513/bin/velocyto run {bamfile_sort_final} /hwfssz5/ST_SUPERCELLS/P21Z10200N0206/xiangrong/database/Drosophila_melanogaster.BDGP6.28.101.gtf --outputfolder {save_folder} --samtools-threads {threads} -M --samtools-memory 100000 -t int -U', shell=True)
    end = time.time()
    print("time:",end - start)