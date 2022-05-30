#!/usr/bin/env python3
"""usage: {} count_intron_exon.py -l <lasso> -b <bam> -o <outdir> -x <offsetx> -y <offsety>
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
    sargs = 'hl:b:o:x:y:'
    largs = ['help','bam=', 'outdir=', 'offsetx=','offsety=']
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
           # self.threads = 4
            if opt == '-h':
                print("count_intron_exon.py -l <lasso> -b <bam> -o <outdir> -x <offsetx> -y <offsety>")
                sys.exit()
            elif opt in ['-l', '--lasso']:
                self.lasso_file = arg
            elif opt in ['-b', '--bam']:
                self.bam_file = arg
            elif opt in ['-o', '--outdir']:
                self.out_dir = arg
            elif opt in ['-x', '--offsetx']:
                self.offsetx = int(arg)
            elif opt in ['-y', '--offsety']:
                self.offsety = int(arg)
            # elif opt in ['-t','--threads']:
            #     self.threads = int(arg)
    @property
    def usage(self):
        return __doc__.format(__file__)

def read_bgi_as_dataframe(path: str) -> pd.DataFrame:

    """Read a BGI read file as a pandas DataFrame.

    Args:
        path: Path to read file.

    Returns:
        Pandas Dataframe with column names `gene`, `x`, `y`, `total` and
        additionally `spliced` and `unspliced` if splicing counts are present.
    """
    return pd.read_csv(
        path,
        sep="\t",
        dtype={
            "geneID": "category",  # geneID
            "x": np.uint32,  # x
            "y": np.uint32,  # y
            "MIDCounts": np.uint16,  # total
        },
        comment="#",
    )



if __name__ == '__main__':
    start = time.time()
    opts = Options(sys.argv[1:])
    
    offsetx = opts.offsetx
    offsety = opts.offsety
    
    save_folder = opts.out_dir
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)
    
    bam_file = opts.bam_file
    lasso_file = opts.lasso_file
    lasso_id = os.path.basename(lasso_file)
    lasso_save = os.path.join(save_folder,f"intron_{lasso_id}")
    ## read bgi 
    df = read_bgi_as_dataframe(lasso_file)
    coords = df['x'].apply(str)+"_"+df['y'].apply(str)
    del df
    coords = set(list(coords))

    print("count======")
    sam = pysam.AlignmentFile(bam_file, 'rb')

    x_coords = []
    y_coords = []
    GE_list = []

    EXONIC_list = []
    INTRONIC_list = []
    MIDCounts_list = []
    #INTERGENIC_list = []
    #i = 0 
    all_df = pd.DataFrame(columns =["geneID","x","y","MIDCounts","EXONIC","INTRONIC"])
    for read in sam.fetch():
        # i = i+1
        # if i>=40000:
        #     break
        ### remove duplicate
        if not read.has_tag('GE'):
            continue

        x = read.get_tag('Cx') - offsetx
        y = read.get_tag('Cy') - offsety

        bin_tag = '{}_{}'.format(x, y)

        if bin_tag not in coords:
            continue
        x_coords.append(x)
        y_coords.append(y)
        GE = read.get_tag('GE')
        GE_list.append(GE)
        GE_last = GE

        XF = read.get_tag('XF')
        if XF == 0:
            EXONIC_list.append(1)
            INTRONIC_list.append(0)
        else:
            EXONIC_list.append(0)
            INTRONIC_list.append(1)
        MIDCounts_list.append(1)

        if len(x_coords)>=10000:
            if GE == GE_last:
                continue
            else:
                df = pd.DataFrame({"geneID":GE_list,"x":x_coords,"y":y_coords,"MIDCounts":MIDCounts_list,"EXONIC":EXONIC_list,"INTRONIC":INTRONIC_list})
                rdf = df.groupby(["geneID","x","y"],as_index =False).sum()
                all_df = pd.concat([all_df,rdf])

                x_coords = []
                y_coords = []
                GE_list = []

                EXONIC_list = []
                INTRONIC_list = []
                MIDCounts_list = []
    sam.close()
    ## count final
    if len(x_coords) > 0:
        df = pd.DataFrame({"geneID":GE_list,"x":x_coords,"y":y_coords,"MIDCounts":MIDCounts_list,"EXONIC":EXONIC_list,"INTRONIC":INTRONIC_list})
        rdf = df.groupby(["geneID","x","y"],as_index =False).sum()
        all_df = pd.concat([all_df,rdf])
    
    all_df.to_csv(lasso_save,index = False, sep ="\t")
    end = time.time()
    print("time:",end - start)