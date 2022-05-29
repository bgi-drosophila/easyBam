# easyBam
This is the flow used to process bam files.

# 1. filter bam and add CB and UB tag 

* -b bam from stereo-seq web
* -o outputdir
* -x offset x 
* -y offset y
* -t threads for samtools (deafult 4)

``` bash
/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/xiangrong/software/miniconda3/envs/spateo0513/bin/python /hwfssz5/ST_SUPERCELLS/P21Z10200N0206/project/drosophila/03.ST/01.rawdata/03.bam_test/addtag.py -b /hwfssz5/ST_SUPERCELLS/P21Z10200N0206/project/drosophila/03.ST/01.rawdata/03.bam_test/SS200000101TR_A3_web_0.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam -o /hwfssz5/ST_SUPERCELLS/P21Z10200N0206/project/drosophila/03.ST/01.rawdata/04.bam_test -x 0 -y 8800
```

# 2. get every slice's barcodes

lasso coords is consistent with stereo-seq web

Example:

```bash
zcat SS200000101TR_A3.gem.gz|awk '{print($2"_"$3)}'|sort|uniq| >SS200000101TR_A3_S01.txt
```

# 3.velocyto 

```
Usage: velocyto run [OPTIONS] BAMFILE... GTFFILE

          Runs the velocity analysis outputting a loom file

          BAMFILE bam file with sorted reads

          GTFFILE genome annotation file

        Options:
          -b, --bcfile FILE               Valid barcodes file, to filter the bam. If --bcfile is not specified all the cell barcodes will be included.
                                          Cell barcodes should be specified in the bcfile as the `CB` tag for each read
          -o, --outputfolder PATH         Output folder, if it does not exist it will be created.
          -e, --sampleid PATH             The sample name that will be used to retrieve informations from metadatatable
          -s, --metadatatable FILE        Table containing metadata of the various samples (csv formatted, rows are samples and cols are entries)
          -m, --mask FILE                 .gtf file containing intervals to mask
          -c, --onefilepercell            If this flag is used every bamfile passed is interpreted as an independent cell, otherwise multiple files are interpreted as batch of different cells to be analyzed together.
                                          Important: cells reads should not be distributed over multiple bamfiles is not supported!!
                                          (default: off)
          -l, --logic TEXT                The logic to use for the filtering (default: Default)
          -U, --without-umi               If this flag is used the data is assumed UMI-less and reads are counted instead of molecules (default: off)
          -u, --umi-extension TEXT        In case UMI is too short to guarantee uniqueness (without information from the ampping) set this parameter to `chr`, `Gene` ro `[N]bp`
                                          If set to `chr` the mapping position (binned to 10Gb intervals) will be appended to `UB` (ideal for InDrops+dropEst). If set to
                                          `Gene` then the `GX` tag will be appended to the `UB` tag.
                                          If set to `[N]bp` the first N bases of the sequence will be used to extend `UB` (ideal for STRT). (Default: `no`)
          -M, --multimap                  Consider not unique mappings (not reccomended)
          -@, --samtools-threads INTEGER  The number of threads to use to sort the bam by cellID file using samtools
          --samtools-memory INTEGER       The number of MB used for every thread by samtools to sort the bam file
          -t, --dtype TEXT                The dtype of the loom file layers - if more than 6000 molecules/reads per gene per cell are expected set uint32 to avoid truncation (default run: uint32)
          -d, --dump TEXT                 For debugging purposes only: it will dump a molecular mapping report to hdf5. --dump N, saves a cell every N cells. If p is prepended a more complete (but huge) pickle report is printed (default: 0)
          -v, --verbose                   Set the vebosity level: -v (only warnings) -vv (warnings and info) -vvv (warnings, info and debug)
          --help                          Show this message and exit.
```

Example

```bash
/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/xiangrong/software/miniconda3/envs/spateo0513/bin/velocyto run /hwfssz5/ST_SUPERCELLS/P21Z10200N0206/project/drosophila/03.ST/01.rawdata/03.bam_test/SS200000101TR_A3_web_0.filtered.sort.addtag_dedup.sorted.bam /hwfssz5/ST_SUPERCELLS/P21Z10200N0206/xiangrong/database/Drosophila_melanogaster.BDGP6.28.101.gtf --outputfolder /hwfssz5/ST_SUPERCELLS/P21Z10200N0206/project/drosophila/03.ST/01.rawdata/03.bam_test/velocyto --samtools-threads 3 --samtools-memory 100000 -t int -u no -v -b /hwfssz5/ST_SUPERCELLS/P21Z10200N0206/project/drosophila/03.ST/01.rawdata/03.bam_test/SS200000101TR_A3_S01.txt -e SS200000101TR_A3_test2l
```




