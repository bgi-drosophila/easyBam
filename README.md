# easyBam

This is the flow used to process bam files.



# 1. filter and sort bam 

* filter duplicate

  * input_bam_file from stereo-seq web
  * bamfile_filter output

  ```bash
  /share/app/samtools/1.11/bin/samtools view -h -@ {threads} -b -q 30 -F 4 -F 256 -F 1024 {input_bam_file} >{bamfile_filter}
  ```

* sort 

  ```bash
  /share/app/sambamba/0.7.1/sambamba sort {bamfile_filter} -t {threads} && rm {bamfile_filter}
  ```

# 2. count

* -l lasso file,lasso coords is consistent with stereo-seq web
* -b bam from first step(sorted)
* -o outputdir
* -x offset x (from stereo-seq web)
* -y offset y (from stereo-seq web)

``` bash
/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/xiangrong/software/miniconda3/envs/spateo0513/bin/python count_intron_exon.py -l /hwfssz5/ST_SUPERCELLS/P21Z10200N0206/project/drosophila/03.ST/01.rawdata/03.bam_test/SS200000101TR_A3.gem.gz -b /hwfssz5/ST_SUPERCELLS/P21Z10200N0206/project/drosophila/03.ST/01.rawdata/03.bam_test/SS200000101TR_A3_web_0.filtered.sorted.bam -o /hwfssz5/ST_SUPERCELLS/P21Z10200N0206/project/drosophila/03.ST/01.rawdata/04.bam_test -x 0 -y 8800
```

# 
