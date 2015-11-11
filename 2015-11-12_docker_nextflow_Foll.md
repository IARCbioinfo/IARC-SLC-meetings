# Writing reproducible and scalable bioinformatics pipelines using [nextflow](http://www.nextflow.io), [docker](https://www.docker.com) and [github](https://github.com)

Author: Matthieu Foll

Date: 12-11-2015

## Get some data to play with

In this repo I have some BAM files with TP53 from 1000 genome project with a fasta ref of chr17 and a bed files containing the coordinates of TP53. Let's download it:
```bash
$ git clone --depth=1 https://github.com/mfoll/NGS_data_test.git
$ cd NGS_data_test/1000G_CEU_TP53/
```

## The magic of all tools together

```bash
$ nextflow run iarcbioinfo/needlestack -with-docker iarcbioinfo/needlestack --bed TP53_all.bed --bam_folder BAM/ --fasta_ref 17.fasta.gz
N E X T F L O W  ~  version 0.16.1
Launching 'iarcbioinfo/needlestack' - revision: 35ee5283b4 [master]
--------------------------------------------------
NEEDLESTACK: A MULTI-SAMPLE SOMATIC VARIANT CALLER
--------------------------------------------------
Copyright (C) 2015 Matthieu Foll and Tiffany Delhomme
This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE.txt
This is free software, and you are welcome to redistribute it
under certain conditions; see LICENSE.txt for details.
--------------------------------------------------
Input BAM folder (--bam_folder)                                 : BAM/
Reference in fasta format (--fasta_ref)                         : 17.fasta.gz
Intervals for calling (--bed)                                   : TP53_all.bed
Number of regions to splot (--nsplit)                           : 1
To consider a site for calling:
     minimum coverage (--min_dp)                                : 50
     minimum of alternative reads (--min_ao)                    : 5
Phred-scale qvalue threshold (--min_qval)                       : 50
Strand bias measure (--sb_type)                                 : SOR
Strand bias threshold for SNVs (--sb_snv)                       : 100
Strand bias threshold for indels (--sb_indel)                   : 100
Samtools minimum mapping quality (--map_qual)                   : 20
Samtools minimum base quality (--base_qual)                     : 20
Samtools maximum coverage before downsampling (--max_DP)        : 30000
Sample names definition (--use_file_name)                       : BAM
Output all SNVs (--all_SNVs)                                    : no
PDF regression plots (--no_plots)                               : yes
Skip indels (--no_indels)                                       : no
output folder (--out_folder)                                    : BAM/


[warm up] executor > local
[ba/b25523] Submitted process > split_bed (1)
[52/f6085b] Submitted process > samtools_mpileup (17_7572814-17_7579962_regions)
[b0/072020] Submitted process > mpileup2table (17_7572814-17_7579962_regions)
[2f/81bcb4] Submitted process > R_regression (17_7572814-17_7579962_regions)
[9f/e437ca] Submitted process > collect_vcf_result (1)
```

## Docker

Install docker as explained [here](https://docs.docker.com/engine/installation/).



```bash
$ docker run -it --rm -v $PWD:$PWD -w $PWD --entrypoint /bin/bash iarcbioinfo/needlestack -c "samtools view -H NA11930.bam | head"
@HD	VN:1.0	SO:coordinate
@SQ	SN:1	LN:249250621	M5:1b22b98cdeb4a9304cb5d48026a85128	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz	AS:NCBI37	SP:Human
```

```bash
$ NGS_docker () {  docker run -it --rm -v $PWD:$PWD -w $PWD --entrypoint /bin/bash iarcbioinfo/needlestack -c "$@"; }
$ NGS_docker "samtools view -H NA11930.bam | head"
@HD	VN:1.0	SO:coordinate
@SQ	SN:1	LN:249250621	M5:1b22b98cdeb4a9304cb5d48026a85128	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz	AS:NCBI37	SP:Human
```

```bash
$ samtools_docker () { docker run -it --rm -v $PWD:$PWD -w $PWD --entrypoint /bin/bash iarcbioinfo/needlestack -c "samtools $@"; }
$ samtools_docker "view -H NA11930.bam | head"
@HD	VN:1.0	SO:coordinate
@SQ	SN:1	LN:249250621	M5:1b22b98cdeb4a9304cb5d48026a85128	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz	AS:NCBI37	SP:Human
```
