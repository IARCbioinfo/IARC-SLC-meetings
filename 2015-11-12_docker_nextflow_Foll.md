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

Installing docker is very system specific (but quite easy in most cases), follow  [docker documentation](https://docs.docker.com/installation/). Also follow the optional configuration step called `Create a Docker group` in their documentation.

Let's first download an image of Ubuntu latest version (from https://hub.docker.com):
```bash
$ docker run -it ubuntu
Unable to find image 'ubuntu:latest' locally
latest: Pulling from library/ubuntu
2332d8973c93: Pull complete 
ea358092da77: Pull complete 
a467a7c6794f: Pull complete 
ca4d7b1b9a51: Pull complete 
library/ubuntu:latest: The image you are pulling has been verified. Important: image verification is a tech preview feature and should not be relied on to provide security.
Digest: sha256:f91f9bab1fe6d0db0bfecc751d127a29d36e85483b1c68e69a246cf1df9b4251
Status: Downloaded newer image for ubuntu:latest
root@7c94b70daa8a:/# exit
```

Let's look at the images:
```bash
$ docker images
REPOSITORY                TAG                 IMAGE ID            CREATED             VIRTUAL SIZE
ubuntu                    latest              ca4d7b1b9a51        39 hours ago        187.9 MB
```

Be careful, our container (the instance of the image) is still here:
```bash
$ docker ps -a
CONTAINER ID        IMAGE               COMMAND             CREATED             STATUS                     PORTS               NAMES
7c94b70daa8a        ubuntu              "/bin/bash"         2 minutes ago       Exited (0) 2 minutes ago                       jovial_hypatia
```

You can kill it with `docker rm 7c94b70daa8a`. You can also kill all running containers using `docker rm $(docker ps -a -q)`. Similarly you can delete our ubuntu image using `docker rmi ca4d7b1b9a51` and all images using `docker rmi $(docker images -q)`.

Now let's install something ([samtools](http://www.htslib.org)) in our container:
```bash
$ docker run -it ubuntu
$ apt-get update
$ apt-get install samtools
$ samtools

Program: samtools (Tools for alignments in the SAM format)
Version: 0.1.19-96b5f2294a

Usage:   samtools <command> [options]
$ exit
```

Be careful, now if you again type `docker run -it ubuntu`, samtools won't be here: the image remained the same when we installed samtools, only the container changed, but now you stopped it when you exited. You can start it again:
```bash
$ docker ps -a
CONTAINER ID        IMAGE               COMMAND             CREATED             STATUS                          PORTS               NAMES
e6f92d19b512        ubuntu              "/bin/bash"         3 minutes ago       Exited (1) About a minute ago                       mad_goodall
$ docker start -ia e6f92d19b512
$ samtools

Program: samtools (Tools for alignments in the SAM format)
Version: 0.1.19-96b5f2294a

Usage:   samtools <command> [options]
$ exit
```

Now if we want to keep our container as an image (to be able to re-use it), you can do it simply with:
```bash
$ docker commit -m "added samtools" e6f92d19b512 samtools_img
e664fa7478a7a74addbef19379412a901bc05c8164d21c2de58e62a32e0117c3
$ docker images
REPOSITORY                TAG                 IMAGE ID            CREATED             VIRTUAL SIZE
samtools_img              latest              e664fa7478a7        4 seconds ago       211.5 MB
ubuntu                    latest              ca4d7b1b9a51        40 hours ago        187.9 MB
```

Now that we saved it as an image, a good practice is to delete containers when we are done using them. You can simply do this by adding `--rm` when you create a container:
```bash
$ docker rm $(docker ps -a -q)
$ docker run --rm -it samtools_img
$ samtools

Program: samtools (Tools for alignments in the SAM format)
Version: 0.1.19-96b5f2294a

Usage:   samtools <command> [options]
$ exit
$ docker ps -a
CONTAINER ID        IMAGE               COMMAND             CREATED             STATUS              PORTS               NAMES
```

Now the final problem is that the docker container is disconnected from our machine, so we have no way of accessing any data... We can "mount" a directory from our machine in the container using `-v` option and also automatically move to this directory when the container starts. For example:
```bash
$ docker run -it --rm -v /Users/follm/DEMO_IARC_SLC/NGS_data_test/1000G_CEU_TP53/:/data -w /data samtools_img
$ samtools view -H BAM/NA11830.bam | head
@HD	VN:1.0	SO:coordinate
@SQ	SN:1	LN:249250621	M5:1b22b98cdeb4a9304cb5d48026a85128	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz	AS:NCBI37	SP:Human
$ exit
```

Note `-v` option only accepts absolute path, so a good trick is to systematically mount your current working directory using `$PWD` into the same target directory. Note that you can have multiple `-v` options when you run a container. So it makes something like this:
```bash
$ docker run -it --rm -v $PWD:$PWD -w $PWD samtools_img
$ samtools view -H BAM/NA11830.bam | head
```

So far we always run interactively docker containers, but you can actually run it automatically. Still the same example becomes:
```bash
$ docker run -it --rm -v $PWD:$PWD -w $PWD --entrypoint /bin/bash samtools_img -c "samtools view -H BAM/NA11830.bam | head"
@HD	VN:1.0	SO:coordinate
@SQ	SN:1	LN:249250621	M5:1b22b98cdeb4a9304cb5d48026a85128	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz	AS:NCBI37	SP:Human
```

One more trick tu make the command simpler to type by creating a bash function:
```bash
$ samtools_docker () { docker run -it --rm -v $PWD:$PWD -w $PWD --entrypoint /bin/bash samtools_img -c "samtools $@"; }
$ samtools_docker "view -H BAM/NA11830.bam | head"
@HD	VN:1.0	SO:coordinate
@SQ	SN:1	LN:249250621	M5:1b22b98cdeb4a9304cb5d48026a85128	UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz	AS:NCBI37	SP:Human
```

And voilà! You can save your image using `docker save 'samtools_img' > samtools_img.tar` and distribute it (check software licences...). You can also host it freely on [Docker Hub](https://hub.docker.com). You can also make the process automatic by using a [Dockerfile](http://docs.docker.com/engine/reference/builder/). In our case it looks like:
```
FROM ubuntu
MAINTAINER Matthieu Foll <follm@iarc.fr>
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install samtools
```

If the file is called Dockerfile you simply need to type: `docker build -t samtools_img2 .` to create an image called `samtools_img2`. The Dockerfile is easier to distribute than the actual image. Note that if you host your Dockerfile on [github](https://github.com), [Docker Hub](https://hub.docker.com) can make automatically build it whenever you modify it. You can browse [Docker Hub](https://hub.docker.com) and you will see that you can find three different classes of images:
- Base images, just containing a given Linux distributio (centos, debian, ubuntu etc.)
- Images to distribute a specific software ([R](https://hub.docker.com/_/r-base/) for example)
- Images bundling several software for a given task ([NGSeasy](https://hub.docker.com/r/compbio/ngseasy-base/) for example)

## Nextflow

Install [nextflow](http://www.nextflow.io/) (you will need [java](https://java.com/download/) JRE if you don't already have it):
```bash
$ curl -fsSL get.nextflow.io | bash
```

This is a non-working simplified version of `needlestack` to understand the logic:
```java
#! /usr/bin/env nextflow

bed = file( params.bed )
fasta_ref = file( params.fasta_ref )

bam = Channel.fromPath( params.bam_folder+'/*.bam' ).toList()   
bai = Channel.fromPath( params.bam_folder+'/*.bam.bai' ).toList()

params.nsplit = 1

params.out_folder = params.bam_folder

process split_bed {
	  
	intput:
	file bed
	   
	output:
	file '*_regions' into split_bed mode flatten 

	shell:
	'''
	bedtools merge -i !{bed} | bed_large_cut.r !{params.nsplit}
	'''
}

process samtools_mpileup {
		
	input:
	file split_bed
	file bam
	file bai  
	file fasta_ref
	
	output:
	file pileup
	   
 	shell:
	'''
	samtools mpileup --fasta-ref !{fasta_ref} --region !{split_bed} !{bam} > pileup
	'''
}

process mpileup2table {
	
	input:
	file pileup
	
	output:
	file 'sample*.txt' into table
	   
 	shell:
 	'''
	pileup2baseindel.pl -i !{pileup}
	'''
}

process R_regression {
		
	input:
	file table
	file fasta_ref
	
	output:
	file vcf
	file '*.pdf' into PDF
	   
 	shell:
 	'''
	needlestack.r --out_file=vcf --fasta_ref=!{fasta_ref} 
	'''
}

process collect_vcf_result {

	storeDir { params.out_folder }

	input:
	file '*.vcf' from vcf.toList()
  
	output:
	file 'all_variants.vcf' 

	shell:
	vcfoverlay *.vcf > all_variants.vcf
	'''
}
```

A good practice is to keep (and publish) the `.nextflow.log` file create during the pipeline process, as it contains useful information for reproducibility (full command line, software versions etc.). You should also add the option `-with-trace` in the `nextflow run` command line that will create an additional `trace.csv` file containing even more information to keep for records. The option `-with-timeline` also creates a nice processes execution timeline file (a web page). You can easily create aliases to avoid having to always type long command lines. For example with `needlestack`:
```bash
alias needlestack='nextflow run iarcbioinfo/needlestack -with-docker iarcbioinfo/needlestack -latest -with-trace --with-timeline'
```

It will allow you to do this:
```bash
needlestack --bed TP53_all.bed --bam_folder BAM/ --fasta_ref 17.fasta.gz
```
