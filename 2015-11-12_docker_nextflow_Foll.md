# Writing reproducible and scalable bioinformatics pipelines using [nextflow](http://www.nextflow.io), [docker](https://www.docker.com) and [github](https://github.com)

Author: Matthieu Foll

Date: 12-11-2015

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
