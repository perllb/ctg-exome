# NGS-tools Singularity Container

Singularity container with NGS related tools.
It is mainly based on conda environment, and the `environment.yml` file will be used to create conda environment inside the container.

### Requirements

- bcl2fastq source code (bcl2fastq2-v2.20.0.422-Source.tar.gz)

### HOW-TO-Build

Make sure the `environment.yml` file is in the current dir when running the command

```
sudo -E singularity build ngs-tools.sif ngs-tools-build
```

### Included tools

- bcl2fastq (v2.20.0.422)
- fastqc (v0.11.9)
- illumina-interop (v1.1.21)
- multiqc (v1.10.1)
- numpy v(1.18.5)
- pandas v(1.0.5)
- python (v3.8.3)
- star (v2.7.1a)
- subread (v2.0.1) 
