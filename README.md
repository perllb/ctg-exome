# ctg-exome 

Nextflow pipeline for exome analysis with Illumina Dragen server

## The following steps are performed by the pipeline:

* `Demultiplexing` (dragen bcl-conversion): Converts raw basecalls to fastq, and demultiplex samples based on index. Adapters are trimmed if added to samplesheet [Settings].
* `Alignment` and `Variant calling` (dragen align + calling): Fastq from each sample is aligned to the reference genome. Variants (SNV) are called in the target region. 
* `Dragen metrics`: Compiling Dragen alignment and coverage metrics to table.
* `FastQC`: FastQC calculates quality metrics on raw sequencing reads (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). 
* `MultiQC`: Summarizes FastQC and Dragen metrics into one document (https://multiqc.info/).

Currently supports the following panels
- `Twist Core Exome` (https://www.twistbioscience.com/resources/bed-file/ngs-human-core-exome-panel-bed-files)   
- `Twist Comprehensive Exome`: (https://www.twistbioscience.com/resources/bed-file/twist-human-comprehensive-exome-panel-bed-files)

Currently supports the following references
- hg38
- mm10
- Custom references can be added in config file (add "custom" to the "panel" column in samplesheet)


### Output:
* CTG-output
    * `fastq`: Contains raw fastq files from demultiplexing.
    * `dragen`: Output from dragen alignment + variant calling. Here you find the .bam files for the aligned sample and .vcf files. The *hard-filtered.vcf contain the variants that passed filters.
    * `qc`: Quality control output. 
        * fastqc output (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
        * dragen metrics: Summarized metrics for each sample.
        * multiqc output: Summarizing FastQC, dragen metrics and demultiplexing (https://multiqc.info/)

## Example Samplesheet
```
[Header]
IEMFileVersion,5
Date,2021-04-29
Workflow,GenerateFASTQ
Application,NovaSeq FASTQ Only
Instrument Type,NovaSeq
Assay,TWIST
"Index Adapters,""IDT- UD Indexes (96 Indexes)"""
Chemistry,Amplicon

[Reads]
151
151

[Settings]
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Index_Plate_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Sample_Species,panel,analyze
34D,,,,D01,UDI0004,GAGACGAT,UDI0004,ACCGGTTA,proj46,hg38,comprehensive,y
49D,,,,E02,UDI0013,TACGCTAC,UDI0013,CGTGTGAT,proj46,hg38,comprehensive,y
P1017_1,,,,C12,UDI0091,CCTTGTAG,UDI0091,CCTTGGAA,proj38,hg38,core,y
P1017_2,,,,D12,UDI0092,GAACATCG,UDI0092,TCGACAAG,proj38,hg38,core,y
```

Note that the format is standard IEM generated sheet, with additional columns:
| Column | Supported values |
| ------ | -------- |
| Sample_Species |  hg38 / mm10 / custom : hg38 and mm10 are embedded in singularity container - if custom, the nextflow.config needs to contain full path to custome reference |
| panel | comprehensive / core : Twist-Comprehensive or Twist-Core panel bed files are embedded in container |
| analyze | y / n : set 'y' if alignment + vc should be performed. If 'n' then only demux and QC |

## Container
- `ngs-tools` Singularity container contain NGS-related tools, embedded in the repo: 
https://github.com/perllb/ctg-wgs/tree/master/container 

## Specs
- Dragen version: 3.8
    - User guide: (https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/dragen-bio-it/Illumina-DRAGEN-Bio-IT-Platform-User-Guide-1000000141465-00.pdf). 
    - Relase notes: (https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/dragen/Illumina-DRAGEN-Bio-IT-Platform-3.7-Release-Notes-1000000142362-v00.pdf)
- Reference genome: hg38
- Target: Twist Core Exome (https://www.twistbioscience.com/resources/bed-file/ngs-human-core-exome-panel-bed-files)   
- Target: Twist Comprehensive Exome: (https://www.twistbioscience.com/resources/bed-file/twist-human-comprehensive-exome-panel-bed-files)
- Padding: 20bp



