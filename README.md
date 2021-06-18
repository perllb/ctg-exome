# ctg-exome 

Nextflow pipeline for exome analysis with Illumina Dragen server

## The following steps are performed by the pipeline:

* `Demultiplexing` (dragen bcl-conversion): Converts raw basecalls to fastq, and demultiplex samples based on index. Adapters are trimmed if added to samplesheet [Settings].
* `Alignment` (dragen): Map and align fastqs to reference genome
* `Variant calling` (dragen): Variants (SNV and SV) are called in the target region (specified with bed). CNV not yet available - will be available with Dragen 3.9, but will require panel of normal.  
* `Annotation` (nirvana): Clinical grade annotation of all variants passing basic filtering in Dragen. 
* `Dragen metrics`: Compiling Dragen alignment and coverage metrics to table.
* `FastQC`: FastQC calculates quality metrics on raw sequencing reads (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). 
* `MultiQC`: Summarizes FastQC and Dragen metrics into one document (https://multiqc.info/).

**Currently supports the following panels**
- `Twist Core Exome` (https://www.twistbioscience.com/resources/bed-file/ngs-human-core-exome-panel-bed-files)   
- `Twist Comprehensive Exome`: (https://www.twistbioscience.com/resources/bed-file/twist-human-comprehensive-exome-panel-bed-files)

**Padding**
Default padding is 20bp. Can be altered in nextflow.config.

**Filters for Coverage Reports**
By default, the following filters are applied for coverage reports:
- `--remove-duplicates true`: Duplicate reads are ignored
- `--qc-coverage-ignore-overlaps true`: Resolve all of the alignments for each fragment and avoid double-counting any overlapping bases
- `--qc-coverage-filters-1 "mapq<20,bq<20"`: Reads with MAPQ<20 and BQ<20 are ignored

**Currently supports the following references**
- hg38
- hg37
- mm10
- Custom references can be added in config file (set "panel" column to "custom" in samplesheet, and add customgenome=/path/to/custom/reference in nextflow.config)


### Output:
* CTG-output
    * `fastq`: Contains raw fastq files from demultiplexing.
    * `dragen`: 
      * `metrics`: Output from dragen alignment + variant calling. Contains metrics from alignment and calling, raw .vcfs and .bam files. 
         * `sv`: Output metrics and candidate vcf from sv calling.       
      * `vcf`: 
         * `*hard-filtered.vcf` contain the SNV variants (and some SVs) that passed filters. 
         * `*DiploidSV.vcf` contain SVs that passes the candidate filters. 
         * Annotated and filtered vcf are also written here.
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
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Index_Plate_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Sample_Species,Sample_Ref,panel,annotate
34D,,,,D01,UDI0004,GAGACGAT,UDI0004,ACCGGTTA,2021_046,human,hg38,comprehensive,y                                                                      
49D,,,,E02,UDI0013,TACGCTAC,UDI0013,CGTGTGAT,2021_046,human,hg38,comprehensive,y                                                                      
P1017_3,,,,F12,UDI0094,AACCGTTC,UDI0094,CTAGCTCA,2021_038,human,hg38,core,y                                                                           
P1017_4,,,,G12,UDI0095,TGGTACAG,UDI0095,TCGAGAGT,2021_038,human,hg38,core,y   
```

Note that the format is standard IEM generated sheet, with additional columns:
| Column | Supported values |
| ------ | -------- |
| Sample_Species | human / mouse |
| Sample_Ref | hg37 / hg38 / mm10 : hg37, hg38 and mm10 are currently set up on dragen |
| panel | comprehensive / core : Twist-Comprehensive or Twist-Core panel bed files are embedded in container |
| annotate | y / n : set 'y' for nirvana annotation |

## Container
- `ngs-tools` Singularity container contain NGS-related tools, embedded in the repo: 
https://github.com/perllb/ctg-wgs/tree/master/container 

## References
- Dragen version: 3.8
    - User guide: (https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/dragen-bio-it/Illumina-DRAGEN-Bio-IT-Platform-User-Guide-1000000141465-00.pdf). 
    - Relase notes: (https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/dragen/Illumina-DRAGEN-Bio-IT-Platform-3.7-Release-Notes-1000000142362-v00.pdf)
- Nirvana annotator: https://illumina.github.io/NirvanaDocumentation/
- Twist Core Exome (https://www.twistbioscience.com/resources/bed-file/ngs-human-core-exome-panel-bed-files)   
- Twist Comprehensive Exome: (https://www.twistbioscience.com/resources/bed-file/twist-human-comprehensive-exome-panel-bed-files)




