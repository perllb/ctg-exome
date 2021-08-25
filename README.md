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


## Input files

1. Samplesheet (`CTG_SampleSheet.exome.csv`)

- The samplesheet format is standard IEM generated sheet, with additional columns added after [Data]:

| Column | Supported values |
| ------ | -------- |
| Sample_Ref | hg38 / mm10 : hg38 and mm10 are currently set up for dragen |
| panel | comprehensive / core : Twist-Comprehensive or Twist-Core panel bed files are embedded in container |
| annotate | y / n : set 'y' for nirvana annotation |

- Also note that no Sample_Name should be added. Leave that column blank!

### Samplesheet template (.csv)

Samplesheet name: `CTG_SampleSheet.exome.csv`

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
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Index_Plate_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Sample_Ref,panel,annotate
34D,,,,D01,UDI0004,GAGACGAT,UDI0004,ACCGGTTA,2021_046,hg38,comprehensive,y                                                                      
49D,,,,E02,UDI0013,TACGCTAC,UDI0013,CGTGTGAT,2021_046,hg38,comprehensive,y                                                                      
P1017_3,,,,F12,UDI0094,AACCGTTC,UDI0094,CTAGCTCA,2021_038,mm10,core,y                                                                           
P1017_4,,,,G12,UDI0095,TGGTACAG,UDI0095,TCGAGAGT,2021_038,mm10,core,y   
```


## USAGE (manual run with nextflow)
Alternative is to run with driver (see below). 

1. Clone and build the Singularity container for this pipeline: (https://github.com/perllb/ctg-exome/tree/master/container). 
2. Edit the nextflow.config file to fit your project and system. Set directory from where you want to run nextflow (containing config and samplesheet) as `basedir`. (from where you execute the pipeline).
3. Set up `slide_area.txt` file and prepare slide .gpr file. See below for info.
4. Edit your samplesheet to match the example samplesheet
5. Run pipeline 
```
nohup nextflow run pipe-exome.nf > log.pipe-exome.txt &
```

## USAGE with driver 
For automated execution of pipeline and workflow.

- Must be started from within runfolder root directory.
- Needs:
 1. Runfolder (from where it is started)
 2. Samplesheet (with format as specified in ***Example Samplesheet*** below). If not specified, driver will take `runfolder/CTG_SampleSheet.csv` from runfolder.
 3. Genome reference on dragen. Must be in /staging/$species/reference/$ref. (Species is extracted from `ref` field in SampleSheet (hg* = human, mm* = mouse)
 4. Nirvana annotation. Must be set up before run, and path to nirvanadir must be added to config.
 5. Target bed files. Defined path in nextflow.config
 
```
Usage: exome-driver [ -i META_ID ] [ -s SAMPLESHEET ] [ -a INDEX-TYPE ] [ -b BCL2FASTQ-ARG ] [ -r RESUME ] [ -c CUSTOM-GENOME ] [ -t CUSTOM-TARGET ] [ -p PADDING ] [ -d DEMUX-OFF ] [ -h HELP ]                                                                                                        


    Optional arguments:                                                                                                                             
    META-ID    -i : Set 'meta-id' for runfolder (e.g. 210330-10x). Default: Takes date of runfolder + run ID in runfolder name and adds exome as suffix. E.g. '210330_A00681_0334_AHWFKTDMXX' becomes 210330_0334-exome                                                                                   
    SAMPLESHEET   -s : Set samplesheet used for run (Default: CTG_SampleSheet.csv)                                                                   
    INDEX-TYPE    -a : Set -a if single index uses. (Default: dual)                                                                     
    BCL2FASTQ-ARG -b : String with bcl2fastq argument. e.g. '--minimum-trimmed-read-length 20 --mask-short-adapter-reads 20' 
    RESUME        -r : Set if to resume nf-pipeline                          
    CUSTOM-GENOME -c : Path to custom reference genome if needed. Skip if human/mouse defined in samplesheet 
    CUSTOM-TARGET -t : Path to custom target bed file if TWIST Comprehensive or Core is not used. Skip if core/comprehensive is defined in samplesheet
    PADDING       -p : Set bp padding around target coordinates (default: 20) 
    DEMUX-OFF     -d : Set flag to skip mkfastq (then fastq must be in FQDIR)
    HELP          -h : print help message
```

***Run driver with default settings***
This requires the current files and directories to be in correct name and location:
- `CTG_SampleSheet.csv` in runfolder

```
cd runfolder 
exome-driver
```

***Run driver with single-index samples and bcl2fastq-arguments location***
```
cd runfolder 
exome-driver -a -b '--minimum-trimmed-read-length 20 --mask-short-adapter-reads 20' 
```

***Run driver with specific samplesheet (not CTG_SampleSheet.csv in runfolder)***
```
cd runfolder 
exome-driver -s /path/to/customsheet.csv
```

## Functions of exome-driver
1. Creates project folder, containing:
   - nextflow.config (copied from ctg-exome pipeline dir, and edited based on default driver params and user-specified parameters)
   - pipe-exome.nf (copied from ctg-exome pipeline dir)
   - samplesheet (copied from the one specified in driver)
2. Creates pipeline output directory
   - default is specified in driver script (/projets/fs1/nas-sync/ctg-delivery/exome/<metaid>)
3. Creates QC log output directory
   - in which qc output of pipeline is copied 
4. Starts pipe-exome.nf

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
- mm10
- Custom references can be added in config file (set "panel" column to "custom" in samplesheet, and add customgenome=/path/to/custom/reference in nextflow.config)


## Output:
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




