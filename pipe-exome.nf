#!/usr/bin/env nextFlow

// set variables
runfolder = params.runfolder
basedir = params.basedir
metaID = params.metaid
OUTDIR = params.outdir
FQDIR = params.fqdir
QCDIR = params.qcdir
FQCDIR = params.fqcdir
CTGQC = params.ctgqc
demux = params.demux
b2farg = params.bcl2fastq_arg
index = params.index
nirvanadir = params.nirvanadir

// Read and process sample sheet
sheet = file(params.sheet)

// samplesheet to be parsed as input channel (take everything below [Data])
channel_sheet = file("$basedir/samplesheet.channel.nf.sc-rna-10x.csv")

// create new samplesheet parsed to fit the format for dragen demux
newsheet = "${basedir}/samplesheet.demux.nf.exome.csv"

// Read and process sample sheet
all_lines = sheet.readLines()
write_b = false // if next lines has sample info
channel_sheet.text=""     

for ( line in all_lines ) {

    if ( write_b ) {
	channel_sheet.append(line + "\n")
    }
    if (line.contains("[Data]")) {
	write_b = true
    }
}


println "============================="
println ">>> exome dragen pipeline "
println ""
println "> INPUT: "
println "> runfolder		: $runfolder "
println "> sample-sheet		: $sheet "
println "> run-meta-id		: $metaID "
println "> basedir		: $basedir "
println "> bcl2fastq args	: $b2farg "
println ""
println "> OUTPUT: "
println "> output-dir		: $OUTDIR "
println "> fastq-dir		: $FQDIR "
println "> qc-dir		: $QCDIR "
println "> fastqc-dir		: $FQCDIR "
println "> ctg-qc-dir		: $CTGQC "
println "============================="

// sample info
Channel
    .fromPath(channel_sheet)
    .splitCsv(header:true)
    .map { row -> tuple( row.Sample_ID, row.Sample_Project, row.Sample_Ref, row.panel, row.annotate ) }
    .unique()
    .tap{infoSamples}
    .into{ move_fastq_csv; analyze_csv  }

// project info
Channel
    .fromPath(channel_sheet)
    .splitCsv(header:true)
    .map { row -> tuple(  row.Sample_Project, row.panel ) }
    .unique()
    .tap{infoProjects}
    .set{ dragen_summary  }

println " > Samples to process: "
infoSamples.subscribe{ println "Sample: $it" }


println " > Projects to process: "
infoProjects.subscribe{ println "Sample: $it" }

// Parse samplesheet
process parsesheet {

	tag "$metaID"

	input:
	val sheet
	val index

	output:
	val newsheet into demux_sheet

	when:
	demux == 'y'

	"""
python $basedir/bin/ctg-parse-samplesheet.dragen-exome.py -s $sheet -o $newsheet -i $index
	"""
}

// dragen demux
process demux {

    tag "$metaID"	
    label 'dragen'

    input:
    val newsheet from demux_sheet

    output:
    val "x" into mv_fastq

    when:
    demux = "y"
        
    """
    export LC_ALL=C        

    mkdir -p ${FQDIR}
    /opt/edico/bin/dragen --force --bcl-conversion-only=true \\
           --bcl-input-directory ${runfolder} \\
	   --output-directory ${FQDIR} \\
	   --sample-sheet ${newsheet} \\
	   --no-lane-splitting true \\
	   ${b2farg}

     """
}

process moveFastq {

    tag "${sid}_${projid}"

    input:
    val x from mv_fastq
    set sid, projid, ref, panel, annotate from move_fastq_csv

    output:
    val "y" into run_analysis
    set sid, projid into fastqc_go

    when:
    demux = 'y'

    """
    mkdir -p ${OUTDIR}/${projid}
    mkdir -p ${OUTDIR}/${projid}/fastq

    mkdir -p ${OUTDIR}/${projid}/fastq/$sid

    # If there is a directory per project
    if [ -d ${FQDIR}/${projid}/ ]; then
        if [ -d ${FQDIR}/${projid}/$sid ]; then
	    mv ${FQDIR}/${projid}/$sid ${OUTDIR}/${projid}/fastq/
	# If there is a flat structure under project dir
	else
	    mv ${FQDIR}/${projid}/$sid* ${OUTDIR}/${projid}/fastq/$sid/
	fi
    # If there is a flat structure with all samples for all projects in one - create a sid folder for each sample
    else
	mv ${FQDIR}/$sid* ${OUTDIR}/${projid}/fastq/$sid/
    fi
    """

}

// Channel to start analysis if demux == 'n'
// Projects
if ( demux == 'n' ) {

   Channel
	 .from("1")
    	 .set{ run_analysis }

}

// dragen run : align, vc + metrics
process dragen_align_vc {

    tag "${sid}_${projid}"
    label 'dragen' 

    input:
    val x from run_analysis
    set sid, projid, ref, panel, annotate from analyze_csv

    output:
    val x into done_analyze
    val projid into dragen_metrics
    set sid, projid, ref, annotate, val("${OUTDIR}/${projid}/dragen/vcf/${sid}/${sid}.hard-filtered.vcf.gz"), val("${OUTDIR}/${projid}/dragen/vcf/${sid}/${sid}.diploidSV.vcf.gz") into annotate_vcf
    
    """
    export LC_ALL=C

    R1=\$(echo ${OUTDIR}/${projid}/fastq/${sid}/${sid}*_R1_*fastq.gz)
    R2=\$(echo ${OUTDIR}/${projid}/fastq/${sid}/${sid}*_R2_*fastq.gz)

    # Get target panel file
    if [ $panel == "comprehensive" ] || [ $panel == "Comprehensive" ]
    then
	if [ $ref == 'hg38' ]; then
	    targetfile=${params.target_twist_comprehensive_hg38}
	elif [ $ref == 'mm10' ]; then
	    targetfile=${params.target_twist_comprehensive_mm10}
        fi
    elif [ $panel == "core" ] || [ $panel == "Core" ]
    then
        if [ $ref == 'hg38' ]; then
	    targetfile=${params.target_twist_core_hg38}
	elif [ $ref == 'mm10' ]; then
	    targetfile=${params.target_twist_core_mm10}
        fi
    elif [ $panel == "custom"  ] || [ $panel == "Custom" ] 
    then
        targetfile=${params.custom_target}
    else
        echo '>PANEL NOT RECOGNIZED!'
	echo 'in samplesheet - only 'comprehensive', 'core' and 'custom' can be specified in 'panel' section'
        targetfile='ERR'
    fi

    # Get species based on ref
    if [[ $ref == hg* ]]
    then 
    	species='human'
    elif [[ $ref == mm* ]]
    then
	species='mouse'
    else
        echo 'species cannot be extracted from reference: $ref'
	exit 1;
    fi

    dragendir=${OUTDIR}/${projid}/dragen/
    mkdir -p \$dragendir
    mkdir -p \${dragendir}/metrics
    mkdir -p \${dragendir}/vcf/
    mkdir -p \${dragendir}/vcf/${sid}

    outdir=\${dragendir}/metrics/${sid}
    mkdir -p \$outdir

    echo "R1: '\${R1}'"
    echo "R2: '\${R2}'"
    echo "sid: '${sid}'"
    echo "padding: '${params.padding}'"
    echo "outdir: '\${outdir}'"
    echo "targetfile: '\${targetfile}'"
    echo "species: '\${species}'"   

    /opt/edico/bin/dragen -f -r /staging/\$species/reference/$ref \\
        -1 \${R1} \\
        -2 \${R2} \\
        --RGID ${projid}_${sid} \\
        --RGSM $sid \\
        --intermediate-results-dir /staging/tmp/ \\
        --enable-map-align true \\
        --enable-map-align-output true \\
        --vc-target-bed \$targetfile \\
        --vc-target-bed-padding ${params.padding} \\
        --output-format bam \\
        --output-directory \$outdir \\
        --enable-variant-caller true \\
        --enable-sv true \\
        --output-file-prefix $sid \\
        --remove-duplicates true \\
        --qc-coverage-region-1 \$targetfile \\
	--qc-coverage-region-padding-1 ${params.padding} \\
        --qc-coverage-ignore-overlaps true \\
        --qc-coverage-filters-1 "mapq<20,bq<20"
	
    # move vcfs to vcf dir     
    # SNV
    mv \${outdir}/${sid}.hard-filtered.vcf.gz \${dragendir}/vcf/${sid}/
    # SV
    mv \${outdir}/sv/results/variants/diploidSV.vcf.gz \${dragendir}/vcf/${sid}/${sid}.diploidSV.vcf.gz
    """
}

// Annotate vcf
process annotate {

	tag "${sid}"
	label 'dragen' 	

	input:
	set sid, projid, ref, annotate, snv, sv from annotate_vcf
	
	output:
	set sid, projid, ref, annotate into filter_vcf

	when:
	annotate == 'y'

	"""
	nirvanaref=$nirvanadir/$ref/

	vcfdir='${OUTDIR}/${projid}/dragen/vcf/$sid'

	outSNV_BN=\$(basename $snv .vcf.gz)
	outSV_BN=\$(basename $sv .vcf.gz)
	
	outSNV=\$vcfdir/\${outSNV_BN}.annotated.nirvana.txt
	outSV=\$vcfdir/\${outSV_BN}.annotated.nirvana.txt

	# Convert hg37/38 to GRCh37/8 
	Gref=\$(echo $ref | sed 's/hg/GRCh/')	

	# SNV annotate
	/opt/edico/share/nirvana/Nirvana -c \$nirvanaref/Cache/\$Gref/Both -r \$nirvanaref/References/Homo_sapiens.\$Gref.Nirvana.dat --sd \$nirvanaref/SupplementaryAnnotation/\$Gref -i $snv -o \$outSNV

	# SNV annotate		 
	/opt/edico/share/nirvana/Nirvana -c \$nirvanaref/Cache/\$Gref/Both -r \$nirvanaref/References/Homo_sapiens.\$Gref.Nirvana.dat --sd \$nirvanaref/SupplementaryAnnotation/\$Gref -i $sv -o \$outSV

	"""
}


// fastqc 
process fastqc {

	tag "${sid}_${projid}"
	
	input:
	set sid, projid from fastqc_go

        output:
        val projid into multiqc_fastqc

	
	"""

        qcdir=${OUTDIR}/${projid}/qc
        fqcdir=${OUTDIR}/${projid}/qc/fastqc
        mkdir -p \$qcdir
        mkdir -p \$fqcdir 

	read1=\$(echo ${OUTDIR}/$projid/fastq/$sid/$sid*R1*fastq.gz)
   	read2=\$(echo ${OUTDIR}/$projid/fastq/$sid/$sid*R2*fastq.gz)

    	# Check if fastq is not found with pattern above (due to sample fastq are not put in sample id folder
    	if [[ \$read1 == *"*R1*"* ]]; then
       	   read1=\$(echo ${FQDIR}/${projid}/${sid}/${sid}*R1*fastq.gz)
       	   read2=\$(echo ${FQDIR}/${projid}/${sid}/${sid}*R2*fastq.gz)
    	fi

	fastqc -t $task.cpus --outdir \$fqcdir \$read1
	fastqc -t $task.cpus --outdir \$fqcdir \$read2
 	
	"""
    
}

process dragen_stats {

        tag "${projid}"

	input: 
	val "y" from dragen_metrics.collect()
	set projid, panel  from dragen_summary.unique()
	
	output:
	val projid into multiqc_dragen
	val "x" into completed

	"""

    if [ $panel == "comprehensive" ] || [ $panel == "Comprehensive" ]
    then
        targetfile='TWIST-comprehensive'
    elif [ $panel == "core" ] || [ $panel == "Core" ]
    then
        targetfile='TWIST-core'
    else
        targetfile='Custom panel'
    fi

	mkdir -p ${OUTDIR}/$projid/qc/dragen

	${basedir}/bin/ctg-dragen-stats-panel -p $projid -i ${OUTDIR}/$projid/dragen/ -o ${OUTDIR}/$projid/qc/dragen/ -t \$targetfile -a ${params.padding}
	
	"""
}



process multiqc {

    tag "${projid}"

    input:
    set projid, projid2 from multiqc_fastqc.unique().phase(multiqc_dragen.unique())


    """
    
    cd $OUTDIR
    multiqc -f ${OUTDIR}/$projid/  --outdir ${OUTDIR}/$projid/qc/multiqc/ -n ${projid}_exome_dragen_report.html

    mkdir -p ${CTGQC}
    mkdir -p ${CTGQC}/$projid

    cp -r ${OUTDIR}/$projid/qc ${CTGQC}/$projid/

    """

}


process multiqc_run {

	tag "$metaID"

	input:
	val x from completed.collect()

	"""

	cd $OUTDIR
	mkdir -p ${OUTDIR}/qc

	multiqc -f ${OUTDIR} $params.runfolder/ctg-interop --outdir ${OUTDIR}/qc/ -n exome_dragen_run_${metaID}_multiqc_report.html

	mkdir -p ${CTGQC}

	cp -r ${OUTDIR}/qc ${CTGQC}
	"""
}
