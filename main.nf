#!/usr/bin/env nextflow
/*
========================================================================================
                         mpozuelo/QC_demux
========================================================================================
mpozuelo/QC_demux Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/mpozuelo/QC_demux
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info mpozueloHeader()
    log.info """

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run mpozuelo/QC_demux --input '*.txt' -profile docker

    Mandatory arguments:
      --input [file]                Samplesheet samples information (only samples that want to show in the same MultiQC)
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity.
      --cluster_path                Cluster path to store data and to search for input data (Default: /datos/ngs/dato-activo)
      --project                     Name of the project to create the folder and run the quality steps

    QC:
      --skipQC                      Skip all QC steps apart from MultiQC
      --skipFastQC                  Skip FastQC


    Other options
      --outdir                      The output directory where the results will be saved
      -w/--work-dir                 The temporary directory where intermediate data will be saved
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    """.stripIndent()
}


// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


/*
 * SET UP CONFIGURATION VARIABLES
 */


 // Has the run name been specified by the user?
 //  this has the bonus effect of catching both -name and --name
 custom_runName = params.name
 if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
   custom_runName = workflow.runName
 }
 else{
   workflow.runName = params.user + " " + params.timestamp
   custom_runName = workflow.runName
 }


// Validate inputs

if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { exit 1, "Input samplesheet file not specified!" }

if (!params.outdir) {
  params.outdir = params.run
}

// Mandatory arguments for the publishDir option (cluster_path has a default but project is mandatory as input)

cluster_path = params.cluster_path
project = params.project


// Stage multiqc config files
ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_image_docs            = file("$baseDir/assets/figures/Logo_IdisNA_CIMA.png", checkIfExists: true)



// Header log info
log.info mpozueloHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name'] = custom_runName ?: workflow.runName
summary['Project'] = params.project
summary['Input'] = params.input
summary['Max Resources'] = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['User'] = workflow.userName

summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"


// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'mpozuelo-QC_demux-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'mpozuelo/QC_demux Workflow Summary'
    section_href: 'https://github.com/mpozuelo/QC_demux'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}




 /*
  * Parse software version numbers
  */
 process get_software_versions {
     publishDir "${params.outdir}/pipeline_info", mode: 'copy',
         saveAs: { filename ->
             if (filename.indexOf(".csv") > 0) filename
             else null
         }

     output:
     file 'software_versions_mqc.yaml' into software_versions_yaml
     file "software_versions.csv"

     script:
     """
     echo $workflow.manifest.version &> v_ngi_QC.txt
     echo $workflow.nextflow.version &> v_nextflow.txt
     fastqc --version &> v_fastqc.txt
     trim_galore --version &> v_trim_galore.txt
     STAR --version &> v_star.txt
     samtools --version &> v_samtools.txt
     multiqc --version &> v_multiqc.txt
     scrape_software_versions.py &> software_versions_mqc.yaml
     """
 }


/*
 * LOAD SAMPLESHEET and assign get the columns we will use for demultiplexing

 /*
  * LOAD SAMPLESHEET and assign get the columns we will use for demultiplexing
  It contains the following columns:
 1- Name given to the sample
 2- Index
 3- Index2
 4- Barcode
 5- Run ID
 6- Lane
 7- Sequencing date
 8- Protocol
 9- Platform
 10- Sample Source (place where the samples come from)
 11- Genome
 12- User
 13- Coverage
 */

 process modify_samplesheet {
    publishDir "${params.outdir}/samplesheet/", mode: params.publish_dir_mode

    input:
    path samplesheet from ch_input

    output:
    path "samplesheet_validated.csv" into ch_samplesheet

    script:
    out = "samplesheet_validated.csv"

    """
    modify_samplesheet.py $samplesheet $out
    """
  }



  def validate_input(LinkedHashMap sample) {

    def sample_id = sample.sampleID
    def index = sample.index
    def index2 = sample.index2
    def barcode = sample.barcode
    def run = sample.run
    def lane = sample.lane
    def date = sample.date
    def protocol = sample.protocol
    def platform = sample.platform
    def source = sample.source
    def genome = sample.genome
    def user = sample.user
    def star_gtf = sample.star_gtf
    def star_index = sample.star_index
    def fastq1 = sample.fastq1
    def fastq2 = sample.fastq2

    def array = []
    array = [ sample_id, [file(fastq1, checkIfExists: true), file(fastq2, checkIfExists: true)], run, lane, date, protocol, platform, source, genome,
      user, [file(star_gtf, checkIfExists: true), file(star_index, checkIfExists: true)] ]

    return array
  }

  /*
  * Create channels for input fastq files
  */
  ch_samplesheet
  .splitCsv(header:true, sep:',')
  .map { validate_input(it) }
  .into { ch_subset
          ch_fastq }



/*
 * STEP 1 - Make subset
 */
//Get a subset of 10% of the reads in the complete file,
if (!params.complete) {
  process subset_10pc {
    tag "$sample"
    label 'process_low'
    publishDir "${cluster_path}/data/05_QC/${project}/subset_fastq/${sample}", mode: 'copy',
    saveAs: { filename ->
      filename.endsWith("QC.fq.gz") ? filename : null
    }

    input:
    set val(sample), path(reads), val(run_id), val(lane), val(date), val(protocol), val(platform), val(source), val(genome), val(user), path(star) from ch_subset

    output:
    set val(sample), path("*QC*.fq.gz"), val(run_id), val(lane), val(date), val(protocol), val(platform), val(source), val(genome), val(user), path(star) into ch_trimming

    script:
    read1 = "${sample}_${run_id}_${lane}_QC_R1.fq"
    read2 = "${sample}_${run_id}_${lane}_QC_R2.fq"

    //First count the number of total reads in one of the input files (R1 or R2, in this case R1) and get the 10%
    //Then get the subset of samples with seqtk
    """
    subset=(\$(echo \$((\$(echo -e `zcat ${reads[0]} | awk 'NR % 4 == 2' - | wc -l`)*10/100))))
    seqtk sample -s100 ${reads[0]} \$subset > $read1
    seqtk sample -s100 ${reads[1]} \$subset > $read2
    pigz -p $task.cpus $read1
    pigz -p $task.cpus $read2
    """
  }

} else {
  ch_subset
    .into { ch_trimming }
}






/*
* STEP 2 - Trimming
*/
// Make quality and trimming of possible adapters


process trimming {
  tag "$sample"
  label 'process_medium'
  publishDir "${cluster_path}/data/05_QC/${project}/trimgalore/${sample}", mode: 'copy',
  saveAs: { filename ->
    if (filename.endsWith(".log") || filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
    else if (filename.endsWith("_fastqc.html")) "fastqc/$filename"
    else if (filename.endsWith(".zip")) "fastqc/zips/$filename"
    else if (filename.indexOf("_fastqc") > 0) filename
    else if (filename.endsWith("UMI.fq.gz")) "UMIs/$filename"
    else null
  }


  input:
  set val(sample), path(subset), val(run_id), val(lane), val(date), val(protocol), val(platform), val(source), val(genome), val(user), path(star) from ch_trimming

  output:
  set val(sample), path("*R*.fq.gz"), val(run_id), val(lane), val(date), val(protocol), val(platform), val(source), val(genome), val(user), path(star) into ch_star
  file "*trimming_report.txt" into trimgalore_trim_mqc
  file "*_fastqc.{zip,html}" into trimgalore_fastqc_mqc
  file "*UMI.fq.gz" optional true

  script:
  //trimmed1 = "${sample}_${run_id}_${lane}_QC_R1.cutadapt.fq.gz"
  //trimmed2 = "${sample}_${run_id}_${lane}_QC_R2.cutadapt.fq.gz"
  umi = "${sample}_${run_id}_${lane}_UMI.fq.gz"
  woumi1 = "${sample}_${run_id}_${lane}_woUMI_R1.fq.gz"
  woumi2 = "${sample}_${run_id}_${lane}_woUMI_R2.fq.gz"

  if (protocol == 'RNAseq_3_S' | protocol == 'RNAseq_3_ULI') {

    //  # First remove UMI and 30nt at 5' of read1 (Generate file with UMIs)
    //  # Remove polyT tail from 5' end and nextera adaptor from 3' end in read1 (lower case options)
    //  # Remove polyA tail in read2 (upper case options) and truseq adaptor (both in 3'), 18N accounts for 8nt BC+10nt UMI (in the truseq adapter)


    /*cutadapt -q 30 -u 40 -g "polyA_Tail=T{100}" --minimum-length 20 -a Nextera=CTGTCTCTTATACACATCT \
    -A "polyA_Tail=A{100}" -A "TruSeq=N{18}AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT;min_overlap=25" -n 2 --pair-filter=any \
    -o $trimmed1 -p $trimmed2 -j 0 \
    ${subset[0]} ${subset[1]} > "${sample}_${run_id}_${lane}_report.txt"

    fastqc --quiet --threads $task.cpus $trimmed1 $trimmed2
    */
    """
    cutadapt -l 10 -j 0 -o $umi ${subset[0]}

    umi_tools extract -I ${subset[0]} -S $woumi1 --read2-in=${subset[1]} --read2-out=$woumi2 --bc-pattern=NNNNNNNNNN

    trim_galore \\
    -q 30 \\
    --paired \\
    -a " T{100} -a CTGTCTCTTATACACATCT" \\
    -a2 " A{100} -a N{18}AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" \\
    --length 20 \\
    -j $task.cpus \\
    --fastqc \\
    $woumi1 $woumi2
    """

  } else {

    //  # Remove nextera adaptor from 3' end in read1 (lower case options)
    //  # Remove truseq adaptor (both in 3') from read2 (upper case options)
/*    cutadapt -q 30 --minimum-length 20 -a "Nextera=CTGTCTCTTATACACATCT" -a "TruSeq=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" \
    -A "TruSeq=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" -A "Nextera=CTGTCTCTTATACACATCT" -n 2 --pair-filter=any \
    -o $trimmed1 -p $trimmed2 -j 0 \
    ${subset[0]} ${subset[1]} > "${sample}_${run_id}_${lane}_report.txt"

    fastqc --quiet --threads $task.cpus $trimmed1 $trimmed2*/

    """
    trim_galore \\
    -q 30 \\
    --paired \\
    -a " CTGTCTCTTATACACATCT -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" \\
    -a2 " CTGTCTCTTATACACATCT -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" \\
    --length 20 \\
    -j $task.cpus \\
    --fastqc \\
    ${subset[0]} ${subset[1]}
    """
  }
}



/*
 * STEP 3 - STAR alignment
 */

/* SEPARATED BY READ

process star {
  tag "$sample"
  label 'process_high'
  publishDir "${cluster_path}/data/05_QC/${project}/STAR", mode: 'copy',
  saveAs: { filename ->
    filename.endsWith(".log") ? "logs/$filename" : filename
  }

  input:
  set val(sample), path(trimmed), val(run_id), val(lane), val(date), val(protocol), val(platform), val(source), val(genome), val(user) from ch_star

  output:
  set val(sample), path("*.sortedByCoord.out.bam"), val(run_id), val(lane), val(protocol), val(platform), val(date), val(user) into ch_bam_qc
  set val(sample), path("*Unmapped*") optional true
  set val(sample), path("*Log.out") into ch_log_out
  set val(sample), path("*.tab") into ch_tab
  set val(sample), path("*Log.final.out") into ch_log_final
  set val(sample), path("*Log.progress.out") into ch_log_progress


  script:
  if (genome == "mm10") {
    starrefgenome = "${cluster_path}/scripts/genomic_reference_data/STAR/STARgenomes/GENCODE/GRCm38_GencodeM21"
    //picardref = "${cluster_path}/scripts/genomic_reference_data/bowtieIndexes/mm10_Bowtie2/mm10.fa"
    //picardrefflat = "${cluster_path}/scripts/genomic_reference_data/mm10/refFlat.txt"
  } else if (genome == "hg19") {
    starrefgenome = "${cluster_path}/scripts/genomic_reference_data/STAR/STARgenomes/GENCODE/hg19.v19"
    //picardref = "${cluster_path}/scripts/genomic_reference_data/bowtieIndexes/hg19_Bowtie2/hg19_no_r.fa"
    //picardrefflat = "${cluster_path}/scripts/genomic_reference_data/hg19/refFlat.txt"
  }

  read1 = "${reads[0]}"
  read2 = "${reads[1]}"
  bam1 = params.complete ? "${sample}_${run_id}_R1.STAR.${genome}." : "${sample}_${run_id}_R1.STAR.${genome}.QC."
  bam2 = params.complete ? "${sample}_${run_id}_R2.STAR.${genome}." : "${sample}_${run_id}_R2.STAR.${genome}.QC."
  unaligned = params.complete ? "--outReadsUnmapped Fastx" : ''


  //First make the alingment for each read separated to obtain later the metrics per file
  """
  STAR --runThreadN 4 \
  --genomeDir $starrefgenome \
  --readFilesIn $trimmed1 \
  --readFilesCommand zcat \
  --outFileNamePrefix $bam1 \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMunmapped Within \
  --limitBAMsortRAM 20000000000 \
  --outSAMattributes NH HI AS NM MD \
  --outSAMattrRGline ID:${sample}.R1.${protocol} SM:${sample} LB:${protocol} PL:${platform} CN:${source} DT:${date}T00:00:00-0400 \
  $unaligned

  STAR --runThreadN 4 \
  --genomeDir $starrefgenome \
  --readFilesIn $trimmed2 \
  --readFilesCommand zcat \
  --outFileNamePrefix $bam2 \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMunmapped Within \
  --limitBAMsortRAM 20000000000 \
  --outSAMattributes NH HI AS NM MD \
  --outSAMattrRGline ID:${sample}.R2.${protocol} SM:${sample} LB:${protocol} PL:${platform} CN:${source} DT:${date}T00:00:00-0400 \
  $unaligned
  """
  }
}

*/



process star {
  tag "$sample"
  label 'process_high'
  publishDir "${cluster_path}/data/05_QC/${project}/STAR/${sample}", mode: 'copy',
  saveAs: { filename ->
    if (filename.endsWith(".out")) "logs/$filename"
    else if (filename.endsWith("*.fq.gz")) "unmapped/$filename"
    //else if (filename.endsWith(".tab")) ""
    else if (filename.endsWith(".bam")) "bam/$filename"
    else if (filename.endsWith(".bai")) "bam_index/$filename"
  }

  input:
  set val(sample), path(trimmed), val(run_id), val(lane), val(date), val(protocol), val(platform), val(source), val(genome), val(user), path(star) from ch_star

  output:
  set val(sample), path("*.sortedByCoord.out.bam"), path("*Aligned.sortedByCoord.out.bam.bai") into ch_bam_samtools,
                                                                                                    ch_bam_stats
  set val(sample), path("*.sortedByCoord.out.bam") into ch_markduplicates
  //set val(sample), path("*Aligned.sortedByCoord.out.bam.bai") into bam_index
  set val(sample), path("*unmapped*") optional true
  path("*.out") into star_logs //multiqc
  set val(sample), path("*.tab") into ch_tab
  path("*.mapped.tsv") into concatenate_mapped



  script:

  prefix = params.complete ? "${sample}_${run_id}.STAR.${genome}." : "${sample}_${run_id}.STAR.${genome}.QC."
  unaligned = params.complete ? "--outReadsUnmapped Fastx" : ''

  def star_mem = task.memory ?: params.star_memory ?: false
  def avail_mem = star_mem ? "--limitBAMsortRAM ${star_mem.toBytes() - 100000000}" : ''

  //First make the alingment for each read separated to obtain later the metrics per file
  """
  STAR \\
  --sjdbGTFfile ${star[0]} \\
  --runThreadN $task.cpus \\
  --genomeDir ${star[1]} \\
  --readFilesIn $trimmed \\
  --readFilesCommand zcat \\
  --outWigType bedGraph \\
  --outSAMunmapped Within \\
  --outFileNamePrefix $prefix \\
  --outSAMtype BAM SortedByCoordinate $avail_mem \\
  --outSAMattributes NH HI AS NM MD \\
  --outSAMattrRGline ID:${sample}.R1.${protocol} SM:${sample} LB:${protocol} PL:${platform} CN:${source} DT:${date}T00:00:00-0400 \\
  $unaligned


  if [ -f ${sample}.Unmapped.out.mate1 ]; then
    mv ${sample}.Unmapped.out.mate1 ${sample}.unmapped_R1.fq
    pigz -p $task.cpus ${sample}.unmapped_R1.fq
  fi
  if [ -f ${sample}.Unmapped.out.mate2 ]; then
    mv ${prefix}.Unmapped.out.mate2 ${sample}.unmapped_R2.fq
    pigz -p $task.cpus ${sample}.unmapped_R2.fq
  fi

  samtools index ${sample}Aligned.sortedByCoord.out.bam

  uniquely=\$(grep "Uniquely mapped reads %" ${sample}.Log.final.out | | grep -Eo '[0-9.]+%')
  multiple=\$(grep "% of reads mapped to multiple loci" ${sample}.Log.final.out | | grep -Eo '[0-9.]+%)
  many=\$(grep "% of reads mapped to too many loci" ${sample}.Log.final.out | | grep -Eo '[0-9.]+%)
  printf "%s\t%s\t%s\t%s" "${sample}" "\$uniquely" "\$multiple" "\$many" > ${sample}.mapped.tsv
  """
  }



process samtools {
  tag "$sample"
  label 'process_medium'
  publishDir "${cluster_path}/data/05_QC/${project}/samtools/${sample}", mode: 'copy'
  /*saveAs: {filename ->
              if (filename.indexOf("flagstat") > 0) "flagstat/$filename"
         else if (filename.indexOf("idxstats") > 0) "idxstats/$filename"
         else if (filename.indexOf(".stats") > 0) "stats/$filename"
       }*/

  input:
  set val(sample), path(bam), path(index) from ch_bam_samtools

  output:
  //path("*.flagstat") into flagstat //multiqc
  //path("*.idxstats") into indxstats //multiqc
  path("*.stats") into bam_stats //multiqc

/*  samtools flagstat $bam > ${bam}.flagstat
  samtools idxstats $bam > ${bam}.idxstats */

  script:
  """
  samtools stats $bam > ${bam}.stats
  """
}



process rseqc {
  tag "$sample"
  label 'process_medium'
  publishDir "${cluster_path}/data/05_QC/${project}/rseqc/${sample}", mode: 'copy',
  saveAs: {filename ->
         if (filename.indexOf("bam_stat.txt") > 0)                      "bam_stat/$filename"
    else if (filename.indexOf("read_duplication.DupRate_plot.pdf") > 0) "read_duplication/$filename"
    else if (filename.indexOf("read_duplication.DupRate_plot.r") > 0)   "read_duplication/rscripts/$filename"
    else if (filename.indexOf("read_duplication.pos.DupRate.xls") > 0)  "read_duplication/dup_pos/$filename"
    else if (filename.indexOf("read_duplication.seq.DupRate.xls") > 0)  "read_duplication/dup_seq/$filename"
  }

  input:
  set val(sample), path(bam), path(bai) from ch_bam_stats

  output:
  path("${bam.baseName}.bam_stat.txt") into rseqc_bam //multiqc
  path("*.{xls,pdf,r}") into rseqc_dup


  script:
  """
  bam_stat.py -i $bam 2> ${bam.baseName}.bam_stat.txt
  read_duplication.py -i $bam -o ${bam.baseName}.read_duplication
  """

}


//Mark duplicates

process picard {
  tag "${bam.baseName - '.sorted'}"
  label 'process_medium'
  publishDir "${cluster_path}/data/05_QC/${project}/markduplicates/", mode: params.publish_dir_mode,
  saveAs: {filename -> filename.indexOf("_metrics.txt") > 0 ? "metrics/$filename" : "$filename"}

  input:
  set val(name), file(bam) from ch_markduplicates

  output:
  file "${bam.baseName}.markDups_metrics.txt" into picard_mrkd_results
  file "${bam.baseName}.sort.qual_score_dist.txt" into picard_distribution_results

  script:
  markdup_java_options = (task.memory.toGiga() > 8) ? params.markdup_java_options : "\"-Xms" +  (task.memory.toGiga() / 2 )+"g "+ "-Xmx" + (task.memory.toGiga() - 1)+ "g\""

  """
  picard ${markdup_java_options} MarkDuplicates \\
  INPUT=$bam \\
  OUTPUT=${bam.baseName}.markDups.bam \\
  METRICS_FILE=${bam.baseName}.markDups_metrics.txt \\
  REMOVE_DUPLICATES=false \\
  ASSUME_SORTED=true \\
  PROGRAM_RECORD_ID='null' \\
  VALIDATION_STRINGENCY=LENIENT

  picard ${markdup_java_options} QualityScoreDistribution \\
  INPUT=$bam \\
  O=${bam.baseName}.sort.qual_score_dist.txt \
  CHART=${bam.baseName}.sort.qual_score_dist.pdf
  """
}


/*
 * STEP 4 - FastQC subset or complete in case we do not perform any subset
 */
 /*if (!params.skipQC || !params.skipFastQC) {

   process fastqc {
     tag "$sample"
     label 'process_medium'
     publishDir "${cluster_path}/data/05_QC/${project}/FastQC_raw/${sample}", mode: 'copy',
     saveAs: { filename ->
       filename.endsWith(".zip") ? "zips/$filename" : filename
     }

     when:
     !params.complete

     input:
     set val(sample), val(run_id), val(lane), val(date), val(protocol), val(platform), val(source), val(genome), val(user), path(reads) into ch_fastqc_subset

     output:
     path("*_fastqc.{zip,html}") into fastqc_results //multiqc

     script:
     """
     fastqc --quiet --threads $task.cpus $reads
     """
   }
 } else {
   fastqc_results = Channel.empty()
 }


  /*
   * STEP 5 - FastQC of the whole sample
   */

  process fastqc {
     tag "$sample"
     label 'process_low'
     publishDir "${cluster_path}/data/05_QC/${project}/fastqc/${sample}", mode: 'copy',
     saveAs: { filename ->
       filename.endsWith(".zip") ? "zips/$filename" : filename
     }

     input:
     set val(sample), path(reads), val(run_id), val(lane), val(date), val(protocol), val(platform), val(source), val(genome), val(user), path(star) from ch_fastq

     output:
     path("*_fastqc.{zip,html}") into fastqc_results //multiqc
     path("*.total.reads.tsv") into total_reads_merge

     script:
     """
     totalReads=\$(echo \$(echo -e `zcat ${reads[0]} | awk 'NR % 4 == 2' - | wc -l`))
     q301=\$(echo \$(q30.py ${reads[0]}))
     q302=\$(echo \$(q30.py ${reads[1]}))
     printf "%s\t%s\t%s\t%s" "${sample}" "\$totalReads" "\$q301" "\$q302" > "${sample}.total.reads.tsv"
     fastqc --quiet --threads $task.cpus $reads
     """
   }


  process merge_files {
    tag "merge"
    label 'process_low'
    publishDir "${cluster_path}/data/05_QC/${project}/QCTable/", mode: 'copy',
    saveAs: { filename ->
      filename.endsWith(".zip") ? "zips/$filename" : filename
    }

    input:
    path("*") from total_reads_merge.collect().ifEmpty([])
    path("*") from concatenate_mapped.collect().ifEmpty([])

    output:
    path("*QC.table.tsv")

    script:
    """
    printf "%s\t%s\t%s\t%s" "sampleID" "TotalReads" "Q30%1" "Q30%2" "UniquelyMapped%" "MultipleMapped%" "TooManyMapped%" > "${project}.QC.table.tsv
    cat "*.total.reads.tsv" >> total.reads.tsv
    sort total.reads.tsv > total.reads.sort.tsv
    cat "*.mapped.tsv" >> mapped.tsv
    sort mapped.tsv > mapped.sort.tsv
    join -t "\t" -j 1 total.reads.sort.tsv mapped.sort.tsv > QC.table.tsv
    cat QC.table.tsv >> "${project}.QC.table.tsv
    """
  }


 process multiqc {

   publishDir "${cluster_path}/data/05_QC/${project}/multiqc/", mode: 'copy',
   saveAs: { filename ->
     if (filename.endsWith(".html")) filename
   }

   input:
   path multiqc_config from ch_multiqc_config
   //path multiqc_custom_config
   path ('software_versions/*') from software_versions_yaml.collect()
   path workflow_summary from create_workflow_summary(summary)
   path('fastqc/*') from fastqc_results.collect().ifEmpty([])
   path('trimgalore/*') from trimgalore_trim_mqc.collect().ifEmpty([])
   path('trimgalore/fastqc/*') from trimgalore_fastqc_mqc.collect().ifEmpty([])
   path('star/*') from star_logs.collect().ifEmpty([])
   path('samtools/stats/*') from bam_stats.collect().ifEmpty([])
   path('picard_mrkd/*') from picard_mrkd_results.collect().ifEmpty([])
   path('picard_dist_score/*') from picard_distribution_results.collect().ifEmpty([])
   //path ('samtools/flagstat/*') from flagstat.collect().ifEmpty([])
   //path ('samtools/idxstats/*') from indxstats.collect().ifEmpty([])
   path ('rseqc/bam_stat/*') from rseqc_bam.collect().ifEmpty([])
   path ('rseqc/read_duplication/*') from rseqc_dup.collect().ifEmpty([])
   path "*" from ch_image_docs

   output:
   path "*multiqc_report.html"
   path "*_data"
   //path "multiqc_plots"

   script:
   rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
   rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
   """
   multiqc . -f $rtitle $rfilename --config $multiqc_config
   """

 }




/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[mpozuelo/QC_demux] Successful: $workflow.runName"

    if (!workflow.success) {
      subject = "[mpozuelo/QC_demux] FAILED: $workflow.runName"
    }



    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";


    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
        log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}"
        log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}"
    }

    if (workflow.success) {
        log.info "${c_purple}[mpozuelo/QC_demux]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[mpozuelo/QC_demux]${c_red} Pipeline completed with errors${c_reset}"
    }

}

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

def mpozueloHeader() {
  // Log colors ANSI codes
  c_blue = params.monochrome_logs ? '' : "\033[0;34m";
  c_dim = params.monochrome_logs ? '' : "\033[2m";
  c_white = params.monochrome_logs ? '' : "\033[0;37m";
  c_reset = params.monochrome_logs ? '' : "\033[0m";


  return """    -${c_dim}--------------------------------------------------${c_reset}-
  ${c_blue}  __  __  __   __  ___         ${c_reset}
  ${c_blue}  | \\/ | |__| |  |  /  |  |     ${c_reset}
  ${c_blue}  |    | |    |__| /__ |__|         ${c_reset}
  ${c_white}  mpozuelo/QC_demux v${workflow.manifest.version}${c_reset}
  -${c_dim}--------------------------------------------------${c_reset}-
  """.stripIndent()
}


def checkHostname() {
  def c_reset = params.monochrome_logs ? '' : "\033[0m"
  def c_white = params.monochrome_logs ? '' : "\033[0;37m"
  def c_red = params.monochrome_logs ? '' : "\033[1;91m"
  def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
  if (params.hostnames) {
    def hostname = "hostname".execute().text.trim()
    params.hostnames.each { prof, hnames ->
      hnames.each { hname ->
        if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
          log.error "====================================================\n" +
          "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
          "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
          "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
          "============================================================"
        }
      }
    }
  }
}
