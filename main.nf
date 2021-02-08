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
     cutadapt --version &> v_cutadapt.txt
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



Channel
  .from( ch_input )
  .splitCsv(header:false, sep:',')
  .map { it = ["${it[0]}", "${it[4]}", "${it[5]}", "${it[6]}", "${it[7]}", "${it[8]}", "${it[9]}", "${it[10]}", "${it[11]}",
  [file("${cluster_path}/data/04_pfastq/${it[8]}/${it[4]}/${it[5]}/${it[11]}/demux_fastq/${it[0]}_${it[4]}_${it[5]}_R1.fq.gz", checkIfExists: true),
  file("${cluster_path}/data/04_pfastq/${it[8]}/${it[4]}/${it[5]}/${it[11]}/demux_fastq/${it[0]}_${it[4]}_${it[5]}_R2.fq.gz", checkIfExists: true)]]}
  .into { ch_subset
          ch_fastqc_original }



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
    set val(sample), val(run_id), val(lane), val(date), val(protocol), val(platform), val(source), val(genome), val(user), path(reads) from ch_subset

    output:
    set val(sample), val(run_id), val(lane), val(date), val(protocol), val(platform), val(source), val(genome), val(user), path("*QC.fq.gz") into ch_trimming

    script:
    read1 = "${sample}_${run_id}_${lane}_R1.QC.fq"
    read2 = "${sample}_${run_id}_${lane}_R2.QC.fq"

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
  publishDir "${cluster_path}/data/05_QC/${project}/Trimming/${sample}", mode: 'copy',
  saveAs: { filename ->
    if (filename.endsWith(".log")) "logs/$filename"
    else if (filename.endsWith("_fastqc.html")) "FastQC/$filename"
    else if (filename.endsWith("_fastqc.zip")) "FastQC/zips/$filename"
  }

  input:
  set val(sample), val(run_id), val(lane), val(date), val(protocol), val(platform), val(source), val(genome), val(user), path(subset) from ch_trimming

  output:
  set val(sample), path("*cutadapt.fq.gz"), val(run_id), val(lane), val(date), val(protocol), val(platform), val(source), val(genome), val(user) into ch_star
  path("*_fastqc.{zip,html}") into fastqc_trimmed_multiqc //multiqc
  path("*report.txt") into trimmed_multiqc

  script:
  trimmed1 = "${sample}_${run_id}_${lane}_R1.QC.cutadapt.fq.gz"
  trimmed2 = "${sample}_${run_id}_${lane}_R2.QC.cutadapt.fq.gz"
  umi = "${sample}_${run_id}_${lane}_UMI.fq.gz"

  if (protocol == 'RNAseq_3_S' | protocol == 'RNAseq_3_ULI') {

    //  # First remove UMI and 30nt at 5' of read1 (Generate file with UMIs)
    //  # Remove polyT tail from 5' end and nextera adaptor from 3' end in read1 (lower case options)
    //  # Remove polyA tail in read2 (upper case options) and truseq adaptor (both in 3'), 18N accounts for 8nt BC+10nt UMI (in the truseq adapter)

    """
    cutadapt -l 10 -j 0 -o $umi ${subset[0]}

    cutadapt -u 40 -g "polyA_Tail=T{100}" --minimum-length 20 -a Nextera=CTGTCTCTTATACACATCT \
    -A "polyA_Tail=A{100}" -A "TruSeq=N{18}AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT;min_overlap=25" -n 2 --pair-filter=any \
    -o $trimmed1 -p $trimmed2 -j 0 \
    ${subset[0]} ${subset[1]} > "${sample}_${run_id}_${lane}_report.txt"

    fastqc --quiet --threads $task.cpus $trimmed1 $trimmed2
    """
  } else {

    //  # Remove nextera adaptor from 3' end in read1 (lower case options)
    //  # Remove truseq adaptor (both in 3') from read2 (upper case options)

    """
    cutadapt --minimum-length 20 -a "Nextera=CTGTCTCTTATACACATCT" -a "TruSeq=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" \
    -A "TruSeq=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" -A "Nextera=CTGTCTCTTATACACATCT" -n 2 --pair-filter=any \
    -o $trimmed1 -p $trimmed2 -j 0 \
    ${subset[0]} ${subset[1]} > "${sample}_${run_id}_${lane}_report.txt"

    fastqc --quiet --threads $task.cpus $trimmed1 $trimmed2
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
  set val(sample), path(trimmed), val(run_id), val(lane), val(date), val(protocol), val(platform), val(source), val(genome), val(user) from ch_star

  output:
  set val(sample), path("*.sortedByCoord.out.bam"), path("*Aligned.sortedByCoord.out.bam.bai") into ch_bam_samtools

  set val(sample), path("*.sortedByCoord.out.bam") into ch_bam_stats
  //set val(sample), path("*Aligned.sortedByCoord.out.bam.bai") into bam_index
  set val(sample), path("*unmapped*") optional true
  path("*.out") into star_logs //multiqc
  set val(sample), path("*.tab") into ch_tab



  script:
  if (genome == "mm38") {
    gtf = file("${cluster_path}/References/iGenomes/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf", checkIfExists: true)
    index = file("${cluster_path}/References/iGenomes/Mus_musculus/Ensembl/GRCm38/Annotation/Sequence/STARIndex/", checkIfExists: true)
    //picardref = "${cluster_path}/scripts/genomic_reference_data/bowtieIndexes/mm10_Bowtie2/mm10.fa"
    //picardrefflat = "${cluster_path}/scripts/genomic_reference_data/mm10/refFlat.txt"
  } else if (genome == "hg38") {
    gtf = file("${cluster_path}/References/iGenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf", checkIfExists: true)
    index = file("${cluster_path}/References/iGenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Sequence/STARIndex/", checkIfExists: true)

    //picardref = "${cluster_path}/scripts/genomic_reference_data/bowtieIndexes/hg19_Bowtie2/hg19_no_r.fa"
    //picardrefflat = "${cluster_path}/scripts/genomic_reference_data/hg19/refFlat.txt"
  }

  prefix = params.complete ? "${sample}_${run_id}.STAR.${genome}." : "${sample}_${run_id}.STAR.${genome}.QC."
  unaligned = params.complete ? "--outReadsUnmapped Fastx" : ''

  def star_mem = task.memory ?: params.star_memory ?: false
  def avail_mem = star_mem ? "--limitBAMsortRAM ${star_mem.toBytes() - 100000000}" : ''

  //First make the alingment for each read separated to obtain later the metrics per file
  """
  STAR \
  --sjdbGTFfile $gtf \
  --runThreadN ${task.cpus} \
  --genomeDir $index \
  --readFilesIn $trimmed \
  --readFilesCommand zcat \
  --outWigType bedGraph \
  --outFileNamePrefix $prefix \
  --outSAMtype BAM SortedByCoordinate $avail_mem \
  --outSAMattributes NH HI AS NM MD \
  --outSAMattrRGline ID:${sample}.R1.${protocol} SM:${sample} LB:${protocol} PL:${platform} CN:${source} DT:${date}T00:00:00-0400 \
  $unaligned


  if [ -f ${prefix}.Unmapped.out.mate1 ]; then
    mv ${prefix}.Unmapped.out.mate1 ${prefix}.unmapped_R1.fq
    pigz -p $task.cpus ${prefix}.unmapped_R1.fq
  fi
  if [ -f ${prefix}.Unmapped.out.mate2 ]; then
    mv ${prefix}.Unmapped.out.mate2 ${prefix}.unmapped_R2.fq
    pigz -p $task.cpus ${prefix}.unmapped_R2.fq
  fi

  samtools index ${prefix}Aligned.sortedByCoord.out.bam

  """


  }



process samtools {
  tag "$sample"
  label 'process_medium'
  publishDir "${cluster_path}/data/05_QC/${project}/QC_samtools/${sample}", mode: 'copy',
  saveAs: {filename ->
              if (filename.indexOf("flagstat") > 0) "flagstat/$filename"
         else if (filename.indexOf("idxstats") > 0) "idxstats/$filename"
         else if (filename.indexOf(".stats") > 0) "stats/$filename"
       }

  input:
  set val(sample), path(bam), path(index) from ch_bam_samtools

  output:
  path("*.flagstat") into flagstat //multiqc
  path("*.idxstats") into indxstats //multiqc
  path("*.stats") into bam_stats //multiqc


  script:
  """
  samtools flagstat $bam > ${bam}.flagstat
  samtools idxstats $bam > ${bam}.idxstats
  samtools stats $bam > ${bam}.stats
  """
}



process rseqc {
  tag "$sample"
  label 'process_medium'
  publishDir "${cluster_path}/data/05_QC/${project}/QC_rseqc/${sample}", mode: 'copy',
  saveAs: {filename ->
         if (filename.indexOf("bam_stat.txt") > 0)                      "bam_stat/$filename"
    else if (filename.indexOf("read_duplication.DupRate_plot.pdf") > 0) "read_duplication/$filename"
    else if (filename.indexOf("read_duplication.DupRate_plot.r") > 0)   "read_duplication/rscripts/$filename"
    else if (filename.indexOf("read_duplication.pos.DupRate.xls") > 0)  "read_duplication/dup_pos/$filename"
    else if (filename.indexOf("read_duplication.seq.DupRate.xls") > 0)  "read_duplication/dup_seq/$filename"
  }

  input:
  set val(sample), path(bam) from ch_bam_stats

  output:
  path("${bam.baseName}.bam_stat.txt") into rseqc_bam //multiqc
  path("*.{xls,pdf,r}") into rseqc_dup


  script:
  """
  bam_stat.py -i $bam 2> ${bam.baseName}.bam_stat.txt
  read_duplication.py -i $bam -o ${bam.baseName}.read_duplication
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
     label 'process_medium'
     publishDir "${cluster_path}/data/05_QC/${project}/FastQC/${sample}", mode: 'copy',
     saveAs: { filename ->
       filename.endsWith(".zip") ? "zips/$filename" : filename
     }

     input:
     set val(sample), val(run_id), val(lane), val(date), val(protocol), val(platform), val(source), val(genome), val(user), path(reads) from ch_fastqc_original

     output:
     path("*_fastqc.{zip,html}") into fastqc_results //multiqc

     script:
     """
     fastqc --quiet --threads $task.cpus $reads
     """
   }




 process multiqc {

   publishDir "${cluster_path}/data/05_QC/${project}/FastQC_raw/${sample}", mode: 'copy',
   saveAs: { filename ->
     if (filename.endsWith(".html")) filename
   }

   input:
   path multiqc_config from ch_multiqc_config
   //path multiqc_custom_config
   path ('software_versions/*') from software_versions_yaml.collect()
   path workflow_summary from create_workflow_summary(summary)
   path('fastqc/*') from fastqc_results.collect().ifEmpty([])
   path('cutadapt/*') from trimmed_multiqc.collect().ifEmpty([])
   path('cutadapt/fastqc/*') from fastqc_trimmed_multiqc.collect().ifEmpty([])
   path('star/*') from star_logs.collect().ifEmpty([])
   path ('samtools/stats/*') from bam_stats.collect().ifEmpty([])
   path ('samtools/flagstat/*') from flagstat.collect().ifEmpty([])
   path ('samtools/idxstats/*') from indxstats.collect().ifEmpty([])
   path ('rseqc/bam_stat/*') from rseqc_bam.collect().ifEmpty([])
   path ('rseqc/read_duplication/*') from rseqc_dup.collect().ifEmpty([])

   output:
   path "*multiqc_report.html"
   path "*_data"
   path "*_plots"

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
