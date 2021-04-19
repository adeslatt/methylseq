#!/usr/bin/env nextflow
/*
========================================================================================
                         adeslatt/methylseq
========================================================================================
*/

def helpMessage() {
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf --reads '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:
      --aligner [str]                   Alignment tool to use (default: bismark)
                                            Available: bismark, bismark_hisat, bwameth
      --reads [file]                    Path to input data (must be surrounded with quotes)
      -profile [str]                    Configuration profile to use. Can use multiple (comma separated)
                                            Available: conda, docker, singularity, test, awsbatch, <institute> and more

    Options:
     --zenodo_reference [int]           Integer for the Zenodo DOI reference for the reference genome
     --genome [str]                     Name of iGenomes reference
     --single_end [bool]                Specifies that the input is single end reads
     --comprehensive [bool]             Output information for all cytosine contexts
     --cytosine_report [bool]           Output stranded cytosine report during Bismark's bismark_methylation_extractor step.
     --ignore_flags [bool]              Run MethylDackel with the flag to ignore SAM flags.
     --meth_cutoff [int]                Specify a minimum read coverage to report a methylation call during Bismark's bismark_methylation_extractor step.
     --min_depth [int]                  Specify a minimum read coverage for MethylDackel to report a methylation call.
     --methyl_kit [bool]                Run MethylDackel with the --methyl_kit flag to produce files suitable for use with the methylKit R package.
     --skip_deduplication [bool]        Skip deduplication step after alignment. This is turned on automatically if --rrbs is specified
     --non_directional [bool]           Run alignment against all four possible strands
     --save_align_intermeds [bool]      Save aligned intermediates to results directory
     --save_trimmed [bool]              Save trimmed reads to results directory
     --unmapped [bool]                  Save unmapped reads to fastq files
     --relax_mismatches [bool]          Turn on to relax stringency for alignment (set allowed penalty with --num_mismatches)
     --num_mismatches [float]           0.6 will allow a penalty of bp * -0.6 - for 100bp reads (bismark default is 0.2)
     --known_splices [file]             Supply a .gtf file containing known splice sites (bismark_hisat only)
     --slamseq [bool]                   Run bismark in SLAM-seq mode
     --local_alignment [bool]           Allow soft-clipping of reads (potentially useful for single-cell experiments)
     --bismark_align_cpu_per_multicore [int] Specify how many CPUs are required per --multicore for bismark align (default = 3)
     --bismark_align_mem_per_multicore [str] Specify how much memory is required per --multicore for bismark align (default = 13.GB)

    References                          If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta [file]                    Path to fasta reference
      --fasta_index [path]              Path to Fasta Index
      --bismark_index [path]            Path to Bismark index
      --bwa_meth_index [path]           Path to bwameth index
      --save_reference [bool]           Save reference(s) to results directory

    Trimming options:
     --skip_trimming [bool]             Skip read trimming
     --clip_r1 [int]                    Trim the specified number of bases from the 5' end of read 1 (or single-end reads).
     --clip_r2 [int]                    Trim the specified number of bases from the 5' end of read 2 (paired-end only).
     --three_prime_clip_r1 [int]        Trim the specified number of bases from the 3' end of read 1 AFTER adapter/quality trimming
     --three_prime_clip_r2 [int]        Trim the specified number of bases from the 3' end of read 2 AFTER adapter/quality trimming
     --rrbs [bool]                      Turn on if dealing with MspI digested material.

    Trimming presets:
     --pbat [bool]
     --single_cell [bool]
     --epignome [bool]
     --accell [bool]
     --zymo [bool]
     --cegx [bool]

    Other options:
     --outdir [file]                    The output directory where the results will be saved
     --email [email]                    Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
     --email_on_fail [email]            Same as --email, except only send mail if the workflow is not successful
     --max_multiqc_email_size [str]     Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
     -name [str]                        Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch options:
      --awsqueue [str]                  The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]                 The AWS Region for your AWS Batch job to run on
      --awscli [str]                    Path to the AWS CLI tool

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
custom_runName           = params.name
zenodo_doi               = params.zenodo_reference


/*--------------------------------------------------
  Channel setup
---------------------------------------------------*/

// set the input channels

ch_multiqc_config        = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
// To me This is overly complex - difficult to read and takes more than a few seconds to sort out
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs                          = file("$baseDir/docs/output.md", checkIfExists: true)

ch_flagstat_results_for_multiqc         = Channel.from(false)
ch_samtools_stats_results_for_multiqc   = Channel.from(false)
ch_markDups_results_for_multiqc         = Channel.from(false)
ch_methyldackel_results_for_multiqc     = Channel.from(false)

ch_bismark_index_for_bismark_align      = Channel.fromPath ( params.bismark_index )
ch_bismark_index_for_bismark_methXtract = Channel.fromPath ( params.bismark_index )

def summary = [:]
summary['Run Name']  = custom_runName ?: workflow.runName
summary['Reads']     = params.reads
summary['Aligner']   = params.aligner
summary['Data Type'] = params.single_end ? 'Single-End' : 'Paired-End'


Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-methylseq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/methylseq Workflow Summary'
    section_href: 'https://github.com/nf-core/methylseq'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

// Trimming presets
clip_r1 = params.clip_r1
clip_r2 = params.clip_r2
three_prime_clip_r1 = params.three_prime_clip_r1
three_prime_clip_r2 = params.three_prime_clip_r2
if(params.pbat){
    clip_r1 = 9
    clip_r2 = 9
    three_prime_clip_r1 = 9
    three_prime_clip_r2 = 9
}
else if( params.single_cell ){
    clip_r1 = 6
    clip_r2 = 6
    three_prime_clip_r1 = 6
    three_prime_clip_r2 = 6
}
else if( params.epignome ){
    clip_r1 = 8
    clip_r2 = 8
    three_prime_clip_r1 = 8
    three_prime_clip_r2 = 8
}
else if( params.accel || params.zymo ){
    clip_r1 = 10
    clip_r2 = 15
    three_prime_clip_r1 = 10
    three_prime_clip_r2 = 10
}
else if( params.cegx ){
    clip_r1 = 6
    clip_r2 = 6
    three_prime_clip_r1 = 2
    three_prime_clip_r2 = 2
}

/*
 * Get Zenodo Reference files
 */
process getZenodoReference {
    publishDir "data", mode: 'copy'

    input:

    output:
	file "done.txt" into ch_run_me_first
    script:
    """
    chmod 775 bin/*.py
    mkdir -p data
    cd data
    mkdir -p BismarkIndex
    cd BismarkIndex
    wget https://zenodo.org/record/$zenodo_doi/files/Bisulfite_Genome.tar.gz
    tar xzvf Bisulfite_Genome.tar.gz
    cd ..
    wget https://zenodo.org/record/$zenodo_doi/files/hg19_lambda.fa.tar.gz
    tar xzvf hg19_lambda.fa.tar.gz
    samtools faidx hg19_lambda.fa
    touch done.txt
    cd ..
    """
}

process getReadsIntoQC {
    publishDir ".", mode: 'copy'

    input:
    file zenodo_reference from ch_run_me_first
    
    output:
    file "${params.reads}" into ch_read_files_for_fastqc, ch_wherearemyfiles_for_trimgalore, ch_wherearemyfiles_for_alignment
    script:
    """
    
    """

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
    input:
    file zenodo_reference from ch_run_me_first
    file fastq_files from ch_read_files_for_fastqc
    
    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml_for_multiqc
    file "software_versions.csv"

    script:
    """
    echo "$workflow.manifest.version" &> v_ngi_methylseq.txt
    echo "$workflow.nextflow.version" &> v_nextflow.txt
    bismark_genome_preparation --version &> v_bismark_genome_preparation.txt
    fastqc --version &> v_fastqc.txt
    cutadapt --version &> v_cutadapt.txt
    trim_galore --version &> v_trim_galore.txt
    bismark --version &> v_bismark.txt
    deduplicate_bismark --version &> v_deduplicate_bismark.txt
    bismark_methylation_extractor --version &> v_bismark_methylation_extractor.txt
    bismark2report --version &> v_bismark2report.txt
    bismark2summary --version &> v_bismark2summary.txt
    samtools --version &> v_samtools.txt
    hisat2 --version &> v_hisat2.txt
    bwa &> v_bwa.txt 2>&1 || true
    bwameth.py --version &> v_bwameth.txt
    picard MarkDuplicates --version &> v_picard_markdups.txt 2>&1 || true
    MethylDackel --version &> v_methyldackel.txt
    qualimap --version &> v_qualimap.txt || true
    preseq &> v_preseq.txt
    multiqc --version &> v_multiqc.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    label 'process_medium'
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: { filename ->
                      filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
                }

    input:
    set val(name), file(reads) from ch_read_files_for_fastqc

    output:
    file '*_fastqc.{zip,html}' into ch_read_files_for_trim_galore, ch_fastqc_results_for_multiqc

    script:
    """
    fastqc --quiet --threads $task.cpus $reads
    """
}

/*
 * STEP 2 - Trim Galore!
 */
process trim_galore {
    tag "$name"
    publishDir "${params.outdir}/trim_galore", mode: 'copy',
        saveAs: {filename ->
            if( filename.indexOf("_fastqc") > 0 ) "FastQC/$filename"
            else if( filename.indexOf("trimming_report.txt" ) > 0) "logs/$filename"
            else if( !params.save_trimmed && filename == "where_are_my_files.txt" ) filename
            else if( params.save_trimmed && filename != "where_are_my_files.txt" ) filename
            else null
        }

    input:
    set val(name), file(reads) from ch_read_files_for_trim_galore
    file wherearemyfiles from ch_wherearemyfiles_for_trimgalore.collect()

    output:
    set val(name), file('*fq.gz') into ch_trimmed_reads_for_alignment
    file "*trimming_report.txt" into ch_trim_galore_results_for_multiqc
    file "*_fastqc.{zip,html}"
    file "where_are_my_files.txt"

    script:
    def c_r1 = clip_r1 > 0 ? "--clip_r1 $clip_r1" : ''
    def c_r2 = clip_r2 > 0 ? "--clip_r2 $clip_r2" : ''
    def tpc_r1 = three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 $three_prime_clip_r1" : ''
    def tpc_r2 = three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 $three_prime_clip_r2" : ''
    def rrbs = params.rrbs ? "--rrbs" : ''
    def cores = 1
    if(task.cpus){
        cores = (task.cpus as int) - 4
        if (params.single_end) cores = (task.cpus as int) - 3
        if (cores < 1) cores = 1
        if (cores > 4) cores = 4
    }
    if( params.single_end ) {
        """
        trim_galore --fastqc --gzip $reads \
          $rrbs $c_r1 $tpc_r1 --cores $cores
        """
    } else {
        """
        trim_galore --fastqc --gzip --paired $reads \
          $rrbs $c_r1 $c_r2 $tpc_r1 $tpc_r2 --cores $cores
        """
    }
}

/*
 * STEP 3.1 - align with Bismark
 */
if( params.aligner =~ /bismark/ ){
    process bismark_align {
        tag "$name"
        publishDir "${params.outdir}/bismark_alignments", mode: 'copy',
            saveAs: {filename ->
                if( filename.indexOf(".fq.gz") > 0 ) "unmapped/$filename"
                else if( filename.indexOf("report.txt") > 0 ) "logs/$filename"
                else if( (!params.save_align_intermeds && !params.skip_deduplication && !params.rrbs).every() && filename == "where_are_my_files.txt" ) filename
                else if( (params.save_align_intermeds || params.skip_deduplication || params.rrbs).any() && filename != "where_are_my_files.txt" ) filename
                else null
            }

        input:
        set val(name), file(reads) from ch_trimmed_reads_for_alignment
        file index from ch_bismark_index_for_bismark_align.collect()

        output:
        set val(name), file("*.bam") into ch_bam_for_bismark_deduplicate, ch_bam_for_bismark_summary, ch_bam_for_preseq
        set val(name), file("*report.txt") into ch_bismark_align_log_for_bismark_report, ch_bismark_align_log_for_bismark_summary, ch_bismark_align_log_for_multiqc
        file "*.fq.gz" optional true
        file "where_are_my_files.txt"

        script:
        // Paired-end or single end input files
        input = params.single_end ? reads : "-1 ${reads[0]} -2 ${reads[1]}"

        // Choice of read aligner
        aligner = params.aligner == "bismark_hisat" ? "--hisat2" : "--bowtie2"

        // Optional extra bismark parameters
        splicesites = params.aligner == "bismark_hisat" && knownsplices.name != 'null' ? "--known-splicesite-infile <(hisat2_extract_splice_sites.py ${knownsplices})" : ''
        pbat = params.pbat ? "--pbat" : ''
        non_directional = params.single_cell || params.zymo || params.non_directional ? "--non_directional" : ''
        unmapped = params.unmapped ? "--unmapped" : ''
        mismatches = params.relax_mismatches ? "--score_min L,0,-${params.num_mismatches}" : ''
        soft_clipping = params.local_alignment ? "--local" : ''

        // Try to assign sensible bismark memory units according to what the task was given
        multicore = ''
        if( task.cpus ){
            // Numbers based on recommendation by Felix for a typical mouse genome
            if( params.single_cell || params.zymo || params.non_directional ){
                cpu_per_multicore = 5
                mem_per_multicore = (18.GB).toBytes()
            } else {
                cpu_per_multicore = 3
                mem_per_multicore = (13.GB).toBytes()
            }
            // Check if the user has specified this and overwrite if so
            if(params.bismark_align_cpu_per_multicore) {
                cpu_per_multicore = (params.bismark_align_cpu_per_multicore as int)
            }
            if(params.bismark_align_mem_per_multicore) {
                mem_per_multicore = (params.bismark_align_mem_per_multicore as nextflow.util.MemoryUnit).toBytes()
            }
            // How many multicore splits can we afford with the cpus we have?
            ccore = ((task.cpus as int) / cpu_per_multicore) as int
            // Check that we have enough memory, assuming 13GB memory per instance (typical for mouse alignment)
            try {
                tmem = (task.memory as nextflow.util.MemoryUnit).toBytes()
                mcore = (tmem / mem_per_multicore) as int
                ccore = Math.min(ccore, mcore)
            } catch (all) {
                log.debug "Warning: Not able to define bismark align multicore based on available memory"
            }
            if( ccore > 1 ){
              multicore = "--multicore $ccore"
            }
        }

        // Main command
        """
        bismark $input \\
            $aligner \\
            --bam $pbat $non_directional $unmapped $mismatches $multicore \\
            --genome $index \\
            $reads \\
            $soft_clipping \\
            $splicesites
        """
    }

    /*
     * STEP 4 - Bismark deduplicate
     */
    if( params.skip_deduplication || params.rrbs ) {
        ch_bam_for_bismark_deduplicate.into { ch_bam_dedup_for_bismark_methXtract; ch_bam_dedup_for_qualimap }
        ch_bismark_dedup_log_for_bismark_report = Channel.from(false)
        ch_bismark_dedup_log_for_bismark_summary = Channel.from(false)
        ch_bismark_dedup_log_for_multiqc  = Channel.from(false)
    } else {
        process bismark_deduplicate {
            tag "$name"
            publishDir "${params.outdir}/bismark_deduplicated", mode: 'copy',
                saveAs: {filename -> filename.indexOf(".bam") == -1 ? "logs/$filename" : "$filename"}

            input:
            set val(name), file(bam) from ch_bam_for_bismark_deduplicate

            output:
            set val(name), file("*.deduplicated.bam") into ch_bam_dedup_for_bismark_methXtract, ch_bam_dedup_for_qualimap
            set val(name), file("*.deduplication_report.txt") into ch_bismark_dedup_log_for_bismark_report, ch_bismark_dedup_log_for_bismark_summary, ch_bismark_dedup_log_for_multiqc

            script:
            fq_type = params.single_end ? '-s' : '-p'
            """
            deduplicate_bismark $fq_type --bam $bam
            """
        }
    }

    /*
     * STEP 5 - Bismark methylation extraction
     */
    process bismark_methXtract {
        tag "$name"
        publishDir "${params.outdir}/bismark_methylation_calls", mode: 'copy',
            saveAs: {filename ->
                if( filename.indexOf("splitting_report.txt" ) > 0 ) "logs/$filename"
                else if( filename.indexOf("M-bias" ) > 0) "m-bias/$filename"
                else if( filename.indexOf(".cov" ) > 0 ) "methylation_coverage/$filename"
                else if( filename.indexOf("bedGraph" ) > 0 ) "bedGraph/$filename"
                else if( filename.indexOf("CpG_report" ) > 0 ) "stranded_CpG_report/$filename"
                else "methylation_calls/$filename"
            }

        input:
        set val(name), file(bam) from ch_bam_dedup_for_bismark_methXtract
        file index from ch_bismark_index_for_bismark_methXtract.collect()

        output:
        set val(name), file("*splitting_report.txt") into ch_bismark_splitting_report_for_bismark_report, ch_bismark_splitting_report_for_bismark_summary, ch_bismark_splitting_report_for_multiqc
        set val(name), file("*.M-bias.txt") into ch_bismark_mbias_for_bismark_report, ch_bismark_mbias_for_bismark_summary, ch_bismark_mbias_for_multiqc
        file '*.{png,gz}'

        script:
        comprehensive = params.comprehensive ? '--comprehensive --merge_non_CpG' : ''
        cytosine_report = params.cytosine_report ? "--cytosine_report --genome_folder ${index} " : ''
        meth_cutoff = params.meth_cutoff ? "--cutoff ${params.meth_cutoff}" : ''
        multicore = ''
        if( task.cpus ){
            // Numbers based on Bismark docs
            ccore = ((task.cpus as int) / 3) as int
            if( ccore > 1 ){
              multicore = "--multicore $ccore"
            }
        }
        buffer = ''
        if( task.memory ){
            mbuffer = (task.memory as nextflow.util.MemoryUnit) - 2.GB
            // only set if we have more than 6GB available
            if( mbuffer.compareTo(4.GB) == 1 ){
              buffer = "--buffer_size ${mbuffer.toGiga()}G"
            }
        }
        if(params.single_end) {
            """
            bismark_methylation_extractor $comprehensive $meth_cutoff \\
                $multicore $buffer $cytosine_report \\
                --bedGraph \\
                --counts \\
                --gzip \\
                -s \\
                --report \\
                $bam
            """
        } else {
            """
            bismark_methylation_extractor $comprehensive $meth_cutoff \\
                $multicore $buffer $cytosine_report \\
                --ignore_r2 2 \\
                --ignore_3prime_r2 2 \\
                --bedGraph \\
                --counts \\
                --gzip \\
                -p \\
                --no_overlap \\
                --report \\
                $bam
            """
        }
    }

    ch_bismark_align_log_for_bismark_report
     .join(ch_bismark_dedup_log_for_bismark_report)
     .join(ch_bismark_splitting_report_for_bismark_report)
     .join(ch_bismark_mbias_for_bismark_report)
     .set{ ch_bismark_logs_for_bismark_report }


    /*
     * STEP 6 - Bismark Sample Report
     */
    process bismark_report {
        tag "$name"
        publishDir "${params.outdir}/bismark_reports", mode: 'copy'

        input:
        set val(name), file(align_log), file(dedup_log), file(splitting_report), file(mbias) from ch_bismark_logs_for_bismark_report

        output:
        file '*{html,txt}' into ch_bismark_reports_results_for_multiqc

        script:
        """
        bismark2report \\
            --alignment_report $align_log \\
            --dedup_report $dedup_log \\
            --splitting_report $splitting_report \\
            --mbias_report $mbias
        """
    }

    /*
     * STEP 7 - Bismark Summary Report
     */
    process bismark_summary {
        publishDir "${params.outdir}/bismark_summary", mode: 'copy'

        input:
        file ('*') from ch_bam_for_bismark_summary.collect()
        file ('*') from ch_bismark_align_log_for_bismark_summary.collect()
        file ('*') from ch_bismark_dedup_log_for_bismark_summary.collect()
        file ('*') from ch_bismark_splitting_report_for_bismark_summary.collect()
        file ('*') from ch_bismark_mbias_for_bismark_summary.collect()

        output:
        file '*{html,txt}' into ch_bismark_summary_results_for_multiqc

        script:
        """
        bismark2summary
        """
    }
} // End of bismark processing block
else {
    ch_bismark_align_log_for_multiqc = Channel.from(false)
    ch_bismark_dedup_log_for_multiqc = Channel.from(false)
    ch_bismark_splitting_report_for_multiqc = Channel.from(false)
    ch_bismark_mbias_for_multiqc = Channel.from(false)
    ch_bismark_reports_results_for_multiqc = Channel.from(false)
    ch_bismark_summary_results_for_multiqc = Channel.from(false)
}


/*
 * STEP 8 - Qualimap
 */
process qualimap {
    tag "$name"
    publishDir "${params.outdir}/qualimap", mode: 'copy'

    input:
    set val(name), file(bam) from ch_bam_dedup_for_qualimap

    output:
    file "${bam.baseName}_qualimap" into ch_qualimap_results_for_multiqc

    script:
    gcref = params.genome.toString().startsWith('GRCh') ? '-gd HUMAN' : ''
    gcref = params.genome.toString().startsWith('GRCm') ? '-gd MOUSE' : ''
    def avail_mem = task.memory ? ((task.memory.toGiga() - 6) / task.cpus).trunc() : false
    def sort_mem = avail_mem && avail_mem > 2 ? "-m ${avail_mem}G" : ''
    """
    samtools sort $bam \\
        -@ ${task.cpus} $sort_mem \\
        -o ${bam.baseName}.sorted.bam
    qualimap bamqc $gcref \\
        -bam ${bam.baseName}.sorted.bam \\
        -outdir ${bam.baseName}_qualimap \\
        --collect-overlap-pairs \\
        --java-mem-size=${task.memory.toGiga()}G \\
        -nt ${task.cpus}
    """
}

/*
 * STEP 9 - preseq
 */
process preseq {
    tag "$name"
    publishDir "${params.outdir}/preseq", mode: 'copy'

    input:
    set val(name), file(bam) from ch_bam_for_preseq

    output:
    file "${bam.baseName}.ccurve.txt" into preseq_results

    script:
    def avail_mem = task.memory ? ((task.memory.toGiga() - 6) / task.cpus).trunc() : false
    def sort_mem = avail_mem && avail_mem > 2 ? "-m ${avail_mem}G" : ''
    """
    samtools sort $bam \\
        -@ ${task.cpus} $sort_mem \\
        -o ${bam.baseName}.sorted.bam
    preseq lc_extrap -v -B ${bam.baseName}.sorted.bam -o ${bam.baseName}.ccurve.txt
    """
}

/*
 * STEP 10 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file (multiqc_config) from ch_multiqc_config
    file (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
    file ('fastqc/*') from ch_fastqc_results_for_multiqc.collect().ifEmpty([])
    file ('trimgalore/*') from ch_trim_galore_results_for_multiqc.collect().ifEmpty([])
    file ('bismark/*') from ch_bismark_align_log_for_multiqc.collect().ifEmpty([])
    file ('bismark/*') from ch_bismark_dedup_log_for_multiqc.collect().ifEmpty([])
    file ('bismark/*') from ch_bismark_splitting_report_for_multiqc.collect().ifEmpty([])
    file ('bismark/*') from ch_bismark_mbias_for_multiqc.collect().ifEmpty([])
    file ('bismark/*') from ch_bismark_reports_results_for_multiqc.collect().ifEmpty([])
    file ('bismark/*') from ch_bismark_summary_results_for_multiqc.collect().ifEmpty([])
    file ('samtools/*') from ch_flagstat_results_for_multiqc.flatten().collect().ifEmpty([])
    file ('samtools/*') from ch_samtools_stats_results_for_multiqc.flatten().collect().ifEmpty([])
    file ('picard/*') from ch_markDups_results_for_multiqc.flatten().collect().ifEmpty([])
    file ('methyldackel/*') from ch_methyldackel_results_for_multiqc.flatten().collect().ifEmpty([])
    file ('qualimap/*') from ch_qualimap_results_for_multiqc.collect().ifEmpty([])
    file ('preseq/*') from preseq_results.collect().ifEmpty([])
    file ('software_versions/*') from ch_software_versions_yaml_for_multiqc.collect()
    file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")

    output:
    file "*multiqc_report.html" into ch_multiqc_report
    file "*_data"
    file "multiqc_plots"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc -f $rtitle $rfilename $custom_config_file . \\
        -m custom_content -m picard -m qualimap -m bismark -m samtools -m preseq -m cutadapt -m fastqc
    """
}

/*
 * STEP 11 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/methylseq] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/methylseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/methylseq] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
                }
        }
    } catch (all) {
        log.warn "[nfcore/methylseq] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/methylseq] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            [ 'mail', '-s', subject, email_address ].execute() << email_txt
            log.info "[nf-core/methylseq] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/methylseq]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/methylseq]${c_red} Pipeline completed with errors${c_reset}-"
    }

}

def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/methylseq v${workflow.manifest.version}${c_reset}
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
