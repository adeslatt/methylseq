
/*
 * This works today on April 19, 2020
 *
 * uses channel setup with the word "fromPath"
 *
 */
def helpMessage() {
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf --reads '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:

      --zenodo_reference [int]          Integer for the Zenodo DOI reference for the reference genome

      --aligner [str]                   Alignment tool to use (default: bismark)
                                            Available: bismark, bismark_hisat, bwameth
      --reads [file]                    Path to input data (must be surrounded with quotes)
      -profile [str]                    Configuration profile to use. Can use multiple (comma separated)
                                            Available: conda, docker, singularity, test, awsbatch, <institute> and more
    Options:
     --single_end [bool]                Specifies that the input is single end reads
     --comprehensive [bool]             Output information for all cytosine contexts
     --cytosine_report [bool]           Output stranded cytosine report during Bismark's bismark_methylation_extractor step.
     --meth_cutoff [int]                Specify a minimum read coverage to report a methylation call during Bismark's bismark_methylation_extractor step.
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
      --bismark_index [path]            Path to Bismark index
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
     --email [email]                    Receive a summary e-mail with details of the run sent to you when the workflow exits
     --email_on_fail [email]            Same as --email, except only send mail if the workflow is not successful
     --max_multiqc_email_size [str]     Threshold size for MultiQC report emailed (Default: 25MB)
     -name [str]                        Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

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
custom_runName       = params.name
zenodo_doi           = params.zenodo_reference
gate_keeper_input    = file("$baseDir/data/done.txt")
chmod_markdown_py    = file("$baseDir/bin/markdown_to_html.py")
chmod_scrape_py      = file("$baseDir/bin/scrape_software_versions.py")
params.bismark_index = file("$baseDir/data/BismarkIndex")

/*--------------------------------------------------
  Channel setup
---------------------------------------------------*/

// set the input channels

ch_multiqc_config        = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs           = file("$baseDir/docs/output.md", checkIfExists: true)

ch_flagstat_results_for_multiqc         = Channel.from(false)
ch_samtools_stats_results_for_multiqc   = Channel.from(false)
ch_markDups_results_for_multiqc         = Channel.from(false)
ch_methyldackel_results_for_multiqc     = Channel.from(false)


//ch_first                                = Channel.fromPath(gate_keeper_input)
//ch_second                               = Channel.fromPath(gate_keeper_input)

ch_bismark_index_for_bismark_align      = Channel.fromPath ( params.bismark_index )
ch_bismark_index_for_bismark_methXtract = Channel.fromPath ( params.bismark_index )

//ch_read_files_for_fastqc                = Channel.fromPath(params.reads)
ch_read_files_for_summary               = Channel.fromPath(params.reads)
//ch_wherearemyfiles_for_trimgalore       = Channel.fromPath(params.reads)
//ch_wherearemyfiles_for_alignment        = Channel.fromPath(params.reads)


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


process getZenodoReference {
    publishDir "data", mode: 'copy'

    input:

    output:
        file "$gate_keeper_input" into ch_first, ch_second
    script:

    """
    chmod 775 $chmod_markdown_py
    chmod 775 $chmod_scrape_py
    if test -f $gate_keeper_input; then
        echo "$gate_keeper_input exists - moving to next step."
    else 
        echo "$gate_keeper_input does not exist - retrieving file from zenodo"
        mkdir -p data
        cd data
        mkdir -p BismarkIndex
        cd BismarkIndex
        wget https://zenodo.org/record/$zenodo_doi/files/Bisulfite_Genome.tar.gz
        tar xzvf Bisulfite_Genome.tar.gz
        cd ..
        touch $gate_keeper_input
    fi
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
    file zenodo_reference from ch_first
    file fastq_files from ch_read_files_for_summary
    
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

process getReadsIntoQC {
    input:
    file zenodo_reference from ch_second
    
    output:
    file "${params.reads}" into ch_read_files_for_fastqc, ch_wherearemyfiles_for_trimgalore, ch_wherearemyfiles_for_alignment

    script:
    """
    echo "${params.reads} being processed by fastqc, trimgalore and alignment"
    """
}
