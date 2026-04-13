#!/usr/bin/env nextflow

nextflow.preview.output = true

/*
 * Raichu + pyHICCUPS Multi-Resolution Pipeline
 * 
 * Normalises Hi-C data with raichu and calls loops at 2kb, 5kb, and 10kb
 * resolutions using pyHICCUPS, then combines and converts to arc format.
 */

// Define parameters
params.input = null
params.outdir = 'results'
params.resolutions = '2000,5000,10000'
params.raichu_options = ''
params.pyhiccups_options = ''


// Validate required inputs
if (!params.input) {
    error "ERROR: --input parameter is required (path or glob pattern to input .cool or .mcool files)"
}

// Log parameters
log.info """\
    RAICHU + pyHICCUPS PIPELINE
    ===========================
    input          : ${params.input}
    outdir         : ${params.outdir}
    resolutions    : ${params.resolutions}
    """
    .stripIndent()

/*
 * Process: Normalise Hi-C data using raichu
 */
 /*
process RAICHU_NORMALISE {
    tag "${sample_id}"
    conda '/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/conda_envs/raichu'
    publishDir "${params.outdir}/${sample_id}/raichu", mode: 'copy'


    input:
    tuple val(sample_id), path(cool_file)

    output:
    tuple val(sample_id), path(cool_file), emit: normalised_cool

    script:
    def extra_args = params.raichu_options ?: ''
    """
    # raichu will add the normalisation weights to the same cool file, so we can just pass it through as output

    # Run raichu normalisation for each resulution (if mcool) or single resolution (if cool)
    if [[ -f "${cool_file}" && "${cool_file}" == *.mcool ]]; then
        # If it's an mcool file, we need to specify the resolution in the URI
        resolutions=\$(echo "${params.resolutions}" | tr ',' '\n')
        for res in \$resolutions; do
            echo "Processing resolution \$res for sample ${sample_id}"
            raichu --cool-uri ${cool_file}::/resolutions/\$res --window-size 200 -t 12 -n raichu_weight -f ${extra_args}
        done
    else
        # If it's a cool file, we can only use one resolution (the one specified in the filename or default)
        raichu --cool-uri ${cool_file} --window-size 200 -t 12 -n raichu_weight -f ${extra_args}
    fi   
    """
}
*/

/*
 * Process: Call loops per chromosome at a specific resolution using pyHICCUPS
 */
process PYHICCUPS_CALL_LOOPS_PER_CHROM {
    tag "${sample_id}_pyhiccups_${resolution}_${chromosome}"
    conda '/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/conda_envs/raichu'

    input:
    tuple val(sample_id), path(normalised_cool), val(resolution), val(chromosome)

    output:
    tuple val(sample_id), val(resolution), path("loops_${resolution}_${chromosome}.bedpe"), emit: loops

    script:
    def extra_args = params.pyhiccups_options ?: ''
    def chr_clean = chromosome.replaceFirst(/^chr/, '')
    """
    # Call loops for this chromosome at specified resolution
    # hard coded in as suggested by Wang et al (2026) for microC data
    
    if (( ${resolution} <= 10000 && ${resolution} > 5000 )); then
        max_apart=4000000
    elif (( ${resolution} <= 5000 && ${resolution} > 2000 )); then
        max_apart=2000000
    elif (( ${resolution} <= 2000 )); then
        max_apart=1000000
    else
        #use the default
        max_apart=10000000
    fi

    if [[ "${normalised_cool}" == *.mcool ]]; then
        # If it's an mcool file, we need to specify the resolution in the URI
        cool_uri=${normalised_cool}::/resolutions/${resolution}
    else
        # If it's a cool file, we can only use one resolution (the one specified in the filename or default)
        cool_uri=${normalised_cool}
    fi

    pyHICCUPS -p \$cool_uri --chroms ${chr_clean} -O loops_${resolution}_${chromosome}.bedpe \
            --pw 2 3 4 --ww 5 6 7 --maxww 7 --only-anchors --min-local-reads 100 --nproc ${task.cpus} \
            --clr-weight-name raichu_weight --maxapart \$max_apart --logFile ${sample_id}_${resolution}_${chromosome}_hicpeaks.raichu.log \
            ${extra_args}    
    """
}

/*
 * Process: Merge per-chromosome loops for a given resolution
 */
process MERGE_CHROMOSOME_LOOPS {
    tag "${sample_id}_resolution_${resolution}"
    
    input:
    tuple val(sample_id), val(resolution), path('chrom_loops_*.bedpe')

    output:
    tuple val(sample_id), path("loops_${resolution}.bedpe"), emit: loops

    script:
    """
    cat chrom_loops_*.bedpe > loops_${resolution}.bedpe
    
    """
}

/*
 * Process: Combine loops from all resolutions for a single sample
 */
process COMBINE_LOOPS {
    tag "${sample_id}"
    conda '/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/conda_envs/raichu'

    input:
    tuple val(sample_id), path('loops_*.bedpe')

    output:
    tuple val(sample_id), path('combined_loops.bedpe'), emit: combined

    script:
    def resolutions = params.resolutions.replace(',', ' ')
    def max_res = params.resolutions.split(',').collect { it as Integer }.max()
    """
    # Hard coded parameter values as suggested by Wang et al (2026) for microC data
    combine-resolutions \\
        -O combined_loops.bedpe \\
        -p loops_*.bedpe \\
        -R ${resolutions} \\
        -G ${max_res} \\
        -M 100000 \\
        --max-res ${max_res} \\
        --minapart 8
    """
}

/*
 * Process: Convert BEDPE to arc format for visualisation
 */
process BEDPE_TO_ARC {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/visualisation", mode: 'copy'
    

    input:
    tuple val(sample_id), path(bedpe)

    output:
    tuple val(sample_id), path('loops.arc'), emit: arc

    script:
    """
    #--output-delimiter=\$'\\t' 
    cut -f 1-6,8 ${bedpe} > loops.arc    
    """
}

/*
 * Main workflow
 */
workflow {
    main:

    // Create input channel from multiple cool files
    // Extract sample ID from filename (remove .cool or .mcool extension)
    cool_files_ch = channel.fromPath(params.input, checkIfExists: true)
        .map { file -> 
            def sample_id = file.baseName.replaceAll(/\.(m)?cool$/, '')
            tuple(sample_id, file)
        }
    
    // Parse resolutions
    resolutions_list = params.resolutions.tokenize(',').collect { v -> v.trim() }
    resolutions_ch = channel.of(resolutions_list).flatten()

    //TODO: filechecking
    // resolutions must match those in the mcool file
    // if cool file is given, can only use one resoulution
    // also check named cool or mcool
    
    // Step 1: Normalise with raichu
    // RAICHU_NORMALISE(cool_files_ch)
    
    // Step 2: Create a channel of chromosome names: chr1-chr22, chrX, chrY
    // Create chromosome channel from file for each sample
    chromosomes = channel.of((1..22).collect { n -> "chr${n}" } + ['chrX', 'chrY'] ).flatten()

    // Create sample-chromosome combinations
    chromosomes_per_sample = cool_files_ch
        .map { sample_id, file -> sample_id }
        .combine(chromosomes)
        
    
    // Step 3: Call loops per chromosome per resolution (parallel execution)
    // Create all combinations: sample_id, resolution, chromosome
    // Then combine with the normalised cool file
    res_chrom_combinations = chromosomes_per_sample
        .combine(resolutions_ch)
        .map { sample_id, chrom, res -> tuple(sample_id, chrom, res) }
        //.combine(RAICHU_NORMALISE.out.normalised_cool, by: 0)
        .combine(cool_files_ch, by: 0)
        .map { sample_id, chrom, res, cool -> tuple(sample_id, cool, res, chrom) }
    
    PYHICCUPS_CALL_LOOPS_PER_CHROM(res_chrom_combinations)
    
    // Step 4: Group by sample_id and resolution, then merge chromosome results
    loops_by_sample_resolution = PYHICCUPS_CALL_LOOPS_PER_CHROM.out.loops
        .groupTuple(by: [0, 1])
    
    MERGE_CHROMOSOME_LOOPS(loops_by_sample_resolution)
    
    // Step 5: Combine all resolutions per sample
    loops_by_sample = MERGE_CHROMOSOME_LOOPS.out.loops
        .groupTuple(by: 0)
    
    COMBINE_LOOPS(loops_by_sample)
    
    // Step 6: Convert to arc format
    BEDPE_TO_ARC(COMBINE_LOOPS.out.combined)
    
}

workflow.onComplete {
    log.info """\
        Pipeline completed!
        Status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
        Results: ${params.outdir}
        """
        .stripIndent()
}
