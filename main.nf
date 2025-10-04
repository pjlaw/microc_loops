#!/usr/bin/env nextflow

nextflow.preview.output = true

/*
 * Raichu + pyHICCUPS Multi-Resolution Pipeline
 * 
 * Normalizes Hi-C data with raichu and calls loops at 2kb, 5kb, and 10kb
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
    error "ERROR: --input parameter is required (path to input .cool or .mcool file)"
}

// Log parameters
log.info """\
    RAICHU + pyHICCUPS PIPELINE
    ===========================
    input       : ${params.input}
    outdir      : ${params.outdir}
    resolutions : ${params.resolutions}
    """
    .stripIndent()

/*
 * Process: Normalize Hi-C data using raichu
 */
process RAICHU_NORMALIZE {
    publishDir "${params.outdir}/raichu", mode: 'copy'
    container 'python:3.10-slim'

    input:
    path cool_file

    output:
    path 'normalized.cool', emit: normalized_cool

    script:
    def extra_args = params.raichu_options ?: ''
    """
    # Install raichu
    pip install --quiet raichu-hic

    # Run raichu normalization
    raichu normalize ${extra_args} ${cool_file} normalized.cool
    """
}

/*
 * Process: Call loops at a specific resolution using pyHICCUPS
 */
process PYHICCUPS_CALL_LOOPS {
    tag "resolution_${resolution}"
    publishDir "${params.outdir}/loops/res_${resolution}", mode: 'copy'
    container 'python:3.10-slim'

    input:
    path normalized_cool
    val resolution

    output:
    path "loops_${resolution}.bedpe", emit: loops

    script:
    def extra_args = params.pyhiccups_options ?: ''
    """
    # Install pyHICCUPS and dependencies
    pip install --quiet pyhiccups cooler numpy

    # Call loops at specified resolution
    pyhiccups call ${extra_args} \
        --resolution ${resolution} \
        --output loops_${resolution}.bedpe \
        ${normalized_cool}
    """
}

/*
 * Process: Combine loops from all resolutions
 */
process COMBINE_LOOPS {
    publishDir "${params.outdir}/combined", mode: 'copy'
    container 'python:3.10-slim'

    input:
    path 'loops_*.bedpe'

    output:
    path 'combined_loops.bedpe', emit: combined

    script:
    """
    #!/usr/bin/env python3
    import sys
    from pathlib import Path

    # Get all input BEDPE files
    bedpe_files = sorted(Path('.').glob('loops_*.bedpe'))
    
    if not bedpe_files:
        print("ERROR: No loop files found", file=sys.stderr)
        sys.exit(1)

    # Combine all loops
    all_loops = []
    header = None
    
    for bedpe_file in bedpe_files:
        with open(bedpe_file) as f:
            lines = f.readlines()
            if lines:
                # Store header from first file
                if header is None and lines[0].startswith('#'):
                    header = lines[0]
                    lines = lines[1:]
                elif lines[0].startswith('#'):
                    lines = lines[1:]
                all_loops.extend(lines)
    
    # Write combined output
    with open('combined_loops.bedpe', 'w') as out:
        if header:
            out.write(header)
        out.writelines(all_loops)
    
    print(f"Combined {len(all_loops)} loops from {len(bedpe_files)} files")
    """
}

/*
 * Process: Convert BEDPE to arc format for visualization
 */
process BEDPE_TO_ARC {
    publishDir "${params.outdir}/visualization", mode: 'copy'
    container 'python:3.10-slim'

    input:
    path bedpe

    output:
    path 'loops.arc', emit: arc

    script:
    """
    #!/usr/bin/env python3
    import sys

    def bedpe_to_arc(bedpe_file, arc_file):
        with open(bedpe_file) as f_in, open(arc_file, 'w') as f_out:
            for line in f_in:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\\t')
                if len(fields) < 6:
                    continue
                
                chrom1, start1, end1 = fields[0], int(fields[1]), int(fields[2])
                chrom2, start2, end2 = fields[3], int(fields[4]), int(fields[5])
                
                # Skip inter-chromosomal loops
                if chrom1 != chrom2:
                    continue
                
                # Calculate midpoints
                mid1 = (start1 + end1) // 2
                mid2 = (start2 + end2) // 2
                
                # Arc format: chr start end score
                # Use the span of the loop as a simple score
                score = abs(mid2 - mid1)
                
                # Write arc (from smaller to larger coordinate)
                if mid1 < mid2:
                    f_out.write(f"{chrom1}\\t{mid1}\\t{mid2}\\t{score}\\n")
                else:
                    f_out.write(f"{chrom1}\\t{mid2}\\t{mid1}\\t{score}\\n")

    bedpe_to_arc('${bedpe}', 'loops.arc')
    print("Conversion to arc format complete")
    """
}

/*
 * Process: Generate summary statistics
 */
process GENERATE_STATS {
    publishDir "${params.outdir}/stats", mode: 'copy'
    container 'python:3.10-slim'

    input:
    path normalized_cool
    path combined_bedpe
    path arc_file

    output:
    path 'pipeline_stats.txt', emit: stats

    script:
    """
    #!/usr/bin/env python3
    import os
    from pathlib import Path

    stats = []
    stats.append("=== Pipeline Statistics ===\\n")
    
    # File sizes
    if Path('${normalized_cool}').exists():
        size = Path('${normalized_cool}').stat().st_size
        stats.append(f"Normalized cool file size: {size:,} bytes\\n")
    
    # Count loops
    if Path('${combined_bedpe}').exists():
        with open('${combined_bedpe}') as f:
            loop_count = sum(1 for line in f if not line.startswith('#'))
        stats.append(f"Total loops called: {loop_count:,}\\n")
    
    # Count arcs
    if Path('${arc_file}').exists():
        with open('${arc_file}') as f:
            arc_count = sum(1 for line in f if line.strip())
        stats.append(f"Intra-chromosomal loops (arcs): {arc_count:,}\\n")
    
    # Write stats
    with open('pipeline_stats.txt', 'w') as f:
        f.writelines(stats)
    
    # Print to stdout
    print(''.join(stats))
    """
}

/*
 * Main workflow
 */
output {
    normalized_data {
        path 'raichu'
    }
    loop_calls {
        path 'loops'
    }
    combined_results {
        path 'combined'
    }
    visualization {
        path 'visualization'
    }
    statistics {
        path 'stats'
    }
}

workflow {
    main:
    // Create input channel
    cool_file_ch = channel.fromPath(params.input, checkIfExists: true)
    
    // Parse resolutions
    resolutions_list = params.resolutions.tokenize(',').collect { it.trim() }
    resolutions_ch = channel.of(resolutions_list).flatten()
    
    // Step 1: Normalize with raichu
    RAICHU_NORMALIZE(cool_file_ch)
    
    // Step 2: Call loops at each resolution
    PYHICCUPS_CALL_LOOPS(
        RAICHU_NORMALIZE.out.normalized_cool,
        resolutions_ch
    )
    
    // Step 3: Combine all loops
    all_loops = PYHICCUPS_CALL_LOOPS.out.loops.collect()
    COMBINE_LOOPS(all_loops)
    
    // Step 4: Convert to arc format
    BEDPE_TO_ARC(COMBINE_LOOPS.out.combined)
    
    // Step 5: Generate statistics
    GENERATE_STATS(
        RAICHU_NORMALIZE.out.normalized_cool,
        COMBINE_LOOPS.out.combined,
        BEDPE_TO_ARC.out.arc
    )

    publish:
    normalized_data = RAICHU_NORMALIZE.out.normalized_cool
    loop_calls = PYHICCUPS_CALL_LOOPS.out.loops
    combined_results = COMBINE_LOOPS.out.combined
    visualization = BEDPE_TO_ARC.out.arc
    statistics = GENERATE_STATS.out.stats
}

workflow.onComplete {
    log.info """\
        Pipeline completed!
        Status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
        Results: ${params.outdir}
        """
        .stripIndent()
}
