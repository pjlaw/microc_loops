# Micro-C Loops Pipeline

A Nextflow pipeline for detecting chromatin loops from Hi-C/Micro-C data using raichu normalisation and pyHICCUPS loop calling.

## Overview

This pipeline performs multi-resolution loop detection from Hi-C/Micro-C contact matrices:

1. **normalisation** - Uses [Raichu](https://github.com/XiaoTaoWang/Raichu) for normalisation
2. **Loop Calling** - Applies [pyHICCUPS as implemented in HiCPeaks](https://github.com/XiaoTaoWang/HiCPeaks) per chromosome at multiple resolutions (default: 2kb, 5kb, 10kb) for maximum parallelisation
3. **Merging** - Combines per-chromosome results for each resolution
4. **Combination** - Merges loops from all resolutions
5. **Format Conversion** - Converts BEDPE to arc format for genome browser visualisation (e.g. pygenometracks)


### Parallelisation

The pipeline parallelizes loop calling by **chromosome and resolution**, meaning if you have 23 chromosomes and 3 resolutions, it will run up to **69 parallel tasks** (23 × 3), significantly reducing runtime on HPC clusters.

## Quick Start

### Prerequisites

- Nextflow (≥23.04.0)
- Docker, Singularity, or Conda
- Input: `.cool` or `.mcool` file with Hi-C/Micro-C contact data

### Installation

```bash
git clone https://github.com/plaw/microc_loops.git
cd microc_loops
```

### Basic Usage

```bash
# Using Docker
nextflow run main.nf \
    --input data/your_data.cool \
    --outdir results \
    -profile docker

# Using Singularity
nextflow run main.nf \
    --input data/your_data.cool \
    --outdir results \
    -profile singularity

# Using Conda
nextflow run main.nf \
    --input data/your_data.cool \
    --outdir results \
    -profile conda
```

## Parameters

### Required

- `--input`: Path to input `.cool` or `.mcool` file

### Optional

- `--outdir`: Output directory (default: `results`)
- `--resolutions`: Comma-separated resolutions in bp (default: `2000,5000,10000`)
- `--raichu_options`: Additional options for raichu (default: `''`)
- `--pyhiccups_options`: Additional options for pyHICCUPS (default: `''`)

### Examples

```bash
# Single resolution
nextflow run main.nf --input data.cool --resolutions 5000

# Multiple custom resolutions
nextflow run main.nf --input data.cool --resolutions 1000,2500,5000,10000

# With custom options
nextflow run main.nf \
    --input data.cool \
    --raichu_options "--min-nnz 10" \
    --pyhiccups_options "--fdr 0.01"
```

## Output Structure

```
results/
├── raichu/
│   └── normalized.cool          # Normalized contact matrix
├── chromosomes/
│   └── chromosomes.txt          # List of chromosomes
├── loops/
│   ├── res_2000/
│   │   ├── per_chr/            # Per-chromosome BEDPE files
│   │   │   ├── chr1_loops_2000.bedpe
│   │   │   ├── chr2_loops_2000.bedpe
│   │   │   └── ...
│   │   └── loops_2000.bedpe    # Merged loops at 2kb resolution
│   ├── res_5000/
│   │   ├── per_chr/
│   │   └── loops_5000.bedpe    # Merged loops at 5kb resolution
│   └── res_10000/
│   │   ├── per_chr/
│       └── loops_10000.bedpe    # Merged loops at 10kb resolution
├── combined/
│   └── combined_loops.bedpe     # All loops merged
├── visualisation/
│   └── loops.arc                # Arc format for genome browsers
├── stats/
│   └── pipeline_stats.txt       # Summary statistics
└── reports/
    ├── execution_report.html    # Nextflow execution report
    ├── timeline.html            # Timeline visualisation
    ├── trace.txt               # Resource usage trace
    └── pipeline_dag.svg        # Pipeline DAG diagram
```

## File Formats

### BEDPE Format

Standard 6-column BEDPE format for loops:

```
chr1  start1  end1  chr2  start2  end2  [additional_columns]
```

### Arc Format

Simple 4-column format for visualisation:

```
chr  start  end  score
```

Where `start` and `end` are the midpoints of the loop anchors, and `score` is the genomic distance.

## Configuration Profiles

### Available Profiles

- `standard`: Local execution (default)
- `docker`: Use Docker containers
- `singularity`: Use Singularity containers
- `conda`: Use Conda environments
- `slurm`: SLURM cluster execution
- `aws`: AWS Batch execution
- `gcp`: Google Cloud Batch execution
- `test`: Run with test data

### Custom Configuration

Edit `nextflow.config` to customize:

- Resource requirements (CPU, memory, time)
- Executor settings
- Container images
- Tool-specific parameters

### Example: SLURM Cluster

```bash
nextflow run main.nf \
    --input data.cool \
    -profile slurm,singularity \
    -c cluster.config
```

## Resource Requirements

### Default Resources

| Process | CPUs | Memory | Time | Parallelisation |
|---------|------|--------|------|----------------|
| RAICHU_NORMALIZE | 4 | 8 GB | 6h | 1 task |
| GET_CHROMOSOMES | 1 | 2 GB | 1h | 1 task |
| PYHICCUPS_CALL_LOOPS | 4 | 8 GB | 8h | N chroms × M resolutions |
| MERGE_CHROMOSOME_LOOPS | 1 | 2 GB | 1h | M resolutions |
| COMBINE_LOOPS | 1 | 2 GB | 1h | 1 task |
| BEDPE_TO_ARC | 1 | 2 GB | 1h | 1 task |
| GENERATE_STATS | 1 | 2 GB | 1h | 1 task |

**Performance Note**: The pipeline parallelizes loop calling per chromosome, so if you have access to an HPC cluster with many cores, you can process all chromosomes simultaneously, reducing total runtime from days to hours.

Adjust resources in `nextflow.config` based on your data size and available compute.

## Acknowledgments

Built assisted by Seqera AI.
