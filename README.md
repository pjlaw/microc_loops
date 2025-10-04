# Micro-C Loops Pipeline

A Nextflow pipeline for detecting chromatin loops from Hi-C/Micro-C data using raichu normalization and pyHICCUPS loop calling.

## Overview

This pipeline performs multi-resolution loop detection from Hi-C/Micro-C contact matrices:

1. **Normalization** - Uses [raichu](https://github.com/open2c/raichu) for ICE-like normalization
2. **Loop Calling** - Applies [pyHICCUPS](https://github.com/ParkerLab/pyHICCUPS) at multiple resolutions (default: 2kb, 5kb, 10kb)
3. **Combination** - Merges loops from all resolutions
4. **Format Conversion** - Converts BEDPE to arc format for genome browser visualization
5. **Statistics** - Generates summary statistics

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
├── loops/
│   ├── res_2000/
│   │   └── loops_2000.bedpe    # Loops at 2kb resolution
│   ├── res_5000/
│   │   └── loops_5000.bedpe    # Loops at 5kb resolution
│   └── res_10000/
│       └── loops_10000.bedpe    # Loops at 10kb resolution
├── combined/
│   └── combined_loops.bedpe     # All loops merged
├── visualization/
│   └── loops.arc                # Arc format for genome browsers
├── stats/
│   └── pipeline_stats.txt       # Summary statistics
└── reports/
    ├── execution_report.html    # Nextflow execution report
    ├── timeline.html            # Timeline visualization
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

Simple 4-column format for visualization:

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

| Process | CPUs | Memory | Time |
|---------|------|--------|------|
| RAICHU_NORMALIZE | 4 | 8 GB | 6h |
| PYHICCUPS_CALL_LOOPS | 4 | 8 GB | 8h |
| COMBINE_LOOPS | 1 | 2 GB | 1h |
| BEDPE_TO_ARC | 1 | 2 GB | 1h |
| GENERATE_STATS | 1 | 2 GB | 1h |

Adjust in `nextflow.config` based on your data size.

## Troubleshooting

### Common Issues

1. **Out of memory errors**
   ```bash
   # Increase memory in nextflow.config
   process.memory = '16.GB'
   ```

2. **Long running times**
   ```bash
   # Reduce resolution count or increase CPUs
   nextflow run main.nf --resolutions 5000,10000
   ```

3. **Container pull failures**
   ```bash
   # Pre-pull containers or use conda profile
   nextflow run main.nf -profile conda
   ```

### Logs and Debugging

```bash
# Run with debug output
nextflow run main.nf --input data.cool -with-trace -with-report

# Check work directory for intermediate files
ls -la work/

# View specific task logs
cat .nextflow.log
```

## Citation

If you use this pipeline, please cite:

- **Nextflow**: Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology 35, 316–319.
- **raichu**: [Citation for raichu tool]
- **pyHICCUPS**: [Citation for pyHICCUPS tool]

## License

MIT License - see LICENSE file for details.

## Contact

For issues and questions:
- GitHub Issues: https://github.com/plaw/microc_loops/issues
- Discussions: https://github.com/plaw/microc_loops/discussions

## Acknowledgments

Built with ❤️ using Nextflow and assisted by Seqera AI.
