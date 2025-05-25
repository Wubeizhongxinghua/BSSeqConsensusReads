# BSSeqConsensusReads

A tool for calling duplex consensus reads generated from BS-seq / EM-seq libraries with 2-side UMIs.

## Prerequisites

- Python 3.6+
- Required tools:
  - `fgbio` (v1.5+)
  - `Picard` (v2.25+)
  - `bwameth` (v0.2+)
  - `samtools` (v1.10+)
- Python packages:
  ```bash
  pip install pyyaml click pysam
  ```

## Installation

```bash
git clone https://github.com/Wubeizhongxinghua/BSSeqConsensusReads.git
cd methylation-pipeline
```

## Configuration

Edit `config.yaml` to customize parameters:
```yaml
global:
  tmp_dir: "tmp"
  picard_path: "/path/to/picard.jar"
  # Tool paths and global settings

rules:
  call_consensus_reads_molecular:
    threads: 20
    xmx: "100g"
    others: [{"extra_param": "value"}]  # Custom arguments
  # Other rule configurations...
```

## Usage

Basic execution:
```bash
python main.py \
  --input-bam input/sample_grouped.bam \
  --output-bam output/sample_final.bam \
  --dataset project_name \
  --sample-id sample123
```

Options:
- `--config-file`: Custom config file (default: config.yaml)
- `--tmp-dir`: Temporary directory for intermediates (auto-cleaned by default)
- `--keep-tmp`: Retain temporary files

## Workflow Overview

1. Generate molecular consensus reads
2. Convert to paired FastQ
3. Initial alignment with bwameth
4. Merge aligned/unaligned reads
5. Duplex consensus processing:
   - B-strain conversion
   - Gap extension
   - Duplex consensus calling
6. Final alignment

## Notes

- **Memory Requirements**: Minimum 100GB RAM recommended for consensus steps
- **Intermediate Files**: Stored in `output/{dataset}/bwameth_results/`
- **Logs**: Available in `output/{dataset}/log/`
- Customize parameters via `others` field in config.yaml for tool-specific flags

