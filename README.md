# BSSeqConsensusReads

A tool for calling duplex consensus reads generated from BS-seq / EM-seq libraries with 2-side UMIs.

## Input & Output Information

The input bam file is the output of `fgbio GroupReadsByUmi` command under parameter `-s Paired`.

The output bam file is like the output of `fgbio CallDuplexConsensusReads --min-reads=0 --consensus-call-overlapping-bases=true`, which is **not filtered**.

## Prerequisites

- Python 3.6+
- Required tools:
  - `fgbio` (v1.5+)
  - `Picard` (v2.25+)
  - `bwameth` (v0.2+)
  - `samtools` (v1.10+)
- Python packages:
  ```bash
  pip install pyyaml click pysam rich rich_click
  ```

## Installation

```bash
git clone https://github.com/Wubeizhongxinghua/BSSeqConsensusReads.git
cd BSSeqDuplexConsensusReads
```

## Configuration

1. Edit `config.yaml`, make sure all the parameters are set.
```yaml
genome_dir: "/path/to/genome_dir"
genome_fasta_file_name: "genome.fa"
tools_dir: "./tools"
tmp: "/path/to/tmp"

fgbio: "/path/to/fgbio"
java: "/path/to/java"
bwameth: "/path/to/bwameth.py"
samtools: "/path/to/samtools"
python3: "/path/to/python3"
picard_path: "/path/to/picard.jar"
```

2. Move or softlink your input bam into ./input directory.

> [!IMPORTANT] The input bam requires:
> - Input bam must be the output of `fgbio GroupReadsByUmi`, all prerequisite steps follow standard fgbio pipeline.
> - The RX tag stands for UMI information.
> - The MI tag stands for UMI group information, requiring `-s Paired` of `fgbio GroupReadsByUmi`.


## Usage

Basic execution:
```bash
snakemake -s main.snake.py -c 20 --rerun-incomplete -pr --rerun-triggers mtime --configfile config_self.yaml --config bam="input/test.bam"
```

Options:
- `--configfile`: Custom config file (default: config.yaml).
- `--config`: Set your bam file.

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
- **Intermediate Files**: Stored in `output/`

