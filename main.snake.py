'''
Copyright (c) 2025-05-26 by LiMingyang, YiLab, Peking University.

Author: Li Mingyang (limingyang200101@gmail.com)

Institute: AAIS, Peking University

File Name: /lustre2/cqyi/myli/DNA_5mC_analysis_pipeline/tools/BSSeqDuplexConsensusReads/main.snake.py

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
'''
from rich import print, pretty
from rich.traceback import install
pretty.install()
install(show_locals=True)

configfile: 'config.yaml'
genome_dir = config['genome_dir']
genome_fasta_file_name = config['genome_fasta_file_name']
tools_dir = config['tools_dir']
tmp = config['tmp']

fgbio = config['fgbio']
java = config['java']
bwameth = config['bwameth']
samtools = config['samtools']
python3 = config['python3']
picard_path = config['picard_path']

config['sample'] = config['bam'].split('/')[-1].replace('.bam','')

rule all:
	input:
		expand("output/{sample_id}_consensus_duplex_unfiltered_bwameth.bam", sample_id = config['sample'])
		#expand("output/{sample_id}_consensus_duplex_unfiltered_bwameth.bam", sample_id = 'test')


rule call_consensus_reads_molecular:
	input:
		bamgrouped = "input/{sample_id}.bam"
	output:
		bam_unaligned_consensus = "output/{sample_id}_unalignedConsensus_molecular.bam"
	threads: 20
	shell:
		"""
		{fgbio} -Xmx100g --tmp-dir={tmp} CallMolecularConsensusReads --input={input.bamgrouped} --output={output.bam_unaligned_consensus} --error-rate-pre-umi=45 --min-consensus-base-quality=0 --consensus-call-overlapping-bases=true --error-rate-post-umi=30 --min-input-base-quality=0 --threads={threads} --min-reads=1
		"""


rule consensus_to_fq_unfiltered:
	input:
		bam_unaligned_consensus = "output/{sample_id}_unalignedConsensus_molecular.bam"
	output:
		fq_unaligned_consensus1 = "output/{sample_id}_unalignedConsensus_unfiltered_1.fq.gz",
		fq_unaligned_consensus2 = "output/{sample_id}_unalignedConsensus_unfiltered_2.fq.gz"
	threads: 5
	shell:
		"""
		{java} -jar {picard_path} SamToFastq I={input.bam_unaligned_consensus} F={output.fq_unaligned_consensus1} F2={output.fq_unaligned_consensus2}
		"""

rule consensus_to_fq:
	input:
		bam_unaligned_consensus = "output/{sample_id}_unalignedConsensus_molecular_filtered.bam"
	output:
		fq_unaligned_consensus1 = "output/{sample_id}_unalignedConsensus_1.fq.gz",
		fq_unaligned_consensus2 = "output/{sample_id}_unalignedConsensus_2.fq.gz"
	threads: 5
	shell:
		"""
		{java} -jar {picard_path} SamToFastq I={input.bam_unaligned_consensus} F={output.fq_unaligned_consensus1} F2={output.fq_unaligned_consensus2}
		"""
#
rule align_consensus_unfiltered:
	input:
		fastq1 = "output/{sample_id}_unalignedConsensus_unfiltered_1.fq.gz",
		fastq2 = "output/{sample_id}_unalignedConsensus_unfiltered_2.fq.gz"
	output:
		bam_consensus = "output/{sample_id}_consensus_unfiltered.bam"
	log:
		"output/log/bwameth_results/{sample_id}_consensus_unfiltered.log"
	threads: 20
	shell:
		r"""
		{bwameth} --reference {genome_dir}/{genome_fasta_file_name} -t {threads} {input.fastq1} {input.fastq2} 2> {log} | {samtools} view -h -b -o {output.bam_consensus}
		"""


rule mergeAunA_consensus:
	input:
		bamaligned = "output/{sample_id}_consensus_unfiltered.bam",
		bamunaligned = "output/{sample_id}_unalignedConsensus_molecular.bam"
	output:
		bammerged = "output/{sample_id}_consensus_unfiltered_aunamerged.bam" 
	threads: 20
	shell:
		"""
		{samtools} sort -@ {threads} --no-PG -n {input.bamaligned} | {fgbio} -Xmx100G --tmp-dir=tmp ZipperBams --unmapped {input.bamunaligned} --ref {genome_dir}/{genome_fasta_file_name} --output {output.bammerged} --sort Coordinate
		"""


rule mergeAunA_consensus_grepaligned:
	input:
		bammerged = "output/{sample_id}_consensus_unfiltered_aunamerged.bam" 
	output:
		bamout = "output/{sample_id}_consensus_unfiltered_aunamerged_aligned.bam"
	threads: 10
	shell:
		"""
		{samtools} view -h -b -F 4 -o {output.bamout} {input.bammerged}
		"""

rule convert_Bstrain:
	input:
		bam_consensus_aligned = "output/{sample_id}_consensus_unfiltered_aunamerged_aligned.bam"
	output:
		bam_consensus_converted = temp("output/{sample_id}_consensus_unfiltered_aunamerged_converted.bam")
	threads: 20
	shell:
		"""
		{python3} {tools_dir}/1.convert_AG_to_CT.py --reference {genome_dir}/{genome_fasta_file_name} {input.bam_consensus_aligned} {output.bam_consensus_converted}
		"""

rule extend:
	input:
		bam_consensus_converted = "output/{sample_id}_consensus_unfiltered_aunamerged_converted.bam"
	output:
		bam_consensus_converted = "output/{sample_id}_consensus_unfiltered_aunamerged_converted_extended.bam"
	threads: 10
	shell:
		"""
		{python3} {tools_dir}/2.extend_gap.py -i {input.bam_consensus_converted} -o {output.bam_consensus_converted}
		"""


rule groupsort_convert:
	input:
		bam_consensus_converted = "output/{sample_id}_consensus_unfiltered_aunamerged_converted_extended.bam"
	output:
		bam_consensus_converted = "output/{sample_id}_consensus_unfiltered_aunamerged_converted_extended_groupsort.bam"
	threads: 20
	shell:
		"""
		{fgbio} -Xmx60G SortBam -s TemplateCoordinate -i {input.bam_consensus_converted} -o {output.bam_consensus_converted}
		"""

rule callduplex:
	input:
		bam_consensus_converted = "output/{sample_id}_consensus_unfiltered_aunamerged_converted_extended_groupsort.bam"
	output:
		bam_consensus_converted = "output/{sample_id}_consensus_unfiltered_aunamerged_converted_extended_duplexconsensus.bam"
	threads: 20
	shell:
		"""
		{fgbio} -Xmx100g --tmp-dir=tmp CallDuplexConsensusReads --input={input.bam_consensus_converted} --output={output.bam_consensus_converted} --consensus-call-overlapping-bases=true --error-rate-pre-umi=45 --error-rate-post-umi=30 --min-input-base-quality=0 --threads={threads} --min-reads=0
		"""


rule consensusduplex_to_fq:
	input:
		bam_unaligned_consensus = "output/{sample_id}_consensus_unfiltered_aunamerged_converted_extended_duplexconsensus.bam"
	output:
		fq_unaligned_consensus1 = "output/{sample_id}_unalignedConsensus_duplex_1.fq.gz",
		fq_unaligned_consensus2 = "output/{sample_id}_unalignedConsensus_duplex_2.fq.gz"
	threads: 5
	shell:
		"""
		{java} -jar {picard_path} SamToFastq I={input.bam_unaligned_consensus} F={output.fq_unaligned_consensus1} F2={output.fq_unaligned_consensus2}
		"""

rule align_consensus_unfiltered_duplex:
	input:
		fastq1 = "output/{sample_id}_unalignedConsensus_duplex_1.fq.gz",
		fastq2 = "output/{sample_id}_unalignedConsensus_duplex_2.fq.gz"
	output:
		bam_consensus = "output/{sample_id}_consensus_duplex_unfiltered_bwameth.bam"
	threads: 20
	shell:
		r"""
		{bwameth} --reference {genome_dir}/{genome_fasta_file_name} -t {threads} {input.fastq1} {input.fastq2} | {samtools} view -h -b -o {output.bam_consensus}
		"""
