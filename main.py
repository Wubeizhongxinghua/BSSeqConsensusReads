'''
Copyright (c) 2025-05-26 by LiMingyang, YiLab, Peking University.

Author: Li Mingyang (limingyang200101@gmail.com)

Institute: AAIS, Peking University

File Name: BSSeqDuplexConsensusReads/main.py

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

#!/usr/bin/env python3
import os
import click
import yaml
import tempfile
import shutil
from pathlib import Path

class PipelineRunner:
	def __init__(self, config_path):
		self.config = self.load_config(config_path)
		self.temp_dir = None
		self.keep_temp = False

	def load_config(self, config_path):
		with open(config_path) as f:
			return yaml.safe_load(f)

	def resolve_paths(self, dataset, sample_id):
		base = Path("output") / dataset / "bwameth_results"
		return {
			"bamgrouped": f"output/{dataset}/bwameth_results/{sample_id}_umigrouped.bam",
			"final_bam": f"output/{dataset}/bwameth_results/bam_rmdup_fgbio/{sample_id}_consensus_duplex_unfiltered_bwameth.bam",
			"log_dir": Path(f"output/{dataset}/log/bwameth_results")
		}

	def build_command(self, rule, context):
		template = self.config['rules'][rule]['params']
		cmd = [self.config['global'].get('fgbio', 'fgbio')]
		
		# 替换模板变量
		formatted = template.format(**context)
		cmd.extend(formatted.split())
		
		# 处理特殊命令
		if rule == 'consensus_to_fq_unfiltered':
			cmd = [self.config['global']['picard']] + cmd
		elif rule == 'align_consensus_unfiltered':
			cmd = [self.config['global']['bwameth']] + cmd
		elif rule == 'mergeAunA_consensus_grepaligned':
			cmd = [self.config['global']['samtools']] + cmd
		elif rule in ['convert_Bstrain', 'extend']:
			cmd = [self.config['global']['python']] + cmd
			
		return ' '.join(cmd)

	def execute_step(self, step_name, context):
		cmd = self.build_command(step_name, context)
		print(f"Executing {step_name}:\n{cmd}")
		os.system(cmd)

	def run(self, params):
		# 创建目录结构
		Path(params['log_dir']).mkdir(parents=True, exist_ok=True)
		
		# 设置临时目录
		if self.config['global']['tmp_dir'] == 'auto_temp' and not self.keep_temp:
			self.temp_dir = tempfile.TemporaryDirectory()
			tmp_path = self.temp_dir.name
		else:
			tmp_path = self.config['global']['tmp_dir']
		
		try:
			# 定义执行上下文
			ctx = {
				'tmp_dir': tmp_path,
				'genome_dir': params['genome_dir'],
				'genome_fasta': params['genome_fasta'],
				'dataset': params['dataset'],
				'sample_id': params['sample_id']
			}

			# 执行所有步骤
			steps = [
				('call_consensus_reads_molecular', 
				 {'input': params['input_bam'], 
				  'output': f"output/{params['dataset']}/bwameth_results/{params['sample_id']}_unalignedConsensus_molecular.bam"}),

				('consensus_to_fq_unfiltered',
				 {'input': f"output/{params['dataset']}/bwameth_results/{params['sample_id']}_unalignedConsensus_molecular.bam",
				  'output1': f"output/{params['dataset']}/bwameth_results/bam_rmdup_fgbio/{params['sample_id']}_unalignedConsensus_unfiltered_1.fq.gz",
				  'output2': f"output/{params['dataset']}/bwameth_results/bam_rmdup_fgbio/{params['sample_id']}_unalignedConsensus_unfiltered_2.fq.gz"}),

				('align_consensus_unfiltered',
				 {'input1': f"output/{params['dataset']}/bwameth_results/bam_rmdup_fgbio/{params['sample_id']}_unalignedConsensus_unfiltered_1.fq.gz",
				  'input2': f"output/{params['dataset']}/bwameth_results/bam_rmdup_fgbio/{params['sample_id']}_unalignedConsensus_unfiltered_2.fq.gz",
				  'output': f"output/{params['dataset']}/bwameth_results/bam_rmdup_fgbio/{params['sample_id']}_consensus_unfiltered.bam",
				  'log': f"output/{params['dataset']}/log/bwameth_results/{params['sample_id']}_consensus_unfiltered.log"}),

				('mergeAunA_consensus',
				 {'input': f"output/{params['dataset']}/bwameth_results/bam_rmdup_fgbio/{params['sample_id']}_consensus_unfiltered.bam",
				  'unmapped': f"output/{params['dataset']}/bwameth_results/{params['sample_id']}_unalignedConsensus_molecular.bam",
				  'output': f"output/{params['dataset']}/bwameth_results/bam_rmdup_fgbio/{params['sample_id']}_consensus_unfiltered_aunamerged.bam"}),

				('mergeAunA_consensus_grepaligned',
				 {'input': f"output/{params['dataset']}/bwameth_results/bam_rmdup_fgbio/{params['sample_id']}_consensus_unfiltered_aunamerged.bam",
				  'output': f"output/{params['dataset']}/bwameth_results/bam_rmdup_fgbio/{params['sample_id']}_consensus_unfiltered_aunamerged_aligned.bam"}),

				('convert_Bstrain',
				 {'input': f"output/{params['dataset']}/bwameth_results/bam_rmdup_fgbio/{params['sample_id']}_consensus_unfiltered_aunamerged_aligned.bam",
				  'output': f"output/{params['dataset']}/bwameth_results/bam_rmdup_fgbio/{params['sample_id']}_consensus_unfiltered_aunamerged_converted.bam"}),

				('extend',
				 {'input': f"output/{params['dataset']}/bwameth_results/bam_rmdup_fgbio/{params['sample_id']}_consensus_unfiltered_aunamerged_converted.bam",
				  'output': f"output/{params['dataset']}/bwameth_results/bam_rmdup_fgbio/{params['sample_id']}_consensus_unfiltered_aunamerged_converted_extended.bam"}),

				('groupsort_convert',
				 {'input': f"output/{params['dataset']}/bwameth_results/bam_rmdup_fgbio/{params['sample_id']}_consensus_unfiltered_aunamerged_converted_extended.bam",
				  'output': f"output/{params['dataset']}/bwameth_results/bam_rmdup_fgbio/{params['sample_id']}_consensus_unfiltered_aunamerged_converted_extended_groupsort.bam"}),

				('callduplex',
				 {'input': f"output/{params['dataset']}/bwameth_results/bam_rmdup_fgbio/{params['sample_id']}_consensus_unfiltered_aunamerged_converted_extended_groupsort.bam",
				  'output': f"output/{params['dataset']}/bwameth_results/bam_rmdup_fgbio/{params['sample_id']}_consensus_unfiltered_aunamerged_converted_extended_duplexconsensus.bam"}),

				('align_consensus_unfiltered_duplex',
				 {'input1': f"output/{params['dataset']}/bwameth_results/bam_rmdup_fgbio/{params['sample_id']}_unalignedConsensus_duplex_1.fq.gz",
				  'input2': f"output/{params['dataset']}/bwameth_results/bam_rmdup_fgbio/{params['sample_id']}_unalignedConsensus_duplex_2.fq.gz",
				  'output': params['output_bam'],
				  'log': f"output/{params['dataset']}/log/bwameth_results/{params['sample_id']}_consensus_duplex_unfiltered.log"})
			]

			for step, context in steps:
				self.execute_step(step, {**ctx, **context})

		finally:
			if self.temp_dir and not self.keep_temp:
				self.temp_dir.cleanup()

@click.command()
@click.option("--input-bam", required=True, help="Input grouped BAM file")
@click.option("--output-bam", required=True, help="Final output BAM file")
@click.option("--dataset", required=True, help="Dataset name")
@click.option("--sample-id", required=True, help="Sample ID")
@click.option("--config-file", default="config.yaml", help="Config file path")
@click.option("--tmp-dir", help="Custom temporary directory")
@click.option("--keep-tmp", is_flag=True, help="Keep temporary files")
def main(input_bam, output_bam, dataset, sample_id, config_file, tmp_dir, keep_tmp):
	runner = PipelineRunner(config_file)
	runner.keep_temp = keep_tmp
	
	if tmp_dir:
		runner.config['global']['tmp_dir'] = tmp_dir
	genome_dir = runner.config['global']['genome_dir'] 
	genome_fasta = runner.config['global']['genome_fasta'] 
	runner.run({
		'input_bam': input_bam,
		'output_bam': output_bam,
		'genome_dir': genome_dir,
		'genome_fasta': genome_fasta,
		'dataset': dataset,
		'sample_id': sample_id
	})

if __name__ == "__main__":
	main()
