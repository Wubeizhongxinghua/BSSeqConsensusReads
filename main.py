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
import tempfile
import shutil
import click
import subprocess
import yaml
from pathlib import Path
from typing import Dict, Any

def load_config(config_path: str) -> Dict[str, Any]:
	with open(config_path) as f:
		return yaml.safe_load(f)

def build_command(base_cmd: str, rule_config: Dict[str, Any], rule_name: str) -> str:
	param_mapping = {
		'call_consensus_reads_molecular': {
			'xmx': '-Xmx{}',
			'error_rate_pre_umi': '--error-rate-pre-umi={}',
			'min_consensus_base_quality': '--min-consensus-base-quality={}',
			'consensus_call_overlapping_bases': '--consensus-call-overlapping-bases={}',
			'error_rate_post_umi': '--error-rate-post-umi={}',
			'min_input_base_quality': '--min-input-base-quality={}',
			'min_reads': '--min-reads={}',
			'threads': '--threads={}'
		},
		'duplex_metrics': {
			'threads': '--threads={}'
		},
		'consensus_to_fq_unfiltered': {
			'compression_level': 'COMPRESSION_LEVEL={}',
			'include_non_pf_reads': 'INCLUDE_NON_PF_READS={}',
			'threads': '--threads={}'
		},
		'align_consensus_unfiltered': {
			'aligner_options': '',
			'samtools_options': '',
			'threads': '-t {}'
		},
		'mergeAunA_consensus': {
			'xmx': '-Xmx{}',
			'sort_memory': '-m {}',
			'validation_stringency': '--validation-stringency={}',
			'threads': '-@ {}'
		},
		'mergeAunA_consensus_grepaligned': {
			'samtools_flags': '',
			'threads': '-@ {}'
		},
		'convert_Bstrain': {
			'threads': '--threads={}'
		},
		'extend': {
			'extend_length': '-e {}',
			'threads': '--threads={}'
		},
		'groupsort_convert': {
			'xmx': '-Xmx{}',
			'sort_order': '-s {}',
			'tmp_dir': '--tmp-dir={}',
			'threads': '-@ {}'
		},
		'callduplex': {
			'xmx': '-Xmx{}',
			'error_rate_pre_umi': '--error-rate-pre-umi={}',
			'error_rate_post_umi': '--error-rate-post-umi={}',
			'min_input_base_quality': '--min-input-base-quality={}',
			'min_reads': '--min-reads={}',
			'consensus_overlap_mode': '--consensus-overlap-mode={}',
			'threads': '--threads={}'
		},
		'align_consensus_unfiltered_duplex': {
			'aligner_options': '',
			'samtools_options': '',
			'threads': '-t {}'
		}
	}

	cmd = [base_cmd]
	mapping = param_mapping.get(rule_name, {})
	
	for param, fmt in mapping.items():
		if param in rule_config:
			value = rule_config[param]
			if isinstance(value, bool):
				value = str(value).lower()
			cmd.append(fmt.format(value))
	
	# 处理自定义参数
	for param in rule_config.get('others', []):
		for k, v in param.items():
			cmd.append(f"--{k}={v}")
	
	return ' '.join(cmd)

@click.command()
@click.option("--input-bam", required=True, help="Input grouped BAM file")
@click.option("--output-bam", required=True, help="Final output BAM file")
@click.option("--dataset", required=True, help="Dataset name")
@click.option("--sample-id", required=True, help="Sample ID")
@click.option("--config-file", default="config.yaml", help="Config file path")
@click.option("--tmp-dir", help="Temporary directory for intermediates")
@click.option("--keep-tmp", is_flag=True, help="Keep temporary files")
def main(input_bam, output_bam, dataset, sample_id, config_file, tmp_dir, keep_tmp):
	config = load_config(config_file)
	global_cfg = config['global']
	rules_cfg = config['rules']
  
	genome_dir = global_cfg['genome_dir']
	genome_fasta = global_cfg['genome_fasta']

	# 创建目录结构
	base_dir = Path("output") / dataset
	bwameth_dir = base_dir / "bwameth_results"
	bam_rmdup_dir = bwameth_dir / "bam_rmdup_fgbio"
	log_dir = base_dir / "log" / "bwameth_results"
	
	for d in [bwameth_dir, bam_rmdup_dir, log_dir]:
		d.mkdir(parents=True, exist_ok=True)

	# 临时目录处理
	temp_dir_obj = None
	if not tmp_dir:
		temp_dir_obj = tempfile.TemporaryDirectory()
		tmp_dir = temp_dir_obj.name
	else:
		Path(tmp_dir).mkdir(parents=True, exist_ok=True)

	try:
		# 定义中间文件路径
		intermediates = {
			'consensus': bwameth_dir / f"{sample_id}_unalignedConsensus_molecular.bam",
			'duplex_metrics': bwameth_dir / f"{sample_id}_duplex_metrics.txt",
			'fq_unfiltered_1': bam_rmdup_dir / f"{sample_id}_unalignedConsensus_unfiltered_1.fq.gz",
			'fq_unfiltered_2': bam_rmdup_dir / f"{sample_id}_unalignedConsensus_unfiltered_2.fq.gz",
			'aligned_unfiltered': bam_rmdup_dir / f"{sample_id}_consensus_unfiltered.bam",
			'merged_auna': bam_rmdup_dir / f"{sample_id}_consensus_unfiltered_aunamerged.bam",
			'grep_aligned': bam_rmdup_dir / f"{sample_id}_consensus_unfiltered_aunamerged_aligned.bam",
			'converted': bam_rmdup_dir / f"{sample_id}_consensus_unfiltered_aunamerged_converted.bam",
			'extended': bam_rmdup_dir / f"{sample_id}_consensus_unfiltered_aunamerged_converted_extended.bam",
			'groupsorted': bam_rmdup_dir / f"{sample_id}_consensus_unfiltered_aunamerged_converted_extended_groupsort.bam",
			'duplex_consensus': bam_rmdup_dir / f"{sample_id}_consensus_unfiltered_aunamerged_converted_extended_duplexconsensus.bam",
			'duplex_fq_1': bam_rmdup_dir / f"{sample_id}_unalignedConsensus_duplex_1.fq.gz",
			'duplex_fq_2': bam_rmdup_dir / f"{sample_id}_unalignedConsensus_duplex_2.fq.gz"
		}

		# 步骤1: CallMolecularConsensusReads
		cmd = build_command(
			base_cmd=f"{global_cfg['fgbio_path']} CallMolecularConsensusReads",
			rule_config=rules_cfg['call_consensus_reads_molecular'],
			rule_name='call_consensus_reads_molecular'
		)
		cmd += f" --input={input_bam} --output={intermediates['consensus']}"
		subprocess.run(cmd, shell=True, check=True)

		# 步骤2: SamToFastq (未过滤)
		cmd = build_command(
			base_cmd=f"java -jar {global_cfg['picard_path']} SamToFastq",
			rule_config=rules_cfg['consensus_to_fq_unfiltered'],
			rule_name='consensus_to_fq_unfiltered'
		)
		cmd += f" I={intermediates['consensus']}"
		cmd += f" F={intermediates['fq_unfiltered_1']}"
		cmd += f" F2={intermediates['fq_unfiltered_2']}"
		subprocess.run(cmd, shell=True, check=True)

		# 步骤3: 比对未过滤的consensus
		log_file = log_dir / f"{sample_id}_consensus_unfiltered.log"
		cmd = build_command(
			base_cmd=f"{global_cfg['bwameth_path']} --reference {genome_dir}/{genome_fasta}",
			rule_config=rules_cfg['align_consensus_unfiltered'],
			rule_name='align_consensus_unfiltered'
		)
		cmd += f" {intermediates['fq_unfiltered_1']} {intermediates['fq_unfiltered_2']}"
		cmd += f" 2> {log_file} | {global_cfg['samtools_path']} view -h -b -o {intermediates['aligned_unfiltered']}"
		subprocess.run(cmd, shell=True, check=True)

		# 步骤4: 合并比对结果
		cmd = build_command(
			base_cmd=f"{global_cfg['samtools_path']} sort",
			rule_config=rules_cfg['mergeAunA_consensus'],
			rule_name='mergeAunA_consensus'
		)
		cmd += f" {intermediates['aligned_unfiltered']} | "
		cmd += build_command(
			base_cmd=f"{global_cfg['fgbio_path']} ZipperBams",
			rule_config=rules_cfg['mergeAunA_consensus'],
			rule_name='mergeAunA_consensus'
		)
		cmd += f" --unmapped {intermediates['consensus']} --ref {genome_dir}/{genome_fasta}"
		cmd += f" --output {intermediates['merged_auna']}"
		subprocess.run(cmd, shell=True, check=True)

		# 步骤5: 过滤已比对的reads
		cmd = build_command(
			base_cmd=f"{global_cfg['samtools_path']} view",
			rule_config=rules_cfg['mergeAunA_consensus_grepaligned'],
			rule_name='mergeAunA_consensus_grepaligned'
		)
		cmd += f" -o {intermediates['grep_aligned']} {intermediates['merged_auna']}"
		subprocess.run(cmd, shell=True, check=True)

		# 步骤6: AG到CT转换
		cmd = build_command(
			base_cmd=f"{global_cfg['python_path']} {global_cfg['tools_dir']}/1.convert_AG_to_CT.py --reference {genome_dir}/{genome_fasta}",
			rule_config=rules_cfg['convert_Bstrain'],
			rule_name='convert_Bstrain'
		)
		cmd += f" {intermediates['grep_aligned']} {intermediates['converted']}"
		subprocess.run(cmd, shell=True, check=True)

		# 步骤7: 扩展间隙
		cmd = build_command(
			base_cmd=f"{global_cfg['python_path']} {global_cfg['tools_dir']}/2.extend_gap.py",
			rule_config=rules_cfg['extend'],
			rule_name='extend'
		)
		cmd += f" -i {intermediates['converted']} -o {intermediates['extended']}"
		subprocess.run(cmd, shell=True, check=True)

		# 步骤8: 分组排序
		cmd = build_command(
			base_cmd=f"{global_cfg['fgbio_path']} SortBam",
			rule_config=rules_cfg['groupsort_convert'],
			rule_name='groupsort_convert'
		)
		cmd += f" -i {intermediates['extended']} -o {intermediates['groupsorted']}"
		subprocess.run(cmd, shell=True, check=True)

		# 步骤9: 双链共识分析
		cmd = build_command(
			base_cmd=f"{global_cfg['fgbio_path']} CallDuplexConsensusReads",
			rule_config=rules_cfg['callduplex'],
			rule_name='callduplex'
		)
		cmd += f" --input={intermediates['groupsorted']} --output={intermediates['duplex_consensus']}"
		subprocess.run(cmd, shell=True, check=True)

		# 步骤10: 转换为FastQ
		cmd = build_command(
			base_cmd=f"java -jar {global_cfg['picard_path']} SamToFastq",
			rule_config=rules_cfg['consensus_to_fq_unfiltered'],
			rule_name='consensus_to_fq_unfiltered'
		)
		cmd += f" I={intermediates['duplex_consensus']}"
		cmd += f" F={intermediates['duplex_fq_1']} F2={intermediates['duplex_fq_2']}"
		subprocess.run(cmd, shell=True, check=True)

		# 步骤11: 最终比对
		log_file = log_dir / f"{sample_id}_consensus_duplex_unfiltered.log"
		cmd = build_command(
			base_cmd=f"{global_cfg['bwameth_path']} --reference {genome_dir}/{genome_fasta}",
			rule_config=rules_cfg['align_consensus_unfiltered_duplex'],
			rule_name='align_consensus_unfiltered_duplex'
		)
		cmd += f" {intermediates['duplex_fq_1']} {intermediates['duplex_fq_2']}"
		cmd += f" 2> {log_file} | {global_cfg['samtools_path']} view -h -b -o {output_bam}"
		subprocess.run(cmd, shell=True, check=True)

	finally:
		if temp_dir_obj and not keep_tmp:
			temp_dir_obj.cleanup()

if __name__ == "__main__":
	main()
