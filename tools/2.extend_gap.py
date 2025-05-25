'''
Copyright (c) 2025-05-25 by LiMingyang, YiLab, Peking University.

Author: Li Mingyang (limingyang200101@gmail.com)

Institute: AAIS, Peking University

File Name: /lustre2/cqyi/myli/DNA_5mC_analysis_pipeline/output/20250307-XML_umi/bwameth_results/bam_rmdup_fgbio/2.extend_gap.py

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

import rich_click as click
import pysam
from collections import defaultdict
from typing import List, Tuple, Dict
from tqdm import tqdm
def remove_softclips(seq, qual, cigartuples):
	"""
	Remove soft-clipped bases from the beginning and end of the read.
	Returns updated (seq, qual, cigartuples).
	"""
	if cigartuples is None:
		return seq, qual, cigartuples
	cigars = list(cigartuples)
	# Remove leading softclip (CIGAR op code 4)
	if cigars and cigars[0][0] == 4:
		sc_len = cigars[0][1]
		seq = seq[sc_len:]
		if qual is not None:
			qual = qual[sc_len:]
		cigars = cigars[1:]
	# Remove trailing softclip (CIGAR op code 4)
	if cigars and cigars[-1][0] == 4:
		sc_len = cigars[-1][1]
		seq = seq[:-sc_len]
		if qual is not None:
			qual = qual[:-sc_len]
		cigars = cigars[:-1]
	return seq, qual, cigars

def has_hardclips(read) -> bool:
	"""Check if read contains hardclips in its CIGAR string."""
	return any(op == 5 for op, _ in read.cigartuples) if read.cigartuples else False

def process_read_pair(read1: pysam.AlignedSegment, read2: pysam.AlignedSegment) -> Tuple[pysam.AlignedSegment, pysam.AlignedSegment]:
	"""Process a pair of reads to align their positions."""
	# Determine which read is leftmost
	if read1.flag in {83, 163}: # start, 83 163 is in left.
		left_read, right_read = read1, read2
	else:
		left_read, right_read = read2, read1
	quals_left = left_read.qual if left_read.qual else None
	quals_right = right_read.qual if right_read.qual else None
	
	# Handle start position adjustment
	# if right_read.reference_start == left_read.reference_start + 1:
	if int(left_read.get_tag('LA')) == 1:
		# Extract the first base and quality from left read
		first_base = left_read.query_sequence[0]
		first_qual = quals_left[0] if quals_left else None
		
		# Update right read
		right_read.query_sequence = first_base + right_read.query_sequence
		if quals_right:
			quals_right = first_qual + quals_right
		right_read.reference_start -= 1
		right_read.cigartuples = [(0, 1)] + list(right_read.cigartuples)

	elif int(left_read.get_tag('LA')) == 0:
		pass
	elif left_read.flag in {163} and right_read.flag in {99} and int(left_read.get_tag('LA')) not in {0,1}: #对于这两个reads来说开头对齐重要
		raise ValueError(f"{right_read.query_name} with flag {right_read.flag} does not match the start position. Right: {right_read.reference_start}, left: {left_read.reference_start}")
	
	# Handle end position adjustment
	left_end = left_read.reference_end if left_read.reference_end else float('inf')
	right_end = right_read.reference_end if right_read.reference_end else float('inf')
	
	# if right_end == left_end + 1:
	if int(left_read.get_tag('RD')) == 1:
		# Extract the last base and quality from left read
		last_base = right_read.query_sequence[-1]
		last_qual = quals_right[-1] if quals_right else None
		
		# Update right read
		left_read.query_sequence = left_read.query_sequence + last_base
		if quals_left:
			quals_left = quals_left + last_qual
		left_read.cigartuples = list(left_read.cigartuples) + [(0, 1)]
		
	elif int(left_read.get_tag('RD')) == 0:
		pass
	elif left_read.flag in {83} and right_read.flag in {147} and int(left_read.get_tag('RD')) not in {0,1}: #对于这两个reads来说开头对齐重要
		raise ValueError(f"{right_read.query_name} with flag {right_read.flag} does not match the end position. Right: {right_end}, left: {left_end}")
		
	right_read.qual = quals_right
	left_read.qual = quals_left
	return left_read, right_read

def process_read_group(reads: List[pysam.AlignedSegment]) -> List[pysam.AlignedSegment]:
	"""Process a group of reads (should be 4 reads if valid)."""
	if len(reads) != 4:
		return reads
	
	# Group reads by their flags
	flag_groups = defaultdict(list)
	for read in reads:
		flag_groups[read.flag].append(read)
	
	# Process first pair (99, 163)
	if 99 in flag_groups and 163 in flag_groups:
		flag_groups[99][0], flag_groups[163][0] = process_read_pair(
			flag_groups[99][0], flag_groups[163][0]
		)
	
	# Process second pair (83, 147)
	if 83 in flag_groups and 147 in flag_groups:
		flag_groups[83][0], flag_groups[147][0] = process_read_pair(
			flag_groups[83][0], flag_groups[147][0]
		)
	
	# Combine all processed reads
	processed_reads = []
	for flag in [99, 163, 83, 147]:
		if flag in flag_groups:
			processed_reads.extend(flag_groups[flag])
	
	return processed_reads

@click.command()
@click.option('-i', '--input-bam', required=True, help='Input BAM file path')
@click.option('-o', '--output-bam', required=True, help='Output BAM file path')
def main(input_bam: str, output_bam: str):
	"""Process BAM file to align read pairs and handle soft/hard clips."""
	
	# Open input BAM file
	inbam = pysam.AlignmentFile(input_bam, "rb")
	
	# Create output BAM file with same header
	outbam = pysam.AlignmentFile(output_bam, "wb", header=inbam.header)
	
	# Group reads by MI tag
	read_groups: Dict[str, List[pysam.AlignedSegment]] = defaultdict(list)
	
	# First pass: group reads and remove those with hardclips
	for read in tqdm(inbam, desc="Grouping input bam"):
		# Skip reads with hardclips
		if has_hardclips(read):
			continue
			
		# Get MI tag
		mi_tag = read.get_tag("MI") if read.has_tag("MI") else None
		if mi_tag:
			mi_tag = mi_tag.split('/')[0]
			# Remove softclips if present
			if read.cigartuples and any(op == 4 for op, _ in read.cigartuples):
				new_seq, new_qual, new_cigars = remove_softclips(
					read.query_sequence,
					read.query_qualities,
					read.cigartuples
				)
				read.query_sequence = new_seq
				read.query_qualities = new_qual
				read.cigartuples = new_cigars
			
			read_groups[mi_tag].append(read)
		else:
			raise ValueError(f"{read.query_name} does not have MI tag.")
	
	# Process each group and write to output
	for mi_tag, reads in tqdm(read_groups.items(), desc="Processing groups"):
		processed_reads = process_read_group(reads)
		for read in processed_reads:
			outbam.write(read)
	
	# Close files
	inbam.close()
	outbam.close()

if __name__ == '__main__':
	main()
