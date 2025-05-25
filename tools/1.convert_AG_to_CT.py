'''
Copyright (c) 2025-02-10 by LiMingyang, YiLab, Peking University.

Author: Li Mingyang (limingyang200101@gmail.com)

Institute: AAIS, Peking University

File Name: /gpfs3/chengqiyi_pkuhpc/limingyang/nipt/output/20240801_placenta/bwameth_results/bam_rmdup/test_convert.py

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
from pysam import bcftools
from tqdm.auto import tqdm

@click.command()
@click.argument('input_bam', type=click.Path(exists=True))
@click.argument('output_bam', type=click.Path())
@click.option('--reference', type=click.Path(exists=True), required=True, help='Reference genome FASTA file')
def main(input_bam, output_bam, reference):
	# Open reference FASTA
	fasta = pysam.FastaFile(reference)

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



	
	with pysam.AlignmentFile(input_bam, 'rb') as in_bam:
		with pysam.AlignmentFile(output_bam, 'wb', template=in_bam) as out_bam:
			for read in tqdm(in_bam):
				if read.flag in {0, 99, 147}:
					out_bam.write(read)
					continue
				elif read.flag in {1, 83, 163}:

					left_add = 0
					right_del = 0
					# process 16, 83, 163
					# Remove reads with insertion or deletion
					if any(op in {1, 2, 5} for op, _ in read.cigartuples):
						continue
					trimmed_seq, trimmed_qual, trimmed_cigar = remove_softclips(
						read.seq, read.qual, read.cigartuples
					)
					# Add 'N' to the front of the sequence
					original_seq = trimmed_seq
					readseqlen_ori = len(original_seq)
					modified_seq = 'N' + original_seq # 前面加了一个N
					left_add = 1
					modified_length = len(modified_seq)
					
					# Adjust POS (0-based)
					new_pos = max(read.pos - 1, 0)
					
					# Adjust CIGAR: prepend a 1M to represent the extra base.
					new_cigar = []
					new_cigar.append((pysam.CMATCH, 1))
					if trimmed_cigar:
						new_cigar.extend(trimmed_cigar)
					else:
						new_cigar.append((pysam.CMATCH, modified_length - 1))
					
					# Fetch reference sequence (new_pos to new_pos + modified_length + 1)
					ref_name = in_bam.get_reference_name(read.reference_id)
					ref_start = new_pos # bam中记录的pos是忽略了softclip的，因此新位置不需要考虑softclip信息。
					ref_end = new_pos + modified_length + 1 
					try:
						ref_seq = fasta.fetch(ref_name, ref_start, ref_end).upper() # ref往前多取了一位
					except:
						ref_seq = 'N' * (modified_length + 1)
					
					# Reverse complement reference if read is reverse strand
					#if read.is_reverse:
					#	ref_seq = reverse_complement(ref_seq)
					
					# Ensure ref_seq is long enough
					if len(ref_seq) < modified_length + 1:
						ref_seq += 'N' * (modified_length + 1 - len(ref_seq))
					
					# Process each base in modified_seq
					modified_list = list(modified_seq)
					modified_list[0] = ref_seq[0]
					i = 0
					while i < modified_length:
						read_base = modified_list[i]
						ref_base = ref_seq[i] if i < len(ref_seq) else 'N'
						
						if read_base == 'A':
							if ref_base == 'A':
								modified_list[i] = 'A'
							elif ref_base == 'G':
								modified_list[i] = 'G'
						elif read_base == 'C':
							if i < len(ref_seq) - 1 and ref_seq[i] == 'C' and ref_seq[i+1] == 'G':
								# In CG context
								if i + 1 < modified_length: # not for the last C.
									next_read_base = modified_list[i+1]
									if next_read_base == 'A':
										modified_list[i] = 'T'
										modified_list[i+1] = 'G'
										i += 1  # Skip next base as processed
							else:
								# Not in CG context
								modified_list[i] = 'T'
						elif read_base == 'G':
							# G not in CG remains G (simplified)
							pass
						elif read_base == 'T':
							# T remains T
							pass
						i += 1
					
					# Update modified_seq
					modified_seq = ''.join(modified_list)
					readseqlen_after = len(modified_list) 
					# Check for trimming
					# 如果read最后一个是C，而下一个base是G，那这个C不知道甲基化状态，就删掉
					extra_ref_base = ref_seq[modified_length] if modified_length < len(ref_seq) else 'N'
					if extra_ref_base == 'G' and modified_seq and modified_seq[-1] == 'C':
						modified_seq = modified_seq[:-1]
						right_del = 1
						# Adjust CIGAR
						if new_cigar:
							last_op, last_len = new_cigar[-1]
							if last_len > 1:
								new_cigar[-1] = (last_op, last_len - 1)
							else:
								new_cigar.pop()
						# Adjust Qual
						if trimmed_qual:
							trimmed_qual = trimmed_qual[:-1]
					
					# Update read properties
					read.seq = modified_seq
					if readseqlen_after - readseqlen_ori == 1: #给第一个新加的碱基加质量值。
						# Here we prepend a dummy quality (e.g. 'I') for the extra base.
						# (Adjust as needed for your quality encoding.)
						read.qual = 'I' + (trimmed_qual if trimmed_qual is not None else '')
					elif readseqlen_after != readseqlen_ori and readseqlen_after - readseqlen_ori != 1:
						raise ValueError(f"Strange lengths, before: {readseqlen_ori}, after: {readseqlen_after}.")
					read.pos = new_pos
					read.cigar = new_cigar
					read.set_tag("RD", right_del, 'i')
					read.set_tag("LA", left_add, 'i')

					
					out_bam.write(read)

if __name__ == '__main__':
	main()
