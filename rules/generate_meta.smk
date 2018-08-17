"""Generate all the meta data files"""


rule create_dict:
	input: reference_file
	output: '{reference_prefix}.dict'
	params:
		picard = config['SOFTWARE']['PICARD'],
		java = config['SOFTWARE']['JAVA'],
		TMPDIR = config['LOCAL']['TMPDIR'],
		sge_opts = "-cwd -j y -l h_data=8G,h_rt=2:00:00 -M "+email+" -m bea -o cluster_logs -e cluster_logs -V"
	shell:
		"""{params.java} -jar -Djava.io.tmpdir={params.TMPDIR} {params.picard} CreateSequenceDictionary\
		REFERENCE={input}\
		OUTPUT={output}
		sleep 500
		"""

rule reduce_gtf:
	input:
		reference_dict=expand('{reference_prefix}.dict', reference_prefix=reference_prefix),
		annotation=annotation_file
	params:
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory = config['LOCAL']['MEMORY'],
		sge_opts = "-cwd -j y -l h_data=8G,h_rt=2:00:00 -M "+email+" -m bea -o cluster_logs -e cluster_logs -V"
	output: '{annotation_prefix}.reduced.gtf'
	shell:
		"""{params.DROPSEQ_wrapper} -m {params.memory} -p ReduceGTF\
		GTF={input.annotation}\
		OUTPUT={output}\
		SEQUENCE_DICTIONARY={input.reference_dict}\
		IGNORE_FUNC_TYPE='null'\
		ENHANCE_GTF='false'
		sleep 500
		"""

rule create_refFlat:
	input:
		annotation=annotation_file,
		reference_dict=expand('{reference_prefix}.dict', reference_prefix=reference_prefix)
	params:
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory = config['LOCAL']['MEMORY'],
		sge_opts = "-cwd -j y -l h_data=8G,h_rt=2:00:00 -M "+email+" -m bea -o cluster_logs -e cluster_logs -V"
	output: '{annotation_prefix}.refFlat'
	shell:
		"""{params.DROPSEQ_wrapper} -m {params.memory} -p ConvertToRefFlat\
		ANNOTATIONS_FILE={input.annotation}\
		OUTPUT={output}\
		SEQUENCE_DICTIONARY={input.reference_dict}
		sleep 500
		"""

rule create_intervals:
	input:
		annotation_reduced = expand('{annotation_prefix}.reduced.gtf', annotation_prefix=annotation_prefix),
		reference_dict = expand('{reference_prefix}.dict', reference_prefix=reference_prefix)
	params:
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory = config['LOCAL']['MEMORY'],
		reference_folder = config['META']['reference_folder'],
		reference_prefix = config['META']['reference_file'].split('.fasta')[0],
		sge_opts = "-cwd -j y -l h_data=8G,h_rt=2:00:00 -M "+email+" -m bea -o cluster_logs -e cluster_logs -V"
	output:
		'{reference_prefix}.rRNA.intervals'
	shell:
		"""{params.DROPSEQ_wrapper} -m {params.memory} -p CreateIntervalsFiles\
		REDUCED_GTF={input.annotation_reduced}\
		SEQUENCE_DICTIONARY={input.reference_dict}\
		O={params.reference_folder}\
		PREFIX={params.reference_prefix}
		sleep 500
		"""
def get_sjdbOverhang(wildcards):
	return(int(wildcards.read_length)-1)

rule create_star_index:
	input:
		reference_file = reference_file,
		annotation_file = annotation_file
	params:
		STAR = config['SOFTWARE']['STAR'],
		STAR_FOLDER = config['META']['reference_folder'],
		star_prefix = star_index_prefix,
		read_len = read_lengths,
		sjdbOverhang = lambda wildcards: get_sjdbOverhang(wildcards),
		sge_opts="-cwd -j y -pe shared 4 -l h_data=10G,h_rt=2:00:00 -M "+email+" -m bea -o cluster_logs -e cluster_logs -V"
	output:
		# '{star_index_prefix}_{read_length}/SA'
		'{star_index_prefix}_{read_length}/'
	shell:
		"""
		{params.STAR}\
		--runThreadN 4\
		--runMode genomeGenerate\
		--genomeDir {output}\
		--genomeFastaFiles {input.reference_file}\
		--sjdbGTFfile {input.annotation_file}\
		--limitGenomeGenerateRAM 30000000000\
		--sjdbOverhang {params.sjdbOverhang}
		sleep 500
		"""