"""Align the data with STAR."""

ruleorder: plot_knee_plot_whitelist > plot_knee_plot

rule STAR_align:
	input: 
		data="data/{sample}_filtered.fastq.gz",
		# index = expand('{star_index_prefix}_{read_length}/SA', star_index_prefix=star_index_prefix, read_length=read_lengths)
		index = expand('{star_index_prefix}_{read_length}/', star_index_prefix=star_index_prefix, read_length=read_lengths)
	output:
		sam = temp('logs/{sample}.Aligned.out.bam'),
		log_out = 'logs/{sample}.Log.final.out'
	params:
		prefix = '{sample}.',
		STAR = config['SOFTWARE']['STAR'],
		outFilterMismatchNmax = config['STAR_PARAMETERS']['outFilterMismatchNmax'],
		outFilterMismatchNoverLmax = config['STAR_PARAMETERS']['outFilterMismatchNoverLmax'],
		outFilterMismatchNoverReadLmax = config['STAR_PARAMETERS']['outFilterMismatchNoverReadLmax'],
		outFilterMatchNmin = config['STAR_PARAMETERS']['outFilterMatchNmin'],
		index = expand('{star_index_prefix}_{read_length}/', star_index_prefix=star_index_prefix, read_length=read_lengths),
		sge_opts = "-cwd -j y -pe shared 6 -l h_data=10G,h_rt=2:00:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs"
	shell:
		"""{params.STAR}\
		--genomeDir {params.index}\
		--readFilesCommand zcat\
		--runThreadN 6\
		--readFilesIn {input.data}\
		--outSAMtype BAM Unsorted\
		--outReadsUnmapped Fatsx\
		--outFileNamePrefix logs/{params.prefix}\
		--outFilterMismatchNmax {params.outFilterMismatchNmax}\
		--outFilterMismatchNoverLmax {params.outFilterMismatchNoverLmax}\
		--outFilterMismatchNoverReadLmax {params.outFilterMismatchNoverReadLmax}\
		--outFilterMatchNmin {params.outFilterMatchNmin}
		sleep 500
		"""

rule sort_sam:
	input: 'logs/{sample}.Aligned.out.bam'
	params:
		TMPDIR = config['LOCAL']['TMPDIR'],
		picard = config['SOFTWARE']['PICARD'],
		java = config['SOFTWARE']['JAVA'],
		memory = config['LOCAL']['MEMORY'],
		sge_opts = "-cwd -j y -l h_data=16G,h_rt=4:00:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs"
	output: temp('data/{sample}.Aligned.sorted.bam')
	shell:
		"""{params.java} -Xmx{params.memory} -jar -Djava.io.tmpdir={params.TMPDIR}	{params.picard} SortSam\
		INPUT={input}\
		OUTPUT={output}\
		SORT_ORDER=queryname\
		TMP_DIR={params.TMPDIR}
		sleep 500
		"""

rule multiqc_star:
	input: expand('logs/{sample}.Log.final.out', sample = samples.index)
	output: 'reports/star.html'
	params:
		multiqc = config["SOFTWARE"]["MULTIQC"],
		sge_opts="-cwd -j y -l h_data=8G,h_rt=2:00:00 -M "+email+" -m bea -N "+job_label+" -o cluster_logs -e cluster_logs"
	shell:
		"""{params.multiqc} logs/ -m star -o reports/ -n star.html -f
		sleep 500
		"""


rule MergeBamAlignment:
	input:
		unmapped = 'data/{sample}_trimmed_unmapped.bam',
		mapped = 'data/{sample}.Aligned.sorted.bam',
		dict_file = expand('{reference_prefix}.dict', reference_prefix=reference_prefix)
	output: temp('data/{sample}.Aligned.merged.bam')
	params:
		picard = config['SOFTWARE']['PICARD'],
		java = config['SOFTWARE']['JAVA'],
		reference = reference_file,
		TMPDIR = config['LOCAL']['TMPDIR'],
		memory = config['LOCAL']['MEMORY'],
		sge_opts = "-cwd -j y -l h_data=16G,h_rt=4:00:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs"
	shell:
		"""{params.java} -Djava.io.tmpdir={params.TMPDIR} -Xmx{params.memory} -jar {params.picard} MergeBamAlignment\
		REFERENCE_SEQUENCE={params.reference}\
		UNMAPPED_BAM={input.unmapped}\
		ALIGNED_BAM={input.mapped}\
		INCLUDE_SECONDARY_ALIGNMENTS=false\
		PAIRED_RUN=false\
		OUTPUT={output}
		sleep 500
		"""

rule TagReadWithGeneExon:
	input:
		data = 'data/{sample}.Aligned.merged.bam',
		refFlat = expand('{annotation_prefix}.refFlat', annotation_prefix=annotation_prefix)
	params:
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory = config['LOCAL']['MEMORY'],
		TMPDIR = config['LOCAL']['TMPDIR'],
		sge_opts = "-cwd -j y -l h_data=16G,h_rt=4:00:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs"
	output:
		'data/{sample}_gene_exon_tagged.bam'
	shell:
		"""{params.DROPSEQ_wrapper} -t {params.TMPDIR} -m {params.memory} -p TagReadWithGeneExon\
		INPUT={input.data}\
		OUTPUT={output}\
		ANNOTATIONS_FILE={input.refFlat}\
		TAG=GE\
		CREATE_INDEX=true
		sleep 500
		"""

rule dropEST:
	input: 'data/{sample}_gene_exon_tagged.bam'
	params:
		DROPEST=config['SOFTWARE']['DROPEST'],
		GTF_File=config['META']['gtf_file'],
		XML=config['SOFTWARE']['DROPESTXML'],
		outputName='data/{sample}_dropEST',
		logsName='dropEST_out/{sample}_dropEST',
		sge_opts = "-cwd -V -j y -l h_data=16G,h_rt=4:00:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs"
	output:
		'dropEST_out/{sample}.rds'
	shell:
		"""{params.DROPEST} -f -M -L eEBA -V -w\
		-g {params.GTF_File}\
		-c {params.XML}\
		-o {output}\
		-l {params.logsName}\
		{input}
		sleep 500
		"""

rule dropESTreport:
	input: 'dropEST_out/{sample}.rds'
	params:
		RPATH=config['SOFTWARE']['R_PATH'],
		DROPREPORT=config['SOFTWARE']['dropReport'],
		MTGENES=config['META']['mt_genes'],
		sge_opts = "-cwd -V -j y -l h_data=16G,h_rt=4:00:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs"
	output:
		'dropEST_out/{sample}.html'
	shell:
		"""{params.RPATH}\
		{params.DROPREPORT}\
		-m {params.MTGENES}\
		{input}\
		-o {output}
		sleep 500
		"""

rule dropESTfilterCells:
	input: 'dropEST_out/{sample}.rds'
	params:
		mitoChromName=config['META']['mt_chrom'],
		RPATH=config['SOFTWARE']['R_PATH'],
		script=config['SOFTWARE']['dropESTfilterCells'],
		sge_opts = "-cwd -V -j y -l h_data=8G,h_rt=2:00:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs"
	output:
		'dropEST_out/{sample}/matrix.mtx'
	shell:
		"""{params.RPATH}\
		{params.script}\
		{input}\
		{params.mitoChromName}\
		{output}
		sleep 500
		"""

rule bead_errors_metrics:
	input: 'data/{sample}_gene_exon_tagged.bam'
	output:
		temp('data/{sample}_final.bam')
	params:
		out_stats = 'logs/{sample}_synthesis_stats.txt',
		summary = 'logs/{sample}_synthesis_stats_summary.txt',
		barcodes = lambda wildcards: int(samples.loc[wildcards.sample,'expected_cells']) * 2,
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory =config['LOCAL']['MEMORY'],
		SmartAdapter = config['FILTER']['5PrimeSmartAdapter'],
		TMPDIR = config['LOCAL']['TMPDIR'],
		sge_opts = "-cwd -j y -l h_data=16G,h_rt=4:00:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs"
	shell:
		"""{params.DROPSEQ_wrapper} -t {params.TMPDIR} -m {params.memory} -p DetectBeadSynthesisErrors\
		INPUT={input}\
		OUTPUT={output}\
		OUTPUT_STATS={params.out_stats}\
		SUMMARY={params.summary}\
		NUM_BARCODES={params.barcodes}\
		PRIMER_SEQUENCE={params.SmartAdapter}
		sleep 500
		"""

rule bam_hist:
	input: 'data/{sample}_final.bam'
	params:
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory = config['LOCAL']['MEMORY'],
		TMPDIR=config['LOCAL']['TMPDIR'],
		sge_opts = "-cwd -j y -l h_data=16G,h_rt=4:00:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs"
	output: 'logs/{sample}_hist_out_cell.txt'
	shell:
		"""{params.DROPSEQ_wrapper} -p BAMTagHistogram -m {params.memory} -t {params.TMPDIR}\
		TAG=XC\
		I={input}\
		O={output}
		sleep 500
		"""

rule plot_knee_plot:
	input: 'logs/{sample}_hist_out_cell.txt'
	output:
		pdf = 'plots/{sample}_knee_plot.pdf',
		png = 'plots/png/{sample}_knee_plot.png'
	params: 
		cells = lambda wildcards: int(samples.loc[wildcards.sample,'expected_cells']),
		sge_opts = "-cwd -j y -l h_data=2G,h_rt=0:10:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs"
	script:
		'../scripts/plot_knee_plot.R'

rule plot_knee_plot_whitelist:
	input:
		data = 'logs/{sample}_hist_out_cell.txt',
		barcodes = 'barcodes.csv'
	output:
		plot = 'plots/{sample}_knee_plot.pdf'
	params: 
		cells = lambda wildcards: int(samples.loc[wildcards.sample,'expected_cells']),
		sge_opts = "-cwd -j y -l h_data=2G,h_rt=0:10:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs"
	script:
		'../scripts/plot_knee_plot.R'